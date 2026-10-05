#include "source_distribution.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "Randomize.hh"

namespace py = pybind11;

namespace g4bridge
{

namespace
{

void RequireFinite(double value, const char* name)
{
  if (!std::isfinite(value)) {
    std::ostringstream msg;
    msg << name << " must be finite.";
    throw std::runtime_error(msg.str());
  }
}

std::vector<double> Copy1D(const py::array_t<double, py::array::c_style | py::array::forcecast>& array,
                           const char* name)
{
  const py::buffer_info info = array.request();
  if (info.ndim != 1) {
    std::ostringstream msg;
    msg << name << " must be a 1D NumPy array.";
    throw std::runtime_error(msg.str());
  }

  const auto view = array.unchecked<1>();
  std::vector<double> values(static_cast<std::size_t>(info.shape[0]));
  for (py::ssize_t i = 0; i < info.shape[0]; ++i) {
    values[static_cast<std::size_t>(i)] = view(i);
    RequireFinite(values[static_cast<std::size_t>(i)], name);
  }
  return values;
}

void RequireIncreasingEdges(const std::vector<double>& edges, const char* name)
{
  if (edges.size() < 2) {
    std::ostringstream msg;
    msg << name << " must contain at least two edges.";
    throw std::runtime_error(msg.str());
  }
  for (std::size_t i = 1; i < edges.size(); ++i) {
    if (!(edges[i] > edges[i - 1])) {
      std::ostringstream msg;
      msg << name << " must be strictly increasing.";
      throw std::runtime_error(msg.str());
    }
  }
}

}  // namespace

void SourceComponent::Load(
  const py::array_t<double, py::array::c_style | py::array::forcecast>& box_bounds_mm,
  const py::array_t<double, py::array::c_style | py::array::forcecast>& mu_edges,
  const py::array_t<double, py::array::c_style | py::array::forcecast>& azi_edges,
  const py::array_t<double, py::array::c_style | py::array::forcecast>& energy_edges_mev,
  const py::array_t<double, py::array::c_style | py::array::forcecast>& weights,
  std::size_t n_u,
  std::size_t n_v,
  std::size_t n_events,
  int particle_id)
{
  if (n_events == 0) {
    throw std::runtime_error("source distribution n_events must be > 0.");
  }
  if (n_u == 0 || n_v == 0) {
    throw std::runtime_error("source distribution surface mesh dimensions must be > 0.");
  }

  // load source box bounds
  const py::buffer_info box_info = box_bounds_mm.request();
  if (box_info.ndim != 2 || box_info.shape[0] != 3 || box_info.shape[1] != 2) {
    throw std::runtime_error("box_bounds_mm must have shape (3, 2).");
  }
  const auto box = box_bounds_mm.unchecked<2>();
  for (py::ssize_t axis = 0; axis < 3; ++axis) {
    const double lower = box(axis, 0);
    const double upper = box(axis, 1);
    RequireFinite(lower, "box_bounds_mm");
    RequireFinite(upper, "box_bounds_mm");
    if (!(upper > lower)) {
      throw std::runtime_error("box_bounds_mm entries must satisfy min < max.");
    }
    box_bounds_mm_[static_cast<std::size_t>(axis)] = {lower, upper};
  }

  // load source distribution bin edges
  mu_edges_ = Copy1D(mu_edges, "mu_edges");
  azi_edges_ = Copy1D(azi_edges, "azi_edges");
  energy_edges_mev_ = Copy1D(energy_edges_mev, "energy_edges_mev");
  RequireIncreasingEdges(mu_edges_, "mu_edges");
  RequireIncreasingEdges(azi_edges_, "azi_edges");
  RequireIncreasingEdges(energy_edges_mev_, "energy_edges_mev");

  n_mu_ = mu_edges_.size() - 1;
  n_azi_ = azi_edges_.size() - 1;
  n_energy_ = energy_edges_mev_.size() - 1;
  n_u_ = n_u;
  n_v_ = n_v;
  const std::size_t expected_weight_count = n_mu_ * n_azi_ * n_energy_ * 6 * n_u_ * n_v_;

  const py::buffer_info weight_info = weights.request();
  if (weight_info.ndim != 1 ||
      static_cast<std::size_t>(weight_info.shape[0]) != expected_weight_count) {
    std::ostringstream msg;
    msg << "weights must be a flat array with " << expected_weight_count << " entries.";
    throw std::runtime_error(msg.str());
  }

  // build cumulative weights for source bin sampling
  const auto weight_view = weights.unchecked<1>();
  cdf_.clear();
  cdf_.reserve(expected_weight_count);
  total_weight_ = 0.0;
  for (py::ssize_t i = 0; i < weight_info.shape[0]; ++i) {
    const double weight = weight_view(i);
    RequireFinite(weight, "weights");
    if (weight < 0.0) {
      throw std::runtime_error("weights must be non-negative.");
    }
    total_weight_ += weight;
    cdf_.push_back(total_weight_);
  }
  if (!(total_weight_ > 0.0)) {
    throw std::runtime_error("weights must contain positive total probability.");
  }

  particle_id_ = particle_id;
  n_events_ = n_events;
}

Primary SourceComponent::Sample() const
{
  // sample source distribution bin
  const double source_pick = G4UniformRand() * total_weight_;
  auto source_it = std::lower_bound(cdf_.begin(), cdf_.end(), source_pick);
  std::size_t source_idx = static_cast<std::size_t>(std::distance(cdf_.begin(), source_it));

  // Python flattens weights from shape (mu, azi, energy, face, u, v) in C order.
  const std::size_t i_v = source_idx % n_v_;
  source_idx /= n_v_;
  const std::size_t i_u = source_idx % n_u_;
  source_idx /= n_u_;
  const std::size_t face = source_idx % 6;
  source_idx /= 6;
  const std::size_t i_energy = source_idx % n_energy_;
  source_idx /= n_energy_;
  const std::size_t i_azi = source_idx % n_azi_;
  const std::size_t i_mu = source_idx / n_azi_;

  // initialize particle values
  Primary p{};
  p.particle_id = particle_id_;
  p.weight = total_weight_ / static_cast<double>(n_events_);
  p.time_ns = 0.0;

  // sample particle energy and angle
  const double mu = mu_edges_[i_mu] + (mu_edges_[i_mu + 1] - mu_edges_[i_mu]) * G4UniformRand();
  const double azi =
    azi_edges_[i_azi] + (azi_edges_[i_azi + 1] - azi_edges_[i_azi]) * G4UniformRand();
  const double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  p.ux = sin_theta * std::cos(azi);
  p.uy = sin_theta * std::sin(azi);
  p.uz = mu;
  p.energy_mev = energy_edges_mev_[i_energy] +
                 (energy_edges_mev_[i_energy + 1] - energy_edges_mev_[i_energy]) *
                   G4UniformRand();

  const auto sample_bin = [](double lower, double upper, std::size_t i, std::size_t n) {
    const double width = (upper - lower) / static_cast<double>(n);
    return lower + (static_cast<double>(i) + G4UniformRand()) * width;
  };

  // face order and local axes match MCDC surface_mesh: x:y/z, y:x/z, z:x/y.
  if (face == 0) {
    p.x_mm = box_bounds_mm_[0][0];
    p.y_mm = sample_bin(box_bounds_mm_[1][0], box_bounds_mm_[1][1], i_u, n_u_);
    p.z_mm = sample_bin(box_bounds_mm_[2][0], box_bounds_mm_[2][1], i_v, n_v_);
  } else if (face == 1) {
    p.x_mm = box_bounds_mm_[0][1];
    p.y_mm = sample_bin(box_bounds_mm_[1][0], box_bounds_mm_[1][1], i_u, n_u_);
    p.z_mm = sample_bin(box_bounds_mm_[2][0], box_bounds_mm_[2][1], i_v, n_v_);
  } else if (face == 2) {
    p.y_mm = box_bounds_mm_[1][0];
    p.x_mm = sample_bin(box_bounds_mm_[0][0], box_bounds_mm_[0][1], i_u, n_u_);
    p.z_mm = sample_bin(box_bounds_mm_[2][0], box_bounds_mm_[2][1], i_v, n_v_);
  } else if (face == 3) {
    p.y_mm = box_bounds_mm_[1][1];
    p.x_mm = sample_bin(box_bounds_mm_[0][0], box_bounds_mm_[0][1], i_u, n_u_);
    p.z_mm = sample_bin(box_bounds_mm_[2][0], box_bounds_mm_[2][1], i_v, n_v_);
  } else if (face == 4) {
    p.z_mm = box_bounds_mm_[2][0];
    p.x_mm = sample_bin(box_bounds_mm_[0][0], box_bounds_mm_[0][1], i_u, n_u_);
    p.y_mm = sample_bin(box_bounds_mm_[1][0], box_bounds_mm_[1][1], i_v, n_v_);
  } else {
    p.z_mm = box_bounds_mm_[2][1];
    p.x_mm = sample_bin(box_bounds_mm_[0][0], box_bounds_mm_[0][1], i_u, n_u_);
    p.y_mm = sample_bin(box_bounds_mm_[1][0], box_bounds_mm_[1][1], i_v, n_v_);
  }

  return p;
}

void SourceDistribution::Load(
  const py::array_t<double, py::array::c_style | py::array::forcecast>& box_bounds_mm,
  const py::array_t<double, py::array::c_style | py::array::forcecast>& mu_edges,
  const py::array_t<double, py::array::c_style | py::array::forcecast>& azi_edges,
  const py::array_t<double, py::array::c_style | py::array::forcecast>& energy_edges_mev,
  const py::array_t<double, py::array::c_style | py::array::forcecast>& weights,
  std::size_t n_u,
  std::size_t n_v,
  std::size_t n_events,
  int particle_id)
{
  SourceComponent component;
  component.Load(
    box_bounds_mm, mu_edges, azi_edges, energy_edges_mev, weights, n_u, n_v, n_events,
    particle_id);

  n_events_total_ += component.Size();
  event_ends_.push_back(n_events_total_);
  components_.push_back(std::move(component));
}

void SourceDistribution::Clear()
{
  components_.clear();
  event_ends_.clear();
  n_events_total_ = 0;
}

Primary SourceDistribution::Sample(std::size_t event_id) const
{
  // event ranges are assigned to components in load order
  const auto it = std::upper_bound(event_ends_.begin(), event_ends_.end(), event_id);
  if (it == event_ends_.end()) {
    std::ostringstream msg;
    msg << "Event " << event_id << " is outside the " << n_events_total_
        << " loaded source-distribution events.";
    throw std::runtime_error(msg.str());
  }
  return components_[static_cast<std::size_t>(std::distance(event_ends_.begin(), it))].Sample();
}

}  // namespace g4bridge
