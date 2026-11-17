#pragma once

#include <Eigen/Core>
#include <cstddef>
#include <memory>

#include "mesh.hpp"
#include "settings.hpp"

struct Particles {
    Eigen::MatrixXf positions;

    Particles() = default;
    Particles(std::size_t num_particles, const std::unique_ptr<Mesh>& mesh_ptr);

    constexpr auto size() const -> std::size_t { return positions.rows(); }
};
