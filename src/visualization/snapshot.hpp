#pragma once

#include <Eigen/Core>
#include <cstddef>
#include <memory>
#include <simulation/grid.hpp>

struct GridSnapshot {
    Grid::Config config;

    Eigen::VectorXf pressure;
    Eigen::VectorXf temperature;
    Eigen::VectorXf viscosity;
    Eigen::VectorXi cell_state;
};

struct VelocitySnapshot {
    Eigen::MatrixXf vectors;  // N x 3, cell-centered
};

struct ParticleSnapshot {
    Eigen::MatrixXf positions;  // N x 3
};

struct MeshSnapshot {
    std::size_t generation{};
    std::size_t step{};
    Eigen::MatrixXf vertices;
    Eigen::MatrixXi faces;
};

struct SimulationSnapshot {
    // Changes on reset/load, so results from the old simulation can be ignored.
    std::size_t generation{};
    std::size_t step{};
    double time{};

    GridSnapshot grid;
    VelocitySnapshot velocity;
    ParticleSnapshot particles;
    // Loaded once per reset; keeps the optional original-mesh display
    // available.
    std::shared_ptr<const MeshSnapshot> source_mesh;
};
