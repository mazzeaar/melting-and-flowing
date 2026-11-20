#pragma once

#include <Eigen/Core>

#include "core/benchmark.hpp"
#include "grid.hpp"

struct SurfaceReconstructor {
    int resolution_upfactor = 2;
    float radius_factor = 1.5f;
    benchmark::BenchmarkBase benchmark;

    void reconstruct(const Eigen::MatrixXf& points, const Grid::Config& config,
                     Eigen::MatrixXf& V, Eigen::MatrixXi& F);

   private:
    Eigen::VectorXf node_positions;
    Eigen::VectorXf phi;
};
