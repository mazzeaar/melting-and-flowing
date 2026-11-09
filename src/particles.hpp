#pragma once
#include <Eigen/Core>

#include "mesh/mesh.hpp"

// forward declaration
namespace grid {
struct Config;
}

struct Particles {
    using PositionMatrix =
        Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

    using PositionView =
        Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>, 0,
                   Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >;

    struct Config {
        std::size_t n_particles = 10000;
    };

   private:
    PositionMatrix positions;

   public:
    Particles() = default;
    Particles(Config conf, grid::Config grid_config,
              const std::unique_ptr<Mesh>& mesh_ptr);

    const PositionMatrix& pos() const { return positions; }
    PositionMatrix& pos() { return positions; }
};
