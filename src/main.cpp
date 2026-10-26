#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>

#include <Eigen/Dense>
#include <array>
#include <core/paths.hpp>
#include <grid/example.hpp>
#include <grid/grid.hpp>
#include <print>
#include <solver/solver.hpp>
#include <vector>

#include "grid/core.hpp"
#include "grid/storage.hpp"

int main() {
    auto grid_extents = std::array{20, 20, 20};
    auto grid_size = std::array{20.0, 20.0, 20.0};
    auto grid_anchor = std::array{0.0, 0.0, 0.0};

    auto grid_storage = grid::Storage{
        grid::Config(grid_extents, grid_size, grid_anchor),
    };

    auto grid_view = grid::GridView{grid_storage};

    auto solver_config = solver::Solver::Config{
        .grid = grid_view,
        .dt = 1.0,
        .max_steps = 1000,
    };

    auto solver = solver::Solver{};

    solver.prepare(solver_config);

    for (int i = 0; i < 10; ++i) {
        solver.step(solver_config);
    }
}
