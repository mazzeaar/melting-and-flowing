#pragma once

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <memory>

#include "core/exception.hpp"
#include "core/paths.hpp"
#include "grid/grid.hpp"
#include "particles.hpp"
#include "polyscope/point_cloud.h"
#include "solver/solver.hpp"

struct Simulation {
    grid::Storage storage;
    solver::Solver solver;
    Particles particles;

    static constexpr auto default_mesh = "bunny.obj";
    std::filesystem::path mesh_path = core::paths::mesh(default_mesh);
    std::size_t step_count = 0;
    double time = 0.0;

    Simulation(auto grid_extents, auto grid_size, auto grid_anchor)
        : storage{grid::Config{grid_extents, grid_size, grid_anchor}},
          solver{} {}

    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) = delete;
    Simulation& operator=(Simulation&&) = delete;

    auto write_mesh(Mesh& mesh) { Log("TODO: implement grid to mesh"); }

    auto load_mesh(const std::unique_ptr<Mesh>& mesh) {
        if (!mesh) {
            throw core::Error("Simulation::load_mesh(): mesh is nullptr :(");
        }

        particles = Particles(Particles::Config{}, storage.config, mesh);

        grid().initialise_cell_states(*mesh);

        // re-initialise grid from mesh
        //  - no rescaling == no reallocations
        // maybe we should add rescaling as well?

        solver.prepare(solver::Solver::Config{
            .grid = grid(),
            .particles = Particles::PositionView{particles.pos()},
            .dt = 0.0,
            .max_steps = 1,
        });
    }

    auto grid() -> grid::GridView { return grid::GridView{storage}; }

    void step(double dt) {
        solver.step(solver::Solver::Config{
            .grid = grid(),
            .particles = Particles::PositionView{particles.pos()},
            .dt = dt,
            .max_steps = 1,
        });

        grid().update_cell_states(particles.pos());

        time += dt;
        ++step_count;
    }
};
