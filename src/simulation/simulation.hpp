#pragma once

#include <cstddef>
#include <filesystem>
#include <memory>

#include "core/exception.hpp"
#include "core/paths.hpp"
#include "grid.hpp"
#include "particles.hpp"
#include "settings.hpp"
#include "solver.hpp"

struct Simulation {
    SimulationSettings settings;
    Grid grid;
    Particles particles;
    Solver solver;

    static constexpr auto default_mesh = "bunny.obj";
    std::filesystem::path mesh_path = core::paths::mesh(default_mesh);
    std::size_t step_count = 0;
    double time = 0.0;

    Simulation(SimulationSettings settings)
        : settings{settings}, grid{settings.grid, settings.physics} {}

    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) = delete;
    Simulation& operator=(Simulation&&) = delete;

    auto write_mesh(Mesh& mesh) { Log("TODO: implement grid to mesh"); }

    auto load_mesh(const std::unique_ptr<Mesh>& mesh) {
        if (!mesh) {
            throw core::Error("Simulation::load_mesh(): mesh is nullptr :(");
        }

        const Grid::Config grid_config{settings.grid};
        mesh->fit_to_box(grid_config.min_corner(), grid_config.max_corner(),
                         settings.initialization.mesh_padding);
        // Reset all fields as well as particles, using the current settings.
        grid.resize(grid_config, settings.initialization.initial_temperature,
                    settings.physics.reference_viscosity);

        particles = Particles(settings.initialization.particle_count, mesh);
        mesh_path = mesh->file_path;

        time = 0.0;
        step_count = 0;

        initialise_cell_states(grid, mesh.get(),
                               settings.initialization.samples_per_cell_axis);

        solver.prepare(Solver::Config{
            .grid = grid,
            .particles = particles.positions,
            .settings = settings.solver,
            .physics = settings.physics,
        });
    }

    void step(double dt) {
        if (settings.solver.max_steps != 0 &&
            step_count >= settings.solver.max_steps) {
            return;
        }
        settings.solver.dt = dt;

        solver.step(Solver::Config{
            .grid = grid,
            .particles = particles.positions,
            .settings = settings.solver,
            .physics = settings.physics,
        });

        update_cell_states(grid, particles.positions);
        // chaddy said to add this, idk if its needed?
        enforce_velocity_boundaries(grid);

        time += dt;
        ++step_count;
    }
};
