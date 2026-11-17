#pragma once

#include <Eigen/Core>
#include <cstddef>

struct SimulationSettings {
    struct GridSettings {
        Eigen::Vector3f origin = Eigen::Vector3f::Zero();
        Eigen::Vector3i cell_counts = Eigen::Vector3i::Constant(32);
        float cell_size = 1.0f / 32.0f;
    } grid;

    struct SolverSettings {
        float dt = 1.0f / 60.0f;
        float cfl = 0.5f;
        std::size_t max_substeps = 1024;  // for eigen solvers
        std::size_t max_steps = 0;

        int max_solver_iterations = 1000;
        float solver_tolerance = 1e-6f;
    } solver;

    struct PhysicsSettings {
        float gravity = -9.81f;
        float density = 1000.0f;
        float diffusivity = 1.4e-7f;

        float reference_temperature = 293.15f;
        float reference_viscosity = 1e-6f;
        float thermal_sensitivity = 0.025f;
        float min_viscosity = 1e-7f;
        float max_viscosity = 1e-2f;
    } physics;

    struct InitializationSettings {
        std::size_t particle_count = 20'000;
        float initial_temperature = 293.15f;
        float mesh_padding = 0.1f;
        int samples_per_cell_axis = 4;
    } initialization;
};
