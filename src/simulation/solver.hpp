#pragma once

#include <Eigen/Sparse>
#include <chrono>
#include <core/log.hpp>
#include <cstddef>
#include <vector>

#include "grid.hpp"
#include "settings.hpp"

using Duration = std::chrono::milliseconds;

struct Solver {
    struct Config {
        Grid& grid;
        Eigen::MatrixXf& particles;
        SimulationSettings::SolverSettings settings;
        const SimulationSettings::PhysicsSettings& physics;
    };

    struct Statistics {
        template <class T> class Statistic {
            std::vector<T> history;

           public:
            void clear() { history.clear(); }
            void record(T value) { history.push_back(value); }

            T latest() const { return history.back(); }
        };

        Statistic<Duration> step_duration;
        Statistic<std::size_t> iterations;
        Statistic<std::size_t> cfl_substeps;
    };

    using float_type = float;

   private:
    friend struct PressureTest;

    Statistics stats;

    Eigen::VectorXf buff_mat_u1, buff_mat_v1, buff_mat_w1;
    Eigen::VectorXf buff_mat_u2, buff_mat_v2, buff_mat_w2;
    Eigen::VectorXf buff_mat_cell1, buff_mat_cell2;
    Eigen::VectorXi buff_mat_cell_state;

    std::vector<Eigen::Triplet<float_type>> trip_u, trip_v, trip_w;
    std::vector<Eigen::Triplet<float_type>> trip_cell;

    // Viscous diffusion buffers
    // Eigen::VectorXf viscosity_rhs_u, viscosity_rhs_v, viscosity_rhs_w;
    // Eigen::VectorXf viscosity_u2, viscosity_v2, viscosity_w2;
    // std::vector<Eigen::Triplet<float_type>> viscosity_triplets_u,
    //     viscosity_triplets_v, viscosity_triplets_w;
    Eigen::SparseMatrix<float_type> viscosity_A_u, viscosity_A_v, viscosity_A_w;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float_type>>
        viscosity_solver_u, viscosity_solver_v, viscosity_solver_w;

    Eigen::SparseMatrix<float_type> temperature_A;

    std::vector<Eigen::Triplet<float_type>> temperature_triplets;

    Eigen::ConjugateGradient<Eigen::SparseMatrix<float_type>,
                             Eigen::Lower | Eigen::Upper>
        temperature_solver;

   public:
    Solver() = default;

    // precomputations/conditioning/allocations go here
    // should only be called on resizing and first call?
    void prepare(Config cfg);

    const auto& get_stats() const { return stats; }

    // one full timestep
    // solves for du/dt
    void step(Config cfg);

   private:
    // one full timestep, including all substeps
    // does not solve for du/dt !!!
    void single_step(Config cfg);

    // Step 1: Advection and external forces
    // u1 = u_old - dt * (u_old · ∇)u_old + dt * g
    void advect_velocity(Config cfg);

    // Step 2: Viscous diffusion
    // u2 - dt * nu * ∇²u2 = u1
    void diffuse_velocity(Config cfg);

    // Step 3: Pressure Poisson equation
    // ∇²p = (rho / dt) * ∇·u2
    void solve_pressure(Config cfg);

    // Step 4: Pressure projection
    // u_new = u2 - (dt / rho) * ∇p
    void project_velocity(Config cfg);

    // Step 4.1: project particles to cell state
    void project_particles(Config cfg);

    // Step 5: Temperature advection
    // T1 = T_old - dt * (u_new · ∇)T_old
    void advect_temperature(Config cfg);

    // Step 6: Temperature diffusion
    // T_new - dt * kappa * ∇²T_new = T1
    void diffuse_temperature(Config cfg);

    // Step 7: Update viscosity from temperature
    // eta_new = eta(T_new)
   public:
    void update_viscosity(Config cfg);

   private:
    // Step 8: update particle positions
    void advect_particles(Config cfg);
};  // namespace detail
