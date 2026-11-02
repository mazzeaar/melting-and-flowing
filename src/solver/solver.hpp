#pragma once

#include <chrono>
#include <core/log.hpp>
#include <cstddef>
#include <grid/grid.hpp>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Sparse"
#include "core/exception.hpp"
#include "grid/core.hpp"

namespace solver {

using Duration = std::chrono::milliseconds;

struct Solver {
    struct Config {
        GridView grid;
        double dt;
        std::size_t max_steps;

        // force vectors, constants later, ...
    };

    /// Discrete pressure Poisson setup
    ///
    /// Assumptions:
    /// - Dimensions (nx, ny, nz) count interior cells only, the ghost cells
    /// are not counted
    /// - The stencil is pure, meaning scaling information needs to be
    /// factorized into the RHS of the possion solve
    struct PressureConfig {
        using float_type = grid::Config::float_type;
        using Matrix = Eigen::SparseMatrix<float_type>;
        using Vector = Eigen::VectorX<float_type>;

        std::size_t nx{}, ny{}, nz{};
        double dx{}, dy{}, dz{};

        void prepare(const grid::Config& cfg) {
            const auto geo = grid::GeometryContext{cfg};
            const auto& ext = cfg.grid_extents;
            nx = ext.extent(0) - 2;
            ny = ext.extent(1) - 2;
            nz = ext.extent(2) - 2;
            dx = geo.cell_size[0];
            dy = geo.cell_size[1];
            dz = geo.cell_size[2];

            if (dx != dy || dx != dz) {
                throw core::Error(
                    "PressureConfig: Grid cells are not uniformly sized - not "
                    "supported right "
                    "now");
            }

            std::vector<Eigen::Triplet<float_type>> triplets;
            triplets.reserve(7 * n_unkowns() + 1);

            for (std::size_t i = 1; i <= nx; i++) {
                for (std::size_t j = 1; j <= ny; j++) {
                    for (std::size_t k = 1; k <= nz; k++) {
                        std::size_t n_fluid = 0;
                        auto r = row(i, j, k);

                        auto couple = [&](std::size_t ni, std::size_t nj,
                                          std::size_t nk) {
                            triplets.emplace_back(r, row(ni, nj, nk), -1);
                            n_fluid++;
                        };

                        if (i > 1) {
                            couple(i - 1, j, k);
                        }
                        if (i < nx) {
                            couple(i + 1, j, k);
                        }
                        if (j > 1) {
                            couple(i, j - 1, k);
                        }
                        if (j < ny) {
                            couple(i, j + 1, k);
                        }
                        if (k > 1) {
                            couple(i, j, k - 1);
                        }
                        if (k < nz) {
                            couple(i, j, k + 1);
                        }

                        // For example when we are on a left edge, we only use 4
                        // pressure values hence we count for this using n_fluid
                        triplets.emplace_back(r, r,
                                              static_cast<double>(n_fluid));
                    }
                }
            }

            // We need to pin a certain value, otherwise all rows of A sum to 0
            // and the trivial solution A * 1 =
            triplets.emplace_back(0, 0, 1.0);

            stencil.resize(n_unkowns(), n_unkowns());
            stencil.setFromTriplets(triplets.begin(), triplets.end());

            chol.compute(stencil);

            if (chol.info() != Eigen::Success) {
                throw core::Error(
                    "PressureConfig: Cholesky factorization failed");
            }
        }

        auto solve(const Vector& rhs) const -> Vector {
            return chol.solve(rhs);
        }

        std::size_t n_unkowns() const { return nx * ny * nz; }

        std::size_t row(std::size_t i, std::size_t j, std::size_t k) const {
            return (i - 1) + nx * ((j - 1) + ny * (k - 1));
        }

       private:
        Matrix stencil;
        Eigen::SimplicialCholesky<Matrix> chol;
    };

    // do we even need this?
    struct Result {
        // solver should own the buffers it needs and only provide spans
        std::span<double> some_result_buffer_maybe_idk;
    };

    struct Statistics {
        template <class T> class Statistic {
            std::vector<T> history;

           public:
            void clear() { history.clear(); }
            void record(T value) { history.push_back(value); }

            T latest() const { return history.back(); }
        };

        Statistic<Duration> duration;
        Statistic<std::size_t> iterations;
    };

    using float_type = grid::Config::float_type;

   private:
    friend struct PressureTest;

    Statistics stats;
    PressureConfig pressure_cfg;
    Eigen::SparseMatrix<float_type> some_solver;
    std::vector<float_type> some_pingpong_buffer;
    std::vector<float_type> advect_velocity_buffer;

   public:
    Solver() = default;

    // precomputations/conditioning/allocations go here
    // should only be called on resizing and first call?
    void prepare(Config cfg);

    const auto& get_stats() const { return stats; }

    // one full timestep
    // solves for du/dt
    auto step(Config cfg) -> Result;

   private:
    // ================
    // everything below is generated by chatgpt because im lazy.
    //
    // one cpp file per step?
    // ================

    // one full timestep, including all substeps
    // does not solve for du/dt !!!
    void single_step(Config cfg);

    // Step 1: Advection and external forces
    // u1 = u_old - dt * (u_old · ∇)u_old + dt * g
    void advect_velocity(Config cfg);

    // Step 2: Viscous diffusion
    // u2 - dt * nu * ∇²u2 = u1
    void diffuse_velocity(Config cfg) {}

    // Step 3: Pressure Poisson equation
    // ∇²p = (rho / dt) * ∇·u2
    void solve_pressure(Config cfg);

    // Step 4: Pressure projection
    // u_new = u2 - (dt / rho) * ∇p
    void project_velocity(Config cfg);

    // Step 5: Temperature advection
    // T1 = T_old - dt * (u_new · ∇)T_old
    void advect_temperature(Config cfg) {}

    // Step 6: Temperature diffusion
    // T_new - dt * kappa * ∇²T_new = T1
    void diffuse_temperature(Config cfg) {}

    // Step 7: Update viscosity from temperature
    // eta_new = eta(T_new)
    void update_viscosity(Config cfg) {}
};  // namespace detail

}  // namespace solver
