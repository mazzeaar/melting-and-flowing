#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include "../solver.hpp"

// Step 6: Temperature diffusion
// T_new - dt * kappa * ∇²T_new = T1
void Solver::diffuse_temperature(Config cfg) {
    using float_t = Solver::float_type;

    auto& grid = cfg.grid;
    const auto gc = grid.config();

    const auto ex = gc.extents();
    const auto h = gc.cell_size();

    const int nx = ex.x();
    const int ny = ex.y();
    const int nz = ex.z();

    // TODO: Replace with material thermal diffusivity.
    const auto dt = cfg.settings.dt;
    const auto thermal_diffusivity = cfg.physics.diffusivity;

    const float_t ax = dt * thermal_diffusivity / (h.x() * h.x());
    const float_t ay = dt * thermal_diffusivity / (h.y() * h.y());
    const float_t az = dt * thermal_diffusivity / (h.z() * h.z());

    temperature_rhs = grid.t;
    temperature_triplets.clear();

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const auto r = gc.cell(i, j, k);

                float_t diag = 1.0;

                const auto neighbour = [&](int ni, int nj, int nk, float_t a) {
                    temperature_triplets.emplace_back(r, gc.cell(ni, nj, nk),
                                                      -a);
                    diag += a;
                };

                // X neighbours
                if (i > 0) {
                    neighbour(i - 1, j, k, ax);
                }
                if (i + 1 < nx) {
                    neighbour(i + 1, j, k, ax);
                }

                // Y neighbours
                if (j > 0) {
                    neighbour(i, j - 1, k, ay);
                }
                if (j + 1 < ny) {
                    neighbour(i, j + 1, k, ay);
                }

                // Z neighbours
                if (k > 0) {
                    neighbour(i, j, k - 1, az);
                }
                if (k + 1 < nz) {
                    neighbour(i, j, k + 1, az);
                }

                temperature_triplets.emplace_back(r, r, diag);
            }
        }
    }

    temperature_A.setFromTriplets(temperature_triplets.begin(),
                                  temperature_triplets.end());

    temperature_solver.setMaxIterations(cfg.settings.max_solver_iterations);
    temperature_solver.setTolerance(cfg.settings.solver_tolerance);
    temperature_solver.compute(temperature_A);

    if (temperature_solver.info() != Eigen::Success) {
        throw std::runtime_error("Temperature diffusion setup failed.");
    }

    temperature_T2 = temperature_solver.solve(temperature_rhs);

    if (temperature_solver.info() != Eigen::Success) {
        throw std::runtime_error("Temperature diffusion failed.");
    }

    grid.t = temperature_T2;
}
