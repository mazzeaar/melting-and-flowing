#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include "../solver.hpp"

// Step 6: Temperature diffusion
// T_new - dt * kappa * ∇²T_new = T1
void solver::Solver::diffuse_temperature(Config cfg) {
    using float_t = Solver::float_type;

    auto temp = cfg.grid.view<grid::Field::T>();
    const auto& gc = cfg.grid.config();

    const auto ex = temp.extents();
    const std::size_t nx = ex.extent(0);
    const std::size_t ny = ex.extent(1);
    const std::size_t nz = ex.extent(2);

    const float_t dx = gc.grid_size[0] / nx;
    const float_t dy = gc.grid_size[1] / ny;
    const float_t dz = gc.grid_size[2] / nz;

    // TODO: replace with material thermal diffusivity
    const float_t kappa = 1e-6;

    const float_t ax = cfg.dt * kappa / (dx * dx);
    const float_t ay = cfg.dt * kappa / (dy * dy);
    const float_t az = cfg.dt * kappa / (dz * dz);

    const auto row = [ny, nz](std::size_t i, std::size_t j,
                              std::size_t k) -> Eigen::Index {
        return (i * ny + j) * nz + k;
    };

    temperature_triplets.clear();

    for (std::size_t i = 0; i < nx; ++i) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t k = 0; k < nz; ++k) {
                const auto r = row(i, j, k);

                temperature_rhs(r) = temp[i, j, k];

                float_t diag = 1.0;

                const auto neighbour = [&](std::size_t ni, std::size_t nj,
                                           std::size_t nk, float_t a) {
                    temperature_triplets.emplace_back(r, row(ni, nj, nk), -a);
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

    temperature_solver.compute(temperature_A);

    if (temperature_solver.info() != Eigen::Success)
        throw std::runtime_error("Temperature diffusion setup failed.");

    temperature_T2 = temperature_solver.solve(temperature_rhs);

    if (temperature_solver.info() != Eigen::Success)
        throw std::runtime_error("Temperature diffusion failed.");

    // Write solutions back to grid.
    for (std::size_t i = 0; i < nx; ++i) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t k = 0; k < nz; ++k) {
                temp[i, j, k] = temperature_T2(row(i, j, k));
            }
        }
    }
}
