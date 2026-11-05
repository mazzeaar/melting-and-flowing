#include <cstddef>

#include "../solver.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCholesky"
#include "Eigen/SparseCore"
#include "core/exception.hpp"
#include "grid/core.hpp"

namespace solver {

/// Determines the pressure gradient to apply to velocity field to get
/// divergence free field Boundary conditions are taken from R. Bridson and M.
/// Müller-Fischer: "Fluid Simulation"
void Solver::solve_pressure(Config cfg) {
    using enum grid::Field;
    using enum grid::CellState;

    auto u = cfg.grid.view<U>();
    auto v = cfg.grid.view<V>();
    auto w = cfg.grid.view<W>();
    auto p = cfg.grid.view<P>();
    auto c = cfg.grid.view<State>();

    const auto& pcfg = pressure_cfg;

    std::vector<Eigen::Triplet<float_type>> triplets;
    triplets.reserve(7 * pcfg.n_unkowns() + 1);
    Eigen::VectorX<float_type> rhs(pcfg.n_unkowns());

    for (std::size_t i = 1; i <= pcfg.nx; i++) {
        for (std::size_t j = 1; j <= pcfg.ny; j++) {
            for (std::size_t k = 1; k <= pcfg.nz; k++) {
                auto r = pcfg.row(i, j, k);

                // This is the u_solid implicitely
                // Furthermore, we keep empty rows for easier mapping back
                if (c[i, j, k] != FULL && c[i, j, k] != SURFACE) {
                    triplets.emplace_back(r, r, 1.0);
                    rhs(r) = 0;
                    continue;
                }

                auto div = (u[i, j, k] - u[i - 1, j, k]) +
                           (v[i, j, k] - v[i, j - 1, k]) +
                           (w[i, j, k] - w[i, j, k - 1]);
                rhs(r) = -div / pcfg.dx;

                std::size_t n_fluid = 0;
                auto couple = [&](std::size_t ni, std::size_t nj,
                                  std::size_t nk) {
                    if (c[ni, nj, nk] == SOLID) return;
                    n_fluid++;
                    if (c[ni, nj, nk] == AIR) return;
                    triplets.emplace_back(r, pcfg.row(ni, nj, nk), -1);
                };

                if (i > 1) {
                    couple(i - 1, j, k);
                }
                if (i < pcfg.nx) {
                    couple(i + 1, j, k);
                }
                if (j > 1) {
                    couple(i, j - 1, k);
                }
                if (j < pcfg.ny) {
                    couple(i, j + 1, k);
                }
                if (k > 1) {
                    couple(i, j, k - 1);
                }
                if (k < pcfg.nz) {
                    couple(i, j, k + 1);
                }

                triplets.emplace_back(r, r, static_cast<double>(n_fluid));
            }
        }
    }

    Eigen::SparseMatrix<float_type> A(pcfg.n_unkowns(), pcfg.n_unkowns());
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<float_type>> chol(A);
    if (chol.info() != Eigen::Success)
        throw core::Error("pressure: factorization failed");

    rhs = rhs * pcfg.dx * pcfg.dx / cfg.dt;
    Eigen::VectorX<float_type> p_new = chol.solve(rhs);

    for (std::size_t i = 1; i <= pcfg.nx; ++i) {
        for (std::size_t j = 1; j <= pcfg.ny; ++j) {
            for (std::size_t k = 1; k <= pcfg.nz; ++k) {
                p[i, j, k] = p_new(pcfg.row(i, j, k));
            }
        }
    }
}
}  // namespace solver
