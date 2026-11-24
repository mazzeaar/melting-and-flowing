
#include <cstddef>

#include "../boundary_conditions.hpp"
#include "../solver.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCholesky"
#include "Eigen/SparseCore"
#include "core/exception.hpp"

/// Determines the pressure gradient to apply to velocity field to get
/// divergence free field Boundary conditions are taken from R. Bridson and M.
/// Müller-Fischer: "Fluid Simulation"
void Solver::solve_pressure(Config cfg) {
    using enum CellState;

    auto& u = cfg.grid.u;
    auto& v = cfg.grid.v;
    auto& w = cfg.grid.w;
    auto& p = cfg.grid.p;
    auto& c = cfg.grid.cell_state;

    const auto gc = cfg.grid.config();
    const auto ex = gc.extents();
    const auto dx = gc.cell_size().x();
    const auto dt = cfg.settings.dt;

    assert(dx == gc.cell_size().y() && dx == gc.cell_size().z());

    const int nx = ex.x();
    const int ny = ex.y();
    const int nz = ex.z();

    // TODO: move this buffer outside
    std::vector<Eigen::Triplet<float_type>> triplets;
    triplets.reserve(7 * gc.size() + 1);
    Eigen::VectorX<float_type> rhs(gc.size());

    std::vector<boundary::PressureCell> pinned(gc.size());
    boundary::pressure_reference_cells(cfg.grid, pinned);

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                auto r = gc.cell(i, j, k);

                // This is the u_solid implicitely
                // Furthermore, we keep empty rows for easier mapping back
                if (!boundary::is_liquid(c[r]) ||
                    pinned[r] == boundary::PressureCell::Pinned) {
                    triplets.emplace_back(r, r, 1.0);
                    rhs(r) = 0;
                    continue;
                }

                auto div = (u[gc.u(i + 2, j + 1, k + 1)] -
                            u[gc.u(i + 1, j + 1, k + 1)]) +
                           (v[gc.v(i + 1, j + 2, k + 1)] -
                            v[gc.v(i + 1, j + 1, k + 1)]) +
                           (w[gc.w(i + 1, j + 1, k + 2)] -
                            w[gc.w(i + 1, j + 1, k + 1)]);

                rhs(r) = -div / dx;

                std::size_t n_fluid = 0;
                auto couple = [&](int ni, int nj, int nk) {
                    const auto n = gc.cell(ni, nj, nk);
                    const auto state =
                        boundary::cell_state(cfg.grid, ni, nj, nk);
                    if (state == CellState::SOLID) {
                        return;
                    }

                    n_fluid++;
                    if (state == CellState::AIR ||
                        pinned[n] == boundary::PressureCell::Pinned) {
                        return;
                    }
                    triplets.emplace_back(r, n, -1);
                };

                if (i > 0) {
                    couple(i - 1, j, k);
                }
                if (i + 1 < nx) {
                    couple(i + 1, j, k);
                }
                if (j > 0) {
                    couple(i, j - 1, k);
                }
                if (j + 1 < ny) {
                    couple(i, j + 1, k);
                }
                if (k > 0) {
                    couple(i, j, k - 1);
                }
                if (k + 1 < nz) {
                    couple(i, j, k + 1);
                }

                triplets.emplace_back(r, r, static_cast<double>(n_fluid));
            }
        }
    }

    Eigen::SparseMatrix<float_type> A(gc.size(), gc.size());
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<float_type>> chol(A);
    if (chol.info() != Eigen::Success) {
        throw core::Error("pressure: factorization failed");
    }

    rhs = rhs * dx * dx / dt;
    p = chol.solve(rhs);
}
