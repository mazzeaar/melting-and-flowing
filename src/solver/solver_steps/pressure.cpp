#include <cstddef>
#include <vector>

#include "../solver.hpp"
#include "Eigen/Core"
#include "Eigen/Sparse"

namespace solver {

void Solver::solve_pressure(Config cfg) {
    using enum grid::Field;
    typedef Eigen::Triplet<float_type> Triplet;

    auto u = cfg.grid.view<U>();
    auto v = cfg.grid.view<V>();
    auto w = cfg.grid.view<W>();
    auto p = cfg.grid.view<P>();

    const auto& pcfg = pressure_cfg;

    // Initialize the rhs vector
    Eigen::VectorX<float_type> rhs(pcfg.n_unkowns());
    Eigen::SparseMatrix<float_type> A(pcfg.n_unkowns(), pcfg.n_unkowns());
    std::vector<Triplet> tripletList;
    tripletList.reserve(7 * pcfg.n_unkowns() + 1);

    // For now we assume that the grid cells are uniform
    assert(pcfg.dx == pcfg.dy);
    assert(pcfg.dx == pcfg.dz);

    // Iterate thorugh u, v, w and store divergence in rhs
    for (std::size_t i = 1; i <= pcfg.nx; i++) {
        for (std::size_t j = 1; j <= pcfg.ny; j++) {
            for (std::size_t k = 1; k <= pcfg.nz; k++) {
                auto r = pcfg.row(i, j, k);

                // RHS
                auto div = (u[i, j, k] - u[i - 1, j, k]) +
                           (v[i, j, k] - v[i, j - 1, k]) +
                           (w[i, j, k] - w[i, j, k - 1]);
                rhs(r) = -div / pcfg.dx;

                // Stencil Matrix
                std::size_t n_fluid = 0;
                auto couple = [&](std::size_t ni, std::size_t nj,
                                  std::size_t nk) {
                    tripletList.emplace_back(r, pcfg.row(ni, nj, nk), -1);
                    n_fluid++;
                };

                // TODO: Whole stencil and cholesky prefactor can be done in
                // Solver::prepare() once, this is just clarity for now

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

                // For example when we are on a left edge, we only use 4
                // pressure values hence we count for this using n_fluid
                tripletList.emplace_back(
                    Triplet(r, r, static_cast<double>(n_fluid)));
            }
        }
    }

    // We need to pin a certain value, otherwise all rows of A sum to 0 and
    // the trivial solution A * 1 =
    tripletList.emplace_back(0, 0, 1.0);

    A.setFromTriplets(tripletList.begin(), tripletList.end());
    rhs = rhs * pcfg.dx * pcfg.dx / cfg.dt;

    // Use cholesky solver to solve for pressures
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<float_type>> chol(A);
    Eigen::VectorX<float_type> p_new = chol.solve(rhs);

    for (std::size_t i = 1; i <= pcfg.nx; ++i) {
        for (std::size_t j = 1; j <= pcfg.ny; ++j) {
            for (std::size_t k = 1; k <= pcfg.nz; ++k) {
                p[i, j, k] = p_new(pcfg.row(i, j, k));

                // Fill ghost values into palceholder so pressure gradient
                // in next step will automatically be 0 at every wall face
                p[0, j, k] = p[1, j, k];
                p[pcfg.nx + 1, j, k] = p[pcfg.nx, j, k];

                p[i, 0, k] = p[i, 1, k];
                p[i, pcfg.ny + 1, k] = p[i, pcfg.ny, k];

                p[i, j, 0] = p[i, j, 1];
                p[i, j, pcfg.nz + 1] = p[i, j, pcfg.nz];
            }
        }
    }
}
}  // namespace solver
