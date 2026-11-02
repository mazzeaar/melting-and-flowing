#include <cstddef>

#include "../solver.hpp"
#include "Eigen/Core"

namespace solver {

void Solver::solve_pressure(Config cfg) {
    using enum grid::Field;

    auto u = cfg.grid.view<U>();
    auto v = cfg.grid.view<V>();
    auto w = cfg.grid.view<W>();
    auto p = cfg.grid.view<P>();

    const auto& pcfg = pressure_cfg;

    // Initialize the rhs vector
    Eigen::VectorX<float_type> rhs(pcfg.n_unkowns());

    // Iterate thorugh u, v, w and store divergence in rhs
    for (std::size_t i = 1; i <= pcfg.nx; i++) {
        for (std::size_t j = 1; j <= pcfg.ny; j++) {
            for (std::size_t k = 1; k <= pcfg.nz; k++) {
                auto r = pcfg.row(i, j, k);
                auto div = (u[i, j, k] - u[i - 1, j, k]) +
                           (v[i, j, k] - v[i, j - 1, k]) +
                           (w[i, j, k] - w[i, j, k - 1]);
                rhs(r) = -div / pcfg.dx;
            }
        }
    }

    rhs = rhs * pcfg.dx * pcfg.dx / cfg.dt;
    Eigen::VectorX<float_type> p_new = pcfg.solve(rhs);

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
