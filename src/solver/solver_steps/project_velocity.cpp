#include "../solver.hpp"

namespace solver {

void Solver::project_velocity(Config cfg) {
    using enum grid::Field;

    auto u = cfg.grid.view<U>();
    auto v = cfg.grid.view<V>();
    auto w = cfg.grid.view<W>();
    auto p = cfg.grid.view<P>();

    const auto& pcfg = pressure_cfg;
    const auto rho = 1.0;  // Density

    for (std::size_t i = 1; i <= pcfg.nx; i++) {
        for (std::size_t j = 1; j <= pcfg.ny; j++) {
            for (std::size_t k = 1; k <= pcfg.nz; k++) {
                u[i, j, k] -=
                    cfg.dt / rho * (p[i + 1, j, k] - p[i, j, k]) / pcfg.dx;
                v[i, j, k] -=
                    cfg.dt / rho * (p[i, j + 1, k] - p[i, j, k]) / pcfg.dy;
                w[i, j, k] -=
                    cfg.dt / rho * (p[i, j, k + 1] - p[i, j, k]) / pcfg.dz;
            }
        }
    }
}
}  // namespace solver
