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

    const auto nx = pcfg.nx, ny = pcfg.ny, nz = pcfg.nz;
    const auto s = cfg.dt / rho;

    // Interior faces: each component skips its own wall face at index n
    for (std::size_t i = 1; i <= nx; i++) {
        for (std::size_t j = 1; j <= ny; j++) {
            for (std::size_t k = 1; k <= nz; k++) {
                const auto pc = p[i, j, k];
                if (i < nx) u[i, j, k] -= s * (p[i + 1, j, k] - pc) / pcfg.dx;
                if (j < ny) v[i, j, k] -= s * (p[i, j + 1, k] - pc) / pcfg.dy;
                if (k < nz) w[i, j, k] -= s * (p[i, j, k + 1] - pc) / pcfg.dz;
            }
        }
    }

    // Wall faces: u_solid = 0
    for (std::size_t j = 1; j <= ny; j++) {
        for (std::size_t k = 1; k <= nz; k++) {
            u[0, j, k] = u[nx, j, k] = 0;
        }
    }

    for (std::size_t i = 1; i <= nx; i++) {
        for (std::size_t k = 1; k <= nz; k++) {
            v[i, 0, k] = v[i, ny, k] = 0;
        }
    }

    for (std::size_t i = 1; i <= nx; i++) {
        for (std::size_t j = 1; j <= ny; j++) {
            w[i, j, 0] = w[i, j, nz] = 0;
        }
    }
}
}  // namespace solver
