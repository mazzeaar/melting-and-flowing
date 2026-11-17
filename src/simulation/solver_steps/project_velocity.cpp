#include "../solver.hpp"

void Solver::project_velocity(Config cfg) {
    auto& u = cfg.grid.u;
    auto& v = cfg.grid.v;
    auto& w = cfg.grid.w;
    auto& p = cfg.grid.p;

    const auto gc = cfg.grid.config();
    const auto h = gc.cell_size();
    const auto rho = 1.0;  // Density

    const auto ex = gc.extents();
    const auto nx = ex.x(), ny = ex.y(), nz = ex.z();
    const auto s = cfg.settings.dt / rho;

    // Interior faces: each component skips its own wall face at index n
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                const auto pc = p[gc.cell(i, j, k)];

                if (i + 1 < nx) {
                    u[gc.u(i + 2, j + 1, k + 1)] -=
                        s * (p[gc.cell(i + 1, j, k)] - pc) / h.x();
                }

                if (j + 1 < ny) {
                    v[gc.v(i + 1, j + 2, k + 1)] -=
                        s * (p[gc.cell(i, j + 1, k)] - pc) / h.y();
                }

                if (k + 1 < nz) {
                    w[gc.w(i + 1, j + 1, k + 2)] -=
                        s * (p[gc.cell(i, j, k + 1)] - pc) / h.z();
                }
            }
        }
    }

    // Wall faces: u_solid = 0
    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            u[gc.u(1, j + 1, k + 1)] = u[gc.u(nx + 1, j + 1, k + 1)] = 0;
        }
    }

    for (int i = 0; i < nx; i++) {
        for (int k = 0; k < nz; k++) {
            v[gc.v(i + 1, 1, k + 1)] = v[gc.v(i + 1, ny + 1, k + 1)] = 0;
        }
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            w[gc.w(i + 1, j + 1, 1)] = w[gc.w(i + 1, j + 1, nz + 1)] = 0;
        }
    }
}
