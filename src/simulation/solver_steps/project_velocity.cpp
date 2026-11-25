#include "../boundary_conditions.hpp"
#include "../solver.hpp"

auto is_face_state_full(auto di, auto dj, auto dk, auto axis) {
    return [=](const auto& grid, auto i, auto j, auto k) {
        return boundary::face_state(grid, axis, i + di, j + dj, k + dk) ==
               CellState::FULL;
    };
}

void Solver::project_velocity(Config cfg) {
    auto& u = cfg.grid.u;
    auto& v = cfg.grid.v;
    auto& w = cfg.grid.w;
    auto& p = cfg.grid.p;

    const auto& grid = cfg.grid;
    const auto gc = grid.config();
    const auto h = gc.cell_size();
    const auto rho = 1.0;  // Density

    const auto ex = gc.extents();
    const auto nx = ex.x(), ny = ex.y(), nz = ex.z();
    const auto s = cfg.settings.dt / rho;

    auto is_u_face_full = is_face_state_full(2, 1, 1, 0);
    auto is_v_face_full = is_face_state_full(1, 2, 1, 1);
    auto is_w_face_full = is_face_state_full(1, 1, 2, 2);

    // Interior faces: each component skips its own wall face at index n
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                const auto pc = p[gc.cell(i, j, k)];

                if (i + 1 < nx && is_u_face_full(grid, i, j, k)) {
                    u[gc.u(i + 2, j + 1, k + 1)] -=
                        s * (p[gc.cell(i + 1, j, k)] - pc) / h.x();
                }

                if (j + 1 < ny && is_v_face_full(grid, i, j, k)) {
                    v[gc.v(i + 1, j + 2, k + 1)] -=
                        s * (p[gc.cell(i, j + 1, k)] - pc) / h.y();
                }

                if (k + 1 < nz && is_w_face_full(grid, i, j, k)) {
                    w[gc.w(i + 1, j + 1, k + 2)] -=
                        s * (p[gc.cell(i, j, k + 1)] - pc) / h.z();
                }
            }
        }
    }

    // Wall faces: u_solid = 0
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            u[gc.u(1, j + 1, k + 1)] = u[gc.u(nx + 1, j + 1, k + 1)] = 0;
        }
    }

    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            v[gc.v(i + 1, 1, k + 1)] = v[gc.v(i + 1, ny + 1, k + 1)] = 0;
        }
    }

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            w[gc.w(i + 1, j + 1, 1)] = w[gc.w(i + 1, j + 1, nz + 1)] = 0;
        }
    }
}
