#include <cassert>

#include "../solver.hpp"

// Step 5: Temperature advection
// T1 = T_old - dt * (u_new · ∇)T_old
void Solver::advect_temperature(Solver::Config cfg) {
    auto& grid = cfg.grid;
    auto& temp = grid.t;

    const auto gc = grid.config();
    const auto ex = gc.extents();
    const auto h = gc.cell_size();
    const auto dt = cfg.settings.dt;
    const auto origin = gc.min_corner();

    assert(buff_mat_cell1.size() == gc.size());

    for (int k = 0; k < ex.z(); ++k) {
        for (int j = 0; j < ex.y(); ++j) {
            for (int i = 0; i < ex.x(); ++i) {
                const auto idx = gc.cell(i, j, k);

                if (grid.cell_state[idx] != CellState::FULL &&
                    grid.cell_state[idx] != CellState::SURFACE) {
                    buff_mat_cell1[idx] = temp[idx];
                    continue;
                }

                const Eigen::Vector3f cell{
                    static_cast<float>(i),
                    static_cast<float>(j),
                    static_cast<float>(k),
                };

                const Eigen::Vector3f pos =
                    origin +
                    (cell + Eigen::Vector3f::Constant(0.5f)).cwiseProduct(h);

                const auto velocity = interpolate_velocity(grid, pos);
                const Eigen::Vector3f prev_pos =
                    pos - static_cast<float>(dt) * velocity;

                const Eigen::Vector3f clamped =
                    prev_pos.cwiseMax(gc.min_corner())
                        .cwiseMin(gc.max_corner());

                const Eigen::Vector3f grid_pos =
                    (clamped - origin).cwiseQuotient(h);

                buff_mat_cell1[idx] = interpolate_field(
                    temp, ex, Eigen::Vector3f::Constant(0.5f), grid_pos);
            }
        }
    }

    temp = buff_mat_cell1;
}
