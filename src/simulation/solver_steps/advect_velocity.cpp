#include <utility>

#include "../boundary_conditions.hpp"
#include "../solver.hpp"
#include "Eigen/Core"

// Step 1: Advection and external forces
// old: (eulerian advection)
// u1 = u_old - dt * (u_old · ∇)u_old + dt * g
// what i did: (semi-lagrangian advection)
// u*(x) = u_old(x - dt * u_old(x)) + dt * g(x)
void Solver::advect_velocity(Config cfg) {
    auto& grid = cfg.grid;
    const auto gcfg = grid.config();

    const auto dt = cfg.settings.dt;
    const auto h = gcfg.cell_size();
    const auto origin = gcfg.min_corner();

    const auto gravity = Eigen::Vector3f{0.0f, cfg.physics.gravity, 0.0f};

    enum class Axis { X, Y, Z };

    // this is slow
    // advect one staggered velocity component.
    auto advect_component = [&](const Eigen::VectorXf& field,
                                Eigen::VectorXf& buffer,
                                Eigen::Vector3i extents, Eigen::Vector3f offset,
                                auto&& index_func, Axis axis) {

#pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < extents.z(); ++k) {
            for (int j = 0; j < extents.y(); ++j) {
                for (int i = 0; i < extents.x(); ++i) {
                    const auto idx = index_func(i, j, k);

                    // only faces touching liquid are advected.
                    // walls & ghosts are enforced by
                    // enforce_velocity_boundaries
                    if (boundary::face_state(grid, std::to_underlying(axis), i,
                                             j, k) != CellState::FULL) {
                        buffer[idx] = field[idx];
                        continue;
                    }

                    const Eigen::Vector3f cell{
                        static_cast<float>(i),
                        static_cast<float>(j),
                        static_cast<float>(k),
                    };

                    // staggered field position in physical coordinates
                    const auto pos = origin + (cell + offset).cwiseProduct(h);

                    // this is slow
                    // backtrace through the old velocity field
                    const auto velocity = interpolate_velocity(grid, pos);

                    const auto prev_pos = pos - dt * velocity;

                    // clamp to the physical domain
                    const Eigen::Vector3f clamped =
                        prev_pos.cwiseMax(gcfg.min_corner())
                            .cwiseMin(gcfg.max_corner());

                    // convert physical coordinates to cell coordinates
                    const Eigen::Vector3f grid_pos =
                        (clamped - origin).cwiseQuotient(h);

                    const float advected =
                        interpolate_field(field, extents, offset, grid_pos);

                    buffer[idx] =
                        advected + dt * gravity[std::to_underlying(axis)];
                }
            }
        }
    };

    const Eigen::Vector3f offset_u{-1.0f, -0.5f, -0.5f};
    const Eigen::Vector3f offset_v{-0.5f, -1.0f, -0.5f};
    const Eigen::Vector3f offset_w{-0.5f, -0.5f, -1.0f};

    // clang-format off
    const auto idx_u = [&](auto&&... is) { return gcfg.u(is...); };
    const auto idx_v = [&](auto&&... is) { return gcfg.v(is...); };
    const auto idx_w = [&](auto&&... is) { return gcfg.w(is...); };

    advect_component(grid.u, buff_mat_u1, gcfg.extents_u(), offset_u, idx_u, Axis::X);
    advect_component(grid.v, buff_mat_v1, gcfg.extents_v(), offset_v, idx_v, Axis::Y);
    advect_component(grid.w, buff_mat_w1, gcfg.extents_w(), offset_w, idx_w, Axis::Z);
    // clang-format on

    // update the grid
    grid.u = buff_mat_u1;
    grid.v = buff_mat_v1;
    grid.w = buff_mat_w1;
}
