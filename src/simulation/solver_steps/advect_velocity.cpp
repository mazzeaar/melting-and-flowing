#include <utility>

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

    const auto nu = gcfg.size_u();
    const auto nv = gcfg.size_v();
    const auto nw = gcfg.size_w();

    const auto dt = cfg.settings.dt;
    const auto h = gcfg.cell_size();
    const auto origin = gcfg.min_corner();

    const auto gravity = Eigen::Vector3f{0.0f, cfg.physics.gravity, 0.0f};

    enum class Axis { X, Y, Z };

    // this is slow
    // advect one staggered velocity component.
    auto advect_component = [&](const Eigen::VectorXf& field,
                                Eigen::Vector3i extents, Eigen::Vector3f offset,
                                Eigen::VectorXf& buffer, auto&& index_func,
                                Axis axis) {
        for (int i = 0; i < extents.x(); ++i) {
            for (int j = 0; j < extents.y(); ++j) {
                for (int k = 0; k < extents.z(); ++k) {
                    const auto idx = index_func(i, j, k);

                    const Eigen::Vector3f cell{
                        static_cast<float>(i),
                        static_cast<float>(j),
                        static_cast<float>(k),
                    };

                    // staggered field position in physical coordinates
                    const auto pos = origin + (cell + offset).cwiseProduct(h);

                    // keep boundary values unchanged
                    if (!grid.contains_point(pos)) {
                        buffer[idx] = field[idx];
                        continue;
                    }

                    // this is slow
                    // backtrace through the old velocity field
                    const auto velocity = interpolate_velocity(grid, pos);

                    const auto prev_pos = pos - dt * velocity;

                    // clamp to the physical domain
                    const Eigen::Vector3f clamped =
                        prev_pos.cwiseMax(gcfg.min_corner())
                            .cwiseMin(gcfg.max_corner());

                    // convert physical coordinates to cell coordinates
                    const auto grid_pos = (clamped - origin).cwiseQuotient(h);

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
    // clang-format on

    // advect everything
    // this is very slow
    advect_component(grid.u, gcfg.extents_u(), offset_u, buff_mat_u1, idx_u,
                     Axis::X);
    advect_component(grid.v, gcfg.extents_v(), offset_v, buff_mat_v1, idx_v,
                     Axis::Y);
    advect_component(grid.w, gcfg.extents_w(), offset_w, buff_mat_w1, idx_w,
                     Axis::Z);

    // update the grid
    grid.u = buff_mat_u1;
    grid.v = buff_mat_v1;
    grid.w = buff_mat_w1;
}
