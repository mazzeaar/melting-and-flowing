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

    assert(advect_velocity_buffer.size() == nu + nv + nw);

    const auto dt = cfg.settings.dt;
    const auto h = gcfg.cell_size();
    const auto origin = gcfg.min_corner();

    const Eigen::Vector3f gravity{0.0f, cfg.physics.gravity, 0.0f};

    // advect one staggered velocity component.
    auto advect_component = [&](const Eigen::VectorXf& field,
                                Eigen::Vector3i extents, Eigen::Vector3f offset,
                                std::size_t buffer_offset, int component) {
        const int nx = extents.x();
        const int ny = extents.y();
        const int nz = extents.z();

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    const auto idx =
                        static_cast<std::size_t>(i + nx * (j + ny * k));

                    const Eigen::Vector3f cell{
                        static_cast<float>(i),
                        static_cast<float>(j),
                        static_cast<float>(k),
                    };

                    // staggered field position in physical coordinates
                    const Eigen::Vector3f pos =
                        origin + (cell + offset).cwiseProduct(h);

                    // keep boundary values unchanged
                    if (!contains_point(grid, pos)) {
                        advect_velocity_buffer[buffer_offset + idx] =
                            field[idx];
                        continue;
                    }

                    // backtrace through the old velocity field
                    const Eigen::Vector3f velocity =
                        interpolate_velocity(grid, pos);

                    const Eigen::Vector3f prev_pos = pos - dt * velocity;

                    // clamp to the physical domain
                    const Eigen::Vector3f clamped =
                        prev_pos.cwiseMax(gcfg.min_corner())
                            .cwiseMin(gcfg.max_corner());

                    // convert physical coordinates to cell coordinates
                    const Eigen::Vector3f grid_pos =
                        (clamped - origin).cwiseQuotient(h);

                    const float advected =
                        interpolate_field(field, extents, offset, grid_pos);

                    advect_velocity_buffer[buffer_offset + idx] =
                        advected + dt * gravity[component];
                }
            }
        }
    };

    const Eigen::Vector3f offset_u{-1.0f, -0.5f, -0.5f};
    const Eigen::Vector3f offset_v{-0.5f, -1.0f, -0.5f};
    const Eigen::Vector3f offset_w{-0.5f, -0.5f, -1.0f};

    advect_component(grid.u, gcfg.extents_u(), offset_u, 0, 0);
    advect_component(grid.v, gcfg.extents_v(), offset_v, nu, 1);
    advect_component(grid.w, gcfg.extents_w(), offset_w, nu + nv, 2);

    // update the grid
    grid.u = advect_velocity_buffer.segment(0, nu);
    grid.v = advect_velocity_buffer.segment(nu, nv);
    grid.w = advect_velocity_buffer.segment(nu + nv, nw);
}
