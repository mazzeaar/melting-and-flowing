#include "../solver.hpp"

// Step 8: update particle positions
void Solver::advect_particles(Config cfg) {
    const auto& grid = cfg.grid;
    const auto& conf = grid.config();

    const float dt = cfg.settings.dt;
    auto& particles = cfg.particles;

    assert(particles.cols() == 3 && particles.rows() > 1);

    const Eigen::Vector3f lower = conf.min_corner();
    const Eigen::Vector3f upper = conf.max_corner();

    const auto n = particles.rows();

#pragma omp parallel for schedule(static)
    for (Eigen::Index row = 0; row < n; ++row) {
        auto point = Eigen::Vector3f{
            particles(row, 0),
            particles(row, 1),
            particles(row, 2),
        };

        const auto velocity = interpolate_velocity(grid, point);

        for (int d = 0; d < 3; ++d) {
            particles(row, d) =
                std::clamp(point[d] + dt * velocity[d], lower[d], upper[d]);
        }
    }
}
