#include "../solver.hpp"

template <grid::Field F>
auto interpolate_field(grid::FieldInterface<F> field,
                       const Eigen::Vector3<grid::Config::float_type>& pos) {
    auto stencil = field.interpolation_stencil(pos);
    return std::ranges::fold_left(
        std::views::zip(stencil.weights(), stencil.neighbours()),
        grid::Config::float_type{0}, [&](auto val, auto&& e) {
            auto&& [w, neigh] = e;
            auto [i, j, k] = neigh;
            return val + (w * field[i, j, k]);
        });
}

auto interpolate_velocity(GridView grid,
                          const Eigen::Vector3<grid::Config::float_type>& pos)
    -> Eigen::Vector3<grid::Config::float_type> {
    return {
        interpolate_field(grid.view<grid::Field::U>(), pos),
        interpolate_field(grid.view<grid::Field::V>(), pos),
        interpolate_field(grid.view<grid::Field::W>(), pos),
    };
}

// Step 8: update particle positions
void solver::Solver::advect_particles(Config cfg) {
    using Vec3 = Eigen::Vector3<grid::Config::float_type>;

    const auto dt = cfg.dt;
    auto grid = cfg.grid;
    auto particles = cfg.particles;

    if (!particles) {
        throw core::Error(
            "Solver::advect_particles received config without particles");
    }

    const auto& conf = grid.config();

    const Vec3 lower{
        conf.lowest_corner[0],
        conf.lowest_corner[1],
        conf.lowest_corner[2],
    };

    const Vec3 upper = lower + Vec3{
                                   conf.grid_size[0],
                                   conf.grid_size[1],
                                   conf.grid_size[2],
                               };

    for (Eigen::Index row = 0; row < particles->rows(); ++row) {
        const Vec3 point = particles->row(row).transpose();

        const Vec3 velocity = interpolate_velocity(grid, point);

        const Vec3 next =
            (point + dt * velocity).cwiseMax(lower).cwiseMin(upper);

        particles->row(row) = next.transpose();
    }
}
