#include "../solver.hpp"
#include "Eigen/Core"
#include "grid/core.hpp"
#include "grid/grid.hpp"
#include "grid/interface.hpp"

namespace solver {

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

constexpr auto offset = []<grid::Field F>(grid::FieldInterface<F>) constexpr {
    if constexpr (F == grid::Field::U) {
        return 0;
    } else if constexpr (F == grid::Field::V) {
        return 1;
    } else {
        return 2;
    }
};

// Step 1: Advection and external forces
// old: (eulerian advection)
// u1 = u_old - dt * (u_old · ∇)u_old + dt * g
// what i did: (semi-lagrangian advection)
// u*(x) = u_old(x - dt * u_old(x)) + dt * g(x)
void Solver::advect_velocity(Config cfg) {
    using float_t = Solver::float_type;
    using namespace grid;

    auto grid = cfg.grid;
    assert(advect_velocity_buffer.size() ==
           grid.view<Field::Velocities>().size());

    const auto body_forces = [](const Eigen::Vector3<float_type>& x) {
        return Eigen::Vector3<float_type>{0.0, -9.81, 0.0};
    };

    auto fields = std::tuple{
        grid.view<Field::U>(),
        grid.view<Field::V>(),
        grid.view<Field::W>(),
    };

    const auto ex = std::get<0>(fields).extents();
    const auto N = std::get<0>(fields).size();

    auto evaluate_component = [&](auto field, auto i, auto j, auto k) {
        if (field.cell_type(i, j, k) != CellType::INTERIOR) {
            return field[i, j, k];
        }

        const auto grid_pos = field.position_at(i, j, k);

        {  // sanity check
            // making sure the interpolated fields actually contain our point
            const auto& [U, V, W] = fields;
            assert(U.contains_point(grid_pos));
            assert(V.contains_point(grid_pos));
            assert(W.contains_point(grid_pos));
        }

        const auto velocity = interpolate_velocity(grid, grid_pos);
        const auto prev_pos = grid_pos - cfg.dt * velocity;

        {  // sanity check
            // making sure this point is also contained
            assert(field.contains_point(prev_pos));
        }

        const auto advected = interpolate_field(field, prev_pos);
        const auto acceleration = body_forces(grid_pos);

        return advected + cfg.dt * acceleration[offset(field)];
    };

    auto advect_component = [&](auto field, auto i, auto j, auto k, auto g) {
        advect_velocity_buffer[g + (N * offset(field))] =
            evaluate_component(field, i, j, k);
    };

    std::size_t glob_id = 0;
    for (std::size_t i = 0; i < ex.extent(0); ++i) {
        for (std::size_t j = 0; j < ex.extent(1); ++j) {
            for (std::size_t k = 0; k < ex.extent(2); ++k) {
                std::apply(
                    [&](auto&&... fields) {
                        // fields dont depend on eachother, so this could be
                        // paralellized later
                        (advect_component(fields, i, j, k, glob_id), ...);
                    },
                    fields);
                ++glob_id;
            }
        }
    }

    auto velocity_view = grid.view<Field::Velocities>();
    for (std::size_t i = 0; i < advect_velocity_buffer.size(); ++i) {
        velocity_view[i] = advect_velocity_buffer[i];
    }
}

}  // namespace solver
