#include "../solver.hpp"

template <grid::Field F>
inline auto interpolate_field(
    grid::FieldInterface<F> field,
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

inline auto interpolate_velocity(
    GridView grid, const Eigen::Vector3<grid::Config::float_type>& pos)
    -> Eigen::Vector3<grid::Config::float_type> {
    return {
        interpolate_field(grid.view<grid::Field::U>(), pos),
        interpolate_field(grid.view<grid::Field::V>(), pos),
        interpolate_field(grid.view<grid::Field::W>(), pos),
    };
}

// Step 5: Temperature advection
// T1 = T_old - dt * (u_new · ∇)T_old
void solver::Solver::advect_temperature(solver::Solver::Config cfg) {
    using namespace grid;

    auto grid = cfg.grid;
    auto temp_field = grid.view<Field::T>();
    auto temp = temp_field.flat_view();

    const auto ex = temp_field.extents();

    assert(advect_temperature_buffer.size() == temp_field.size());

    std::size_t idx = 0;

    for (std::size_t i = 0; i < ex.extent(0); ++i) {
        for (std::size_t j = 0; j < ex.extent(1); ++j) {
            for (std::size_t k = 0; k < ex.extent(2); ++k) {
                if (temp_field.cell_type(i, j, k) != CellType::INTERIOR) {
                    advect_temperature_buffer[idx++] = temp_field[i, j, k];
                    continue;
                }

                const auto pos = temp_field.position_at(i, j, k);
                const auto velocity = interpolate_velocity(grid, pos);
                const auto prev_pos = pos - cfg.dt * velocity;

                assert(temp_field.contains_point(prev_pos));

                advect_temperature_buffer[idx++] =
                    interpolate_field(temp_field, prev_pos);
            }
        }
    }

    for (std::size_t i = 0; i < temp.size(); ++i) {
        temp[i] = advect_temperature_buffer[i];
    }
}
