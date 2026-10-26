#pragma once

#include <grid/grid.hpp>

namespace grid::example {

void func() {
    using namespace grid;

    // define resolution, position and size of grid
    auto extents       = Config::Extents{50, 50, 50};
    auto grid_size     = std::array{12.0, 13.0, 14.0};
    auto lowest_corner = std::array{12.0, 13.0, 14.0};

    auto conf         = Config(extents, grid_size, lowest_corner);
    auto grid_storage = GridStorage{conf};

    // free to move around, copy and destroy
    auto grid_view = GridView{grid_storage};
    static_assert(std::is_trivially_copyable_v<GridView>);
    static_assert(std::is_trivially_copy_constructible_v<GridView>);
    static_assert(std::is_trivially_move_constructible_v<GridView>);

    auto u = grid_view.view<Field::U>();
    static_assert(std::is_trivially_copyable_v<decltype(u)>);
    static_assert(std::is_trivially_copy_constructible_v<decltype(u)>);
    static_assert(std::is_trivially_move_constructible_v<decltype(u)>);

    {  // easy indexing
        auto indexing_in_3d = u[1, 2, 3];
        auto flat_indexing  = u[1];
    }

    {  // global and 3d indexing offsets
        // probably only usefull for the implicit stencils
        auto i = 5;
        auto j = 4;
        auto k = 3;

        // get directional offsets for flat indexing
        auto [di, dj, dk] = u.flat_offsets();
        assert(dj == u.flat_offset<Axis::Y>());

        auto glob_idx = u.flat_index(i, j, k);

        u[i, j, k] = 42;
        assert((u[glob_idx]) == (u[i, j, k]));

        u[i - 1, j + 1, k + 2] = u[i, j, k];
        auto new_idx           = glob_idx - di + dj + (2 * dk);
        assert((u[new_idx]) == (u[i, j, k]));
    }

    {  // classify cell types
        // this obviously only passes if its actually the interior
        auto type = u.cell_type(1, 2, 3);
        assert(type == CellType::INTERIOR);

        auto bounds = CellType::DIRECT_BOUNDARY;
        assert(detail::matches(bounds, CellType::BOUNDARY));
    }

    {  // different views of the view lol
        auto u_3d_view  = u.md_view();
        auto u_1d_view  = u.flat_view();
        auto eigen_view = u.eigen_view();
    }

    {  // use buffers without headaches
        std::vector<double> my_buffer;
        grid_view.resize_like<Field::P>(my_buffer);
    }
}
}  // namespace grid::example
