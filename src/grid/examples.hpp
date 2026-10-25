#pragma once

#include <grid/view.hpp>  // include grid
#include <print>
#include <type_traits>
#include <unordered_map>

namespace grid::examples {
auto init_grid(auto Nx, auto Ny, auto Nz) -> decltype(auto) {
    // owns the data
    auto storage = std::make_unique<Storage>(Nx, Ny, Nz);

    // exposes:
    // size of a flattened buffer
    // -> all buffers have the same size
    auto buff_size = storage->buffer_size();

    // number of cells in each direction
    auto ext = storage->extents();

    // (i, j, k) -> linear index
    auto map = storage->mapping();

    // linear view of grid
    // val = flat_view[i]
    auto flat_view = storage->template flat_view<Field::W>();

    // 3d grid view
    // val = grid_view[i, j, k]
    auto grid_view = storage->template md_view<Field::U>();

    return storage;
}

auto use_view() {
    // real/physical grid size
    constexpr auto Nx = 4;
    constexpr auto Ny = 4;
    constexpr auto Nz = 5;

    auto data_ptr = init_grid(Nx, Ny, Nz);  // std::unique_ptr
    auto grid     = View{data_ptr.get()};   // non-owning view

    // costs nothing to pass around
    static_assert(std::is_trivially_copyable_v<View>);
    static_assert(std::is_trivially_destructible_v<View>);

    // get a 3d grid of some field
    constexpr auto FT       = Field::U;  // u-velocity
    auto           my_field = grid.get<FT>();

    auto grad = std::vector(grid.buffer_size(), 0.0);

    auto dx = 0.2;
    for (auto i : grid.iextent<Axis::X>()) {
        for (auto j : grid.iextent<Axis::Y>()) {
            for (auto k : grid.iextent<Axis::Z>()) {
                auto flat_index = grid.flat_index<FT>(i, j, k);

                if (grid.cell_type<FT>(i, j, k) == CellType::PLACEHOLDER) {
                    continue;
                }

                if (grid.cell_type<FT>(i, j, k) == CellType::BOUNDARY) {
                    // handle bc
                    continue;
                }

                // d^2 my_field/dx^2
                if (grid.cell_type<FT>(i, j, k) == CellType::INTERIOR) {
                    auto r_val       = my_field[i, j, k];
                    auto l_val       = my_field[i - 1, j, k];
                    grad[flat_index] = (r_val - l_val) / (2.0 * dx);
                }
            }
        }
    }

    // missing:
    // - a way to wrap grad in a 3d view as well?
    //      basically: have some way of passing an array and turning it into a view
    //      grid_view.make_md_view<Field::Cell?>(my_range);
    // - a way to swap buffers is probably also needed? maybe only allow this in owning_grid?
    //      - then solver takes buffers + grid view and updates buffers, then after solver we update
    //      the state
}

// void grid_example() {
//     using namespace grid;
//     constexpr auto Nx = 4;
//     constexpr auto Ny = 4;
//     constexpr auto Nz = 5;
//
//     std::println("starting grid example");
//
//     constexpr auto FT = Field::U;
//
//     // owns data
//     auto storage = std::make_shared<Storage>(Nx, Ny, Nz);
//
//     // non-owning view
//     auto grid_view = View{storage};
//
//     // non owning 3d view of field FT
//     auto field = grid_view.get<FT>();
//
//     auto m = std::unordered_map<CellType, std::pair<double, char>>{
//         {
//             CellType::INTERIOR,
//             {0.0, 'I'},
//         },
//         {
//             CellType::PLACEHOLDER,
//             {1.0, 'P'},
//         },
//         {
//             CellType::DIRECT_BOUNDARY,
//             {2.0, 'D'},
//         },
//         {
//             CellType::GHOST_BOUNDARY,
//             {3.0, 'G'},
//         },
//     };
//
//     static constexpr auto init_val    = -1.0;
//     auto                  init_lambda = [&](auto i, auto j, auto k) { field[i, j, k] = init_val;
//     };
//
//     // iterate over all indices associated with FT and execute init_lambda
//     grid_view.for_each<FT>(init_lambda);
//
//     auto setter_lambda = [&](auto v) {
//         return [&, v](auto i, auto j, auto k) {
//             if (field[i, j, k] != init_val) {
//                 std::println("error: field overlap at ({}, {}, {})", i, j, k);
//             }
//             field[i, j, k] = v;
//         };
//     };
//
//     // iterate over given cell type and set the value accordingly
//     grid_view.for_each<FT, CellType::INTERIOR>(setter_lambda(0.0));
//     grid_view.for_each<FT, CellType::PLACEHOLDER>(setter_lambda(1.0));
//     grid_view.for_each<FT, CellType::DIRECT_BOUNDARY>(setter_lambda(2.0));
//     grid_view.for_each<FT, CellType::GHOST_BOUNDARY>(setter_lambda(3.0));
//
//     // print the grid layer by layer
//     // check for overwritten cells
//     for (auto i : grid_view.iextent<Axis::X>()) {
//         std::println("x-slice {}", i);
//         for (auto j : grid_view.iextent<Axis::Y>()) {
//             for (auto k : grid_view.iextent<Axis::Z>()) {
//                 auto cell_type = grid_view.cell_type<FT>(i, j, k);
//                 auto [v, c]    = m.at(cell_type);
//                 if (field[i, j, k] != v) {
//                     std::println("\nerror: cell ({}, {}, {}) was overwritten\n", i, j, k);
//                 } else {
//                     std::print(" {} ", c);
//                 }
//             }
//             std::println();
//         }
//         std::println();
//     }
// }

}  // namespace grid::examples
