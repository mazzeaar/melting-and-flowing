#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>

#include <Eigen/Dense>
#include <array>
#include <core/paths.hpp>
#include <grid/view.hpp>
#include <print>
#include <unordered_map>
#include <vector>

void grid_example() {
    using namespace grid;
    constexpr auto Nx = 4;
    constexpr auto Ny = 4;
    constexpr auto Nz = 5;

    std::println("starting grid example");

    constexpr auto FT = Field::U;

    // owns data
    auto storage = std::make_shared<Storage>(Nx, Ny, Nz);

    // non-owning view
    auto grid_view = View{storage};

    // non owning 3d view of field FT
    auto field = grid_view.get<FT>();

    auto m = std::unordered_map<CellType, std::pair<double, char>>{
        {
            CellType::INTERIOR,
            {0.0, 'I'},
        },
        {
            CellType::PLACEHOLDER,
            {1.0, 'P'},
        },
        {
            CellType::DIRECT_BOUNDARY,
            {2.0, 'D'},
        },
        {
            CellType::GHOST_BOUNDARY,
            {3.0, 'G'},
        },
    };

    static constexpr auto init_val    = -1.0;
    auto                  init_lambda = [&](auto i, auto j, auto k) { field[i, j, k] = init_val; };

    // iterate over all indices associated with FT and execute init_lambda
    grid_view.for_each<FT>(init_lambda);

    auto setter_lambda = [&](auto v) {
        return [&, v](auto i, auto j, auto k) {
            if (field[i, j, k] != init_val) {
                std::println("error: field overlap at ({}, {}, {})", i, j, k);
            }
            field[i, j, k] = v;
        };
    };

    // iterate over given cell type and set the value accordingly
    grid_view.for_each<FT, CellType::INTERIOR>(setter_lambda(0.0));
    grid_view.for_each<FT, CellType::PLACEHOLDER>(setter_lambda(1.0));
    grid_view.for_each<FT, CellType::DIRECT_BOUNDARY>(setter_lambda(2.0));
    grid_view.for_each<FT, CellType::GHOST_BOUNDARY>(setter_lambda(3.0));

    // print the grid layer by layer
    // check for overwritten cells
    for (auto i : grid_view.iextent<Axis::X>()) {
        std::println("x-slice {}", i);
        for (auto j : grid_view.iextent<Axis::Y>()) {
            for (auto k : grid_view.iextent<Axis::Z>()) {
                auto cell_type = grid_view.cell_type<FT>(i, j, k);
                auto [v, c]    = m.at(cell_type);
                if (field[i, j, k] != v) {
                    std::println("\nerror: cell ({}, {}, {}) was overwritten\n", i, j, k);
                } else {
                    std::print(" {} ", c);
                }
            }
            std::println();
        }
        std::println();
    }
}

int main() {
    // Eigen: basic vector operations
    Eigen::Vector3d position{1.0, 2.0, 3.0};
    Eigen::Vector3d direction{0.0, 1.0, 0.0};

    std::println("Position: {}", position.transpose());
    std::println("Length: {}", position.norm());
    std::println("Dot product: {}", position.dot(direction));

    // Polyscope: visualize a point cloud
    std::vector<std::array<double, 3>> points{
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 1.0, 1.0},
    };

    grid_example();

    polyscope::init();
    polyscope::registerPointCloud("Points", points);
    polyscope::show();
}
