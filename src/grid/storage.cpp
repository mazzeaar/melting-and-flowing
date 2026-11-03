#include "grid/storage.hpp"

#include "../mesh/mesh.hpp"
#include "core.hpp"

namespace grid {

void GridView::initialise_cell_states(const Mesh& mesh) {
    const auto cfg = config();
    const auto& ext = cfg.grid_extents;

    const Eigen::Vector3i extents{static_cast<int>(ext.extent(0)),
                                  static_cast<int>(ext.extent(1)),
                                  static_cast<int>(ext.extent(2))};

    const Eigen::Vector3d anchor{cfg.lowest_corner.data()};
    const Eigen::Vector3d grid_size{cfg.grid_size.data()};

    constexpr int substep_count = 4;
    constexpr int samples_per_cell =
        substep_count * substep_count * substep_count;

    const Eigen::Vector3d cell_size =
        grid_size.cwiseQuotient(extents.cast<double>());

    const Eigen::Vector3d substep = cell_size / substep_count;

    const Eigen::Index n_cells = extents.cast<Eigen::Index>().prod();

    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> points(
        n_cells * samples_per_cell, 3);

    Eigen::Index index = 0;

    // Generate all sample points
    for (int i = 0; i < extents.x(); ++i) {
        for (int j = 0; j < extents.y(); ++j) {
            for (int k = 0; k < extents.z(); ++k) {
                const Eigen::Vector3d cell_origin =
                    anchor +
                    Eigen::Vector3d{i * 1.0, j * 1.0, k * 1.0}.cwiseProduct(
                        cell_size);

                for (int ii = 0; ii < substep_count; ++ii) {
                    for (int jj = 0; jj < substep_count; ++jj) {
                        for (int kk = 0; kk < substep_count; ++kk) {
                            const Eigen::Vector3d offset{ii + 0.5, jj + 0.5,
                                                         kk + 0.5};

                            points.row(index++) =
                                (cell_origin + offset.cwiseProduct(substep))
                                    .transpose();
                        }
                    }
                }
            }
        }
    }

    // One batched mesh query
    const auto inside = mesh.contains_many(points);

    auto state = view<Field::State>();

    index = 0;
    for (int i = 0; i < extents.x(); ++i) {
        for (int j = 0; j < extents.y(); ++j) {
            for (int k = 0; k < extents.z(); ++k) {
                const bool is_boundary = i == 0 || i == extents.x() - 1 ||
                                         j == 0 || j == extents.y() - 1 ||
                                         k == 0 || k == extents.z() - 1;

                int inside_count = 0;
                for (int s = 0; s < samples_per_cell; ++s) {
                    inside_count += static_cast<bool>(inside[index++]);
                }

                if (is_boundary) {
                    state[i, j, k] = CellState::SOLID;
                } else if (inside_count == 0) {
                    state[i, j, k] = CellState::EMPTY;
                } else if (inside_count == samples_per_cell) {
                    state[i, j, k] = CellState::SOLID;
                } else {
                    state[i, j, k] = CellState::SURFACE;
                }
            }
        }
    }
}

void GridView::advect_particles(Eigen::MatrixXd& positions) const {
    // solver/advect (step 0) provides helpers for velocity interpolation
    // already,
    // - i should move those to gridview or smth.
    // - i think velocity interpolation is also needed by temperature advection
    Log("particle advection is not implemented yet");
}

}  // namespace grid
