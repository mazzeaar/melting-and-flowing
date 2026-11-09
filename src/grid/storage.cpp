#include "grid/storage.hpp"

#include <cstddef>

#include "../mesh/mesh.hpp"
#include "core.hpp"

namespace grid {

constexpr auto cell_index_helper(GridView gv) {
    return [field = gv.view<Field::P>()](const auto& vec) {
        auto v = vec;
        return field.index_at(vec);
    };
}

void GridView::update_cell_states(
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& pos) {
    const auto ex = config().grid_extents;
    const Eigen::Vector3i dims = Eigen::Vector3i{
        static_cast<int>(ex.extent(0)), static_cast<int>(ex.extent(1)),
        static_cast<int>(ex.extent(2))};

    const int nx = dims.x(), ny = dims.y(), nz = dims.z();

    auto index = [=](int i, int j, int k) { return i + nx * (j + ny * k); };

    auto valid = [=](int i, int j, int k) {
        return i >= 0 && j >= 0 && k >= 0 && i < nx && j < ny && k < nz;
    };

    auto cell_index = cell_index_helper(*this);

    auto state = view<Field::State>();
    std::vector<CellState> region(static_cast<std::size_t>(nx * ny * nz),
                                  CellState::UNKNOWN);

    // 1. Voxelize the surface particles.
    for (Eigen::Index p = 0; p < pos.rows(); ++p) {
        Eigen::Vector3i cell = cell_index(pos.row(p).transpose());

        if (valid(cell.x(), cell.y(), cell.z())) {
            region[index(cell.x(), cell.y(), cell.z())] = CellState::SURFACE;
        }
    }

    // 2. Flood-fill air from domain boundaries.
    std::queue<Eigen::Vector3i> queue;

    auto add_air = [&](int i, int j, int k) {
        if (!valid(i, j, k)) {
            return;
        }

        auto& state = region[index(i, j, k)];
        if (state != CellState::UNKNOWN) {
            return;
        }

        state = CellState::AIR;
        queue.emplace(i, j, k);
    };

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            add_air(0, j, k);
            add_air(nx - 1, j, k);
        }
    }

    for (int k = 0; k < nz; ++k) {
        for (int i = 0; i < nx; ++i) {
            add_air(i, 0, k);
            add_air(i, ny - 1, k);
        }
    }

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            add_air(i, j, 0);
            add_air(i, j, nz - 1);
        }
    }

    constexpr std::array<Eigen::Vector3i, 6> neighbors = {
        {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}};

    while (!queue.empty()) {
        auto cell = queue.front();
        queue.pop();

        for (const auto& d : neighbors) {
            auto n = cell + d;
            add_air(n.x(), n.y(), n.z());
        }
    }

    // 3. Update the actual grid state.
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const bool is_boundary = i == 0 || i == nx - 1 || j == 0 ||
                                         j == ny - 1 || k == 0 || k == nz - 1;

                if (is_boundary) {
                    state[i, j, k] = CellState::SOLID;
                    continue;
                }

                auto r = region[index(i, j, k)];

                if (r == CellState::UNKNOWN) {
                    r = CellState::FULL;
                }

                state[i, j, k] = r;
            }
        }
    }
}

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
                    state[i, j, k] = CellState::AIR;
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
