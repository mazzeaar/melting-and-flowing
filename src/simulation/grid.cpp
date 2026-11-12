#include "grid.hpp"

#include <Eigen/Core>

#include "core/exception.hpp"

Grid::Config::Config(Eigen::Vector3f min_corner, Eigen::Vector3f physical_size,
                     Eigen::Vector3i extents)
    : _min_corner{min_corner},
      _physical_size{physical_size},
      _cell_extents{extents} {
    const auto all_finite =
        _min_corner.allFinite() && _physical_size.allFinite();
    const auto all_positive = (_physical_size.array() > 0).all() &&
                              (_cell_extents.array() >= 2).all();

    if (!all_finite || !all_positive) {
        throw core::Error(
            "Grid requires positive finite dimensions and at least 2 "
            "cells per axis");
    }
}

auto Grid::Config::Uniform(Eigen::Vector3f min_corner, float dx, int extents)
    -> Grid::Config {
    return Config{
        min_corner,
        Eigen::Vector3f::Constant(dx * extents),
        Eigen::Vector3i::Constant(extents),
    };
}

// ============================================================
// Grid::Config implementation
// ============================================================

auto Grid::Config::min_corner() const -> Eigen::Vector3f { return _min_corner; }

auto Grid::Config::physical_size() const -> Eigen::Vector3f {
    return _physical_size;
}

auto Grid::Config::max_corner() const -> Eigen::Vector3f {
    return min_corner() + physical_size();
}

auto Grid::Config::cell_size() const -> Eigen::Vector3f {
    return _physical_size.cwiseQuotient(_cell_extents.cast<float>());
}

// ============================================================
// Flat sizes
// ============================================================

auto Grid::Config::size() const -> std::size_t {
    const auto e = extents();
    return std::size_t(e.x()) * e.y() * e.z();
}

auto Grid::Config::size_u() const -> std::size_t {
    const auto e = extents_u();
    return std::size_t(e.x()) * e.y() * e.z();
}

auto Grid::Config::size_v() const -> std::size_t {
    const auto e = extents_v();
    return std::size_t(e.x()) * e.y() * e.z();
}

auto Grid::Config::size_w() const -> std::size_t {
    const auto e = extents_w();
    return std::size_t(e.x()) * e.y() * e.z();
}

// ============================================================
// Grid implementation
// ============================================================

Grid::Grid(Config grid_config) : _config(grid_config) { resize(grid_config); }

void Grid::resize(Config cfg) {
    _config = cfg;

    u.setZero(static_cast<Eigen::Index>(cfg.size_u()));
    v.setZero(static_cast<Eigen::Index>(cfg.size_v()));
    w.setZero(static_cast<Eigen::Index>(cfg.size_w()));

    const auto n = static_cast<Eigen::Index>(cfg.size());

    p.setZero(n);
    t.setConstant(n, 300.0f);  // K, matching update_viscosity's T_ref
    eta.setConstant(n, 1e5f);
    cell_state.setConstant(n, std::to_underlying(CellState::AIR));
}

// ============================================================
// Geometry
// ============================================================

auto contains_point(const Grid& grid, const Eigen::Vector3f& point) -> bool {
    const auto min = grid.config().min_corner();
    const Eigen::Vector3f max = grid.config().max_corner();

    return (point.array() >= min.array()).all() &&
           (point.array() <= max.array()).all();
}

auto index_of(const Grid& grid, const Eigen::Vector3f& point) -> std::size_t {
    assert(contains_point(grid, point));

    const auto& cfg = grid.config();

    const Eigen::Vector3f local =
        (point - cfg.min_corner()).cwiseQuotient(cfg.cell_size());

    const Eigen::Vector3i ijk =
        local.array().floor().cast<int>().matrix().cwiseMax(0).cwiseMin(
            (cfg.extents() - Eigen::Vector3i::Ones()));

    return cfg.cell(ijk.x(), ijk.y(), ijk.z());
}

// ============================================================
// Trilinear interpolation
// ============================================================

auto interpolate_field(const Eigen::VectorXf& field,
                       const Eigen::Vector3i& field_extents,
                       const Eigen::Vector3f& field_offset,
                       const Eigen::Vector3f& point) -> float {
    assert((field_extents.array() >= 2).all());
    assert(field.size() == field_extents.cast<Eigen::Index>().prod());
    assert(point.allFinite());

    const Eigen::Vector3f grid_pos = point - field_offset;

    const Eigen::Vector3i max_index =
        field_extents - Eigen::Vector3i::Constant(2);

    const Eigen::Vector3i i0 =
        grid_pos.array().floor().cast<int>().matrix().cwiseMax(0).cwiseMin(
            max_index);

    const Eigen::Vector3f t =
        (grid_pos - i0.cast<float>()).cwiseMax(0.0f).cwiseMin(1.0f);

    const auto nx = std::size_t(field_extents.x());
    const auto ny = std::size_t(field_extents.y());

    const auto index = [&](int i, int j, int k) -> std::size_t {
        return std::size_t(i) + nx * (std::size_t(j) + ny * std::size_t(k));
    };

    const int i = i0.x();
    const int j = i0.y();
    const int k = i0.z();

    const float c000 = field[index(i, j, k)];
    const float c100 = field[index(i + 1, j, k)];
    const float c010 = field[index(i, j + 1, k)];
    const float c110 = field[index(i + 1, j + 1, k)];

    const float c001 = field[index(i, j, k + 1)];
    const float c101 = field[index(i + 1, j, k + 1)];
    const float c011 = field[index(i, j + 1, k + 1)];
    const float c111 = field[index(i + 1, j + 1, k + 1)];

    const float x = t.x();
    const float y = t.y();
    const float z = t.z();

    const float c00 = c000 * (1 - x) + c100 * x;
    const float c10 = c010 * (1 - x) + c110 * x;
    const float c01 = c001 * (1 - x) + c101 * x;
    const float c11 = c011 * (1 - x) + c111 * x;

    const float c0 = c00 * (1 - y) + c10 * y;
    const float c1 = c01 * (1 - y) + c11 * y;

    return c0 * (1 - z) + c1 * z;
}

// ============================================================
// Velocity interpolation
// ============================================================

auto interpolate_velocity(const Grid& grid, const Eigen::Vector3f& point)
    -> Eigen::Vector3f {
    const auto& cfg = grid.config();

    // Convert physical position to cell coordinates.
    const Eigen::Vector3f pos =
        (point - cfg.min_corner()).cwiseQuotient(cfg.cell_size());

    // First stored sample locations, including ghost cells.
    const Eigen::Vector3f offset_u(-1.0f, -0.5f, -0.5f);
    const Eigen::Vector3f offset_v(-0.5f, -1.0f, -0.5f);
    const Eigen::Vector3f offset_w(-0.5f, -0.5f, -1.0f);

    return {
        interpolate_field(grid.u, cfg.extents_u(), offset_u, pos),
        interpolate_field(grid.v, cfg.extents_v(), offset_v, pos),
        interpolate_field(grid.w, cfg.extents_w(), offset_w, pos),
    };
}

// The scalar grid has no ghost cells. All velocity indices above are STORAGE
// indices: physical face (i,j,k) is stored at (i+1,j+1,k+1).
// Domain walls are faces, rather than an extra ring of solid scalar cells.
void enforce_velocity_boundaries(Grid& grid) {
    const auto gc = grid.config();
    const Eigen::Vector3i n = gc.extents();
    const auto solid = [&](int i, int j, int k) {
        return i < 0 || j < 0 || k < 0 || i >= n.x() || j >= n.y() ||
               k >= n.z() ||
               grid.cell_state[gc.cell(i, j, k)] == CellState::SOLID;
    };

    for (int k = 0; k < n.z(); ++k) {
        for (int j = 0; j < n.y(); ++j) {
            for (int i = 0; i <= n.x(); ++i) {
                if (solid(i - 1, j, k) || solid(i, j, k)) {
                    grid.u[gc.u(i + 1, j + 1, k + 1)] = 0;
                }
            }
        }
    }

    for (int k = 0; k < n.z(); ++k) {
        for (int j = 0; j <= n.y(); ++j) {
            for (int i = 0; i < n.x(); ++i) {
                if (solid(i, j - 1, k) || solid(i, j, k)) {
                    grid.v[gc.v(i + 1, j + 1, k + 1)] = 0;
                }
            }
        }
    }

    for (int k = 0; k <= n.z(); ++k) {
        for (int j = 0; j < n.y(); ++j) {
            for (int i = 0; i < n.x(); ++i) {
                if (solid(i, j, k - 1) || solid(i, j, k)) {
                    grid.w[gc.w(i + 1, j + 1, k + 1)] = 0;
                }
            }
        }
    }

    // Free-slip exterior ghosts: extend tangential velocities constantly;
    // reflect normal velocities about their zero wall face.
    const auto fill = [&](auto& field, Eigen::Vector3i ex, const auto& row,
                          int axis) {
        for (int k = 0; k < ex.z(); ++k) {
            for (int j = 0; j < ex.y(); ++j) {
                for (int i = 0; i < ex.x(); ++i) {
                    Eigen::Vector3i q{i, j, k};
                    if ((q.array() > 0).all() &&
                        (q.array() < ex.array() - 1).all()) {
                        continue;
                    }
                    Eigen::Vector3i src = q.cwiseMax(1).cwiseMin(
                        ex - Eigen::Vector3i::Constant(2));
                    float sign = 1;
                    if (q[axis] == 0) {
                        src[axis] = 2;
                        sign = -1;
                    }
                    if (q[axis] == ex[axis] - 1) {
                        src[axis] = ex[axis] - 3;
                        sign = -1;
                    }
                    field[row(i, j, k)] =
                        sign * field[row(src.x(), src.y(), src.z())];
                }
            }
        }
    };
    fill(
        grid.u, gc.extents_u(),
        [&](int i, int j, int k) { return gc.u(i, j, k); }, 0);
    fill(
        grid.v, gc.extents_v(),
        [&](int i, int j, int k) { return gc.v(i, j, k); }, 1);
    fill(
        grid.w, gc.extents_w(),
        [&](int i, int j, int k) { return gc.w(i, j, k); }, 2);
}

void mark_surface_cells(Grid& grid) {
    const auto gc = grid.config();
    const auto ex = gc.extents();

    static constexpr std::array<Eigen::Vector3i, 6> dirs = {
        {
            {1, 0, 0},
            {-1, 0, 0},
            {0, 1, 0},
            {0, -1, 0},
            {0, 0, 1},
            {0, 0, -1},
        },
    };

    for (int k = 0; k < ex.z(); ++k) {
        for (int j = 0; j < ex.y(); ++j) {
            for (int i = 0; i < ex.x(); ++i) {
                auto& state = grid.cell_state[gc.cell(i, j, k)];

                if (state != CellState::FULL) {
                    continue;
                }

                for (const auto& d : dirs) {
                    auto q = Eigen::Vector3i{i, j, k} + d;

                    auto a = (q.array() < ex.array()).all();
                    auto b = (q.array() < 0).any();
                    auto valid_grid_index = a && b;
                    if (!valid_grid_index) {
                        continue;
                    }

                    auto is_air =
                        grid.cell_state[gc.cell(q.x(), q.y(), q.z())] ==
                        CellState::AIR;

                    if (is_air) {
                        state = std::to_underlying(CellState::SURFACE);
                        break;
                    }
                }
            }
        }
    }
}

void initialise_cell_states(Grid& grid, const Mesh* mesh) {
    if (!mesh) {
        throw core::Error("Cannot initialize grid without a mesh");
    }

    const auto gc = grid.config();
    const auto ex = gc.extents();

    constexpr int num_samples = 4;
    constexpr int per_cell = num_samples * num_samples * num_samples;

    Eigen::MatrixXf points(static_cast<Eigen::Index>(gc.size()) * per_cell, 3);

    auto sample_cell = [&](auto i, auto j, auto k) {
        const auto cell_anchor = Eigen::Vector3i{i, j, k}.cast<float>();
        const auto cell_center = Eigen::Vector3f::Constant(0.5);

        int s = 0;
        for (int z = 0; z < num_samples; ++z) {
            for (int y = 0; y < num_samples; ++y) {
                for (int x = 0; x < num_samples; ++x) {
                    auto cell_offset =
                        Eigen::Vector3i{x, y, z}.cast<float>() + cell_center;
                    auto point = cell_anchor + cell_offset / num_samples;
                    points.row(gc.cell(i, j, k) * per_cell + s++) =
                        (gc.min_corner() + point.cwiseProduct(gc.cell_size()))
                            .transpose();
                }
            }
        }
    };

    for (int k = 0; k < ex.z(); ++k) {
        for (int j = 0; j < ex.y(); ++j) {
            for (int i = 0; i < ex.x(); ++i) {
                sample_cell(i, j, k);
            }
        }
    }

    const auto inside = mesh->contains_many(points);
    for (Eigen::Index r = 0; r < grid.cell_state.size(); ++r) {
        const bool occupied = inside.segment(r * per_cell, per_cell).any();
        grid.cell_state[r] =
            std::to_underlying(occupied ? CellState::FULL : CellState::AIR);
    }
    mark_surface_cells(grid);
}

// These particles sample the SURFACE, not the fluid volume. Mark their cells,
// flood-fill exterior air, and treat the enclosed region as fluid. Sparse or
// disconnected surface samples can leak; this is an approximate voxelization.
void update_cell_states(Grid& grid, const Eigen::MatrixXf& surface_particles) {
    const auto gc = grid.config();
    const auto ex = gc.extents();
    Eigen::VectorXi region = Eigen::VectorXi::Constant(
        gc.size(), std::to_underlying(CellState::UNKNOWN));
    for (Eigen::Index r = 0; r < region.size(); ++r) {
        if (grid.cell_state[r] == CellState::SOLID) {
            region[r] = std::to_underlying(CellState::SOLID);
        }
    }
    for (Eigen::Index r = 0; r < surface_particles.rows(); ++r) {
        const Eigen::Vector3f point = surface_particles.row(r).transpose();
        if (!point.allFinite() || !contains_point(grid, point)) {
            continue;
        }
        const auto cell = index_of(grid, point);
        if (region[cell] != CellState::SOLID) {
            region[cell] = std::to_underlying(CellState::SURFACE);
        }
    }
    std::queue<Eigen::Vector3i> queue;
    const auto add_air = [&](Eigen::Vector3i q) {
        if ((q.array() < 0).any() || (q.array() >= ex.array()).any()) {
            return;
        }
        auto& state = region[gc.cell(q.x(), q.y(), q.z())];
        if (state != CellState::UNKNOWN) {
            return;
        }
        state = std::to_underlying(CellState::AIR);
        queue.push(q);
    };
    for (int k = 0; k < ex.z(); ++k) {
        for (int j = 0; j < ex.y(); ++j) {
            add_air({0, j, k});
            add_air({ex.x() - 1, j, k});
        }
    }
    for (int k = 0; k < ex.z(); ++k) {
        for (int i = 0; i < ex.x(); ++i) {
            add_air({i, 0, k});
            add_air({i, ex.y() - 1, k});
        }
    }
    for (int j = 0; j < ex.y(); ++j) {
        for (int i = 0; i < ex.x(); ++i) {
            add_air({i, j, 0});
            add_air({i, j, ex.z() - 1});
        }
    }
    const std::array<Eigen::Vector3i, 6> dirs = {
        {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}};
    while (!queue.empty()) {
        const auto q = queue.front();
        queue.pop();
        for (const auto& d : dirs) {
            add_air(q + d);
        }
    }
    for (Eigen::Index r = 0; r < region.size(); ++r) {
        if (region[r] == CellState::UNKNOWN) {
            region[r] = std::to_underlying(CellState::FULL);
        }
    }
    grid.cell_state = std::move(region);
    mark_surface_cells(grid);
}
