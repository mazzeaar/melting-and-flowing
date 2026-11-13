#include "grid.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <queue>
#include <ranges>
#include <utility>

#include "core/exception.hpp"
#include "simulation/settings.hpp"

// ============================================================
// Grid::Config implementation

// ------------------------------------------------------------
// constructors
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

Grid::Config::Config(SimulationSettings::GridSettings settings)
    : Config{settings.origin,
             settings.cell_size * settings.cell_counts.cast<float>(),
             settings.cell_counts} {}

// ------------------------------------------------------------
// helpers

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

// ------------------------------------------------------------
// flat sizes

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

Grid::Grid(Config grid_config,
           const SimulationSettings::PhysicsSettings& physics)
    : _config(grid_config) {
    resize(grid_config, physics.reference_temperature,
           physics.reference_viscosity);
}

void Grid::resize(Config cfg,
                  const SimulationSettings::PhysicsSettings& physics) {
    resize(cfg, physics.reference_temperature, physics.reference_viscosity);
}
void Grid::resize(Config cfg, float initial_temp, float initial_eta) {
    _config = cfg;

    u.setZero(static_cast<Eigen::Index>(cfg.size_u()));
    v.setZero(static_cast<Eigen::Index>(cfg.size_v()));
    w.setZero(static_cast<Eigen::Index>(cfg.size_w()));

    const auto n = static_cast<Eigen::Index>(cfg.size());

    p.setZero(n);
    t.setConstant(n, initial_temp);
    eta.setConstant(n, initial_eta);
    cell_state.setConstant(n, std::to_underlying(CellState::AIR));
}

// ============================================================
// Grid Helpers

// ------------------------------------------------------------
// geometry checks
auto Grid::contains_point(const Eigen::Vector3f& point) const -> bool {
    const auto min = config().min_corner();
    const Eigen::Vector3f max = config().max_corner();

    return (point.array() >= min.array()).all() &&
           (point.array() <= max.array()).all();
}

auto Grid::index_of(const Eigen::Vector3f& point) const -> std::size_t {
    assert(grid.contains_point(point));

    const auto& cfg = config();

    const Eigen::Vector3f local =
        (point - cfg.min_corner()).cwiseQuotient(cfg.cell_size());

    const Eigen::Vector3i ijk =
        local.array().floor().cast<int>().matrix().cwiseMax(0).cwiseMin(
            (cfg.extents() - Eigen::Vector3i::Ones()));

    return cfg.cell(ijk.x(), ijk.y(), ijk.z());
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

                    // make sure index is in bounds
                    if ((q.array() > 0).all() &&
                        (q.array() < ex.array() - 1).all()) {
                        continue;
                    }

                    Eigen::Vector3i src = q.cwiseMax(1).cwiseMin(
                        ex - Eigen::Vector3i::Constant(2));

                    auto sign = 1.f;

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

// this function walks through all cells.
// if a cell is FULL and any (Neumann) neighbour is AIR then its a surface
// particle
// PRE: cell states are already initialised through initialise_cell_states()
void mark_surface_cells(Grid& grid) {
    const auto gc = grid.config();
    const auto ex = gc.extents();

    static constexpr const std::array<Eigen::Vector3i, 6> neighbour_offsets = {
        {
            {1, 0, 0},
            {-1, 0, 0},
            {0, 1, 0},
            {0, -1, 0},
            {0, 0, 1},
            {0, 0, -1},
        },
    };

    const auto cell_is_air = [&](const auto& q) {
        // check bounds
        if ((q.array() >= ex.array()).any() || (q.array() < 0).any()) {
            return false;
        }

        return grid.cell_state[gc.cell(q)] == CellState::AIR;
    };

    for (int k = 0; k < ex.z(); ++k) {
        for (int j = 0; j < ex.y(); ++j) {
            for (int i = 0; i < ex.x(); ++i) {
                auto& state = grid.cell_state[gc.cell(i, j, k)];

                // only need to decide for full cells if they are on a boundary
                // or not
                if (state != CellState::FULL) {
                    continue;
                }

                const auto cell_index = Eigen::Vector3i{i, j, k};
                auto has_empty_neighbours = std::ranges::any_of(
                    neighbour_offsets, cell_is_air,
                    [&](const auto& d) { return d + cell_index; });

                if (has_empty_neighbours) {
                    state = std::to_underlying(CellState::SURFACE);
                }
            }
        }
    }
}

// generates evenly distributed points inside of a unit cube
auto get_cell_samples(int N) {
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> samples(N * N * N,
                                                                     3);

    int row = 0;
    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                samples.row(row++) =
                    (Eigen::Vector3f(x, y, z).array() + 0.5f) / N;
            }
        }
    }

    return samples;
}

void initialise_cell_states(Grid& grid, const Mesh* mesh,
                            int samples_per_cell_axis) {
    if (!mesh) {
        throw core::Error("Cannot initialize grid without a mesh");
    }

    const auto gc = grid.config();
    const auto ex = gc.extents();
    const auto h = gc.cell_size();

    const int N = std::max(1, samples_per_cell_axis);  // samples per axis
    const int spc = N * N * N;                         // samples per cell

    const auto cell_offsets = get_cell_samples(N);

    // sample points for the whole grid
    auto sample_points = Eigen::MatrixXf{gc.size() * spc, 3};

    // scale and shift cell_offsets to fit to a cell ijk and store the
    // coordinates in sample_points
    auto sample_cell = [&](int i, int j, int k) {
        const auto cell_origin =
            gc.min_corner() + Eigen::Vector3f(i, j, k).cwiseProduct(h);

        auto start_id = gc.cell(i, j, k) * spc;

        sample_points.middleRows(start_id, spc) =
            (cell_offsets.array().rowwise() * h.transpose().array())
                .matrix()
                .rowwise() +
            cell_origin.transpose();
    };

    // generate sample points for each cell
    for (int k = 0; k < ex.z(); ++k) {
        for (int j = 0; j < ex.y(); ++j) {
            for (int i = 0; i < ex.x(); ++i) {
                sample_cell(i, j, k);
            }
        }
    }

    const auto inside = mesh->contains_many(sample_points);
    for (Eigen::Index r = 0; r < grid.cell_state.size(); ++r) {
        const auto& segment = inside.segment(r * spc, spc);
        grid.cell_state[r] = std::to_underlying(segment.any() ? CellState::FULL
                                                              : CellState::AIR);
    }

    // this is dumb... can probably be done in the loop above directly?
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
        if (!point.allFinite() || !grid.contains_point(point)) {
            continue;
        }

        const auto cell = grid.index_of(point);
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

    // #####################################################################
    // #####################################################################
    // #####################################################################
    // #####################################################################
    // # TODO: MOVE THIS SOMEWHERE OUTSIDE, ITS USED BY MULTIPLE FUNCTIONS #
    // #####################################################################
    // #####################################################################
    // #####################################################################
    // #####################################################################
    static constexpr std::array<Eigen::Vector3i, 6> dirs = {
        {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}},
    };

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
