#pragma once

#include <Eigen/Core>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "mesh.hpp"

// ============================================================
// CellState enum and comparison helpers.
// von Neumann neighbourhood
enum class CellState {
    SOLID,    // boundaries/walls
    AIR,      // no particles
    FULL,     // #particles > 0 && all neighbours not empty
    SURFACE,  // #particles > 0 && (#empty neighbours > 0)
    UNKNOWN,
};

// enum class is not converted automatically, hence the helper.
// since enums can not be put into Eigen matrices we need this
constexpr auto operator==(CellState s, std::integral auto i) -> bool {
    return std::to_underlying(s) ==
           static_cast<std::underlying_type_t<CellState>>(i);
}

// ============================================================
// grid specific helpers

struct Grid;  // forward declaration

// TODO: we should maybe combine those functions into one.
// in a perfect world the grid does not know about the existence of the mesh.
// -> maybe we defien two sets of particles:
//      1. one for the surface
//      2. one for the liquid itself
//  -> this way we can initialise like:
//      - no particles          -> empty
//      - only interior         -> full
//      - at least one surface  -> surface
//      - grid boundaries       -> solid (for now? maybe more later idk)
void initialise_cell_states(Grid& grid, const Mesh* mesh);
void update_cell_states(Grid& grid, const Eigen::MatrixXf& surface_particles);

void mark_surface_cells(Grid& grid);
void enforce_velocity_boundaries(Grid& grid);

// point contained in physical domain?
auto contains_point(const Grid& grid, const Eigen::Vector3f& point) -> bool;

// get (flat) cell index containing point
auto index_of(const Grid& grid, const Eigen::Vector3f& point) -> std::size_t;

// trilinear interpolation.
// point and field_offset are expressed in cell coordinates.
auto interpolate_field(const Eigen::VectorXf& field,
                       const Eigen::Vector3i& field_extents,
                       const Eigen::Vector3f& field_offset,
                       const Eigen::Vector3f& point) -> float;

// velocity interpolation in physical coordinates
auto interpolate_velocity(const Grid& grid, const Eigen::Vector3f& point)
    -> Eigen::Vector3f;

struct Grid {
    class Config {
        Eigen::Vector3f _min_corner = Eigen::Vector3f::Zero();
        Eigen::Vector3f _physical_size = Eigen::Vector3f::Ones();
        Eigen::Vector3i _cell_extents = Eigen::Vector3i::Constant(10);

       public:
        Config(Eigen::Vector3f min_corner, Eigen::Vector3f physical_size,
               Eigen::Vector3i extents);

        static auto Uniform(Eigen::Vector3f min_corner, float dx, int extents)
            -> Config;

        // geometry stuff
        auto min_corner() const -> Eigen::Vector3f;
        auto max_corner() const -> Eigen::Vector3f;
        auto physical_size() const -> Eigen::Vector3f;
        auto cell_size() const -> Eigen::Vector3f;

        // field extents
        auto extents() const -> Eigen::Vector3i;
        auto extents_u() const -> Eigen::Vector3i;
        auto extents_v() const -> Eigen::Vector3i;
        auto extents_w() const -> Eigen::Vector3i;

        // flat sizes
        auto size() const -> std::size_t;
        auto size_u() const -> std::size_t;
        auto size_v() const -> std::size_t;
        auto size_w() const -> std::size_t;

        // indexing
        auto cell(auto i, auto j, auto k) const -> std::size_t;
        auto u(auto i, auto j, auto k) const -> std::size_t;
        auto v(auto i, auto j, auto k) const -> std::size_t;
        auto w(auto i, auto j, auto k) const -> std::size_t;
    };

   private:
    Config _config;

   public:
    explicit Grid(Config grid_config);

    auto config() const -> const Config& { return _config; }
    void resize(Config cfg);

    Eigen::VectorXf u, v, w, p, t, eta;
    Eigen::VectorXi cell_state;
};

// ============================================================
// Grid::Config helper functions
inline auto Grid::Config::extents() const -> Eigen::Vector3i {
    return _cell_extents;
}

inline auto Grid::Config::extents_u() const -> Eigen::Vector3i {
    return extents() + Eigen::Vector3i(3, 2, 2);
}

inline auto Grid::Config::extents_v() const -> Eigen::Vector3i {
    return extents() + Eigen::Vector3i(2, 3, 2);
}

inline auto Grid::Config::extents_w() const -> Eigen::Vector3i {
    return extents() + Eigen::Vector3i(2, 2, 3);
}

// ============================================================
// indexing
inline auto Grid::Config::cell(auto i, auto j, auto k) const -> std::size_t {
    const auto e = extents();

    assert(i >= 0 && i < e.x());
    assert(j >= 0 && j < e.y());
    assert(k >= 0 && k < e.z());

    return std::size_t(i) +
           std::size_t(e.x()) *
               (std::size_t(j) + std::size_t(e.y()) * std::size_t(k));
}

inline auto Grid::Config::u(auto i, auto j, auto k) const -> std::size_t {
    const auto e = extents_u();

    assert(i >= 0 && i < e.x());
    assert(j >= 0 && j < e.y());
    assert(k >= 0 && k < e.z());

    return std::size_t(i) +
           std::size_t(e.x()) *
               (std::size_t(j) + std::size_t(e.y()) * std::size_t(k));
}

inline auto Grid::Config::v(auto i, auto j, auto k) const -> std::size_t {
    const auto e = extents_v();

    assert(i >= 0 && i < e.x());
    assert(j >= 0 && j < e.y());
    assert(k >= 0 && k < e.z());

    return std::size_t(i) +
           std::size_t(e.x()) *
               (std::size_t(j) + std::size_t(e.y()) * std::size_t(k));
}

inline auto Grid::Config::w(auto i, auto j, auto k) const -> std::size_t {
    const auto e = extents_w();

    assert(i >= 0 && i < e.x());
    assert(j >= 0 && j < e.y());
    assert(k >= 0 && k < e.z());

    return std::size_t(i) +
           std::size_t(e.x()) *
               (std::size_t(j) + std::size_t(e.y()) * std::size_t(k));
}
