#pragma once

#include <Eigen/Core>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "mesh.hpp"
#include "settings.hpp"

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
void initialise_cell_states(Grid& grid, const Mesh* mesh,
                            int samples_per_cell_axis = 4);

void update_cell_states(Grid& grid, const Eigen::MatrixXf& surface_particles);

void mark_surface_cells(Grid& grid);
void enforce_velocity_boundaries(Grid& grid);

// trilinear interpolation.
// point and field_offset are expressed in cell coordinates.
auto interpolate_field(const Eigen::VectorXf& field,
                       const Eigen::Vector3i& field_extents,
                       const Eigen::Vector3f& field_offset,
                       const Eigen::Vector3f& sample_point) -> float;

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

        Config(SimulationSettings::GridSettings settings);

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

        auto cell(const Eigen::Vector3i& id) const -> std::size_t;
        auto u(const Eigen::Vector3i& id) const -> std::size_t;
        auto v(const Eigen::Vector3i& id) const -> std::size_t;
        auto w(const Eigen::Vector3i& id) const -> std::size_t;

       private:
        template <std::integral I, std::integral J, std::integral K>
        static auto index(const Eigen::Vector3i& e, I i, J j, K k)
            -> std::size_t;

        static auto index(const Eigen::Vector3i& e, const Eigen::Vector3i& id)
            -> std::size_t;
    };

   private:
    Config _config;

   public:
    explicit Grid(Config grid_config,
                  const SimulationSettings::PhysicsSettings& physics);

    auto config() const -> const Config& { return _config; }

    void resize(Config cfg, const SimulationSettings::PhysicsSettings& physics);
    void resize(Config cfg, float initial_temp, float initial_eta);

    // point contained in physical domain?
    auto contains_point(const Eigen::Vector3f& point) const -> bool;

    // get (flat) cell index containing point
    auto index_of(const Eigen::Vector3f& point) const -> std::size_t;

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

template <std::integral I, std::integral J, std::integral K>
inline auto Grid::Config::index(const Eigen::Vector3i& e, I i, J j, K k)
    -> std::size_t {
    return static_cast<std::size_t>(i) +
           static_cast<std::size_t>(e.x()) *
               (static_cast<std::size_t>(j) +
                static_cast<std::size_t>(e.y()) * static_cast<std::size_t>(k));
}

inline auto Grid::Config::index(const Eigen::Vector3i& e,
                                const Eigen::Vector3i& id) -> std::size_t {
    return index(e, id.x(), id.y(), id.z());
}

inline auto Grid::Config::cell(auto i, auto j, auto k) const -> std::size_t {
    return index(extents(), i, j, k);
}

inline auto Grid::Config::u(auto i, auto j, auto k) const -> std::size_t {
    return index(extents_u(), i, j, k);
}

inline auto Grid::Config::v(auto i, auto j, auto k) const -> std::size_t {
    return index(extents_v(), i, j, k);
}

inline auto Grid::Config::w(auto i, auto j, auto k) const -> std::size_t {
    return index(extents_w(), i, j, k);
}

inline auto Grid::Config::cell(const Eigen::Vector3i& id) const -> std::size_t {
    return index(extents(), id);
}

inline auto Grid::Config::u(const Eigen::Vector3i& id) const -> std::size_t {
    return index(extents_u(), id);
}

inline auto Grid::Config::v(const Eigen::Vector3i& id) const -> std::size_t {
    return index(extents_v(), id);
}

inline auto Grid::Config::w(const Eigen::Vector3i& id) const -> std::size_t {
    return index(extents_w(), id);
}
