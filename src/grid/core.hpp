#pragma once

#include <Eigen/Core>
#include <cstdint>
#include <mdspan>
#include <utility>

namespace grid {

struct Config {
    using float_type = double;

    template <class T> using VecView = std::span<T, 3>;
    using Extents = std::dextents<std::size_t, 3>;
    using Layout = std::layout_right;
    using Mapping = Layout::mapping<Extents>;

    template <class T> using MdSpan = std::mdspan<T, Extents, Layout>;
    template <class T> using Span = std::span<T>;

    Extents grid_extents;
    Mapping grid_mapping;
    std::array<double, 3> grid_size;
    std::array<double, 3> lowest_corner = {0.0, 0.0, 0.0};

    Config(Extents extents, std::array<double, 3> grid_size,
           std::array<double, 3> lowest_corner)
        : grid_extents{extents},
          grid_mapping{extents},
          grid_size{grid_size},
          lowest_corner{lowest_corner} {}
};

struct GeometryContext {
    using Extents = Config::Extents;

    std::array<double, 3> lowest_corner;
    std::array<double, 3> cell_size;
    Extents extents;

    constexpr GeometryContext(const Config& cfg)
        : lowest_corner{cfg.lowest_corner},
          cell_size{
              cfg.grid_size[0] / static_cast<Config::float_type>(
                                     cfg.grid_extents.extent(0) - 2),
              cfg.grid_size[1] / static_cast<Config::float_type>(
                                     cfg.grid_extents.extent(1) - 2),
              cfg.grid_size[2] / static_cast<Config::float_type>(
                                     cfg.grid_extents.extent(2) - 2),
          },
          extents{cfg.grid_extents} {}
};

// all fields which we store/represent
enum class Field : std::uint8_t {
    U,
    V,
    W,
    P,
    T,
    COUNT,
    Velocities,  // combines all 3 velocities into one, for solvers
    State,
};

// von Neumann neighbourhood
enum class CellState {
    SOLID,    // boundaries/walls
    AIR,      // no particles
    FULL,     // #particles > 0 && all neighbours not empty
    SURFACE,  // #particles > 0 && (#empty neighbours > 0)
    UNKNOWN,
};

namespace detail {
template <Field F> struct FieldTypeHelper {
    static_assert(F != F, "Unhandled field");
};

template <Field F>
    requires(F != Field::State)
struct FieldTypeHelper<F> {
    using type = Config::float_type;
};

template <> struct FieldTypeHelper<Field::State> {
    using type = CellState;
};
}  // namespace detail

template <Field F> using FieldType = typename detail::FieldTypeHelper<F>::type;

enum class Axis : std::uint8_t { X, Y, Z };

enum class CellType : std::uint8_t {
    INTERIOR = 0b0001,
    PLACEHOLDER = 0b0010,
    DIRECT_BOUNDARY = 0b0100,
    GHOST_BOUNDARY = 0b1000,
    BOUNDARY = DIRECT_BOUNDARY | GHOST_BOUNDARY,
    ALL = INTERIOR | PLACEHOLDER | BOUNDARY,
};

namespace detail {
template <class From, class To>
using copy_const_t =
    std::conditional_t<std::is_const_v<std::remove_reference_t<From>>, const To,
                       To>;

template <Field F>
concept GeometricField = F == Field::U || F == Field::V || F == Field::W ||
                         F == Field::P || F == Field::T || F == Field::State;

constexpr auto matches(CellType actual, CellType requested) -> bool {
    return (std::to_underlying(actual) & std::to_underlying(requested)) != 0;
}

enum class Location : std::uint8_t {
    CELL,
    X_FACE,
    Y_FACE,
    Z_FACE,
    XY_EDGE,
    XZ_EDGE,
    YZ_EDGE,
    NODE
};

template <Location L>
consteval auto cell_position() -> std::array<Config::float_type, 3> {
    if constexpr (L == Location::CELL) {
        return {0.5, 0.5, 0.5};
    } else if constexpr (L == Location::X_FACE) {
        return {1.0, 0.5, 0.5};
    } else if constexpr (L == Location::Y_FACE) {
        return {0.5, 1.0, 0.5};
    } else if constexpr (L == Location::Z_FACE) {
        return {0.5, 0.5, 1.0};
    } else if constexpr (L == Location::XY_EDGE) {
        return {1.0, 1.0, 0.5};
    } else if constexpr (L == Location::XZ_EDGE) {
        return {1.0, 0.5, 1.0};
    } else if constexpr (L == Location::YZ_EDGE) {
        return {0.5, 1.0, 1.0};
    } else if constexpr (L == Location::NODE) {
        return {1.0, 1.0, 1.0};
    } else {
        static_assert(false,
                      "Did you add a new grid::Location enum entry and not "
                      "handle it in "
                      "cell_position()?:)");
    }
}

template <Field F>
consteval auto field_location() -> Location
    requires(GeometricField<F>)
{
    if constexpr (F == Field::U) {
        return Location::X_FACE;
    } else if constexpr (F == Field::V) {
        return Location::Y_FACE;
    } else if constexpr (F == Field::W) {
        return Location::Z_FACE;
    } else if constexpr (F == Field::P || F == Field::T || F == Field::State) {
        return Location::CELL;
    } else {
        static_assert(false,
                      "You forgot to handle the new Field enum value in "
                      "field_location()!");
    }
}

}  // namespace detail
}  // namespace grid
