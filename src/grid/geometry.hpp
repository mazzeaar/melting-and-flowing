#pragma once

#include <cassert>
#include <cstddef>
#include <mdspan>
#include <ranges>
#include <utility>

#ifdef NDEBUG
#include <exception>
#include <format>
#endif

namespace grid {

using float_type   = double;
using extents_type = std::dextents<std::size_t, 3>;

enum class Field : std::uint8_t { U, V, W, P, COUNT };

enum class CellType : std::uint8_t {
    INTERIOR        = 0b0001,
    PLACEHOLDER     = 0b0010,
    DIRECT_BOUNDARY = 0b0100,
    GHOST_BOUNDARY  = 0b1000,
    BOUNDARY        = DIRECT_BOUNDARY | GHOST_BOUNDARY,
    ALL             = INTERIOR | PLACEHOLDER | BOUNDARY,
};

constexpr auto matches(CellType actual, CellType requested) -> bool {
    return (std::to_underlying(actual) & std::to_underlying(requested)) != 0;
}

/*
 Default CG coordinate axis according to chatgpt
      +Y, +j
       ^
       |
       |
       |
       O------------> +X, +i
      /
     /
    v
   +Z, +k (toward viewer)
*/
enum class Axis : std::uint8_t { X, Y, Z };
constexpr auto to_index(Axis ax) { return std::to_underlying(ax); }

namespace detail {

    using Index    = std::array<std::size_t, 3>;
    using Position = std::array<double, 3>;

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

    template <Location L> static constexpr auto offset() -> Position {
        if constexpr (L == Location::CELL) {
            return {0.5, 0.5, 0.5};
        } else if constexpr (L == Location::X_FACE) {
            return {0.0, 0.5, 0.5};
        } else if constexpr (L == Location::Y_FACE) {
            return {0.5, 0.0, 0.5};
        } else if constexpr (L == Location::Z_FACE) {
            return {0.5, 0.5, 0.0};
        } else if constexpr (L == Location::XY_EDGE) {
            return {0.0, 0.0, 0.5};
        } else if constexpr (L == Location::XZ_EDGE) {
            return {0.0, 0.5, 0.0};
        } else if constexpr (L == Location::YZ_EDGE) {
            return {0.5, 0.0, 0.0};
        } else {
            return {0.0, 0.0, 0.0};
        }
    }

    template <Location L>
    concept IsFaceValue = L == Location::X_FACE || L == Location::Y_FACE || L == Location::Z_FACE;

    template <Location L>
        requires(IsFaceValue<L>)
    consteval auto normal_axis() -> Axis {
        if constexpr (L == Location::X_FACE) {
            return Axis::X;
        } else if constexpr (L == Location::Y_FACE) {
            return Axis::Y;
        } else {
            return Axis::Z;
        }
    }

    template <Field F> struct FieldTraits;

#define SET_TRAIT(F, L)                         \
    template <> struct FieldTraits<F> {         \
        static constexpr Location location = L; \
    };

    SET_TRAIT(Field::U, Location::X_FACE);
    SET_TRAIT(Field::V, Location::Y_FACE);
    SET_TRAIT(Field::W, Location::Z_FACE);
    SET_TRAIT(Field::P, Location::CELL);

#undef SET_TRAIT

    template <Field F> struct FieldGeometry {
        static constexpr auto location = FieldTraits<F>::location;

        static constexpr auto classify(extents_type extents, Index index) -> CellType {
            static_assert(
                F == Field::P || IsFaceValue<location>,
                "if you can read this, then cell classification for this field type is not "
                "implemented yet in geometry.hpp classify()");

            static constexpr auto axis = std::ranges::iota_view(std::size_t{0}, std::size_t{3});

#ifdef NDEBUG
            for (auto ax : axis) {
                if (index[ax] >= extents.extent(ax)) {
                    throw std::logic_error(
                        std::format("Out of bounds cell classification for index "
                                    "({}, {}, {}), with extents ({}, {}, {})",
                                    index[0], index[1], index[2], extents.extent(0),
                                    extents.extent(1), extents.extent(2)));
                }
            }
#endif

            // pressure: outer layers are placeholders
            if constexpr (F == Field::P) {
                for (auto ax : axis) {
                    if (index[ax] == 0 || index[ax] == extents.extent(ax) - 1) {
                        return CellType::PLACEHOLDER;
                    }
                }

                return CellType::INTERIOR;
            }

            // velocities on cell faces
            if constexpr (IsFaceValue<location>) {
                static constexpr auto normal = std::to_underlying(normal_axis<location>());

                // unused outer position on the positive normal axis
                if (index[normal] == extents.extent(normal) - 1) {
                    return CellType::PLACEHOLDER;
                }

                // ghost layers along tangential axes
                for (auto ax : axis) {
                    if (ax == normal) {
                        continue;
                    }

                    if (index[ax] == 0 || index[ax] == extents.extent(ax) - 1) {
                        return CellType::GHOST_BOUNDARY;
                    }
                }

                // direct boundary values on normal walls
                if (index[normal] == 0 || index[normal] == extents.extent(normal) - 2) {
                    return CellType::DIRECT_BOUNDARY;
                }

                return CellType::INTERIOR;
            }
        }
    };

}  // namespace detail
}  // namespace grid
