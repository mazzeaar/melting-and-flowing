#pragma once

#include <Eigen/Core>
#include <cassert>
#include <cstddef>

#include "core.hpp"

namespace grid::detail {
template <Location L> struct Geometry {
  using float_type = Config::float_type;

  static consteval auto offset() { return detail::cell_position<L>(); };

  [[nodiscard]] static constexpr auto position(std::array<std::size_t, 3> idx,
                                               const GeometryContext& ctx) {
    static constexpr auto off = offset();
    std::array<float_type, 3> out;
    for (auto i : {0, 1, 2}) {
      // cell (1,1,1) has its lower physical corner at ctx.lowest_corner
      out[i] =
          ctx.lowest_corner[i] +
          (static_cast<float_type>(idx[i]) - 1.0 + off[i]) * ctx.cell_size[i];
    }
    return out;
  }

  // DOES NOT DEPEND ON LOCATION AT ALL!!
  // index of the cell containing the point
  [[nodiscard]] static constexpr auto cell(std::span<const double, 3> point,
                                           const GeometryContext& ctx)
      -> std::array<std::size_t, 3> {
    std::array<std::size_t, 3> out;
    for (auto i : {0, 1, 2}) {
      const auto index =
          std::floor((point[i] - ctx.lowest_corner[i]) / ctx.cell_size[i]) +
          1.0;

      assert(std::isfinite(index));
      assert(index >= 0.0);
      assert(index < static_cast<float_type>(ctx.extents.extent(i)));

      out[i] = static_cast<std::size_t>(index);
    }
    return out;
  };

  // nearest lowest neighbour
  // meaning: for any point inside the grid there are 8 values closest to it.
  //          -> the point is contained in a cube formed by those 8 values
  // this function returns the lowest index of this cube
  [[nodiscard]] static constexpr auto min_corner(
      std::span<const double, 3> point, const GeometryContext& ctx)
      -> std::array<std::size_t, 3> {
    constexpr auto off = offset();
    std::array<std::size_t, 3> out;
    for (auto i : {0, 1, 2}) {
      const auto index =
          std::floor((point[i] - ctx.lowest_corner[i]) / ctx.cell_size[i] -
                     off[i]) +
          1.0;

      assert(std::isfinite(index));
      assert(index >= 0.0);
      assert(index + 1.0 < static_cast<float_type>(ctx.extents.extent(i)));

      out[i] = static_cast<std::size_t>(index);
    }
    return out;
  }

  static constexpr auto cell_type(const std::array<std::size_t, 3>& idx,
                                  const GeometryContext& ctx) -> CellType {
    constexpr auto off = offset();

    for (std::size_t axis = 0; axis < 3; ++axis) {
      assert(ctx.extents.extent(axis) >= 3);
      assert(idx.at(axis) < ctx.extents.extent(axis));
    }

    // no early returns, because if we first hit a boundary cell, but later we
    // notice its actually a placeholder we return the wrong thing. order is:
    // placeholder > ghost > direct > interior
    // -> V(Nx, Ny, Nz)
    //      - first axis: ghost
    //      - second: placeholder
    bool is_ghost = false;
    bool is_direct = false;

    for (std::size_t axis = 0; axis < 3; ++axis) {
      const auto index = idx.at(axis);
      const auto last = ctx.extents.extent(axis) - 1;

      if constexpr (L == Location::CELL) {
        // Cell-centered fields have no samples on physical walls.
        // Their surrounding padding is classified as placeholder.
        if (index == 0 || index == last) {
          return CellType::PLACEHOLDER;
        }
      } else {
        if (off.at(axis) ==
            1.0) {  // location lies on cell boundary along this axis
          // the final plane in any direction is unused
          if (index == last) {
            return CellType::PLACEHOLDER;
          }

          // lies on the physical boundary in this direction
          if (index == 0 || index == last - 1) {
            is_direct = true;
          }

        } else {
          // location is centered along this axis
          // -> out planes provide ghost values (indirect boundary)
          if (index == 0 || index == last) {
            is_ghost = true;
          }
        }
      }
    }

    if (is_ghost) {
      return CellType::GHOST_BOUNDARY;
    }

    if (is_direct) {
      return CellType::DIRECT_BOUNDARY;
    }

    return CellType::INTERIOR;
  }
};
}  // namespace grid::detail
