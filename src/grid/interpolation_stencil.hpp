#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <ranges>

#include "grid/core.hpp"

namespace grid {

struct InterpolationStencil {
    static constexpr std::size_t n_neighbours = 8;

    using Index = std::array<std::size_t, 3>;
    using Neighbours = std::array<Index, n_neighbours>;

   private:
    // we could do this smarter
    std::array<Index, n_neighbours> m_neighbours{};
    std::array<double, n_neighbours> m_weights{};

    static constexpr auto neighbour_offsets = std::array{
        // bottom 4 neighbours
        Index{0, 0, 0},
        Index{1, 0, 0},
        Index{0, 0, 1},
        Index{1, 0, 1},

        // top 4 neighbours
        Index{0, 1, 0},
        Index{1, 1, 0},
        Index{0, 1, 1},
        Index{1, 1, 1},
    };

   public:
    // if this does not work with Eigen stuff add a new constructor and convert
    // to std::span
    constexpr InterpolationStencil(
        std::span<const Config::float_type, 3> target_point,
        const Index& min_corner,
        std::span<const Config::float_type, 3> min_corner_pos,
        std::span<const Config::float_type, 3> cell_size) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            assert(cell_size[axis] > 0.0);
        }

        auto t = std::array<double, 3>{};
        for (std::size_t axis = 0; axis < 3; ++axis) {
            t[axis] =
                (target_point[axis] - min_corner_pos[axis]) / cell_size[axis];
            assert(t[axis] >= 0.0 && t[axis] <= 1.0);
        }

        for (std::size_t i = 0; i < n_neighbours; ++i) {
            const auto& curr_offset = neighbour_offsets[i];

            m_weights[i] = 1.0;

            for (std::size_t axis = 0; axis < 3; ++axis) {
                m_neighbours[i][axis] = min_corner[axis] + curr_offset[axis];

                m_weights[i] *=
                    curr_offset[axis] == 0 ? 1.0 - t[axis] : t[axis];
            }
        }
    }

    constexpr InterpolationStencil(
        const Eigen::Vector3<Config::float_type>& target_point,
        const Index& min_corner,
        const Eigen::Vector3<Config::float_type>& min_corner_pos,
        const Eigen::Vector3<Config::float_type>& cell_size)
        : InterpolationStencil{
              std::span<const Config::float_type, 3>{target_point.data(), 3},
              min_corner,
              std::span<const Config::float_type, 3>{min_corner_pos.data(), 3},
              std::span<const Config::float_type, 3>{cell_size.data(), 3},
          } {}

    static constexpr auto size() -> std::size_t { return n_neighbours; }
    constexpr const auto& weights() const { return m_weights; }
    constexpr const auto& neighbours() const { return m_neighbours; }
    constexpr auto zip() const {
        return std::ranges::zip_view(weights(), neighbours());
    }
};

}  // namespace grid
