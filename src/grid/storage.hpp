#pragma once

#include <cstddef>
#include <mdspan>
#include <span>
#include <type_traits>
#include <vector>

#include "geometry.hpp"

namespace grid {

namespace detail {
    template <class From, class To>
    using copy_const_t =
        std::conditional_t<std::is_const_v<std::remove_reference_t<From>>, const To, To>;
}

// owns state/allocates buffers
struct Storage {
   private:
    extents_type                             grid_extents;
    std::layout_right::mapping<extents_type> grid_mapping;
    std::size_t                              number_of_cells;

    using storage_type = std::vector<float_type>;
    storage_type u_field;
    storage_type v_field;
    storage_type w_field;
    storage_type p_field;
    // add more grids here

    // only access fields through this (inside of Storage)
    template <Field F> constexpr decltype(auto) field(this auto&& self) {
        using enum Field;

        if constexpr (F == U) {
            return (self.u_field);
        } else if constexpr (F == V) {
            return (self.v_field);
        } else if constexpr (F == W) {
            return (self.w_field);
        } else if constexpr (F == P) {
            return (self.p_field);
        } else {
            static_assert(F == U || F == V || F == W || F == P, "Unsupported field");
        }
    }

   public:
    Storage(std::size_t nx, std::size_t ny, std::size_t nz)
        : grid_extents{nx + 1, ny + 1, nz + 1}, grid_mapping{grid_extents},
          number_of_cells{buffer_size()}, u_field(number_of_cells, 0.0),
          v_field(number_of_cells, 0.0), w_field(number_of_cells, 0.0),
          p_field(number_of_cells, 0.0) {}

    constexpr auto extents() const { return grid_extents; }
    constexpr auto mapping() const { return grid_mapping; }
    constexpr auto buffer_size() const -> std::size_t { return mapping().required_span_size(); }

    // flat view
    template <Field F> constexpr auto flat_view(this auto&& self) {
        auto& values = self.template field<F>();
        return std::span{values};
    }

    // 3d view, including boundary and placeholder values
    template <Field F> constexpr auto md_view(this auto&& self) {
        auto flat = self.template flat_view<F>();
        return std::mdspan{flat.data(), self.grid_extents};
    }
};

}  // namespace grid
