#pragma once

#include <concepts>
#include <cstddef>
#include <mdspan>
#include <memory>
#include <ranges>
#include <tuple>
#include <utility>

#include "geometry.hpp"
#include "storage.hpp"

namespace grid {

template <class F>
concept CellFunctionLike = std::invocable<F, std::size_t, std::size_t, std::size_t>;

// cheap to copy/move
class View {
    std::shared_ptr<Storage> data;

   public:
    View(std::shared_ptr<Storage>& grid) : data{grid} {}

    template <Field... Fields> constexpr auto get() {
        if constexpr (sizeof...(Fields) == 1) {
            return data->md_view<Fields...>();
        } else {
            return std::tuple{data->md_view<Fields>()...};
        }
    }

    template <Field F>
    constexpr auto cell_type(std::size_t i, std::size_t j, std::size_t k) const -> CellType {
        return detail::FieldGeometry<F>::classify(data->extents(), {i, j, k});
    }

    template <Field F>
    constexpr auto flat_index(std::size_t i, std::size_t j, std::size_t k) const -> std::size_t {
        return data->mapping()(i, j, k);
    }

    template <Axis Ax> constexpr auto offset() const {
        return data->mapping().stride(std::to_underlying(Ax));
    }

    template <Axis Ax> constexpr auto extent() const {
        return data->extents().extent(std::to_underlying(Ax));
    }

    template <Axis Ax> constexpr auto iextent() const {
        return std::ranges::iota_view{std::size_t{0},
                                      data->extents().extent(std::to_underlying(Ax))};
    }

    template <Field F, CellType T = CellType::ALL>
    constexpr void for_each(CellFunctionLike auto&& func) {
        for (auto i : iextent<Axis::X>()) {
            for (auto j : iextent<Axis::Y>()) {
                for (auto k : iextent<Axis::Z>()) {
                    if constexpr (T != CellType::ALL) {
                        if (!matches(cell_type<F>(i, j, k), T)) {
                            continue;
                        }
                    }

                    std::invoke(func, i, j, k);
                }
            }
        }
    }

    constexpr void for_each(CellFunctionLike auto&& func) {
        for (auto i : iextent<Axis::X>()) {
            for (auto j : iextent<Axis::Y>()) {
                for (auto k : iextent<Axis::Z>()) {
                    std::invoke(func, i, j, k);
                }
            }
        }
    }
};
}  // namespace grid
