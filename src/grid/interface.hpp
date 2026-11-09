#pragma once

#include <type_traits>

#include "core.hpp"
#include "geometry.hpp"
#include "interpolation_stencil.hpp"

namespace grid {

// non owning view of a field
template <Field FieldTag> struct FieldInterface {
    using value_type = FieldType<FieldTag>;
    static constexpr auto field_tag = FieldTag;
    static constexpr auto is_geometric = detail::GeometricField<field_tag>;

   private:
    std::span<value_type> data_view;
    // removed Config::Mapping and replacing it with just Config
    Config cfg;
    GeometryContext geo_context;

   public:
    constexpr FieldInterface(std::span<value_type> data_view, Config config)
        : data_view{data_view}, cfg{config}, geo_context{cfg} {}

    constexpr auto wrap(std::span<value_type> other_data) const
        -> FieldInterface {
        return {other_data, cfg};
    }

    // flat indexing
    constexpr auto operator[](this auto&& self, std::size_t g)
        -> decltype(auto) {
        return self.data_view[g];
    }

    // 3d indexing
    constexpr auto operator[](this auto&& self, std::size_t i, std::size_t j,
                              std::size_t k) -> decltype(auto)
        requires(is_geometric)
    {
        return self[self.flat_index(i, j, k)];
    }

    constexpr auto operator[](this auto&& self, std::size_t i, std::size_t j,
                              std::size_t k) -> decltype(auto)
        requires(field_tag == Field::Velocities)
    {
        return self[self.flat_index(i, j, k)];
    }

    constexpr auto size() const -> std::size_t { return data_view.size(); }
    constexpr auto mapping() const -> Config::Mapping {
        return cfg.grid_mapping;
    };
    constexpr auto extents() const -> Config::Extents {
        return mapping().extents();
    }

    // is equivalent for all fields associated with a cell (i,j,k)
    constexpr auto flat_index(std::size_t i, std::size_t j, std::size_t k) const
        -> std::size_t {
        return mapping()(i, j, k);
    };

    // flat stepsize in +A direction
    template <Axis A> constexpr auto flat_offset() const -> std::size_t {
        return mapping().stride(std::to_underlying(A));
    };

    // all flat stepsizes in order
    constexpr auto flat_offsets() const {
        return std::array{flat_offset<Axis::X>(), flat_offset<Axis::Y>(),
                          flat_offset<Axis::Z>()};
    }

    template <typename Derived>
    constexpr auto index_at(const Eigen::MatrixBase<Derived>& p) const
        -> Eigen::Vector3i
        requires(detail::GeometricField<FieldTag>)
    {
        static_assert(Derived::SizeAtCompileTime == 3);
        static constexpr auto Loc = detail::field_location<FieldTag>();
        const Eigen::Vector3d point = p;
        const auto idx = detail::Geometry<Loc>::cell(
            std::span<const double, 3>{point.data(), 3}, geo_context);
        return {static_cast<int>(idx[0]), static_cast<int>(idx[1]),
                static_cast<int>(idx[2])};
    }

    constexpr auto position_at(std::size_t i, std::size_t j,
                               std::size_t k) const
        -> Eigen::Vector3<Config::float_type>
        requires(detail::GeometricField<FieldTag>)
    {
        static constexpr auto Loc = detail::field_location<FieldTag>();
        auto pos = detail::Geometry<Loc>::position({i, j, k}, geo_context);
        return {pos[0], pos[1], pos[2]};
    }

    // todo, use Geometry<F> for this?
    // should return boundary, or interior etc
    // maybe bad name? cell_type could also be liquid, air etc
    constexpr auto cell_type(std::size_t i, std::size_t j, std::size_t k)
        -> CellType {
        static constexpr auto Loc = detail::field_location<FieldTag>();
        return detail::Geometry<Loc>::cell_type({i, j, k}, {cfg});
    }

    constexpr auto contains_point(
        const Eigen::Vector3<grid::Config::float_type>& point) const -> bool
        requires(detail::GeometricField<FieldTag>)
    {
        constexpr auto Loc = detail::field_location<FieldTag>();
        return detail::Geometry<Loc>::contains_point(
            std::span<const Config::float_type, 3>{point.data(), 3}, {cfg});
    }

    constexpr auto interpolation_stencil(
        std::span<const Config::float_type, 3> point) const
        -> InterpolationStencil
        requires(detail::GeometricField<FieldTag>)
    {
        constexpr auto Loc = detail::field_location<FieldTag>();
        using Geo = detail::Geometry<Loc>;

        const auto min_corner = Geo::min_corner(point, geo_context);
        const auto min_corner_pos = Geo::position(min_corner, geo_context);

        return InterpolationStencil{point, min_corner, min_corner_pos,
                                    geo_context.cell_size};
    }

    constexpr auto interpolation_stencil(
        const Eigen::Vector3<Config::float_type>& point) const
        -> InterpolationStencil
        requires(detail::GeometricField<FieldTag>)
    {
        return interpolation_stencil(
            std::span<const Config::float_type, 3>{point.data(), 3});
    }

    // --------------------------------
    // data views
    constexpr auto flat_view(this auto&& self)
        -> std::span<detail::copy_const_t<decltype(self), value_type>> {
        return self.data_view;
    }

    [[nodiscard]]
    constexpr auto md_view(this auto&& self)
        -> Config::MdSpan<detail::copy_const_t<decltype(self), value_type>>
        requires(detail::GeometricField<FieldTag>)
    {
        return {
            self.data_view.data(),
            self.mapping(),
        };
    }

    auto eigen_view(this auto&& self) {
        using Scalar = std::remove_const_t<value_type>;
        using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

        if constexpr (std::is_const_v<
                          std::remove_reference_t<decltype(self)>> ||
                      std::is_const_v<value_type>) {
            return Eigen::Map<const Vector>{
                self.data_view.data(),
                static_cast<Eigen::Index>(self.data_view.size())};
        } else {
            return Eigen::Map<Vector>{
                self.data_view.data(),
                static_cast<Eigen::Index>(self.data_view.size())};
        }
    }
};

}  // namespace grid
