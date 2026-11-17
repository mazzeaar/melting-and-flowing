// ################################################################################
// interpolation is done for single values right now. its by far the slowest
// part of our code. speedup might be achieved by doing batched processing but
// im lazy:)
// ################################################################################

#include "grid.hpp"

// ============================================================
// Trilinear interpolation
// each point in a (3D) field is surrounded by 8 values.
// interpolation is done by finding the lexicographically lowest index of those
// 8 values (i0).
// the weight (t) is the (3D) distance of the sample_point to i0.
// the value at the sample point is interpolated using the weight t and the
// anchor i0.

// ------------------------------------------------------------
// basic interpolation helper
// - field:     the field to be interpolated
// - base:      is linear base index (assumes i is contiguous)
// - nx, nxy:   offsets for y and z respectively
inline auto interpolate(const Eigen::VectorXf& field, std::size_t base,
                        std::size_t nx, std::size_t nxy,
                        const Eigen::Vector3f& t) -> float {
    // lower 4 values
    auto c000 = field[base];
    auto c100 = field[base + 1];
    auto c010 = field[base + nx];
    auto c110 = field[base + nx + 1];

    // upper 4 values
    auto c001 = field[base + nxy];
    auto c101 = field[base + nxy + 1];
    auto c011 = field[base + nxy + nx];
    auto c111 = field[base + nxy + nx + 1];

    // interpolate x-axis
    auto x00 = std::lerp(c000, c100, t.x());
    auto x10 = std::lerp(c010, c110, t.x());
    auto x01 = std::lerp(c001, c101, t.x());
    auto x11 = std::lerp(c011, c111, t.x());

    // interpolate y-axis
    auto y0 = std::lerp(x00, x10, t.y());
    auto y1 = std::lerp(x01, x11, t.y());

    // interpolate z-axis
    return std::lerp(y0, y1, t.z());
}

// ------------------------------------------------------------
// field interpolation
// - field:     the field to be interpolated
// - extents:   the extents (Nx, Ny, Nz) of the field
// - offset:    position of the value to be interpolated on the unit cube
// - grid_pos:  position where we want to interpolate
auto interpolate_field(const Eigen::VectorXf& field,
                       const Eigen::Vector3i& extents,
                       const Eigen::Vector3f& offset,
                       const Eigen::Vector3f& grid_pos) -> float {
    const Eigen::Vector3f p = grid_pos - offset;
    const Eigen::Vector3i i0 =
        p.array().floor().cast<int>().matrix().cwiseMax(0).cwiseMin(
            extents.array().matrix() - Eigen::Vector3i::Constant(2));

    const Eigen::Vector3f t =
        (p - i0.cast<float>()).cwiseMax(0.f).cwiseMin(1.f);

    const std::size_t nx = extents.x();
    const std::size_t nxy = nx * extents.y();

    const std::size_t base = static_cast<std::size_t>(i0.x()) +
                             nx * static_cast<std::size_t>(i0.y()) +
                             nxy * static_cast<std::size_t>(i0.z());

    return interpolate(field, base, nx, nxy, t);
}

auto interpolate_velocity(const Grid& grid, const Eigen::Vector3f& point)
    -> Eigen::Vector3f {
    const auto& cfg = grid.config();

    const Eigen::Vector3f pos =
        (point - cfg.min_corner()).cwiseQuotient(cfg.cell_size());

    return {
        interpolate_field(grid.u, cfg.extents_u(), {-1.f, -0.5f, -0.5f}, pos),
        interpolate_field(grid.v, cfg.extents_v(), {-0.5f, -1.f, -0.5f}, pos),
        interpolate_field(grid.w, cfg.extents_w(), {-0.5f, -0.5f, -1.f}, pos),
    };
}
