#include "surface_reconstruction.hpp"

#include <igl/marching_cubes.h>

#include <algorithm>

void SurfaceReconstructor::reconstruct(const Eigen::MatrixXf& points,
                                       const Grid::Config& config,
                                       Eigen::MatrixXf& V, Eigen::MatrixXi& F) {
    benchmark.tick("mesh_total");

    // Surface reconstruction after Zhu & Bridson 2005, "Animating Sand as a
    // Fluid", eq. 7-10:

    const Eigen::Vector3f min_corner = config.min_corner();

    const Eigen::Vector3f node_spacing =
        config.cell_size() / static_cast<float>(resolution_upfactor);
    const Eigen::Vector3i node_counts =
        config.extents() * resolution_upfactor + Eigen::Vector3i::Ones();

    const int nx = node_counts.x();
    const int ny = node_counts.y();
    const int nz = node_counts.z();
    const int n_nodes = nx * ny * nz;

    auto kernel = [](float s2) {
        auto a = 1.f - s2;
        return std::max(0.f, a * a * a);
    };

    auto node_index = [&](int a, int b, int c) {
        return a + nx * (b + ny * c);
    };

    auto node_position = [&](int a, int b, int c) -> Eigen::Vector3f {
        return min_corner + node_spacing.cwiseProduct(
                                Eigen::Vector3i(a, b, c).cast<float>());
    };

    const float particle_radius = radius_factor * node_spacing.maxCoeff();
    const float kernel_radius = 2.f * particle_radius;
    const float kernel_radius_sq = kernel_radius * kernel_radius;
    const float inv_kernel_radius_sq = 1.f / kernel_radius_sq;

    // Splat every particle into the fine grid nodes within kernel_radius
    benchmark.tick("mesh_splat");
    Eigen::VectorXf weight_sum = Eigen::VectorXf::Zero(n_nodes);
    Eigen::MatrixXf weighted_pos = Eigen::MatrixXf::Zero(n_nodes, 3);

    const Eigen::Index stride = std::max(1, marker_stride);

    for (Eigen::Index i = 0; i < points.rows(); i += stride) {
        const Eigen::Vector3f point = points.row(i).transpose();

        // Range of node indices within kernel_radius, clamped to the grid
        const Eigen::Vector3f local = point - min_corner;
        Eigen::Vector3i lo =
            ((local.array() - kernel_radius) / node_spacing.array())
                .ceil()
                .cast<int>()
                .matrix();
        Eigen::Vector3i hi =
            ((local.array() + kernel_radius) / node_spacing.array())
                .floor()
                .cast<int>()
                .matrix();
        lo = lo.cwiseMax(Eigen::Vector3i::Zero());
        hi = hi.cwiseMin(node_counts - Eigen::Vector3i::Ones());

        const Eigen::Vector3f offset = min_corner - point;

        // Mesh cells around particles
        for (int c = lo.z(); c <= hi.z(); ++c) {
            const float dz = offset.z() + node_spacing.z() * c;
            const float dz_sq = dz * dz;
            if (dz_sq >= kernel_radius_sq) {
                continue;
            }

            for (int b = lo.y(); b <= hi.y(); ++b) {
                const float dy = offset.y() + node_spacing.y() * b;
                const float dyz_sq = dz_sq + dy * dy;
                if (dyz_sq >= kernel_radius_sq) {
                    continue;
                }

                const auto row = node_index(0, b, c);

                for (int a = lo.x(); a <= hi.x(); ++a) {
                    const float dx = offset.x() + node_spacing.x() * a;
                    const float dist_sq = dyz_sq + dx * dx;
                    if (dist_sq >= kernel_radius_sq) {
                        continue;
                    }

                    const float k = kernel(dist_sq * inv_kernel_radius_sq);

                    const auto n = row + a;
                    weight_sum(n) += k;
                    weighted_pos.row(n) += k * point.transpose();
                }
            }
        }
    }
    benchmark.tock("mesh_splat");

    benchmark.tick("mesh_phi");
    Eigen::VectorXf phi(n_nodes);
    Eigen::MatrixXf node_positions(n_nodes, 3);

    for (int c = 0; c < nz; ++c) {
        for (int b = 0; b < ny; ++b) {
            for (int a = 0; a < nx; ++a) {
                const auto n = node_index(a, b, c);
                const Eigen::Vector3f x = node_position(a, b, c);
                node_positions.row(n) = x.transpose();

                // Outermost layer is forced outside to close mesh
                const bool on_boundary = a == 0 || b == 0 || c == 0 ||
                                         a == nx - 1 || b == ny - 1 ||
                                         c == nz - 1;

                // No particle within kernel_radius -> outside
                if (on_boundary || weight_sum(n) == 0.f) {
                    phi(n) = kernel_radius;
                    continue;
                }

                const Eigen::Vector3f x_bar =
                    weighted_pos.row(n).transpose() / weight_sum(n);

                phi(n) = (x - x_bar).norm() - particle_radius;
            }
        }
    }
    benchmark.tock("mesh_phi");

    benchmark.tick("mesh_marching_cubes");
    igl::marching_cubes(phi, node_positions, static_cast<unsigned>(nx),
                        static_cast<unsigned>(ny), static_cast<unsigned>(nz),
                        0.f, V, F);
    benchmark.tock("mesh_marching_cubes");

    benchmark.tock("mesh_total");
}
