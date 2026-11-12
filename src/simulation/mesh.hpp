#pragma once

#include <igl/fast_winding_number.h>
#include <igl/readOBJ.h>

#include <Eigen/Core>
#include <cmath>
#include <filesystem>

#include "core/exception.hpp"

struct Mesh {
    static constexpr auto wn_accuracy = 2.0f;
    static constexpr int taylor_order = 2;
    static constexpr float inside_threshold = 0.5;

    std::filesystem::path file_path;

    Eigen::MatrixXf vertices;
    Eigen::MatrixXi faces;

    igl::FastWindingNumberBVH bvh;

    Mesh() = default;

    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;

    Mesh(Mesh&&) = delete;
    Mesh& operator=(Mesh&&) = delete;

    explicit Mesh(const std::filesystem::path& path) : file_path{path} {
        if (!igl::readOBJ(path.string(), vertices, faces)) {
            throw core::Error("Failed to load mesh: " + path.string());
        }

        igl::fast_winding_number(vertices, faces, taylor_order, bvh);
    }

    void fit_to_box(const Eigen::Vector3f& box_min,
                    const Eigen::Vector3f& box_max, float padding = 0.1) {
        if (vertices.rows() == 0) {
            throw core::Error("Cannot rescale an empty mesh");
        }

        if (padding < 0.0 || padding >= 0.5) {
            throw core::Error("Invalid mesh padding");
        }

        const Eigen::Vector3f mesh_min = vertices.colwise().minCoeff();
        const Eigen::Vector3f mesh_max = vertices.colwise().maxCoeff();

        const Eigen::Vector3f mesh_size = mesh_max - mesh_min;
        const Eigen::Vector3f box_size = box_max - box_min;

        if ((mesh_size.array() <= 0.0).any() ||
            (box_size.array() <= 0.0).any()) {
            throw core::Error("Invalid mesh or box dimensions");
        }

        const Eigen::Vector3f available_size = box_size * (1.0 - 2.0 * padding);

        const auto scale =
            (available_size.array() / mesh_size.array()).minCoeff();

        const Eigen::Vector3f mesh_center = (mesh_min + mesh_max) * 0.5;
        const Eigen::Vector3f box_center = (box_min + box_max) * 0.5;

        vertices.rowwise() -= mesh_center.transpose();
        vertices *= scale;
        vertices.rowwise() += box_center.transpose();

        igl::fast_winding_number(vertices, faces, taylor_order, bvh);
    }

    auto contains(const Eigen::Vector3f& point) const -> bool {
        Eigen::VectorXf winding;
        igl::fast_winding_number(vertices, faces, point, winding);
        return std::abs(winding(0)) > 0.5;
    }

    [[nodiscard]] auto contains_many(const Eigen::MatrixXf& points) const
        -> Eigen::Array<bool, Eigen::Dynamic, 1> {
        Eigen::VectorXf winding;
        igl::fast_winding_number(bvh, wn_accuracy, points, winding);

        return (winding.array().abs() > inside_threshold).eval();
    }
};
