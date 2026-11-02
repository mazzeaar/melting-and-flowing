#pragma once

#include <igl/fast_winding_number.h>
#include <igl/readOBJ.h>

#include <Eigen/Core>
#include <cmath>
#include <filesystem>
#include <memory>
#include <stdexcept>

#pragma once

struct Mesh {
    static constexpr auto wn_accuracy = 2.0f;
    static constexpr int taylor_order = 2;
    static constexpr double inside_threshold = 0.5;

    std::filesystem::path file_path;

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;

    igl::FastWindingNumberBVH bvh;

    Mesh() = default;

    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;

    Mesh(Mesh&&) = delete;
    Mesh& operator=(Mesh&&) = delete;

    explicit Mesh(const std::filesystem::path& path) : file_path{path} {
        if (!igl::readOBJ(path.string(), vertices, faces)) {
            // TODO: replace with custom exception
            throw std::runtime_error("Failed to load mesh: " + path.string());
        }

        igl::fast_winding_number(vertices, faces, taylor_order, bvh);
    }

    template <class Derived>
    auto contains(const Eigen::MatrixBase<Derived>& point) const -> bool {
        const auto w = igl::fast_winding_number(bvh, wn_accuracy, point.eval());

        return std::abs(w) > inside_threshold;
    }

    auto contains_many(const Eigen::MatrixXd& points) const
        -> Eigen::Array<bool, Eigen::Dynamic, 1> {
        Eigen::VectorXd winding;
        igl::fast_winding_number(bvh, wn_accuracy, points, winding);

        return (winding.array().abs() > inside_threshold).eval();
    }
};
