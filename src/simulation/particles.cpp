#include "particles.hpp"

#include <igl/random_points_on_mesh.h>

#include "core/log.hpp"
#include "mesh.hpp"

Particles::Particles(std::size_t num_particles,
                     const std::unique_ptr<Mesh>& mesh_ptr)
    : positions{Eigen::MatrixXf::Zero(num_particles, 3)} {
    if (!mesh_ptr) {
        throw core::Error("Cannot initialize particles without a mesh");
    }

    if (num_particles <= 0) {
        throw core::Error("Number of particles must be positive");
    }

    Eigen::MatrixXf barycentric;
    Eigen::VectorXi face_indices;

    igl::random_points_on_mesh(num_particles, mesh_ptr->vertices,
                               mesh_ptr->faces, barycentric, face_indices,
                               positions);

    Log("Generated {} surface particles", positions.rows());
}
