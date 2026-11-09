#include "particles.hpp"

#include <igl/random_points_on_mesh.h>

#include "grid/grid.hpp"
#include "mesh/mesh.hpp"

Particles::Particles(Config conf, grid::Config grid_config,
                     const std::unique_ptr<Mesh>& mesh_ptr) {
    if (!mesh_ptr) {
        throw core::Error("Cannot initialize particles without a mesh");
    }

    if (conf.n_particles <= 0) {
        throw core::Error("Number of particles must be positive");
    }

    Eigen::MatrixXd barycentric;
    Eigen::VectorXi face_indices;

    igl::random_points_on_mesh(conf.n_particles, mesh_ptr->vertices,
                               mesh_ptr->faces, barycentric, face_indices,
                               positions);

    std::println("Generated {} surface particles", positions.rows());
}
