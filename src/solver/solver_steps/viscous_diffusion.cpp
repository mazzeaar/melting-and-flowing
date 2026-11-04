#include "../solver.hpp"
#include <Eigen/IterativeLinearSolvers>

namespace solver {

// Step 2: Viscous diffusion
// u2 - dt * nu * ∇²u2   = u1 (scalar nu)
// (I - dt * nu * ∇²) u2 = u1
// or
// u2 - ∇ · (nu * ∇u2)   = u1 (vector nu)
// (I - ∇ · (nu * ∇)) u2 = u1

void Solver::diffuse_velocity(Solver::Config cfg) {

    auto u1 = cfg.grid.view<grid::Field::U>(); // X
    auto v1 = cfg.grid.view<grid::Field::V>(); // Y
    auto w1 = cfg.grid.view<grid::Field::W>(); // Z

    const auto& pcfg = pressure_cfg;

    const auto nx = pcfg.nx;
    const auto ny = pcfg.ny;
    const auto nz = pcfg.nz;
    const auto n = pcfg.n_unkowns();

    const auto dx = pcfg.dx;
    const auto dy = pcfg.dy;    
    const auto dz = pcfg.dz;
    const auto dt = cfg.dt;

    const double idx2 = 1.0 / (dx * dx);
    const double idy2 = 1.0 / (dy * dy);
    const double idz2 = 1.0 / (dz * dz);

    // Placeholder...
    const auto nu = 0.3;
    //Eigen::VectorXd nu = Eigen::VectorXd::Constant(n, 0.3);

    // Construct the RHS Vectors (for the interior cells)
    Eigen::VectorXd rhs_u = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd rhs_v = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd rhs_w = Eigen::VectorXd::Zero(n);

    // Construct the matrix A (for scalar nu)
    std::vector<Eigen::Triplet<float_type>> triplets;
    triplets.reserve(7 * n); 
    
    // Eigen is column major but cpp is row major, so... probably swap order?
    for (std::size_t i = 1; i <= nx; ++i) {
        for (std::size_t j = 1; j <= ny; ++j) {
            for (std::size_t k = 1; k <= nz; ++k) {
                const auto r = pcfg.row(i, j, k);

                // Populate RHS Vectors
                rhs_u(r) = u1[i, j, k];
                rhs_v(r) = v1[i, j, k];
                rhs_w(r) = w1[i, j, k];

                // Diagonal Element
                triplets.emplace_back(r, r, 1.0 + dt * nu * 2.0 * (idx2 + idy2 + idz2));

                // X-Axis Neighbors
                if (i > 1 ) { triplets.emplace_back(r, pcfg.row(i - 1, j, k), -dt * nu * idx2); }
                else { } 
                if (i < nx) { triplets.emplace_back(r, pcfg.row(i + 1, j, k), -dt * nu * idx2); } 
                else { } // Boundary Conditions

                // Y-Axis Neighbors
                if (j > 1)  { triplets.emplace_back(r, pcfg.row(i, j - 1, k), -dt * nu * idy2); }
                else { } 
                if (j < ny) { triplets.emplace_back(r, pcfg.row(i, j + 1, k), -dt * nu * idy2); }
                else { } 

                // Z-Axis Neighbors
                if (k > 1)  { triplets.emplace_back(r, pcfg.row(i, j, k - 1), -dt * nu * idz2); }
                else { } 
                if (k < nz) { triplets.emplace_back(r, pcfg.row(i, j, k + 1), -dt * nu * idz2); }
                else { } 
            }
        }
    }

    Eigen::SparseMatrix<float_type> A(n, n);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    // Solve for matrix A
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float_type>> eigen_solver;
    eigen_solver.compute(A);  

    auto u2 = eigen_solver.solve(rhs_u);
    auto v2 = eigen_solver.solve(rhs_v);
    auto w2 = eigen_solver.solve(rhs_w);
    
    // Ugly copying...
    for (std::size_t i = 1; i <= nx; ++i) {
        for (std::size_t j = 1; j <= ny; ++j) {
            for (std::size_t k = 1; k <= nz; ++k) {
                const auto r = pcfg.row(i, j, k);

                u1[i, j, k] = u2(r);
                v1[i, j, k] = v2(r);
                w1[i, j, k] = w2(r);
            }
        }
    }
}

}  // namespace solver 
