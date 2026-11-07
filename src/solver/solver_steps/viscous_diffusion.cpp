#include "../solver.hpp"
#include <Eigen/IterativeLinearSolvers>

namespace solver {

// Step 2: Viscous Diffusion
// u2 - dt * nu * ∇²u2   = u1 (scalar nu)
// (I - dt * nu * ∇²) u2 = u1
// or
// u2 - ∇ · (nu * ∇u2)   = u1 (vector nu)
// (I - ∇ · (nu * ∇)) u2 = u1

enum class WallCondition : std::uint8_t {
    FreeSlip,
    NoSlip
};

enum class Dimension : std::uint8_t {
    X,
    Y,
    Z
};

void Solver::diffuse_velocity(Solver::Config cfg) {

    auto u1 = cfg.grid.view<grid::Field::U>(); // X
    auto v1 = cfg.grid.view<grid::Field::V>(); // Y
    auto w1 = cfg.grid.view<grid::Field::W>(); // Z
    auto state = 0; //cfg.grid.view<grid::Field::State>(); or sthm

    const auto& pcfg = pressure_cfg;

    const auto nx = pcfg.nx;
    const auto ny = pcfg.ny;
    const auto nz = pcfg.nz;
    const auto n = pcfg.n_unkowns();

    // Placeholder...
    const auto nu = 0.3;
    //Eigen::VectorXd nu = Eigen::VectorXd::Constant(n, 0.3);
    WallCondition bc = WallCondition::NoSlip;

    // Construct the RHS Vectors (just for the interior cells)
    Eigen::VectorXd rhs_u = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd rhs_v = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd rhs_w = Eigen::VectorXd::Zero(n);

    // Prepare Triplets for Matrix A
    std::vector<Eigen::Triplet<float_type>> triplets_u;
    std::vector<Eigen::Triplet<float_type>> triplets_v;
    std::vector<Eigen::Triplet<float_type>> triplets_w;
    triplets_u.reserve(7 * n); 
    triplets_v.reserve(7 * n); 
    triplets_w.reserve(7 * n); 
    
    // Eigen is column major but cpp is row major, so...?
    for (std::size_t i = 1; i <= nx; ++i) {
        for (std::size_t j = 1; j <= ny; ++j) {
            for (std::size_t k = 1; k <= nz; ++k) {
                const auto r = pcfg.row(i, j, k);
 
                // Populate RHS Vectors
                rhs_u(r) = u1[i, j, k];
                rhs_v(r) = v1[i, j, k];
                rhs_w(r) = w1[i, j, k];

                // Populate Triplets for A
                triplets(triplets_u, pcfg, state, nu, bc, Dimension::X, i, j, k);
                triplets(triplets_v, pcfg, state, nu, bc, Dimension::Y, i, j, k);
                triplets(triplets_w, pcfg, state, nu, bc, Dimension::Z, i, j, k);
            }
        }
    }

    // Construct Matrix A
    Eigen::SparseMatrix<float_type> A_u(n, n);
    Eigen::SparseMatrix<float_type> A_v(n, n);
    Eigen::SparseMatrix<float_type> A_w(n, n);
    A_u.setFromTriplets(triplets_u.begin(), triplets_u.end());
    A_v.setFromTriplets(triplets_v.begin(), triplets_v.end());
    A_w.setFromTriplets(triplets_w.begin(), triplets_w.end());
    A_u.makeCompressed();
    A_v.makeCompressed();
    A_w.makeCompressed();

    // Solve with Jacobi preconditioner (default)
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float_type>> solver_u;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float_type>> solver_v;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float_type>> solver_w;
    solver_u.compute(A_u);  
    solver_v.compute(A_v); 
    solver_w.compute(A_w);  
    auto u2 = solver_u.solve(rhs_u);
    auto v2 = solver_v.solve(rhs_v);      
    auto w2 = solver_w.solve(rhs_w);

    if (solver_u.info() != Eigen::Success) { throw std::runtime_error("Viscous Diffusion: Solver failed for u-component."); }
    if (solver_v.info() != Eigen::Success) { throw std::runtime_error("Viscous Diffusion: Solver failed for v-component."); }
    if (solver_w.info() != Eigen::Success) { throw std::runtime_error("Viscous Diffusion: Solver failed for w-component."); }

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

void triplets(auto& triplets, const auto& pcfg, const auto& state,
        const auto nu, const WallCondition bc, const Dimension dim,
        const std::size_t i, const std::size_t j, const std::size_t k) {

    const float_type dt = cfg.dt;
    const float_type dx = pcfg.dx;
    const float_type dy = pcfg.dy;    
    const float_type dz = pcfg.dz;

    const float_type idx2 = float_type{1} / (dx * dx);
    const float_type idy2 = float_type{1} / (dy * dy);
    const float_type idz2 = float_type{1} / (dz * dz);

    const size_t r = pcfg.row(i, j, k);

    // Diagonal Element 
    triplets.emplace_back(r, r, 1 + dt * nu * 2 * (idx2 + idy2 + idz2));

    // X-Axis Neighbors
    if (i > 1 ) { triplets.emplace_back(r, pcfg.row(i - 1, j, k), -dt * nu * idx2); }
    else { } 
    if (i < nx) { triplets.emplace_back(r, pcfg.row(i + 1, j, k), -dt * nu * idx2); } 
    else { } // Boundary Conditions...

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

}  // namespace solver 
