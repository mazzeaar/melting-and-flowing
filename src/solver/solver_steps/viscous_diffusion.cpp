#include <Eigen/IterativeLinearSolvers>

#include "../solver.hpp"

namespace solver {
using float_type = Solver::float_type;

// Step 2: Viscous Diffusion
// u2 - dt * nu * ∇²u2   = u1 (scalar nu)
// (I - dt * nu * ∇²) u2 = u1
// or
// u2 - ∇ · (nu * ∇u2)   = u1 (vector nu)
// (I - ∇ · (nu * ∇)) u2 = u1

enum class WallCondition : std::uint8_t { FreeSlip, NoSlip };

enum class Dimension : std::uint8_t { X, Y, Z };

void triplets(auto& triplets, const auto& pcfg, const auto& state,
              const auto nu, const WallCondition bc, const Dimension dim,
              const float_type dt, const std::size_t i, const std::size_t j,
              const std::size_t k) {
    const size_t n = pcfg.n_unkowns();
    const size_t nx = pcfg.nx;
    const size_t ny = pcfg.ny;
    const size_t nz = pcfg.nz;

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
    if (i > 1) {
        triplets.emplace_back(r, pcfg.row(i - 1, j, k), -dt * nu * idx2);
    } else {
    }
    if (i < nx) {
        triplets.emplace_back(r, pcfg.row(i + 1, j, k), -dt * nu * idx2);
    } else {
    }  // Boundary Conditions...

    // Y-Axis Neighbors
    if (j > 1) {
        triplets.emplace_back(r, pcfg.row(i, j - 1, k), -dt * nu * idy2);
    } else {
    }
    if (j < ny) {
        triplets.emplace_back(r, pcfg.row(i, j + 1, k), -dt * nu * idy2);
    } else {
    }

    // Z-Axis Neighbors
    if (k > 1) {
        triplets.emplace_back(r, pcfg.row(i, j, k - 1), -dt * nu * idz2);
    } else {
    }
    if (k < nz) {
        triplets.emplace_back(r, pcfg.row(i, j, k + 1), -dt * nu * idz2);
    } else {
    }
}

void Solver::diffuse_velocity(Solver::Config cfg) {
    const float_type dt = cfg.dt;

    auto u1 = cfg.grid.view<grid::Field::U>();
    auto v1 = cfg.grid.view<grid::Field::V>();
    auto w1 = cfg.grid.view<grid::Field::W>();
    auto state = 0;  // TODO: cell states

    const auto& pcfg = pressure_cfg;

    const std::size_t nx = pcfg.nx;
    const std::size_t ny = pcfg.ny;
    const std::size_t nz = pcfg.nz;

    // TODO: replace with actual nu
    const float_type nu = 0.3;
    // Eigen::VectorXd nu = Eigen::VectorXd::Constant(n, 0.3);
    const WallCondition bc = WallCondition::NoSlip;

    // prep buffers
    viscosity_triplets_u.clear();
    viscosity_triplets_v.clear();
    viscosity_triplets_w.clear();

    // assemble RHS and diffusion matrices
    // Eigen is column major but cpp is row major, so...?
    for (std::size_t i = 1; i <= nx; ++i) {
        for (std::size_t j = 1; j <= ny; ++j) {
            for (std::size_t k = 1; k <= nz; ++k) {
                const auto r = pcfg.row(i, j, k);

                // populate RHS Vectors
                viscosity_rhs_u(r) = u1[i, j, k];
                viscosity_rhs_v(r) = v1[i, j, k];
                viscosity_rhs_w(r) = w1[i, j, k];

                // populate triplets for A
                triplets(viscosity_triplets_u, pcfg, state, nu, bc,
                         Dimension::X, dt, i, j, k);
                triplets(viscosity_triplets_v, pcfg, state, nu, bc,
                         Dimension::Y, dt, i, j, k);
                triplets(viscosity_triplets_w, pcfg, state, nu, bc,
                         Dimension::Z, dt, i, j, k);
            }
        }
    }

    // assemble sparse matrices
    viscosity_A_u.setFromTriplets(viscosity_triplets_u.begin(),
                                  viscosity_triplets_u.end());
    viscosity_A_v.setFromTriplets(viscosity_triplets_v.begin(),
                                  viscosity_triplets_v.end());
    viscosity_A_w.setFromTriplets(viscosity_triplets_w.begin(),
                                  viscosity_triplets_w.end());

    // factorize / prepare iterative solvers
    viscosity_solver_u.compute(viscosity_A_u);
    viscosity_solver_v.compute(viscosity_A_v);
    viscosity_solver_w.compute(viscosity_A_w);

    // solve diffusion systems
    viscosity_u2 = viscosity_solver_u.solve(viscosity_rhs_u);
    viscosity_v2 = viscosity_solver_v.solve(viscosity_rhs_v);
    viscosity_w2 = viscosity_solver_w.solve(viscosity_rhs_w);

    if (viscosity_solver_u.info() != Eigen::Success) {
        throw std::runtime_error("Viscous diffusion failed for U.");
    }
    if (viscosity_solver_v.info() != Eigen::Success) {
        throw std::runtime_error("Viscous diffusion failed for V.");
    }
    if (viscosity_solver_w.info() != Eigen::Success) {
        throw std::runtime_error("Viscous diffusion failed for W.");
    }

    // Ugly copying...
    for (std::size_t i = 1; i <= nx; ++i) {
        for (std::size_t j = 1; j <= ny; ++j) {
            for (std::size_t k = 1; k <= nz; ++k) {
                const auto r = pcfg.row(i, j, k);

                u1[i, j, k] = viscosity_u2(r);
                v1[i, j, k] = viscosity_v2(r);
                w1[i, j, k] = viscosity_w2(r);
            }
        }
    }
}

}  // namespace solver
