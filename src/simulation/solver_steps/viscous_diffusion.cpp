#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "../solver.hpp"

// Step 2: Viscous Diffusion
// u2 - dt * nu * ∇²u2   = u1 (scalar nu)
// (I - dt * nu * ∇²) u2 = u1
// or
// u2 - dt * ∇ · (nu * ∇u2)   = u1 (spatially varying nu)
// (I - dt * ∇ · (nu * ∇)) u2 = u1
// Component-wise approximation, not the full coupled viscous stress operator.

enum class WallCondition : std::uint8_t { FreeSlip, NoSlip };

enum class Dimension : std::uint8_t { X, Y, Z };

void triplets(auto& triplets, const auto& gc, const WallCondition bc,
              const Dimension dim, const float dt, const auto& nu,
              const auto& face_state, const auto& row, const int i, const int j,
              const int k) {
    const auto h = gc.cell_size();

    const float dx = h.x();
    const float dy = h.y();
    const float dz = h.z();

    const float idx2 = 1.0f / (dx * dx);
    const float idy2 = 1.0f / (dy * dy);
    const float idz2 = 1.0f / (dz * dz);

    const auto r = row(i, j, k);

    const int axis = static_cast<int>(dim);
    const Eigen::Vector3i q{i, j, k};

    // Identity rows for ghosts and prescribed wall faces (also air-only faces).
    if (face_state(q, axis) != CellState::FULL) {
        triplets.emplace_back(r, r, 1.0f);
        return;
    }

    // Diagonal Element
    float diag = 1.0f;
    const float nu_center = nu(q, axis);

    const auto neighbour = [&](int ni, int nj, int nk, int direction,
                               float inv_h2) {
        const Eigen::Vector3i n{ni, nj, nk};
        const auto state = face_state(n, axis);

        if (state == CellState::FULL) {
            // Shared coefficients keep the matrix symmetric for CG.
            const float a =
                dt * (0.5f * nu_center + 0.5f * nu(n, axis)) * inv_h2;
            diag += a;
            triplets.emplace_back(r, row(ni, nj, nk), -a);
        } else if (state == CellState::SOLID) {
            // Normal wall faces are zero. Tangential no-slip uses odd
            // reflection; free-slip and air neighbours have zero normal
            // diffusive flux.
            if (direction == axis) {
                diag += dt * nu_center * inv_h2;
            } else if (bc == WallCondition::NoSlip) {
                diag += 2.0f * dt * nu_center * inv_h2;
            }
        }
    };

    // X-Axis Neighbors
    neighbour(i - 1, j, k, 0, idx2);
    neighbour(i + 1, j, k, 0, idx2);  // Boundary Conditions...

    // Y-Axis Neighbors
    neighbour(i, j - 1, k, 1, idy2);
    neighbour(i, j + 1, k, 1, idy2);

    // Z-Axis Neighbors
    neighbour(i, j, k - 1, 2, idz2);
    neighbour(i, j, k + 1, 2, idz2);

    triplets.emplace_back(r, r, diag);
}

void Solver::diffuse_velocity(Solver::Config cfg) {
    const float_type dt = cfg.settings.dt;
    auto& u1 = cfg.grid.u;  // X
    auto& v1 = cfg.grid.v;  // Y
    auto& w1 = cfg.grid.w;  // Z
    const auto& state = cfg.grid.cell_state;

    const auto gc = cfg.grid.config();

    const auto ex = gc.extents();
    const float density = cfg.physics.density;
    if (!std::isfinite(density) || density <= 0.0f ||
        !cfg.grid.eta.allFinite() || (cfg.grid.eta.array() < 0.0f).any()) {
        throw std::runtime_error(
            "Viscous Diffusion: Invalid density or viscosity.");
    }

    // Match enforce_velocity_boundaries()'s free-slip exterior ghosts.
    WallCondition bc = WallCondition::FreeSlip;

    const auto cell_state = [&](const Eigen::Vector3i& q) {
        if ((q.array() < 0).any() || (q.array() >= ex.array()).any()) {
            return CellState::SOLID;
        }
        return static_cast<CellState>(state[gc.cell(q)]);
    };
    const auto fluid = [](CellState s) {
        return s == CellState::FULL || s == CellState::SURFACE;
    };
    const auto face_state = [&](const Eigen::Vector3i& q, int axis) {
        const Eigen::Vector3i upper = q - Eigen::Vector3i::Ones();
        Eigen::Vector3i lower = upper;
        --lower[axis];
        const auto a = cell_state(lower), b = cell_state(upper);
        if (a == CellState::SOLID || b == CellState::SOLID) {
            return CellState::SOLID;
        }
        return fluid(a) || fluid(b) ? CellState::FULL : CellState::AIR;
    };
    const auto nu = [&](const Eigen::Vector3i& q, int axis) {
        Eigen::Vector3i cell = q - Eigen::Vector3i::Ones();
        float eta = 0.0f;
        int count = 0;
        for (int side = 0; side < 2; ++side, --cell[axis]) {
            if (fluid(cell_state(cell))) {
                eta += cfg.grid.eta[gc.cell(cell)];
                ++count;
            }
        }
        return count ? (eta / static_cast<float>(count)) / density : 0.0f;
    };

    // Construct the RHS Vectors (just for the interior cells)
    auto& rhs_u = viscosity_rhs_u;
    auto& rhs_v = viscosity_rhs_v;
    auto& rhs_w = viscosity_rhs_w;

    // Prepare Triplets for Matrix A
    auto& triplets_u = viscosity_triplets_u;
    auto& triplets_v = viscosity_triplets_v;
    auto& triplets_w = viscosity_triplets_w;

    triplets_u.clear();
    triplets_v.clear();
    triplets_w.clear();

    triplets_u.reserve(7 * gc.size_u());
    triplets_v.reserve(7 * gc.size_v());
    triplets_w.reserve(7 * gc.size_w());

    const auto prepare = [&](const auto& field, auto& rhs, auto& entries,
                             Eigen::Vector3i ex, const auto& row,
                             Dimension dim) {
        const int axis = static_cast<int>(dim);

        // Eigen is column major but cpp is row major, so...?
        for (int i = 0; i < ex.x(); ++i) {
            for (int j = 0; j < ex.y(); ++j) {
                for (int k = 0; k < ex.z(); ++k) {
                    const auto r = row(i, j, k);

                    const Eigen::Vector3i q{i, j, k};
                    const auto face = face_state(q, axis);

                    // Populate RHS Vectors
                    rhs(r) = face == CellState::SOLID ? 0.0f : field[r];

                    // Populate Triplets for A
                    triplets(entries, gc, bc, dim, dt, nu, face_state, row, i,
                             j, k);
                }
            }
        }
    };

    prepare(
        u1, rhs_u, triplets_u, gc.extents_u(),
        [&](int i, int j, int k) { return gc.u(i, j, k); }, Dimension::X);

    prepare(
        v1, rhs_v, triplets_v, gc.extents_v(),
        [&](int i, int j, int k) { return gc.v(i, j, k); }, Dimension::Y);

    prepare(
        w1, rhs_w, triplets_w, gc.extents_w(),
        [&](int i, int j, int k) { return gc.w(i, j, k); }, Dimension::Z);

    // Construct Matrix A
    auto& A_u = viscosity_A_u;
    auto& A_v = viscosity_A_v;
    auto& A_w = viscosity_A_w;

    A_u.setFromTriplets(triplets_u.begin(), triplets_u.end());
    A_v.setFromTriplets(triplets_v.begin(), triplets_v.end());
    A_w.setFromTriplets(triplets_w.begin(), triplets_w.end());

    A_u.makeCompressed();
    A_v.makeCompressed();
    A_w.makeCompressed();

    // Solve with Jacobi preconditioner (default)
    auto& solver_u = viscosity_solver_u;
    auto& solver_v = viscosity_solver_v;
    auto& solver_w = viscosity_solver_w;

    for (auto* solver : {&solver_u, &solver_v, &solver_w}) {
        solver->setMaxIterations(cfg.settings.max_solver_iterations);
        solver->setTolerance(cfg.settings.solver_tolerance);
    }

    solver_u.compute(A_u);
    solver_v.compute(A_v);
    solver_w.compute(A_w);

    auto& u2 = viscosity_u2;
    auto& v2 = viscosity_v2;
    auto& w2 = viscosity_w2;

    u2 = solver_u.solve(rhs_u);
    v2 = solver_v.solve(rhs_v);
    w2 = solver_w.solve(rhs_w);

    if (solver_u.info() != Eigen::Success || !u2.allFinite()) {
        throw std::runtime_error(
            "Viscous Diffusion: Solver failed for u-component.");
    }

    if (solver_v.info() != Eigen::Success || !v2.allFinite()) {
        throw std::runtime_error(
            "Viscous Diffusion: Solver failed for v-component.");
    }

    if (solver_w.info() != Eigen::Success || !w2.allFinite()) {
        throw std::runtime_error(
            "Viscous Diffusion: Solver failed for w-component.");
    }

    u1 = u2;
    v1 = v2;
    w1 = w2;
}
