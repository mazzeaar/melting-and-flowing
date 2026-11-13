#include <utility>

#include "../solver.hpp"

// Step 2: Viscous Diffusion
// u2 - dt * nu * ∇²u2   = u1 (scalar nu)
// (I - dt * nu * ∇²) u2 = u1
// or
// u2 - ∇ · (nu * ∇u2)   = u1 (vector nu)
// (I - ∇ · (nu * ∇)) u2 = u1

enum class WallCondition : std::uint8_t { FreeSlip, NoSlip };

enum class Dimension : std::uint8_t { X, Y, Z };

void triplets(auto& triplets, const auto& gc, const WallCondition bc,
              const Dimension dim, const float dt, const float nu,
              const Eigen::Vector3i& ex, const auto& row, const int i,
              const int j, const int k) {
    const auto h = gc.cell_size();

    const float dx = h.x();
    const float dy = h.y();
    const float dz = h.z();

    const float idx2 = 1.0f / (dx * dx);
    const float idy2 = 1.0f / (dy * dy);
    const float idz2 = 1.0f / (dz * dz);

    const auto r = row(i, j, k);

    const int axis = std::to_underlying(dim);

    const auto active = [&](int x, int y, int z) {
        const Eigen::Vector3i q{x, y, z};

        return (q.array() >= 1).all() && (q.array() < ex.array() - 1).all() &&
               q[axis] > 1 && q[axis] < ex[axis] - 2;
    };

    // Identity rows for ghosts and prescribed wall faces.
    if (!active(i, j, k)) {
        triplets.emplace_back(r, r, 1.0f);
        return;
    }

    // Diagonal Element
    triplets.emplace_back(r, r, 1 + dt * nu * 2 * (idx2 + idy2 + idz2));

    // X-Axis Neighbors
    if (active(i - 1, j, k)) {
        triplets.emplace_back(r, row(i - 1, j, k), -dt * nu * idx2);
    } else {
    }
    if (active(i + 1, j, k)) {
        triplets.emplace_back(r, row(i + 1, j, k), -dt * nu * idx2);
    } else {
    }  // Boundary Conditions...

    // Y-Axis Neighbors
    if (active(i, j - 1, k)) {
        triplets.emplace_back(r, row(i, j - 1, k), -dt * nu * idy2);
    } else {
    }
    if (active(i, j + 1, k)) {
        triplets.emplace_back(r, row(i, j + 1, k), -dt * nu * idy2);
    } else {
    }

    // Z-Axis Neighbors
    if (active(i, j, k - 1)) {
        triplets.emplace_back(r, row(i, j, k - 1), -dt * nu * idz2);
    } else {
    }
    if (active(i, j, k + 1)) {
        triplets.emplace_back(r, row(i, j, k + 1), -dt * nu * idz2);
    } else {
    }
}

void Solver::diffuse_velocity(Solver::Config cfg) {
    auto& u1 = cfg.grid.u;  // X
    auto& v1 = cfg.grid.v;  // Y
    auto& w1 = cfg.grid.w;  // Z
    auto& state = cfg.grid.cell_state;

    const auto dt = cfg.settings.dt;
    const auto gc = cfg.grid.config();

    const auto nu = 0.3f;
    // Eigen::VectorXd nu = Eigen::VectorXd::Constant(n, 0.3);
    WallCondition bc = WallCondition::NoSlip;

    // Construct the RHS Vectors (just for the interior cells)
    auto& rhs_u = buff_mat_u1;
    auto& rhs_v = buff_mat_v1;
    auto& rhs_w = buff_mat_w1;

    // Prepare Triplets for Matrix A
    trip_u.clear();
    trip_v.clear();
    trip_w.clear();

    trip_u.reserve(7 * gc.size_u());
    trip_v.reserve(7 * gc.size_v());
    trip_w.reserve(7 * gc.size_w());

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
                    const bool interior = (q.array() >= 1).all() &&
                                          (q.array() < ex.array() - 1).all() &&
                                          q[axis] > 1 && q[axis] < ex[axis] - 2;

                    // Populate RHS Vectors
                    rhs(r) = interior ? field[r] : 0.0f;

                    // Populate Triplets for A
                    triplets(entries, gc, bc, dim, dt, nu, ex, row, i, j, k);
                }
            }
        }
    };

    prepare(
        u1, rhs_u, trip_u, gc.extents_u(),
        [&](int i, int j, int k) { return gc.u(i, j, k); }, Dimension::X);

    prepare(
        v1, rhs_v, trip_v, gc.extents_v(),
        [&](int i, int j, int k) { return gc.v(i, j, k); }, Dimension::Y);

    prepare(
        w1, rhs_w, trip_w, gc.extents_w(),
        [&](int i, int j, int k) { return gc.w(i, j, k); }, Dimension::Z);

    // Construct Matrix A
    auto& A_u = viscosity_A_u;
    auto& A_v = viscosity_A_v;
    auto& A_w = viscosity_A_w;

    A_u.setFromTriplets(trip_u.begin(), trip_u.end());
    A_v.setFromTriplets(trip_v.begin(), trip_v.end());
    A_w.setFromTriplets(trip_w.begin(), trip_w.end());

    A_u.makeCompressed();
    A_v.makeCompressed();
    A_w.makeCompressed();

    // Solve with Jacobi preconditioner (default)
    auto& solver_u = viscosity_solver_u;
    auto& solver_v = viscosity_solver_v;
    auto& solver_w = viscosity_solver_w;

    solver_u.compute(A_u);
    solver_v.compute(A_v);
    solver_w.compute(A_w);

    auto& u2 = buff_mat_u2;
    auto& v2 = buff_mat_v2;
    auto& w2 = buff_mat_w2;

    u2 = solver_u.solve(rhs_u);
    v2 = solver_v.solve(rhs_v);
    w2 = solver_w.solve(rhs_w);

    if (solver_u.info() != Eigen::Success) {
        throw std::runtime_error(
            "Viscous Diffusion: Solver failed for u-component.");
    }

    if (solver_v.info() != Eigen::Success) {
        throw std::runtime_error(
            "Viscous Diffusion: Solver failed for v-component.");
    }

    if (solver_w.info() != Eigen::Success) {
        throw std::runtime_error(
            "Viscous Diffusion: Solver failed for w-component.");
    }

    u1 = u2;
    v1 = v2;
    w1 = w2;
}
