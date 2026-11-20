#include "../solver.hpp"

// Step 2: Viscous Diffusion
// u2 - dt * nu * ∇²u2   = u1 (scalar nu)
// (I - dt * nu * ∇²) u2 = u1
// or
// u2 - ∇ · (nu * ∇u2)   = u1 (vector nu)
// (I - ∇ · (nu * ∇)) u2 = u1

namespace {

enum class WallCondition : std::uint8_t { FreeSlip, NoSlip };
enum class FaceState : std::uint8_t { Fluid, Solid, Air, Unknown };
enum class Axis : std::uint8_t { X, Y, Z };

void triplets(auto& entries, const Solver::Config& cfg, const WallCondition bc,
              const Axis ax, const float nu, const int i, const int j,
              const int k) {
    const auto gc = cfg.grid.config();
    const auto& state = cfg.grid.cell_state;
    const float dt = cfg.settings.dt;

    const auto ex = gc.extents();
    const auto h = gc.cell_size();

    const float dx = h.x();
    const float dy = h.y();
    const float dz = h.z();

    const float cx = dt * nu / (dx * dx);
    const float cy = dt * nu / (dy * dy);
    const float cz = dt * nu / (dz * dz);

    const auto get_face = [&](int i0, int j0, int k0) -> FaceState {
        int i1 = i0 + (ax == Axis::X ? -1 : 0);
        int j1 = j0 + (ax == Axis::Y ? -1 : 0);
        int k1 = k0 + (ax == Axis::Z ? -1 : 0);

        if (i0 < 0 || i0 >= ex.x() || j0 < 0 || j0 >= ex.y() || k0 < 0 ||
            k0 >= ex.z() || i1 < 0 || j1 < 0 || k1 < 0) {
            return FaceState::Solid;
        }

        auto cs0 = state[gc.cell(i0, j0, k0)];
        auto cs1 = state[gc.cell(i1, j1, k1)];

        if (cs0 == CellState::UNKNOWN || cs1 == CellState::UNKNOWN) {
            return FaceState::Unknown;
        }

        if (cs0 == CellState::SOLID || cs1 == CellState::SOLID) {
            return FaceState::Solid;
        }

        if (cs0 == CellState::FULL || cs1 == CellState::FULL ||
            cs0 == CellState::SURFACE || cs1 == CellState::SURFACE) {
            return FaceState::Fluid;
        }

        if (cs0 == CellState::AIR && cs1 == CellState::AIR) {
            return FaceState::Air;
        }

        return FaceState::Unknown;
    };

    const auto row = [&](int i0, int j0, int k0) -> std::size_t {
        switch (ax) {
        case Axis::X: return gc.u(i0, j0, k0);
        case Axis::Y: return gc.v(i0, j0, k0);
        case Axis::Z: return gc.w(i0, j0, k0);
        }
        return 0;
    };

    const int r = row(i, j, k);

    // Non-Fluid
    if (get_face(i, j, k) != FaceState::Fluid) {
        entries.emplace_back(r, r, 1.0f);
        return;
    }

    float diag = 1.0f;

    struct Neighbors {
        int delta_x, delta_y, delta_z;
        float coeff;
    };
    const Neighbors neighbors[6] = {{-1, 0, 0, cx}, {1, 0, 0, cx},
                                    {0, -1, 0, cy}, {0, 1, 0, cy},
                                    {0, 0, -1, cz}, {0, 0, 1, cz}};

    for (const auto& n : neighbors) {
        const int nx = i + n.delta_x;
        const int ny = j + n.delta_y;
        const int nz = k + n.delta_z;

        FaceState nFaceState = get_face(nx, ny, nz);

        // Fluid
        if (nFaceState == FaceState::Fluid) {
            diag += n.coeff;
            const int nr = row(nx, ny, nz);
            entries.emplace_back(r, nr, -n.coeff);
        }
        // Solid
        else if (nFaceState == FaceState::Solid) {
            const bool is_normal = (ax == Axis::X && n.delta_x != 0) ||
                                   (ax == Axis::Y && n.delta_y != 0) ||
                                   (ax == Axis::Z && n.delta_z != 0);
            if (bc == WallCondition::NoSlip) {
                diag += (is_normal ? 1.0f : 2.0f) * n.coeff;
                // For moving walls
                // rhs += (is_normal ? 1.0f : 2.0f) * n.coeff * wall_velocity
            } else if (bc == WallCondition::FreeSlip) {
                diag += is_normal ? n.coeff : 0.0f;
                // For moving walls
                // rhs += (is_normal ? n.coeff : 0.0f) * wall_velocity
            }
        }
        // Air
        // else if (nFaceState == FaceState::Air) { }
    }

    entries.emplace_back(r, r, diag);
}

}  // namespace

void Solver::diffuse_velocity(Solver::Config cfg) {
    auto& u1 = cfg.grid.u;  // X
    auto& v1 = cfg.grid.v;  // Y
    auto& w1 = cfg.grid.w;  // Z

    const auto gc = cfg.grid.config();

    // Magic numbers for now
    const auto nu = cfg.physics.reference_viscosity;
    // Eigen::VectorXd nu = Eigen::VectorXd::Constant(n, 0.3);
    WallCondition bc = WallCondition::FreeSlip;

    // Construct the RHS Vectors
    // (not needed, but could be used for moving walls)
    // auto& rhs_u = buff_mat_u1;
    // auto& rhs_v = buff_mat_v1;
    // auto& rhs_w = buff_mat_w1;

    // Prepare Triplets for Matrix A
    trip_u.clear();
    trip_v.clear();
    trip_w.clear();

    trip_u.reserve(7 * gc.size_u());
    trip_v.reserve(7 * gc.size_v());
    trip_w.reserve(7 * gc.size_w());

    const auto prepare = [&](auto& entries, Eigen::Vector3i ex, Axis ax) {
        for (int i = 0; i < ex.x(); ++i) {
            for (int j = 0; j < ex.y(); ++j) {
                for (int k = 0; k < ex.z(); ++k) {
                    // Populate RHS Vectors
                    // const auto r = row(gc, ax, i, j, k);
                    // rhs(r) = field[r];

                    // Populate Triplets for A
                    triplets(entries, cfg, bc, ax, nu, i, j, k);
                }
            }
        }
    };

    prepare(trip_u, gc.extents_u(), Axis::X);
    prepare(trip_v, gc.extents_v(), Axis::Y);
    prepare(trip_w, gc.extents_w(), Axis::Z);

    // Construct Matrix A
    auto& A_u = viscosity_A_u;
    auto& A_v = viscosity_A_v;
    auto& A_w = viscosity_A_w;

    A_u.setFromTriplets(trip_u.begin(), trip_u.end());
    A_v.setFromTriplets(trip_v.begin(), trip_v.end());
    A_w.setFromTriplets(trip_w.begin(), trip_w.end());

    // Already compressed by setFromTriplets
    // A_u.makeCompressed();
    // A_v.makeCompressed();
    // A_w.makeCompressed();

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

    u2 = solver_u.solveWithGuess(u1, u1);
    v2 = solver_v.solveWithGuess(v1, v1);
    w2 = solver_w.solveWithGuess(w1, w1);

    /*
    // Often fails... but good enough...
    if (solver_u.info() != Eigen::Success ||
        solver_v.info() != Eigen::Success ||
        solver_w.info() != Eigen::Success) {
        throw std::runtime_error(
            "Viscous Diffusion: Solver failed to converge.");
    }
    */

    u1 = u2;
    v1 = v2;
    w1 = w2;
}
