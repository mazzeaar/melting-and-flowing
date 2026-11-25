#include "../boundary_conditions.hpp"
#include "../solver.hpp"

// Step 2: Viscous Diffusion
// u2 - dt * nu * ∇²u2   = u1 (scalar nu)
// (I - dt * nu * ∇²) u2 = u1
// or
// u2 - dt * ∇ · (nu * ∇u2)   = u1 (non-scalar nu)
// (I - dt * ∇ · (nu * ∇)) u2 = u1
// (accidental simplification by author)

namespace {

using float_type = Solver::float_type;
using boundary::cell_state;
using boundary::face_state;
using boundary::is_liquid;

enum class WallCondition : std::uint8_t { FreeSlip, NoSlip };
enum class Axis : std::uint8_t { X, Y, Z };

void triplets(auto& entries, const Solver::Config& cfg, const WallCondition bc, const Axis ax, 
              const int i, const int j, const int k) {
    const auto gc = cfg.grid.config();
    const auto& nu = cfg.grid.eta;

    const auto h = gc.cell_size();
    const float_type dt = cfg.settings.dt;
    const float_type dt_dx2 = dt / (h.x() * h.x());
    const float_type dt_dy2 = dt / (h.y() * h.y());
    const float_type dt_dz2 = dt / (h.z() * h.z());

    const auto row = [&](int i0, int j0, int k0) -> std::size_t {
        switch (ax) {
        case Axis::X: return gc.u(i0, j0, k0);
        case Axis::Y: return gc.v(i0, j0, k0);
        case Axis::Z: return gc.w(i0, j0, k0);
        }
        return 0;
    };

    const auto get_face = [&](int i0, int j0, int k0) -> CellState {
        return face_state(cfg.grid, std::to_underlying(ax), i0, j0, k0);
    };

    const auto get_nu = [&](int i0, int j0, int k0, bool isNormal) -> float_type {
        const Eigen::Vector3i cell{std::min(i, i0) - 1, 
                                   std::min(j, j0) - 1,
                                   std::min(k, k0) - 1};

        if (isNormal) {
            const CellState state = cell_state(cfg.grid, cell.x(), cell.y(), cell.z());
            return is_liquid(state) ? nu[gc.cell(cell)] : 0; 
        }

        Eigen::Vector3i ax_off = Eigen::Vector3i::Zero();
        ax_off[std::to_underlying(ax)] = 1;
        const Eigen::Vector3i tan_off{i0 != i ? 1 : 0,
                                      j0 != j ? 1 : 0,
                                      k0 != k ? 1 : 0};

        const std::array<Eigen::Vector3i, 4> cells{
            cell,
            cell - ax_off,
            cell           + tan_off,
            cell - ax_off  + tan_off,
        };

        double prod  = 1.0;
        int count = 0;
        for (const auto& c : cells) {
            const CellState state = cell_state(cfg.grid, c.x(), c.y(), c.z());
            if (is_liquid(state)) {
                prod *= nu[gc.cell(c)];
                ++count;
            }
        }
        switch (count) {
            case 1: return prod;
            case 2: return std::sqrt(prod);
            case 3: return std::cbrt(prod);
            case 4: return std::sqrt(std::sqrt(prod));
            default: return 0;
        }
    };

    const int r = row(i, j, k);

    // Non-Fluid
    if (get_face(i, j, k) != CellState::FULL) {
        entries.emplace_back(r, r, 1);
        return;
    }

    float_type diag = 1;

    struct Neighbors {
        int delta_x, delta_y, delta_z;
        float_type dt_dh2;
    };
    const Neighbors neighbors[6] = {{-1, 0, 0, dt_dx2}, {1, 0, 0, dt_dx2},
                                    {0, -1, 0, dt_dy2}, {0, 1, 0, dt_dy2},
                                    {0, 0, -1, dt_dz2}, {0, 0, 1, dt_dz2}};

    for (const auto& n : neighbors) {
        const int nx = i + n.delta_x;
        const int ny = j + n.delta_y;
        const int nz = k + n.delta_z;

        CellState nFaceState = get_face(nx, ny, nz);

        // Air
        if (nFaceState == CellState::AIR) { continue; }

        const bool isNormal = (ax == Axis::X && n.delta_x != 0) ||
                              (ax == Axis::Y && n.delta_y != 0) ||
                              (ax == Axis::Z && n.delta_z != 0);

        const float_type coeff = get_nu(nx, ny, nz, isNormal) * n.dt_dh2;
        if (coeff <= 0) { continue; } // skip zero coefficients

        // Fluid
        if (nFaceState == CellState::FULL) {
            diag += coeff;
            const int nr = row(nx, ny, nz);
            entries.emplace_back(r, nr, -coeff);
        }
        // Solid
        else if (nFaceState == CellState::SOLID) {
            if (bc == WallCondition::NoSlip) {
                diag += (isNormal ? 1 : 2) * coeff;
            } else if (bc == WallCondition::FreeSlip) {
                diag += isNormal ? coeff : 0;
            }
        }
    }

    entries.emplace_back(r, r, diag);
}

}  // namespace

void Solver::diffuse_velocity(Solver::Config cfg) {
    const auto gc = cfg.grid.config();
    WallCondition bc = WallCondition::FreeSlip; // Fixed bc, I guess?

    // Prepare Triplets for Matrix A
    trip_u.clear();
    trip_v.clear();
    trip_w.clear();

    trip_u.reserve(7 * gc.size_u());
    trip_v.reserve(7 * gc.size_v());
    trip_w.reserve(7 * gc.size_w());

    const auto prepare = [&](auto& entries, Eigen::Vector3i ex, Axis ax) {
        for (int k = 0; k < ex.z(); ++k) {
            for (int j = 0; j < ex.y(); ++j) {
                for (int i = 0; i < ex.x(); ++i) {

                    // Populate Triplets for A
                    triplets(entries, cfg, bc, ax, i, j, k);
                }
            }
        }
    };

    // Prepare Triplets for Matrix A
    // Note: Most entires are identity!
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

    // Solve with Jacobi preconditioner (default)
    auto& solver_u = viscosity_solver_u;
    auto& solver_v = viscosity_solver_v;
    auto& solver_w = viscosity_solver_w;

    solver_u.compute(A_u);
    solver_v.compute(A_v);
    solver_w.compute(A_w);

    auto& u1 = cfg.grid.u;
    auto& v1 = cfg.grid.v;
    auto& w1 = cfg.grid.w;

    auto& u2 = buff_mat_u2;
    auto& v2 = buff_mat_v2;
    auto& w2 = buff_mat_w2;

    u2 = solver_u.solveWithGuess(u1, u1);
    v2 = solver_v.solveWithGuess(v1, v1);
    w2 = solver_w.solveWithGuess(w1, w1);

    if (solver_u.info() != Eigen::Success ||
        solver_v.info() != Eigen::Success ||
        solver_w.info() != Eigen::Success) {
        throw std::runtime_error(
            "Viscous Diffusion: Solver failed to converge.");
    }

    u1 = u2;
    v1 = v2;
    w1 = w2;
}
