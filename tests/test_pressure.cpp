// checks for the pressure solve + projection (steps 3 and 4)
// run via `ctest --test-dir <build-dir>` or the test_pressure binary directly

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <print>
#include <vector>

#include "core/exception.hpp"
#include "grid/grid.hpp"
#include "solver/solver.hpp"

namespace solver {

// friend of Solver, reaches the private substeps
struct PressureTest {
    static void solve_and_project(Solver& solver, Solver::Config cfg) {
        solver.solve_pressure(cfg);
        solver.project_velocity(cfg);
    }
};

}  // namespace solver

namespace {

using grid::Field;
using PressureConfig = solver::Solver::PressureConfig;

// 4 interior cells per axis with dx = 0.5, so a missing dx or dt factor shows
// up in the result
auto make_config() -> grid::Config {
    return {std::array{6, 6, 6}, std::array{2.0, 2.0, 2.0},
            std::array{0.0, 0.0, 0.0}};
}

auto max_abs_divergence(GridView grid) -> double {
    auto pcfg = PressureConfig{};
    pcfg.prepare(grid.config());
    auto u = grid.view<Field::U>();
    auto v = grid.view<Field::V>();
    auto w = grid.view<Field::W>();

    auto max_div = 0.0;
    for (std::size_t i = 1; i <= pcfg.nx; ++i) {
        for (std::size_t j = 1; j <= pcfg.ny; ++j) {
            for (std::size_t k = 1; k <= pcfg.nz; ++k) {
                const auto div = (u[i, j, k] - u[i - 1, j, k]) / pcfg.dx +
                                 (v[i, j, k] - v[i, j - 1, k]) / pcfg.dy +
                                 (w[i, j, k] - w[i, j, k - 1]) / pcfg.dz;
                max_div = std::max(max_div, std::abs(div));
            }
        }
    }
    return max_div;
}

void run_pressure_steps(GridView grid) {
    auto cfg = solver::Solver::Config{.grid = grid, .dt = 0.1, .max_steps = 1};
    auto solver = solver::Solver{};
    solver.prepare(cfg);
    solver::PressureTest::solve_and_project(solver, cfg);
}

// arbitrary velocities on all interior faces -> divergence must vanish
auto test_removes_divergence() -> bool {
    auto storage = grid::Storage{make_config()};
    auto grid = GridView{storage};
    auto pcfg = PressureConfig{};
    pcfg.prepare(grid.config());
    auto u = grid.view<Field::U>();
    auto v = grid.view<Field::V>();
    auto w = grid.view<Field::W>();

    // wall faces (index 0 and n along their own axis) stay 0
    for (std::size_t i = 1; i <= pcfg.nx; ++i) {
        for (std::size_t j = 1; j <= pcfg.ny; ++j) {
            for (std::size_t k = 1; k <= pcfg.nz; ++k) {
                const auto x = static_cast<double>(i);
                const auto y = static_cast<double>(j);
                const auto z = static_cast<double>(k);
                if (i < pcfg.nx) {
                    u[i, j, k] = std::sin(x + (2.0 * y) + (3.0 * z));
                }
                if (j < pcfg.ny) {
                    v[i, j, k] = std::cos((3.0 * x) + y + (2.0 * z));
                }
                if (k < pcfg.nz) {
                    w[i, j, k] = std::sin((2.0 * x) - y + (1.5 * z));
                }
            }
        }
    }

    const auto before = max_abs_divergence(grid);
    run_pressure_steps(grid);
    const auto after = max_abs_divergence(grid);

    std::println("  max |div| before: {:.3e}, after: {:.3e}", before, after);
    return before > 1e-2 && after < 1e-8;
}

// a small vortex is divergence free -> projection must not change it
auto test_keeps_divergence_free_field() -> bool {
    auto storage = grid::Storage{make_config()};
    auto grid = GridView{storage};
    auto pcfg = PressureConfig{};
    pcfg.prepare(grid.config());
    auto u = grid.view<Field::U>();
    auto v = grid.view<Field::V>();

    // circulation through the four cells around the edge between (1,1) and
    // (2,2) in every xy-slice
    for (std::size_t k = 1; k <= pcfg.nz; ++k) {
        u[1, 1, k] = 1.0;
        v[2, 1, k] = 1.0;
        u[1, 2, k] = -1.0;
        v[1, 1, k] = -1.0;
    }

    const auto before = max_abs_divergence(grid);
    const auto u_before =
        std::vector<double>(u.flat_view().begin(), u.flat_view().end());
    const auto v_before =
        std::vector<double>(v.flat_view().begin(), v.flat_view().end());

    run_pressure_steps(grid);

    auto max_change = 0.0;
    for (std::size_t g = 0; g < u.size(); ++g) {
        max_change = std::max(max_change, std::abs(u[g] - u_before[g]));
        max_change = std::max(max_change, std::abs(v[g] - v_before[g]));
    }

    std::println("  max |div| before: {:.3e}, max velocity change: {:.3e}",
                 before, max_change);
    return before < 1e-12 && max_change < 1e-10;
}

// dz = 0.75 while dx = dy = 0.5 -> prepare must reject the grid
auto test_rejects_non_uniform_cells() -> bool {
    const auto cfg =
        grid::Config{std::array{6, 6, 6}, std::array{2.0, 2.0, 3.0},
                     std::array{0.0, 0.0, 0.0}};
    auto pcfg = PressureConfig{};
    try {
        pcfg.prepare(cfg);
    } catch (const core::Error& e) {
        std::println("  threw: {}", e.what());
        return true;
    }
    return false;
}

}  // namespace

int main() {
    struct Case {
        const char* name;
        bool (*run)();
    };

    const auto cases = std::array{
        Case{"removes divergence", test_removes_divergence},
        Case{"keeps divergence-free field", test_keeps_divergence_free_field},
        Case{"rejects non-uniform cells", test_rejects_non_uniform_cells},
    };

    auto failed = 0;
    for (const auto& c : cases) {
        std::println("[ RUN  ] {}", c.name);
        const auto ok = c.run();
        std::println("[ {} ] {}", ok ? " OK " : "FAIL", c.name);
        failed += ok ? 0 : 1;
    }
    return failed == 0 ? 0 : 1;
}
