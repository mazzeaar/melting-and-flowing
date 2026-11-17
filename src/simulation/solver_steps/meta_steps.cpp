#include <chrono>

#include "../solver.hpp"
#include "core/log.hpp"

void Solver::prepare(Solver::Config cfg) {
    Log("[Solver] Preparing solver");
    update_viscosity(cfg);
    enforce_velocity_boundaries(cfg.grid);

    const auto& gconf = cfg.grid.config();

    const auto n = gconf.size();
    const auto nu = gconf.size_u();
    const auto nv = gconf.size_v();
    const auto nw = gconf.size_w();

    {  // new buffers
        // eigen vector buffers
        buff_mat_u1.resize(nu);
        buff_mat_v1.resize(nv);
        buff_mat_w1.resize(nw);
        buff_mat_cell1.resize(n);
        buff_mat_cell2.resize(n);

        buff_mat_u2.resize(nu);
        buff_mat_v2.resize(nv);
        buff_mat_w2.resize(nw);

        // triplet buffers
        trip_u.resize(nu);
        trip_v.resize(nu);
        trip_w.resize(nu);
        trip_cell.resize(n);

        // cell state buffer
        buff_mat_cell_state.resize(n);
    }

    {  // each staggered component has its own number of unknowns
        // viscosity_rhs_u.resize(nu);
        // viscosity_rhs_v.resize(nv);
        // viscosity_rhs_w.resize(nw);
        // viscosity_u2.resize(nu);
        // viscosity_v2.resize(nv);
        // viscosity_w2.resize(nw);
        // viscosity_triplets_u.reserve(7 * nu);
        // viscosity_triplets_v.reserve(7 * nv);
        // viscosity_triplets_w.reserve(7 * nw);
        viscosity_A_u.resize(nu, nu);
        viscosity_A_v.resize(nv, nv);
        viscosity_A_w.resize(nw, nw);
    }

    {  // temperature advection buffer
        temperature_triplets.reserve(7 * n);
        temperature_A.resize(n, n);
    }
}

auto cfl_dt(const Grid& grid, float cfl) -> float {
    const auto h = grid.config().cell_size();

    const auto rate = grid.u.cwiseAbs().maxCoeff() / h.x() +
                      grid.v.cwiseAbs().maxCoeff() / h.y() +
                      grid.w.cwiseAbs().maxCoeff() / h.z();

    if (rate == 0.0) {
        return std::numeric_limits<double>::infinity();
    }

    return cfl / rate;
}

// one full timestep
// solves for du/dt
void Solver::step(Solver::Config cfg) {
    auto begin = std::chrono::high_resolution_clock::now();
    const auto cfl = cfg.settings.cfl;
    const auto target_dt = cfg.settings.dt;

    if (!std::isfinite(target_dt) || target_dt <= 0) {
        throw core::Error("Invalid timestep");
    }

    auto remaining = target_dt;
    std::size_t substeps = 0;

    // apply material slider changes before the first viscous solve
    update_viscosity(cfg);

    while (remaining > 0.0) {
        if (substeps >= cfg.settings.max_substeps) {
            throw core::Error(
                "CFL substep limit reached:( increase Max substeps or reduce "
                "dt");
        }

        const auto dt = std::min(remaining, cfl_dt(cfg.grid, cfl));
        if (!std::isfinite(dt) &&
            dt != std::numeric_limits<double>::infinity()) {
            throw core::Error("Invalid CFL timestep");
        }

        if (dt <= 0.0) {
            throw core::Error("Non-positive CFL timestep");
        }

        // avoid floating-point stagnation (idk what that means, chaddy told me
        // to do it)
        if (remaining - dt == remaining) {
            throw core::Error("CFL timestep too small");
        }

        cfg.settings.dt = dt;
        single_step(cfg);

        remaining -= dt;
        ++substeps;
    }

    auto end = std::chrono::high_resolution_clock::now();
    stats.step_duration.record(
        std::chrono::duration_cast<Duration>(end - begin));
    stats.cfl_substeps.record(substeps);
}

void Solver::single_step(Solver::Config cfg) {
    advect_velocity(cfg);

    // Step 2: Viscous diffusion
    // u2 - dt * nu * ∇²u2 = u1
    diffuse_velocity(cfg);

    // Step 3: Pressure Poisson equation
    // ∇²p = (rho / dt) * ∇·u2
    solve_pressure(cfg);

    // Step 4: Pressure projection
    // u_new = u2 - (dt / rho) * ∇p
    project_velocity(cfg);

    // Step 5: Temperature advection
    // T1 = T_old - dt * (u_new · ∇)T_old
    advect_temperature(cfg);

    // Step 6: Temperature diffusion
    // T_new - dt * kappa * ∇²T_new = T1
    diffuse_temperature(cfg);

    // Step 7: Update viscosity from temperature
    // eta_new = eta(T_new)
    update_viscosity(cfg);

    advect_particles(cfg);
}
