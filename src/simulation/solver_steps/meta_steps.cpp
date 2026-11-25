#include "../solver.hpp"
#include "core/benchmark.hpp"
#include "core/exception.hpp"
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
    benchmark.tick("total");

    const auto cfl = cfg.settings.cfl;
    const auto target_dt = cfg.settings.dt;

    if (!std::isfinite(target_dt) || target_dt <= 0) {
        throw core::Error("Invalid timestep");
    }

    auto remaining = target_dt;
    std::size_t substeps = 0;

    // apply material slider changes before the first viscous solve
    benchmark.tick("initial_viscosity");
    update_viscosity(cfg);
    benchmark.tock("initial_viscosity");

    while (remaining > 0.0) {
        if (substeps >= cfg.settings.max_substeps) {
            throw core::Error(
                "CFL substep limit reached:( increase Max substeps or reduce "
                "dt");
        }

        benchmark.tick_cfl("cfl_dt");
        const auto dt = std::min(remaining, cfl_dt(cfg.grid, cfl));
        benchmark.tock_cfl("cfl_dt");

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

        // chaddy said to add this, i dont see a difference with/without it
        enforce_velocity_boundaries(cfg.grid);
    }

    benchmark.end_step(substeps);
    benchmark.tock("total");
}

void Solver::single_step(Solver::Config cfg) {
    // Step 1: Advection and external forces
    benchmark.tick_cfl("advect_velocity");
    advect_velocity(cfg);
    enforce_velocity_boundaries(cfg.grid);
    benchmark.tock_cfl("advect_velocity");

    // Step 2: Viscous diffusion
    // u2 - dt * nu * ∇²u2 = u1
    benchmark.tick_cfl("diffuse_velocity");
    diffuse_velocity(cfg);
    enforce_velocity_boundaries(cfg.grid);
    benchmark.tock_cfl("diffuse_velocity");

    // Step 3: Pressure Poisson equation
    // ∇²p = (rho / dt) * ∇·u2
    benchmark.tick_cfl("solve_pressure");
    solve_pressure(cfg);
    benchmark.tock_cfl("solve_pressure");

    // Step 4: Pressure projection
    // u_new = u2 - (dt / rho) * ∇p
    benchmark.tick_cfl("project_velocity");
    project_velocity(cfg);
    enforce_velocity_boundaries(cfg.grid);
    benchmark.tock_cfl("project_velocity");

    // Step 5: Temperature advection
    // T1 = T_old - dt * (u_new · ∇)T_old
    benchmark.tick_cfl("advect_temperature");
    advect_temperature(cfg);
    benchmark.tock_cfl("advect_temperature");

    // Step 6: Temperature diffusion
    // T_new - dt * kappa * ∇²T_new = T1
    benchmark.tick_cfl("diffuse_temperature");
    diffuse_temperature(cfg);
    benchmark.tock_cfl("diffuse_temperature");

    // Step 7: Update viscosity from temperature
    // eta_new = eta(T_new)
    benchmark.tick_cfl("update_viscosity");
    update_viscosity(cfg);
    benchmark.tock_cfl("update_viscosity");

    // Step 8: Update particle positions
    benchmark.tick_cfl("advect_particles");
    advect_particles(cfg);
    benchmark.tock_cfl("advect_particles");
}
