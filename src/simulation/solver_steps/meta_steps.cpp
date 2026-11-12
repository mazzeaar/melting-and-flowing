#include "../solver.hpp"
#include "core/log.hpp"

void Solver::prepare(Solver::Config cfg) {
    Log("[Solver] Preparing solver");
    update_viscosity(cfg);
    enforce_velocity_boundaries(cfg.grid);

    const auto& gconf = cfg.grid.config();

    {  // velocity buffer
        // has size: size(U) + size(V) + size(W)
        advect_velocity_buffer.resize(gconf.size_u() + gconf.size_v() +
                                      gconf.size_w());
    }

    {  // Each staggered component has its own number of unknowns.
        const auto nu = gconf.size_u(), nv = gconf.size_v(),
                   nw = gconf.size_w();
        viscosity_rhs_u.resize(nu);
        viscosity_rhs_v.resize(nv);
        viscosity_rhs_w.resize(nw);
        viscosity_u2.resize(nu);
        viscosity_v2.resize(nv);
        viscosity_w2.resize(nw);
        viscosity_triplets_u.reserve(7 * nu);
        viscosity_triplets_v.reserve(7 * nv);
        viscosity_triplets_w.reserve(7 * nw);
        viscosity_A_u.resize(nu, nu);
        viscosity_A_v.resize(nv, nv);
        viscosity_A_w.resize(nw, nw);
    }

    {  // temperature advection buffer
        const auto n = gconf.size();
        advect_temperature_buffer.resize(n);

        // Temperature diffusion buffers
        temperature_rhs.resize(n);
        temperature_T2.resize(n);

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
    const auto cfl = cfg.settings.cfl;
    const auto target_dt = cfg.settings.dt;

    if (!std::isfinite(target_dt) || target_dt <= 0) {
        throw core::Error("Invalid timestep");
    }

    auto remaining = target_dt;
    std::size_t substeps = 0;

    // Apply material slider changes before the first viscous solve.
    update_viscosity(cfg);

    while (remaining > 0.0) {
        if (substeps >= cfg.settings.max_substeps) {
            throw core::Error(
                "CFL substep limit reached; increase Max substeps or reduce "
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

        // Avoid floating-point stagnation.
        if (remaining - dt == remaining) {
            throw core::Error("CFL timestep too small");
        }

        cfg.settings.dt = dt;
        std::println("[Solver] substep {}", substeps);
        single_step(cfg);

        remaining -= dt;
        ++substeps;
    }
}

// one full timestep, including all substeps
// does not solve for du/dt !!!
// void Solver::single_step(Solver::Config cfg) {
//     enforce_velocity_boundaries(cfg.grid);
//     advect_velocity(cfg);
//     enforce_velocity_boundaries(cfg.grid);
//
//     // Step 2: Viscous diffusion
//     // u2 - dt * nu * ∇²u2 = u1
//     diffuse_velocity(cfg);
//     enforce_velocity_boundaries(cfg.grid);
//
//     // Step 3: Pressure Poisson equation
//     // ∇²p = (rho / dt) * ∇·u2
//     solve_pressure(cfg);
//
//     // Step 4: Pressure projection
//     // u_new = u2 - (dt / rho) * ∇p
//     project_velocity(cfg);
//
//     // Step 5: Temperature advection
//     // T1 = T_old - dt * (u_new · ∇)T_old
//     advect_temperature(cfg);
//
//     // Step 6: Temperature diffusion
//     // T_new - dt * kappa * ∇²T_new = T1
//     diffuse_temperature(cfg);
//
//     // Step 7: Update viscosity from temperature
//     // eta_new = eta(T_new)
//     update_viscosity(cfg);
//
//     advect_particles(cfg);
// }

void Solver::single_step(Solver::Config cfg) {
    using Clock = std::chrono::steady_clock;

    std::string timings;

    const auto timed = [&](std::string_view name, auto&& fn) {
        const auto start = Clock::now();
        fn();
        const auto end = Clock::now();

        const double ms =
            std::chrono::duration<double, std::milli>(end - start).count();

        timings += std::format("  {:<24} {:>10.3f} ms\n", name, ms);
    };

    const auto total_start = Clock::now();

    timed("Enforce boundaries 1",
          [&] { enforce_velocity_boundaries(cfg.grid); });

    timed("Advect velocity", [&] { advect_velocity(cfg); });

    timed("Enforce boundaries 2",
          [&] { enforce_velocity_boundaries(cfg.grid); });

    timed("Diffuse velocity", [&] { diffuse_velocity(cfg); });

    timed("Enforce boundaries 3",
          [&] { enforce_velocity_boundaries(cfg.grid); });

    timed("Solve pressure", [&] { solve_pressure(cfg); });

    timed("Project velocity", [&] { project_velocity(cfg); });

    timed("Advect temperature", [&] { advect_temperature(cfg); });

    timed("Diffuse temperature", [&] { diffuse_temperature(cfg); });

    timed("Update viscosity", [&] { update_viscosity(cfg); });

    timed("Advect particles", [&] { advect_particles(cfg); });

    const double total_ms =
        std::chrono::duration<double, std::milli>(Clock::now() - total_start)
            .count();

    timings += std::format("  {:<24} {:>10.3f} ms", "TOTAL", total_ms);

    std::println("[Solver] Step timings:\n{}", timings);
}
