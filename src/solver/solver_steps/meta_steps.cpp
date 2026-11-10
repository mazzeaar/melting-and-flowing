#include "../solver.hpp"
#include "core/log.hpp"

namespace solver {

void Solver::prepare(Solver::Config cfg) {
    Log("[Solver] Preparing solver");
    pressure_cfg.from(cfg.grid.config());
    Log("[Solver] Finished preparing pressure config");

    {  // velocity buffer
        // has size: size(U) + size(V) + size(W)
        cfg.grid.resize_like<grid::Field::Velocities>(advect_velocity_buffer);
    }

    {  // viscosity buffers
        const auto n = static_cast<Eigen::Index>(pressure_cfg.n_unkowns());

        // RHS buffers
        viscosity_rhs_u.resize(n);
        viscosity_rhs_v.resize(n);
        viscosity_rhs_w.resize(n);

        // Solution buffers
        viscosity_u2.resize(n);
        viscosity_v2.resize(n);
        viscosity_w2.resize(n);

        // Triplet buffers
        viscosity_triplets_u.reserve(7 * n);
        viscosity_triplets_v.reserve(7 * n);
        viscosity_triplets_w.reserve(7 * n);

        // Sparse matrices
        viscosity_A_u.resize(n, n);
        viscosity_A_v.resize(n, n);
        viscosity_A_w.resize(n, n);
    }

    {  // temperature advection buffer
        const auto n = static_cast<Eigen::Index>(
            cfg.grid.config().grid_mapping.required_span_size());
        advect_temperature_buffer.resize(n);

        // Temperature diffusion buffers
        temperature_rhs.resize(n);
        temperature_T2.resize(n);

        temperature_triplets.reserve(7 * n);
        temperature_A.resize(n, n);
    }
}

// one full timestep
// solves for du/dt
auto Solver::step(Solver::Config cfg) -> Solver::Result {
    Log("timestepping start");
    using Clock = std::chrono::high_resolution_clock;

    for (std::size_t step = 0; step < cfg.max_steps; ++step) {
        auto begin = Clock::now();
        single_step(cfg);

        stats.duration.record(
            std::chrono::duration_cast<Duration>(Clock::now() - begin));
    }

    // do we need this?
    return Solver::Result{};
};

// one full timestep, including all substeps
// does not solve for du/dt !!!
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
}  // namespace solver
