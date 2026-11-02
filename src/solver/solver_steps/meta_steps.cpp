#include "../solver.hpp"
#include "core/log.hpp"

namespace solver {

void Solver::prepare(Solver::Config cfg) {
    Log("[Solver] Preparing solver");
    pressure_cfg.prepare(cfg.grid.config());
    Log("[Solver] Finished preparing pressure config");
    // stuff that we only have to do once per grid
    // - memory allocations
    // - set up stencils and then:
    // https://libeigen.gitlab.io/eigen/docs-5.0.1/group__TopicSparseSystems.html#:~:text=The%20Compute%20Step,with%20other%20matrices.
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
}
}  // namespace solver
