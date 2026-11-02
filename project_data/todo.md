1. Particle initialization
2. Cell classification (solver step)
    - In prepare step
    - After each particle advection step
3. Adjust pressure solve
    - Pressure solve must be adjusted to account for air etc -> pressure 0
    - Correct air cells (right now we are assuming its all liquid)
4. Update pressure projection
5. Particle advection
