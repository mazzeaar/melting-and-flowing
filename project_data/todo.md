1. Particle initialization (Aaron)
2. Cell classification (solver step) (Aaron)
    - In prepare step
    - After each particle advection step
3. Adjust pressure solve (Martin)
    - Pressure solve must be adjusted to account for air etc -> pressure 0
    - Correct air cells (right now we are assuming its all liquid)
4. Update pressure projection (Martin)
5. Particle advection (Aaron)
