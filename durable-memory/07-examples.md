# Examples and independent oracles

examples/poisson.lms solves -Delta u=-4 on the unit square with quadratic Dirichlet data and exports a 41x41 sample grid. examples/channel.lms exports a low-Mach plane Poiseuille profile and checks its centerline magnitude. Every example asserts known values before CSV output.

Tests cover quadratic Poisson exactness, P1/P2 sine convergence, the Laplace eigenmode equation -Delta u=2*pi^2*u, heat decay exp(-2*pi^2*t), affine elasticity, 3D box sampling and LBM rest invariants. Native regressions add nonzero elasticity body forces and higher-order patches, mixed boundary conditions, time-dependent data and Ghia cavity comparison.

The channel walls are halfway outside the fluid nodes at y=0 and y=H. Therefore the exact profile is evaluated at y=j+0.5, not at j or at nodes on the wall. For H=12, g=1e-5, nu=.1 the continuum centerline maximum is .0018; the sampled cell center is slightly below the true midpoint. The exported velocity is .001780999998 near the center after 3000 steps.

No upstream solver code or test data is copied. Manufactured formulas and a small factual benchmark datum are independent scientific oracles, with primary citations in 04-numerics.md. DSL literals use decimal notation, while expression strings accept scientific notation and pi. All bytecode and scratch meshes use build/; examples write build/examples/*.csv.
