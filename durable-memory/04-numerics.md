# Numerical methods and scientific references

MFEM H1 spaces on quad/hex meshes implement Poisson diffusion and isotropic linear elasticity. Boundary elimination uses MFEM FormLinearSystem/RecoverFEMSolution. Systems are solved by serial preconditioned conjugate gradients with an explicit convergence check. Heat uses backward Euler, assembling mass/dt + diffusion and the old-field mass RHS at every step. Coefficients, sources and Dirichlet data may depend on time. Diffusion is checked positive at evaluated quadrature points; expressions must remain finite.

P1/P2 expected L2 convergence is order 2/3 for smooth solutions. The unit-square sine mode satisfies -Delta u=2*pi^2*u, providing the Laplace eigenvalue/manufactured-source oracle. Exact polynomial and affine vector patch tests verify nonzero essential boundaries and elasticity. Heat sine mode uses exp(-2*pi^2*t); spatial and temporal errors are separated with appropriate tolerances.

The original LBM code uses D2Q9 BGK, tau=.5+3nu, Guo forcing F=rho*g, periodic streaming and halfway bounceback, including moving-wall momentum correction. Initialization removes the half-force impulse so requested velocities are the reported initial velocities. Cavity corners follow the documented implementation's halfway-wall rule. No solver code is copied from OpenFOAM or GPL implementations.

Quantitative checks use u(y)=g*y*(H-y)/(2nu) for steady plane Poiseuille flow, exp(-nu*k^2*t) for a periodic shear wave, and the Re=100 centerline value u/U=-0.20581 from Ghia, Ghia and Shin (1982). The 32x32 low-Mach cavity check has a finite-grid tolerance; it is a benchmark spot check, not a high-Re CFD accuracy claim.

Primary references: [Ghia et al., JCP 48 (1982)](https://doi.org/10.1016/0021-9991(82)90058-4); [Guo, Zheng and Shi, PRE 65 (2002)](https://doi.org/10.1103/PhysRevE.65.046308); [MFEM v4.8 examples](https://github.com/mfem/mfem/tree/v4.8/examples).
