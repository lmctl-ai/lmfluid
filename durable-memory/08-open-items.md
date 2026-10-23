# Beyond v0.1

- Coupled FEM incompressible Navier-Stokes/Stokes pressure blocks, advection and nonlinear constitutive models.
- Unstructured LBM geometry, obstacle masks, open inlet/outlet conditions, turbulence and multiphase models.
- Adaptive mesh refinement, mixed spaces, Neumann/Robin data builders, arbitrary initial vector fields and eigenvalue spectra.
- Higher-order time integration, GPU/MPI/OpenMP and external sparse solvers.
- Streaming large samples and explicit host object release for very long scripts.

The LBM backend is low-Mach and weakly compressible. FEM conditioning and corner/boundary regularity affect convergence. Imported meshes are accepted within MFEM's serial support and sampling requires points to be inside the domain. These limits are documented rather than silently approximated as unsupported physics.
