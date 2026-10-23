# Open items (ranked)

1. A release tag v0.2.0 is needed before the examples that use the v0.2.0 API (obstacle,
   gradient, matrix, mach) run from `lmctl autorun` without a checkout; CMakeLists already
   says 0.2.0.
2. No inlet/outlet LBM boundaries: channel flow is periodic and body-force driven, so the
   cylinder wake is a periodic row of cylinders. Zou-He velocity/pressure boundaries would
   allow the Schaefer-Turek benchmark proper (C_D, C_L, St).
3. Obstacles are staircase cells (first order on curved walls); interpolated bounceback
   (Bouzidi) would make them second order.
4. No roller / single-component Dirichlet condition in elasticity, so symmetry planes cannot
   be used (plate_hole meshes the whole plate).
5. No eigenvalue solver (REVIEW.md item 3), no Navier-Stokes FEM, no Neumann/Robin data
   builders (a convective heat-transfer fin needs a Robin condition).
6. Mesh formats other than MFEM v1.0/1.2 are rejected (Gmsh would help real geometries).
7. Byte-identical CSVs are checked on one platform at a time.
8. Higher-order time integration, MPI/GPU, adaptive refinement, explicit object release for
   very long scripts.

Closed in v0.2.0: LBM runs above Mach 0.3 and tau close to 0.5 now warn on stderr (REVIEW
item 1); examples (with their drawing hosts) run in test.sh (REVIEW item 4).
