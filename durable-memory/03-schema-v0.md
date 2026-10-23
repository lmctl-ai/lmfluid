# Schema and units

Schema is generated from literal check registrations in cpp/main.cpp. Test.sh compares generated JSON, embedded bytes and registry/compiler signatures. Typed method arguments are checked at the native boundary. Host objects live in stable arenas until process exit; Fields retain mesh/FEC/FES/GridFunction lifetimes through shared data.

Context fluid creates rectangular/box meshes, imported MFEM meshes, Poisson/heat/elasticity problems and channel/cavity/periodic LBM simulations. Problems use builder methods with expressions in x,y,z,t. FEM orders 1–3 are supported. Attribute zero means every boundary. The default is zero Dirichlet on every boundary; choosing a particular attribute replaces that default, and unselected faces are natural homogeneous boundaries. Setting all-boundary data replaces previous declarations.

FEM fields expose point values (component zero in the 2D convenience overload), integrated scalar/vector L2 errors, DOF/component/time metadata and bounding-box sampling. Matrices are zero-based. Scalar CSV columns x,y,z,u; elasticity x,y,z,ux,uy[,uz]; LBM x,y,rho,ux,uy,speed. There is no header. Samples and CSV use deterministic iteration order and 17-digit precision. Imported nonrectangular meshes can have outside-domain bounding-box sample points, which produce an error.

LBM quantities use lattice units and half-cell centers. viscosity nu sets tau=.5+3nu, g is channel x acceleration, and cavity lid_speed is the moving wall velocity. run advances the current state; init resets lattice density to one and time to zero. Physical units require a caller-selected lattice mapping. No automatic unit conversion is done.

load_mesh accepts the conforming MFEM mesh v1.0/v1.2 text formats (what Mesh.save writes) and validates element/boundary vertex indices, geometry codes, attributes (1..100000), counts and the vertex space dimension before MFEM parses the file: MFEM trusts those fields and segfaults on a dangling vertex index. Other formats (Gmsh, VTK, NC, NURBS) are rejected. MFEM stdout diagnostics (mfem::out) are routed to stderr so stdout stays the JSON result.

Only standard version/notices/schema/signature/bytecode CLI queries exist. Dry run loads/verifies bytecode without running or writing files. Usage exits 2; load/run errors exit 3 with SDK diagnostics. Grid/sample/step/DOF limits and finite/positivity checks guard numerical misuse.

## v0.2.0 additions (2026-10-23)

- `Field.gradient(x,y,direction)` and `Field.gradient(x,y,z,component,direction)`: MFEM
  `GetGradient` / `GetVectorGradient` at the located integration point. Used for heat flux
  (`k dT/dx`), strains and stresses (`sigma_xx = (lambda+2mu) du_x/dx + lambda du_y/dy`).
- `LBM.obstacle(expr)`: cells whose centre gives expr > 0 become solid (union over calls).
  Solid cells skip collision; every fluid link that streams into a solid is reflected with
  stationary halfway bounceback (same rule as the domain walls), and `2 f_post e_q` is added to
  the momentum given to the solids. Solid cells report rho = 1, u = 0 and are excluded from
  `mass`, `validate_state` and `init`. An expression leaving no fluid cell, or non-finite, fails.
- `LBM.solid_force(c)`: that momentum exchange averaged over the steps of the most recent
  `run` call (one-step values carry the acoustic noise of a weakly compressible solver);
  `init` resets it. `LBM.is_solid`, `LBM.solid_cells`.
- `LBM.vorticity()`: matrix `x,y,omega,solid`; central differences, one-sided with the wall
  at half a cell next to walls/obstacles (lid velocity for the moving lid). Solid rows are 0.
- `LBM.tau()`, `LBM.mach()` (= max_speed * sqrt 3).
- `fluid.matrix(rows,cols)` (<= 10 million entries) and `Matrix.set(row,col,value)` (finite
  values only) so scripts can record time series and sweeps and save them as CSV; the DSL has
  no other way to emit numbers.
- Warnings on stderr (never stdout, never failing the run): tau < 0.51 at construction, and
  once per flow when the Mach number after a `run` exceeds 0.3.
- `step()` no longer recomputes macroscopic fields at its start; they are current after
  construction, `init`, `obstacle` and every step (results are bit-identical, ~30% faster).
- Boundary attributes (verified by tests/boundary_attributes.lms): rect 1 bottom, 2 right,
  3 top, 4 left; box 1 z=0, 2 y=0, 3 x=lx, 4 y=ly, 5 x=0, 6 z=lz.
