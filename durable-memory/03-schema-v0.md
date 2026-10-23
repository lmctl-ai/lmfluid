# Schema and units

Schema is generated from literal check registrations in cpp/main.cpp. Test.sh compares generated JSON, embedded bytes and registry/compiler signatures. Typed method arguments are checked at the native boundary. Host objects live in stable arenas until process exit; Fields retain mesh/FEC/FES/GridFunction lifetimes through shared data.

Context fluid creates rectangular/box meshes, imported MFEM meshes, Poisson/heat/elasticity problems and channel/cavity/periodic LBM simulations. Problems use builder methods with expressions in x,y,z,t. FEM orders 1–3 are supported. Attribute zero means every boundary. The default is zero Dirichlet on every boundary; choosing a particular attribute replaces that default, and unselected faces are natural homogeneous boundaries. Setting all-boundary data replaces previous declarations.

FEM fields expose point values (component zero in the 2D convenience overload), integrated scalar/vector L2 errors, DOF/component/time metadata and bounding-box sampling. Matrices are zero-based. Scalar CSV columns x,y,z,u; elasticity x,y,z,ux,uy[,uz]; LBM x,y,rho,ux,uy,speed. There is no header. Samples and CSV use deterministic iteration order and 17-digit precision. Imported nonrectangular meshes can have outside-domain bounding-box sample points, which produce an error.

LBM quantities use lattice units and half-cell centers. viscosity nu sets tau=.5+3nu, g is channel x acceleration, and cavity lid_speed is the moving wall velocity. run advances the current state; init resets lattice density to one and time to zero. Physical units require a caller-selected lattice mapping. No automatic unit conversion is done.

Only standard version/notices/schema/signature/bytecode CLI queries exist. Dry run loads/verifies bytecode without running or writing files. Usage exits 2; load/run errors exit 3 with SDK diagnostics. Grid/sample/step/DOF limits and finite/positivity checks guard numerical misuse.
