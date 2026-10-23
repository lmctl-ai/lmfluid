# lmfluid

Finite-element PDEs and lattice-Boltzmann flow for lmscript in one static executable. MFEM 4.8 supplies serial H1 finite elements for Poisson, heat and linear elasticity; an original MIT D2Q9 solver supplies periodic flow, body-force-driven channels and lid-driven cavities.

[![CI](https://github.com/lmctl-ai/lmfluid/actions/workflows/ci.yml/badge.svg)](https://github.com/lmctl-ai/lmfluid/actions/workflows/ci.yml)

```lmscript
p = _dsl_ cli="lmfluid" install="github:lmctl-ai/lmfluid@latest";
let mesh = fluid.rect(12, 12, 1.0, 1.0);
let field = fluid.poisson(mesh, 2).source("-4").dirichlet(0, "x*x+y*y").solve();
fluid.assert_close(field.value(0.5, 0.5), 0.5, 0.000001);
fluid.assert_close(field.l2_error("x*x+y*y"), 0.0, 0.000001);
field.sample(41, 41).save_csv("poisson.csv");
end _dsl_
savedsl(p, "poisson.dslb");
system ${dsl_cli(p)} --dsl poisson.dslb;
```

Run `lmscript poisson.lms`. The compiler selects lmfluid from PATH or installs the matching release. The host CLI accepts only `--dsl FILE...`, `--dry-run --dsl FILE...`, `--version`, `--notices`, `--dsl-schema [--json]` and `--dsl-signatures`. Dry run validates bytecode without solving or writing files. Usage errors exit 2; load/run errors exit 3.

| API | Purpose |
| --- | --- |
| `fluid.rect(nx,ny,lx,ly)`, `fluid.box(nx,ny,nz,lx,ly,lz)` | Cartesian quad/hex mesh |
| `fluid.load_mesh(path)`, `mesh.save(path)` | MFEM mesh file input/output |
| `mesh.dimension()`, `mesh.elements()` | mesh metadata |
| `fluid.poisson(mesh,order)`, `fluid.heat(...)`, `fluid.elasticity(...)` | finite-element problem (orders 1–3) |
| `problem.diffusion("expr")`, `.source("expr")`, `.initial("expr")` | scalar coefficients, source and heat initial condition |
| `problem.dirichlet(attribute,"expr")` | scalar essential boundary condition; attribute 0 selects all boundaries |
| `problem.lame(lambda,mu)`, `.force(fx,fy,fz)` | isotropic elastic material and expression body force |
| `problem.displacement(attribute,ux,uy,uz)` | vector essential boundary condition |
| `problem.tol(relative_tolerance,max_iterations)` | iterative linear solver controls |
| `problem.solve()` | steady Poisson or elasticity field |
| `problem.evolve(end_time,dt)` | heat field at end time, starting at t=0 |
| `field.value(x,y)` / `.value(x,y,z,component)` | value at a mesh point, zero-based component |
| `field.l2_error("exact")`, `.l2_error_vector(ex,ey,ez)` | integrated L2 error |
| `field.sample(nx,ny)` / `.sample(nx,ny,nz)` | sample bounding-box grid into a matrix |
| `field.dofs()`, `.components()`, `.time()` | solution metadata |
| `fluid.channel(nx,ny,nu,g)`, `.cavity(nx,ny,nu,lid_speed)`, `.periodic(nx,ny,nu)` | D2Q9 flow simulation |
| `flow.init(ux_expr,uy_expr)`, `.run(steps)` | initialize velocities and advance integer lattice steps |
| `flow.velocity(ix,iy,component)`, `.density(ix,iy)` | cell values, zero-based indices |
| `flow.mass()`, `.max_speed()`, `.kinetic_energy()`, `.steps()` | flow diagnostics |
| `flow.sample()` | lattice data matrix |
| `matrix.rows()`, `.cols()`, `.get(row,col)`, `.save_csv(path)` | matrix access and deterministic CSV |
| `fluid.eval("expr")`, `.assert_close(a,b,absolute_tolerance)`, `.assert_true(bool)` | mathematical helper and assertions |

Poisson solves `-div(k grad u)=f`. Heat solves `u_t-div(k grad u)=f` using backward Euler; diffusion and source expressions may depend on time. Elasticity solves `-div(lambda*tr(epsilon) I + 2*mu*epsilon)=force`. Material parameters must satisfy the supported coercivity conditions. The default essential boundary data is zero on all faces. Selecting a specific boundary replaces that default, with unselected faces having homogeneous natural conditions. Attribute numbers come from the mesh; `0` avoids depending on a particular face numbering.

FEM expressions use `x,y,z,t`, arithmetic, powers and common functions such as `sin`, `cos`, `exp`, `log` and `sqrt`; use `pi` inside expressions. For 2D, z=0; the third force/displacement expression is unused. A solved field owns its mesh and approximation space. Altering the originating problem does not alter the stored solution. Point queries outside the mesh fail. Sampling an imported nonrectangular mesh fails when a bounding-box sample is outside the domain. CSV is the exchange format for [lmplot](https://github.com/lmctl-ai/lmplot).

FEM CSV columns are `x,y,z,u` or `x,y,z,ux,uy[,uz]`. LBM CSV columns are `x,y,rho,ux,uy,speed`, with cell centers `(ix+0.5,iy+0.5)`. Files contain numeric rows without headers, in deterministic order and 17-digit precision.

## Lattice units and limits

D2Q9 uses BGK relaxation `tau=0.5+3*nu`, Guo forcing and halfway bounceback. Channel flow is periodic in x with walls at y=0 and y=ny. Cavity flow has stationary side/bottom walls and a moving top wall. All-periodic flow supports shear-wave experiments. Each call to `run` continues the existing state; `init` resets the step counter and unit density. Velocity expressions use lattice x,y,t coordinates.

LBM is a weakly compressible, low-Mach approximation in lattice units, not a general high-Reynolds-number CFD package. Keep velocities well below the lattice sound speed `1/sqrt(3)` and resolve gradients. The host validates finite density/velocity and stops clearly when stability checks fail. No turbulence models, free surfaces, complex LBM obstacles or open inlet/outlet boundaries are included.

## Build and verification

```sh
scripts/sync-deps.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
./test.sh
```

Requires CMake 3.24+, C++17, git, Python 3 and lmscript for testing (`npm i -g @lmctl-ai/lmctl@0.2.369`). MFEM is built statically with exceptions and no MPI, METIS, OpenMP, GPU, LAPACK, zlib or other external numerical libraries. No runtime data files are needed. Linux binaries are fully static; macOS links only system libraries; Windows statically links MinGW runtimes.

The suite checks manufactured solutions, P1/P2 spatial convergence, heat eigenmode decay, an elasticity affine patch, 3D Poisson and mesh round trips, periodic mass/momentum, shear-wave decay, the Poiseuille profile and Ghia's Re=100 cavity centerline datum. Test programs and both example CSV exports run with assertions, deterministic reruns, schema/signature parity, pin/notices/static-link checks, dry run, and error cases. See [durable memory](durable-memory/index.md) for measurements and exact scope.

## License

Our host, wrappers, examples and D2Q9 implementation are [MIT](LICENSE). MFEM is BSD-3-Clause and is fetched at the commit in `deps.lock` into ignored `third_party/`; its bundled permissive notices are reproduced in [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES) and embedded in `--notices`. OpenFOAM and GPL/LGPL scientific libraries are excluded. The optional reduced Eigen backend is unnecessary when the serial MFEM build passes the release matrix.
