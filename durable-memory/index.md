# lmfluid durable memory

`lmfluid` is a static, dependency-free lmscript `_dsl_` host for 2D/3D computational fluid dynamics, Navier-Stokes, lattice Boltzmann (LBM), and particle hydrodynamics.
Repository: https://github.com/lmctl-ai/lmfluid.

## 60-second takeover

```sh
cmake -S . -B build                   # Configure CMake build
cmake --build build --parallel 4      # Build binary at build/lmfluid
./test.sh                             # Run full test suite (stability checks, fluid benchmarks)
```

## Document map

| File | Contents |
| --- | --- |
| [01-purpose-and-place.md](01-purpose-and-place.md) | Scope, fluid solver architecture, and relation to sibling hosts |
| [02-libraries-and-licenses.md](02-libraries-and-licenses.md) | Clean-room implementations, third-party libraries, and licenses |
| [03-schema-v0.md](03-schema-v0.md) | Fluid domain models, boundary conditions, and schema definitions |
| [04-numerics.md](04-numerics.md) | Discretization schemes, CFL stability conditions, pressure Poisson solvers |
| [05-build-run-test.md](05-build-run-test.md) | Build targets, determinism guarantees, and verification suite |
| [06-ci.md](06-ci.md) | GitHub Actions CI matrix (5 platforms) and build workflows |
| [07-examples.md](07-examples.md) | Cavity flow, Karman vortex shedding, and dam break benchmarks |
| [08-open-items.md](08-open-items.md) | Roadmap, turbulence models, and future extensions |
| [09-visual-gallery.md](09-visual-gallery.md) | Flow field visualization, streamlines, and gallery artifacts |
| [10-resource-limits.md](10-resource-limits.md) | Memory budget allocation (`LMFLUID_MAX_MEMORY_MB`) and grid scaling |
| [next-release-version.md](next-release-version.md) | Version synchronization protocol (VERSION, CMake, THIRD_PARTY_NOTICES) |

## SDK pin and sync

The vendored SDK header is pinned via `sdk.lock` pointing to `https://github.com/lmctl-ai/lmcode-cpp`:
```sh
scripts/sync-sdk.sh                   # Refetches pinned cpp/dslvm.hpp and verifies sha256
```

## Adding an example or function

### Adding an example
1. Create `examples/<name>.lms` using `_dsl_ cli="lmfluid" install="github:lmctl-ai/lmfluid@latest"`.
2. Run script: `build/lmfluid run examples/<name>.lms` or via `lmscript`.
3. Add the example to `./test.sh` test run list.
4. Save visual output to `docs/images/<name>.png` and update `README.md` demo table with `lmctl autorun` command.

### Adding a host function
1. Declare function/method in `schema/lmfluid.dsl.json`.
2. Register and implement the function in C++ source under `src/` using `dslvm::Registry`.
3. Rebuild and verify schema synchronization with `cmake --build build` and `./test.sh`.

## Releases

Releases are triggered automatically on push of a `v*` tag. The release workflow builds across Linux x64/arm64, macOS x64/arm64, and Windows x64, publishing platform binaries with `SHA256SUMS`, `LICENSE`, and `THIRD_PARTY_NOTICES`. See `next-release-version.md` before tagging.
