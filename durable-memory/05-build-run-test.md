# Build and measured verification

Measured 2026-10-23 on x86_64 Linux, GCC 11 and CMake/Ninja. The clean MFEM-only Release static target compiled 221 objects in 317.04 seconds with four build jobs. Its static archive is about 33 MB; the fully static stripped host is about 13.4 MB. MFEM reports no external TPL includes or dependencies in this configuration. The reduced Eigen fallback has not been needed.

Build: scripts/sync-deps.sh, cmake -S . -B build -DCMAKE_BUILD_TYPE=Release, cmake --build build --parallel 4, then ./test.sh. C++17, git, Python 3 and the lmscript v5 compiler are required; CI pins @lmctl-ai/lmctl@0.2.369. Dependency and SDK pins are verified before building.

The integration suite has 18 grouped checks, two native suites (29 FEM checks and 34 LBM checks), DSL test and example programs, deterministic CSV files, and 25 runtime error snippets (including a dangling-vertex mesh and CG non-convergence, whose stdout must stay one JSON line). It verifies exact version, registry/generated/embedded/compiler schema parity, notices, platform linkage (dynamic glibc >= 2.35 NEEDED allowlist on Linux via `readelf -d`, `otool -L` on macOS, `objdump -p` on Windows), numeric assertions, byte-identical reruns, dry run without output writes and usage/load/run exit codes. Duplicate program basenames are rejected to prevent bytecode outputs masking another test.

Native FEM regressions include manufactured Poisson, spatial convergence, exact polynomial nonzero BC, mixed scalar/vector boundary attributes, time-varying diffusion/source/boundary heat solutions, tiny and zero time intervals, elasticity, mesh round trips, retained field ownership and invalid numerical inputs. LBM checks include independent-run/split-run determinism, mass/momentum, shear decay, full Poiseuille profile and Ghia Re=100 cavity data.

Measured P1 L2 errors on 4/8 grids: 0.0301896 and 0.00758781, rate 1.9923. P2 errors on 2/4 grids: 0.0144453 and 0.00193277, rate 2.90186. The cavity center ratio is -0.197112 compared with -0.20581 (steady after 10k steps). The native test compares all 30 tabulated Ghia Re=100 centerline points (15 u, 15 v) with linear interpolation between cell centers: max deviation 0.0093 (u) / 0.0074 (v) at 32x32 and 0.0056 / 0.0034 at 64x64 (center -0.2031), so the gap is second-order discretization error; tolerance 0.015. DSL convergence (8 -> 16, sine mode) asserts error ratios 3.994 / 7.973 / 15.96 for P1/P2/P3 within 5-6 percent of 4 / 8 / 16; heat at t=0.01 is 0.82212 vs exp(-2 pi^2 t) = 0.82087 (tolerance 0.003). Native Release suites run in roughly 1.5 seconds locally.

GCC and Clang compile our wrappers; static analysis and shellcheck/Python checks are run. ASan/UBSan instrument our FEM wrapper, matrix code, LBM and native tests while linking the already-built unmodified MFEM archive. MFEM itself is not sanitizer-instrumented. No assertion is disabled by NDEBUG. Cross-platform results use numerical tolerances; deterministic CSV comparisons are within one platform.

Successful main CI build-and-test phase times (seconds, including tests; restored build trees may be reused):

| Platform | Seconds |
| --- | --- |
| macos-x64 | 493 |
| linux-arm64 | 202 |
| linux-x64 | 309 |
| macos-arm64 | 213 |
| windows-x64 | 278 |
