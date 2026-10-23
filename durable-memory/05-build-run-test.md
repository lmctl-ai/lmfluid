# Build and measured verification

Measured 2026-10-23 on x86_64 Linux, GCC 11 and CMake/Ninja. The clean MFEM-only Release static target compiled 221 objects in 317.04 seconds with four build jobs. Its static archive is about 33 MB; the fully static stripped host is about 13.4 MB. MFEM reports no external TPL includes or dependencies in this configuration. The reduced Eigen fallback has not been needed.

Build: scripts/sync-deps.sh, cmake -S . -B build -DCMAKE_BUILD_TYPE=Release, cmake --build build --parallel 4, then ./test.sh. C++17, git, Python 3 and the lmscript v5 compiler are required; CI pins @lmctl-ai/lmctl@0.2.369. Dependency and SDK pins are verified before building.

The integration suite has 17 grouped checks, two native suites (29 FEM checks and 34 LBM checks), eight DSL programs (six tests and two examples), two deterministic CSV files, and fifteen runtime error snippets. It verifies exact version, registry/generated/embedded/compiler schema parity, notices, static linking, numeric assertions, byte-identical reruns, dry run without output writes and usage/load/run exit codes. Duplicate program basenames are rejected to prevent bytecode outputs masking another test.

Native FEM regressions include manufactured Poisson, spatial convergence, exact polynomial nonzero BC, mixed scalar/vector boundary attributes, time-varying diffusion/source/boundary heat solutions, tiny and zero time intervals, elasticity, mesh round trips, retained field ownership and invalid numerical inputs. LBM checks include independent-run/split-run determinism, mass/momentum, shear decay, full Poiseuille profile and Ghia Re=100 cavity data.

Measured P1 L2 errors on 4/8 grids: 0.0301896 and 0.00758781, rate 1.9923. P2 errors on 2/4 grids: 0.0144453 and 0.00193277, rate 2.90186. The cavity center ratio is -0.197112 compared with -0.20581; the coarse-grid test allows 0.04 absolute normalized error. Native Release suites run in roughly 1.5 seconds locally.

GCC and Clang compile our wrappers; static analysis and shellcheck/Python checks are run. ASan/UBSan instrument our FEM wrapper, matrix code, LBM and native tests while linking the already-built unmodified MFEM archive. MFEM itself is not sanitizer-instrumented. No assertion is disabled by NDEBUG. Cross-platform results use numerical tolerances; deterministic CSV comparisons are within one platform.

Successful main CI build-and-test phase times (seconds, including tests; restored build trees may be reused):

| Platform | Seconds |
| --- | --- |
| macos-x64 | 493 |
| linux-arm64 | 202 |
| linux-x64 | 309 |
| macos-arm64 | 213 |
| windows-x64 | 278 |
