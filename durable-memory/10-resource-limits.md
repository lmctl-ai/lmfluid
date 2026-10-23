# Resource limits

Why: host objects (meshes, problems, fields, LBM grids, matrices) live in arenas until exit, so a
DSL loop that keeps creating them grows without bound. On 2026-10-23 a runaway lm* host (lmnum)
reached 16 GB RSS and filled the 24 GB WSL machine; every host now has a total budget.

- **Total memory budget**, default 4096 MiB, set with `LMFLUID_MAX_MEMORY_MB` (positive whole MiB;
  anything else exits 2 with the usage text). Shown in the usage message.
- `cpp/memory_budget.{hpp,cpp}` (identical copy in every lm* host, linked into the `lmfluid`
  executable only, not into `fluid_core` or the native tests) replaces the global
  `operator new`/`delete` and counts block sizes from the C allocator.
- Everything lmfluid keeps goes through `operator new`, so nothing needs `charge()`:
  - MFEM's `Memory<T>` uses `new T[]` for the default host memory type `MemoryType::HOST`
    (`general/mem_manager.hpp`, `NewHOST`); the `std::malloc`/`posix_memalign` paths in
    `mem_manager.cpp` are only for the aligned/managed/device memory types, which this serial
    CPU build never selects. Meshes, FE spaces, grid functions, sparse matrices and CG vectors
    are therefore all counted, including temporaries inside one assembly or solve.
  - LBM grids (`std::vector`: `f`, `next` = 2*9*nx*ny doubles plus rho/ux/uy, 168 MB at 1000x1000),
    sampled matrices and the VM's own values are `std::vector`/`std::string`.
- Exceeding it throws `membudget::Exceeded` (a `std::bad_alloc`, fixed-buffer message). Inside MFEM
  it is thrown by `operator new` itself, so it never reaches `mfem_error`/`abort`: MFEM has no
  `catch` blocks or allocating `noexcept` functions on these paths and unwinds normally (coefficient
  callbacks already rely on this). `translate()` turns it into a `dslvm::Error` with the same
  message; main also catches `std::exception` around load/run (budget inside the VM). The first
  failure prints `lmfluid: memory budget exceeded: ... (set LMFLUID_MAX_MEMORY_MB ...)` on stderr,
  the JSON line carries the message, exit 3. 64 MiB of headroom is granted once for the error path.
- Per-object caps stay as before (mesh <= 100000 elements, FE space <= 1e6 dofs, LBM <= 1e6 cells,
  sample <= 5e6 points, 100000 objects per arena); the budget bounds their sum.
- Not covered: stack and static data, and the 64 MiB reserve after the first failure.
- No `setrlimit`: not portable to macOS static builds, and it would crash instead of reporting.
- Measured (x86_64 Linux, WSL): a loop of `fluid.periodic(1000, 1000, 0.1)` stops at 64 MiB
  immediately and at the default 4096 MiB after 14 s with 4079 MiB peak RSS (4074 MiB tracked at
  the failure). A single `poisson(rect(300,300), 3)` solve stops inside MFEM at 64 MiB in 0.3 s
  (73 MB RSS). `test.sh` runs both at 64 MiB everywhere and the default-budget loop on Linux under
  `ulimit -v 6291456`.
