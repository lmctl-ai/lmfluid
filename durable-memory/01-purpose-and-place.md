# Purpose and place

lmfluid follows simulation/PROCEDURE.md and the lmfluid recommendation in FEASIBILITY.md. It is a scripting-only host with context fluid:Fluid, host mesh/problem/field/matrix objects, expression strings and CSV interchange. A single static executable embeds its schema, legal notices, serial MFEM and original D2Q9 code.

MFEM is the first choice because the serial library has no required third-party numerical dependency, supports CMake/MinGW, and provides mature H1 spaces, assembly and iterative solvers. The reduced own FEM/Eigen fallback is reserved for demonstrated five-platform static-build failure, not used merely because MFEM is larger. LBM provides channel/cavity flow while a coupled incompressible FEM Navier-Stokes solver remains beyond v0.1.
