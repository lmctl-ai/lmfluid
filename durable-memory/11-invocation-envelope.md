# Invocation envelope (0.2.2)

SDK origin/main 93bddf1; lmscript 0.2.460. Based on the lmparticle port.
InvocationResults validates metadata and loads all parts before any execution.
Select/run/complete each part, stop on first error, emit one terminal SDK frame.
Exit 0 success, 2 CLI/metadata, 3 load/run, 4 application/output contract failure.
Dry metadata and standalone --dry-run validate without running FEM/LBM or file
operations; contradictory real metadata is refused. Keep memory budgets and
release a serialization reserve before reporting allocation failures.

fluid.print and MFEM diagnostics go to stderr. Scripts use run_lmcode and r.ok;
CSV files carry numerical data. Existing numerical/error/warning/memory assertions
remain; private-input capture replaces exported bytecode tests. Added protocol
controls and offline example prefixes including CSV postprocessing. Renderers
require compatible local binaries, with no network override. Full images can be
regenerated only when lmplot, lmdraw and lmcad answer --lmcode-schema --json.
