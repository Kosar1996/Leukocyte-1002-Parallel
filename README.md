# Leukocyte-Endothelium Transendothelial Migration Solver

Axisymmetric MAC finite-volume Stokes fluid solver, coupled to solid mechanics for a leukocyte and an endothelium, via a monolithic
Newton solve. Models the thin lubricating fluid film between a leukocyte
squeezing through an endothelial gap.

## How to run

**Entry point: `run_input.m`**, the main reference input deck (referenced
by more of the codebase than any other script). For the specific hybrid
gap-1D / exterior-2D pressure formulation used in the most recent
diagnostic runs, use **`run_input_full2D_pressure2.m`** instead — same
solver, different fluid discretization choice. Both call into
`softlube_run_case_global_coupled.m`, the actual solver.

```matlab
run('run_input_full2D_pressure2.m')
```

### Required input data

Both entry points need two prestress files, expected in the same
directory as the script (`softlubeDir = ''` in both):

- `solid_endo_P300_wide.mat` — endothelium prestress. **Present in this
  repo.** Can be regenerated with `generate_prestress_endothelium_wide.m`
  if needed.
- `solid_leu_P600.mat` — leukocyte prestress. **Present in this repo.**
  **No generator script for this one exists anywhere in the codebase
  history** (checked, including deleted files) — if this file is ever
  lost or needs to be regenerated with different parameters, that
  process needs to be recovered or rebuilt from scratch; it isn't
  currently reproducible from code alone.

## Codebase structure

After a cleanup pass (removed ~161 one-off test/debug/investigation
scripts that had zero references from anywhere else in the codebase —
see git log for the itemized reasoning behind each removal), what
remains falls into these groups:

- **Entry points**: `run_input.m`, `run_input_full2D_pressure2.m`, and a
  couple of named variants (`run_pure2dmac_dtlarge_7steps.m`,
  `run_pure2dmac_dtsmall_14steps.m`) for specific dt/step-count studies.
- **Core solver**: `softlube_run_case_global_coupled.m` (the main
  driver), `solve_monolithic_two_solids_fsolve_timestep.m` /
  `solve_monolithic_analytical_timestep.m` (the two monolithic-solve
  paths — the second is used automatically when the leukocyte is set
  rigid), `monolithic_two_solids_residual_jacobian_scaled.m`,
  `solve_finite_def_solid.m` (FEM solid solve), and the fluid solvers
  (`solve_fluid_2D_bodyfitted_MAC.m`,
  `solve_fluid_hybrid_gap1d_exterior2d*.m`).
- **Traction correction**: `apply_bodyfitted_MAC_traction_correction.m`
  is the one actually used in production (tested, converges).
  `apply_bodyfitted_MAC_traction_correction_feedback.m` is a documented
  non-convergent alternative, kept for reference only — see the comment
  at its call site in `softlube_run_case_global_coupled.m` before ever
  enabling it.
- **Plotting/analysis utilities**: `plot_select_native2d_*.m` (pressure,
  stress, velocity heatmaps — these are what generated the report
  figures), plus various `recover_*`/`check_*`/`compare_*` helper
  functions still in active use (i.e., referenced by something else).

## Known gotchas

**`cfg.fluid.*` vs. `cfg.parOverrides.*` duplication.** Several fluid
settings are set twice in the entry-point scripts — once under
`cfg.fluid.*`, once again under `cfg.parOverrides.*`. `parOverrides`
silently wins when both are present, so editing one without the other
will not do what you expect. If you change a fluid setting, change it
in both places, or the change may appear to have no effect.

**No parallelization.** Confirmed via search — no `parfor`, `parpool`,
`spmd`, or `gpuArray` anywhere in the codebase. A single run cannot be
sped up by requesting more CPUs/nodes; only running independent cases
concurrently as separate jobs makes use of extra compute resources.

**Leukocyte elasticity overrides.** If you need to override the
leukocyte's Young's modulus (`cfg.solid.leukocyte.EL`) or Poisson ratio
(`nuL`) away from the prestress file's own baked-in values, verify the
override actually took effect by checking out.par.EL after the run, not just out.cfg.solid.leukocyte.EL.

## Data files

Only 5 `.mat` files remain, all genuinely needed:

- `solid_endo_P300_wide.mat`, `solid_leu_P600.mat` — required prestress
  input data (see "Required input data" above).
- `prestress_IC.mat` — required input, referenced directly by
  `ini_solid_endothelium.m` and `run_input.m`.
- `out_pure2dmac_dtlarge_7steps.mat`, `out_pure2dmac_dtsmall_14steps.mat`
  — output of the two still-live `run_pure2dmac_*` entry-point variants.

29 orphaned old output `.mat` files (leftovers from already-deleted
one-off scripts) have been removed.
