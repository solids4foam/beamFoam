---
name: beamFoam
description: >
  Use this skill for any task involving the beamFoam repository: writing or
  reviewing C++ changes to the beamFoam solver, beam models, cross-section
  models, beamMomentumContribution classes, beam patch fields, function
  objects or numerics, maintaining createBeamMesh or setInitialPositionBeam,
  adjusting OpenFOAM build files, or changing tutorial beam cases. Trigger
  whenever the user mentions beamFoam, createBeamMesh, setInitialPositionBeam,
  beamProperties, coupledTotalLagNewtonRaphsonBeam, crossSectionModel,
  beamMomentumContribution, Simo-Reissner beams, wmake or a beam tutorial.
---

# beamFoam repository guide

This guide is for agents that write or review changes in beamFoam. The rules
in sections 1-5 apply to both. When reviewing, use them as review criteria;
section 6 applies only when making changes.

beamFoam implements geometrically exact Simo-Reissner beams in OpenFOAM with
a cell-centred finite-volume formulation. It is used on its own and as a
library by coupled solvers (moorFV, solids4foam), so library behaviour must
not assume that `beamFoam` is the caller.

## 1) Compatibility, numerics and known pitfalls

- Supported releases are the OpenFOAM.com versions accepted by the top-level
  `Allwmake`. CI builds OpenFOAM v2312 (`.github/workflows/build.yml`). There
  is no foam-extend or OpenFOAM.org support; do not add variant guards for
  them.
- Eigen comes from `ThirdParty/eigen3`. Do not add other dependencies unless
  required and documented.
- Keep `-D__restrict__=` in `src/wireBunchingModels/Make/options`. OpenFOAM
  field loops mark pointers `__restrict__`, but expressions such as
  `tmp(A) & B` reuse `A`'s storage for the result; GCC 11+ at `-O3` then
  miscompiles tensor products and the tutorials diverge on v2312. In new
  code, still prefer named concrete fields over unnamed `tmp` operands of the
  result type.
- Rotations: use the helpers in `numerics/pseudoVector` (`rotationMatrix`,
  `pseudoVector`, `rotationTangentOperator`) and `numerics/spinTensor`
  rather than inline formulas. Check edge cases (zero rotation, half-turn,
  rotations mostly about one axis) and guard divisions, which raise SIGFPE
  when floating-point trapping (`FOAM_SIGFPE`) is on.
- `evolve()` may be called several times per time step by a coupled solver
  (e.g. once per PIMPLE outer corrector). Anything that must happen once per
  time step, such as the Newmark predictor, is guarded with
  `solverCounters().nEvolveCalls(timeIndex)`. Do not tie per-step work to the
  first Newton iteration of an `evolve()` call.
- A case can hold several beam regions. Find a beam model with
  `lookupBeamModel` (by mesh region), not by the `"beamProperties"` object
  name. `beamModel` reuses a region mesh that another model has already
  registered; do not construct a second mesh for the same region.
- Parallel: reduce processor-local contributions (`gSum`, `gMax`, `reduce`)
  when computing global quantities. All ranks must call collectives in the
  same order and number; never reduce inside per-cell loops. See
  `numerics/beamParallelChecks` for the checks used.
- A new dictionary entry needs a default that leaves existing cases
  unchanged, unless the change is intended and documented. New state fields
  must be written and read back exactly on restart.
- Preserve dimensions, the meaning of `W` (mean-line displacement), `Theta`
  (rotation vector), `pointW`, and the total-Lagrangian reference
  configuration. A change to the residual must keep the Newton Jacobian
  consistent, or state that the Jacobian is an intentional approximation.

## 2) Architecture

- Solver: `applications/solvers/beamFoam/beamFoam.C` is thin. It builds a
  `beamModel` from `constant/beamProperties` and loops over
  `beam().evolve()`, `beam().updateTotalFields()` and `beam().writeFields()`.
  Time-step control is in `readBeamFoamControls.H` and `setDeltaT.H`. Keep
  beam mechanics in the library.
- Library: `src/wireBunchingModels`, one build list in `Make/files`.
  - `beamModels/beamModel`: base class, selection, `lookupBeamModel`,
    cross-section models (`circle`, `rectangle`, `annularSector`, ...).
  - `beamModels/coupledTotalLagNewtonRaphsonBeam`: the beam model in use,
    split into `...Evolve.C`, `...Matrix.C` and `...BoundaryConditions.C`.
  - `beamModels/fvPatchFields`: beam boundary conditions
    (`fixedBeamDisplacement`, `fixedBeamRotation`, `forceBeamDisplacementNR`,
    `followerForceBeamDisplacementNR`, `momentBeamRotationNR`,
    `linearSpringForceBeamDisplacementNR`). They are part of the user
    interface; renaming one breaks existing cases and must be documented.
  - `beamMomentumContribution`: runtime-selectable external loads
    (`morisonDrag`, `morisonAddedMass`, `groundContact`, `almDrag`).
  - `functionObjects`: `beamConvergenceData`, `beamDisplacements`,
    `beamEnergyData`, `beamForcesMoments`. Each takes a `region`.
  - `numerics`: rotations, block Eigen solvers, Hermite splines, helpers.
- Utilities: `createBeamMesh` is the canonical mesh generator (no
  `blockMesh`); `setInitialPositionBeam` translates and rotates a beam
  cellZone into its reference configuration.
- Runtime-selectable classes need `TypeName("...")` in the header,
  `defineTypeNameAndDebug` and `addToRunTimeSelectionTable(...)` in the
  source, and the source in `Make/files`. The dictionary `type` must match
  `TypeName` exactly.

## 3) C++ style

- Added and changed C++ lines follow the OpenFOAM Coding Style Guide
  (<https://openfoam.org/dev/coding-style-guide>). Much existing code does
  not fully conform; do not reformat lines the change does not touch. Key
  rules:
  - at most 80 characters per line; indent 4 spaces, never tabs; no trailing
    whitespace;
  - braces on their own lines, and always around control-statement bodies;
    `if (`, `for (` and `while (` take a space, `forAll(` does not;
  - spaces around `+`, `-`, `&`, `^`, `=`, `==`, `!=`, `<`, `>`, `<=`, `>=`,
    `||` and `&&`; none around `*` and `/`;
  - split long lines after `=` (continuation indented 4 spaces) with
    operators at the start of continuation lines;
  - stream output as `Info<< "..."`;
  - pass primitives by value and larger types by reference; use `const`
    wherever it applies;
  - British spelling in comments and text, `-ize` in identifiers.
- Keep the standard structural comments (section separators, class banners,
  `//-` documentation comments). For a new `.H`/`.C` pair, copy the layout of
  the closest class in the same family.
- New and changed files use the beamFoam licence block. CI runs
  `bin/beamFoamUpdateHeaders -check` on every tracked file; run
  `bin/beamFoamUpdateHeaders <files>` to fix headers. Keep existing
  copyright lines.
- Prefer explicit, readable, local code over clever expressions or
  abstraction-heavy rewrites. Comment only non-obvious behaviour.

## 4) Cases, tutorials, tests and scripts

- Case dictionaries use the beamFoam banner (see any tutorial
  `system/controlDict`); `bin/beamFoamUpdateHeaders` writes it.
- Tutorials run `createBeamMesh`, optionally `setInitialPositionBeam`, then
  `beamFoam`, using `RunFunctions`. Plotting (`gnuplot`, `open`) must stay
  optional so the case still runs where those tools are missing.
- `beamTunnel` needs moorFV and is not a standalone beamFoam case.
- Every tutorial with an `Allrun` should have a `README.md`; copy the layout
  of a neighbouring one.
- A numerical change must say which tutorials change and by how much. A
  refactor should leave tutorial output byte-identical; say so in the PR
  and how it was checked.
- Each standalone tutorial has a `regressionTest.sh`, run by
  `tutorials/Alltest-regression` and by CI after the build. It checks
  results against analytical or published references where they exist, and
  against reference values from `main` to a relative tolerance of `1e-5`.
  Helpers are in `tutorials/RegressionFunctions`; use them rather than
  writing new parsing code.
  - A change meant to alter results updates the reference values in the
    affected `regressionTest.sh` files in the same pull request, and says
    why. Never widen a tolerance, or remove a check, just to make a test
    pass.
  - A new standalone tutorial gets a `regressionTest.sh` that copies the
    layout of a neighbouring one, with a physical or published check where
    one exists. Cases that need moorFV (`beamTunnel`) are not tested here.
- Use bash and standard Unix tools for `Allrun`, `Allclean` and test
  scripts. Python is acceptable for post-processing and verification, with
  the standard library in anything run by default.
- Markdown must pass `markdownlint` (CI runs `markdownlint-cli2`). READMEs
  are also published on the solids4foam website.

## 5) Scope of changes

- Make the smallest correct change; preserve existing architecture, naming
  and dictionary-driven configurability.
- Every changed line should trace to the request, a required build fix, or a
  test or documentation update for it. Do not reformat, rename or clean up
  unrelated code; mention it separately instead.
- Reuse before you write. Before adding a helper, search for an existing one:
  `src/wireBunchingModels/numerics` (listed in its `README.md`: rotations,
  spin tensors, Hermite splines, block solver, parallel checks, fluid
  sampling), `beamModels/beamModel` (`lookupBeamModel`,
  `beamSolverCounters`, cross-sections), the `beamMomentumContribution` base
  class (`checkSingleBeam`, `checkCircularSection`, `unitTangents`), and
  OpenFOAM's own `src`. Extend an existing function rather than copying it,
  and say in the pull request why a new helper was needed. A new helper in
  `numerics` is added to its `README.md`.
- Do not add dictionary options, switches, fallbacks or dependencies unless
  required.
- Behaviour changes must be stated, not silent.
- Split unrelated concerns into separate commits or pull requests.

## 6) Working on a change

- State assumptions when the intended behaviour or expected numerical result
  is unclear. If several interpretations are plausible, ask or present the
  trade-off; push back on unnecessary redesign.
- Read nearby code first and follow its patterns. When adding a class,
  update registration, `Make/files`, tutorials and documentation.
- Build with `./Allwmake -j`, or the narrowest target (`wmake libso` in
  `src/wireBunchingModels`, `wmake` in `applications/solvers/beamFoam`).
- Run the regression tests before opening a pull request, with the
  OpenFOAM environment sourced and the build installed:
  `tutorials/Alltest-regression -j 4` (all cases, about a minute) or
  `tutorials/Alltest-regression <case> ...`. Failures are listed with the
  failing checks; each case's log is `tutorialsTest-regression/<case>.log`.
- Verify with the most relevant tutorials, and compare against `main` when
  results should not change. Report what was and was not verified,
  including the regression test result.
- Pull requests target `main` and use the template in
  `.github/pull_request_template.md`. Keep commit messages short and
  descriptive.

Common reference files: `beamModel.H`, `coupledTotalLagNewtonRaphsonBeam.H`,
`crossSectionModel.H`, `beamMomentumContribution.H` and `Make/files` under
`src/wireBunchingModels`, and `README.md` at the top level.
