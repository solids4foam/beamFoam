# Phase 5 plan: several mooring lines on one monolithic rigid body

Status: Stage 0 built, not run (13 October 2026); Stages 1-5 not started.
Follows Phases 0-4 and the
welded joint (see `monolithicSolverTestProgressReport.pdf` and
`source/monolithicCouplingPlan.html`, Phase 5).

## Why it is single-line today

- **moorFV side:** `beamFoamCoupled::solve`
  (`sDoFRGBFvBeam/sixDoFFvBeamSolvers/beamFoamCoupled/beamFoamCoupled.C:76-110`)
  stops with an error if there is more than one `finiteVolumeBeam` restraint.
- **Who owns the body:** each `finiteVolumeBeam` owns its own `beamModel`, and
  each `beamModel` holds its own copy of the body (`rigidBodyEnd`). That class
  mixes two things:
  - the body itself: mass, inertia, Newmark state, external load, constraints;
  - one attachment: patch name, joint type, moment arm and beam force.
- **Who runs Newton:** `beam.evolve()` runs the whole Newton loop internally,
  so nothing outside a line can get in between its iterations.
- **The coupling data:** `RigidBodyMonolithicCoupling`
  (`RigidBodyCouplingData.H`) joins the body to a single attachment cell, and
  BlockEigen solves it with a direct SparseLU.
- **Test case:** the validation case is the floatRigidBody box (one line,
  10 s wave runs in `tutorials/floatRigidBody`) with a second, mirror-image
  line on the other side of the box. It keeps the October fixes (line
  gravity, `refWf` placement, written Lambda, ground contact, region-aware
  function objects) and the existing one-line results to compare against.
  `tutorial/mooringCase-H12T20` (4 lines, last run Aug 2025) is not used.

## Method: exact static condensation inside one Newton loop

For line *i*, write the Newton system with the beam increments ΔUᵢ and the
body increments Δq (6 unknowns):

```
[ Kᵢ  Cᵢ ] [ΔUᵢ]   [rᵢ]      Cᵢ: nonzero only in the attachment-cell rows
[ Bᵢ  Dᵢ ] [Δq ] = [sᵢ]      Bᵢ: nonzero only in the attachment-cell columns
```

Eliminating each line and summing over the lines gives a 6×6 system for the
body:

```
[ K_body + Σ(Dᵢ − Bᵢ Kᵢ⁻¹ Cᵢ) ] Δq = s_body + Σ(sᵢ − Bᵢ Kᵢ⁻¹ rᵢ)
then  ΔUᵢ = Kᵢ⁻¹ (rᵢ − Cᵢ Δq)   for each line
```

- **Exact:** this is algebraically the same as one large matrix, so it keeps
  the quadratic Newton convergence and the stiffness-independent iteration
  counts.
- **Cheap:** Kᵢ⁻¹Cᵢ needs only 6 extra back-substitutions with the line's
  existing LU factorisation.
- **Constraints:** the projection `P A P + C` must be applied once, to the
  summed 6×6 system. Applying it per line would add C once per line.
- **Parallel:** each processor can condense its own lines, and only a 6×6
  matrix plus a 6-vector needs to be summed across processors. This is the
  main reason to prefer condensation over merging all lines into one
  multi-beam mesh.
- Nothing here depends on the number of lines; the two-line case is the
  first target, and the method is written for N.

## Stages

### Stage 0: two-line reference case (no code changes) - built, not run

`tutorials/floatRigidBody_twoLine/floatRigidBody_partitioned`, built from the
tracked files of `tutorials/floatRigidBody/floatRigidBody_partioned`:

- `beamInit/makeBeams` builds `beamone` (anchor (0.25 0.025 -0.15) to the
  box's lower-left edge (0.45 0.025 -0.025), as before) and `beamtwo`, its
  mirror image about x = 0.5 (anchor (0.75 0.025 -0.15) to the lower-right
  edge (0.55 0.025 -0.025)). Same line properties, strain 0.
- `constant/dynamicMeshDict`: restraint `beamtwo` added; box `mass`
  0.124319 -> 0.123863 to balance two line weights (check against the z of
  `Q` on the `right` patch of both lines after `makeBeams`). The two
  horizontal line pulls cancel, so the one-line case's slow drift should go.
- `makeBeams` already copies `g` and the `ref*` fields for each line, and the
  ground contact in `constant/beamMomentumContributionProperties` applies to
  both lines.
- Beam function objects stay on beamone only: their output files are named by
  patch, not region, so beamtwo's would overwrite beamone's. The restraint's
  own files are per line (`axialForcebeamtwo.dat` etc.). Per-region file
  names are a small code change, to do before Stage 5.
- Still to do: run it partitioned with 3 outer correctors (as set) and an
  8-corrector copy over a short window to give the reference. Ask before
  starting any long run.

### Stage 1: refactor one line into condense and back-substitute steps, with no change in behaviour

- Split `rigidBodyEnd` into a body part and an attachment part (patch, joint,
  arm r₀ᵢ, beam force and moment, and their derivative terms). The body keeps
  a list of attachments.
- Break the monolithic Newton iteration (`evolveMonolithicRigidBodyEnd`) into
  four calls:
  - `beginStep` (predictors);
  - `assembleCondensed(q)`, which returns Sᵢ, gᵢ and the line's residual;
  - `applyBodyIncrement(Δq)`, which back-substitutes and returns the line's
    increment norm;
  - `acceptStep`.
- Add a BlockEigen mode that solves extra right-hand-side columns and keeps the
  factorisation. The legacy staged BlockEigen branches stay untouched.
- Fold in report priority 6 here: include the body state in the X-norm
  denominator.
- **Check:** the 54/54 rigidBodyEnd suite and the twist suite pass, the
  histories match the current code to about 1e-12, and the Newton iteration
  counts are identical.

### Stage 2: an N-line coordinator in the beamFoam library

- A new class holds the body and references to N `beamModel`s. Each iteration
  it:
  1. condenses every line;
  2. sums and constrains the 6×6 system;
  3. solves it;
  4. back-substitutes into every line.
- It stops when the body and every line have converged. All lines must share
  Δt.
- **Finite-difference check:** perturb q, re-solve each line with the body
  fixed at that perturbed state, and compare dF/dq against ΣSᵢ.
- Keeping the coordinator in the library means standalone tests can use it
  without interFoam, through a small driver app that reads several beam
  regions.

### Stage 3: standalone verification suite (`tutorials/rigidBodyLines`)

- **Twin-line equivalence:** two identical lines at one point should match one
  line with twice the EA, EI, GJ and ρA. This is exact, and reuses the
  existing analytic cases.
- **Two (and four) splayed lines holding a hanging mass:** each line carries
  mg/2 (mg/4) vertically, the horizontal resultant is zero, and the body
  doesn't rotate.
- **Two opposed taut lines, small oscillation:** surge, heave and pitch
  frequencies compared with linearised line stiffness (the same motions as
  the CFD case).
- **Asymmetric pretension:** the body tilts until Σ rᵢ × Fᵢ balances.
- **Cross-method:** agrees with a tightly iterated partitioned run to within
  1e-4 of peak.
- **Stiffness sweep:** with S_eff = βΔt² Σ EAᵢ/(Lᵢ m), more lines should make
  the partitioned solver fail sooner while the monolithic iteration count
  stays flat.
- The `Allrun` and check script return a non-zero exit code on failure, like
  `checkRigidBodyEnd.py`.

### Stage 4: moorFV integration

- `beamFoamCoupled` collects every `finiteVolumeBeam` and passes its beam
  model to the coordinator.
- Per-line history files stay as they are, and each line can be pin or welded.
- **Check:** the single-line floatRigidBody run is unchanged and
  `checkCoupledState.py` still passes.

### Stage 5: two-line CFD validation (the plan's pass criterion)

- Give the beam function objects per-region file names and add them for
  beamtwo.
- Make `floatRigidBody_twoLine/floatRigidBody_monolithic` from the Stage 0
  case (`solver { type beamFoamCoupled; }`, as in the one-line monolithic
  case).
- Run monolithic at 1, 3 and 8 outer correctors against the Stage 0
  reference, then one long (10 s) run.
- Compare the six body motions, the axial force at each line's fairlead and
  anchor, Newton iterations per step and runtime.
- Optionally run an E×100 variant to show the partitioned solver breaking
  down.
- Later, a 4-line case (for example a 3-D version of the box with lines at
  its corners) once two lines pass.

### After Phase 5

Parallel support (the condensed-sum design makes it a small step) and restarts
from a displaced body. Both are prepared for by this design but are outside
Phase 5's acceptance criteria.

## Risks

- The monolithic and legacy branches are interleaved in
  `coupledTotalLagNewtonRaphsonBeamEvolve.C`. The single-line round-off check
  in Stage 1 is the guard against breaking either.
- The derivative of the moment term rᵢ × Fᵢ with respect to body rotation
  currently uses the single attachment arm. It has to move into each line's
  Dᵢ.
- Lines may have different Newton settings, so the coordinator needs one
  convergence rule for all of them.
- The centre of rotation must still equal the centre of mass, and restarts are
  still rejected.

## Open decision: how the lines' Newton loops are driven

- **One Newton loop for all lines (recommended):** this is true monolithic
  coupling and keeps the 3-4 iterations per step. It needs the Stage 1
  refactor.
- **Nested Newton (the alternative):** for each body iterate, each line solves
  to convergence with its end prescribed, then reports its exact condensed
  stiffness. It barely touches `evolve()` and is still stiffness-robust, but
  costs roughly 2-3× more line iterations.
