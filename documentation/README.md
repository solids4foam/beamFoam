# documentation/

Notes on the beam + rigid body coupling in beamFoam and moorFV: how the
monolithic solver works, how the code is laid out, and how it was developed.

The files fall into four groups. The first group describes the code as it
is now; the last group is history and describes code that has since been
removed.

## Suggested reading order

1. `pimpleToRigidBodySolverCallChain.txt`: how a run reaches the solver
2. `monolithicRigidBodyCodePath.txt`: what the solver does
3. `monolithicMatrixLayout.txt`: the matrix it builds
4. `twoCellMonolithicMatrix.txt`: the same matrix, small enough to read whole
5. `newmarkAndMonolithicBodyRow.txt`: why the body row is the Newmark scheme

## How the current code works

- `pimpleToRigidBodySolverCallChain.txt`: How a moorFV run gets from interFoam's
  PIMPLE loop to `beamFoamCoupled::solve`. Starts with a primer on classes, base
  and derived classes, virtual functions and run-time selection, then lists
  where each class lives and follows the set-up and the per-time-step call
  chain, down to `eigenSolver.solve`.
- `monolithicRigidBodyCodePath.txt`: How the monolithic solver handles the rigid
  body, compared with the partitioned solver: the call chain from
  `beamFoamCoupled::solve` through `finiteVolumeBeam::solveMonolithic`,
  `evolveBeam` and its Newton loop, `setMonolithicCoupling` and
  `BlockEigenSolverOF::solve`, with full paths and line numbers. Ends with what
  was removed with the legacy BlockEigen coupling (2026-10-20).
- `monolithicMatrixLayout.txt`: Where every coefficient of the coupled matrix
  sits: 6 x 6 blocks, the naming of the coupling blocks, the index expressions
  (`beamRow`, `rbRow`), a three-cell example, the diagonal, the attachment face,
  rotation terms, a table of all body blocks. Section 9 explains what each
  coefficient is physically and where it is computed. Section 10 explains why
  the Newmark predictor has no beta dt^2 a_{n+1} term.
- `twoCellMonolithicMatrix.txt`: The full matrix for a two-cell beam with a
  pinned body, written in physical terms. Defines every symbol, a notation
  section (superscripts and subscripts, Greek letters, B_b, K^M_b, A, J,
  OpenFOAM operators such as `&`, `^` and unary `*`), the meaning of every entry
  row by row, where the action-reaction balance sits and why the matrix is not
  symmetric.
- `newmarkAndMonolithicBodyRow.txt`: How the body's matrix row relates to the
  Newmark-beta update: what x_pred is, the residual, its Newton linearisation,
  and a step-by-step proof that a converged Newton iterate satisfies the Newmark
  equation exactly. Also compares with the partitioned `newmarkDisplacement()`.
- `outerCorrectorsAndNewmark.pdf`: PIMPLE outer correctors versus beamFoam
  Newton iterations, and how each relates to the Newmark scheme in the
  partitioned and monolithic frameworks. Includes a linear stability model with
  figures.

## Plans, reports and write-ups

- `monolithicCouplingPlan.pdf`: The development plan for the monolithic coupling
  (Phases 0-5), the target formulation, and progress through Phase 4.
- `monolithicSolverTestProgressReport.pdf`: Test progress and findings for the
  monolithic solver.
- `phase5MultiLinePlan.md`: Plan for Phase 5: several mooring lines on one
  monolithic rigid body. Why it is single-line today, and the stages to get
  there.
- `rigidBodyCouplingPRSummary.pdf`: Pull request summary and code review of the
  beam / rigid body coupling.
- `pullRequestSummary.pdf`: Commit-by-commit summary of both repositories'
  changes, generated from the git histories.
- `paperOutline.txt`: A section-by-section layout for writing the monolithic
  coupling up as a numerical methods paper, with which file in this folder holds
  the material for each section, figures to prepare, and gaps to fill first.

## Sources of the PDFs (`source/`)

- `source/README.md`: How to rebuild `monolithicCouplingPlan.pdf` and
  `outerCorrectorsAndNewmark.pdf`.
- `source/monolithicCouplingPlan.html`, `source/outerCorrectorsAndNewmark.html`:
  HTML sources of those PDFs.
- `source/model.py`: Linear body-spring model of the partitioned loop and
  monolithic update (amplification matrices).
- `source/order.py`: Order-of-accuracy check of that model.
- `source/figs.py`, `source/growth.png`, `source/newton.png`: Figures for
  `outerCorrectorsAndNewmark.pdf`.
- `source/pullRequestSummary/`: Scripts and hand-written commit explanations
  that generate `pullRequestSummary.pdf`. See its own `README.md`.

## History: the legacy BlockEigen coupling (removed 2026-10-20)

These notes describe earlier coupling approaches. The code they describe
(`blockEigen*` switches, `RigidBodyStepData`, `RigidBodySolution`,
`solvedByBlockEigen`, `blockEigenAuthoritativeMotion` and others) was
removed in beamFoam `ca153c8` and moorFV `2f203c7`. Those keywords now stop
the run with an error. Keep these files only as a record of how the design
developed.

- `rigidBodyAndBeamCoupledSolverDevelopmentLog.txt`: Development log of the
  coupled rigid-body / beam BlockEigen path (July 2026).
- `BLOCK_EIGEN_SOLVER_COUPLING_CHANGES.md`: Changes to `BlockEigenSolverOF` to
  add rigid body rows to the beam's matrix.
- `BLOCK_EIGEN_RIGID_BODY_COUPLING_PROGRESS.md`: Progress moving the
  `rigidBodyAndBeam` case from the sixDoF restraint-force path to the
  BlockEigen-solved body motion.
- `BEAM_ATTACHMENT_KINEMATICS_OWNERSHIP.md`: Which component owned the beam
  attachment kinematics at the time.
- `rigidBodyKinematicCouplingDevelopment.txt`: The two competing owners of the
  attachment kinematics, and how that was resolved.
- `rigidBodyImplicitAttachmentCoefficientChanges.txt`: Change of the body
  translational rows to a Newmark / mass-scaled attachment stiffness.

Today a run is switched between partitioned and monolithic coupling in one
place: the moorFV motion solver type (`FvBeamNewmark` or `beamFoamCoupled`),
or `rigidBodyCoupling partitioned;` / `monolithic;` for a beamFoam-only
`rigidBodyEnd` case.
