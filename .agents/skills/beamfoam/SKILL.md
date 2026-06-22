---
name: beamFoam
description: >
  Use this skill for any task involving the beamFoam repository: editing the
  beamFoam solver, changing beam models or cross-section models, modifying
  beamMomentumContribution classes, updating custom beam patch fields,
  maintaining createBeamMesh or setInitialPositionBeam, adjusting OpenFOAM
  build files, or changing tutorial beam cases. Trigger whenever the user
  mentions beamFoam, createBeamMesh, setInitialPositionBeam, beamProperties,
  coupledTotalLagNewtonRaphsonBeam, crossSectionModel, beamMomentumContribution,
  Simo-Reissner beams, or requests C++ or OpenFOAM changes in this codebase.
---

# Codex Guidelines for beamFoam

This document defines how automated coding changes should be made in this
repository.

## 1) Core Principles

- Make the smallest correct change.
- Preserve existing OpenFOAM-style architecture and naming patterns.
- Prefer consistency with nearby solver and library code over introducing a new
  style.
- Avoid broad refactors unless explicitly requested.
- Keep supported OpenFOAM version compatibility in mind.

## 2) Repository Scope

This repository implements geometrically exact Simo-Reissner beam models in
OpenFOAM using a cell-centred finite-volume formulation.

Key entry points:

- Solver: `applications/solvers/beamFoam/beamFoam.C`
- Beam model library: `src/wireBunchingModels`
- Wave and Morison forcing support: `src/morisonForce`
- Beam mesh utility:
  `applications/utilities/createBeamMesh/createBeamMesh.C`
- Initial-configuration utility:
  `applications/utilities/setInitialPositionBeam/setInitialPositionBeam.C`

The top-level build script is `Allwmake`. It checks for OpenFOAM versions
`v2106`, `v2112`, `v2206`, `v2212`, `v2306`, `v2312`, `v2406`, `v2412`, and
`v2506`. The repository also depends on Eigen via `ThirdParty/`.

## 3) Coding Style Rules

### C++ style

- Follow existing OpenFOAM / repository style in surrounding files.
- Use the same indentation, brace style, comment style, and naming conventions
  as nearby code.
- Keep expressions readable; avoid compact or clever rewrites.
- Prefer explicit local changes over abstraction-heavy redesigns.
- Add comments only when behaviour is non-obvious.

### Headers and includes

- Preserve the existing license/header block format in C++ files.
- Keep include ordering aligned with nearby files.
- Do not change copyright or license text unless explicitly requested.

### Scripts and docs

- Match existing shell style in `Allwmake`, `Allwclean`, tutorial `Allrun`, and
  helper scripts.
- Keep Markdown concise, practical, and repository-specific.
- Do not rewrite documentation just to restyle it.

## 4) OpenFOAM Conventions to Follow

- Respect runtime type and selection table conventions:
  - `TypeName("...")` in headers
  - `defineTypeNameAndDebug(...)` in source files
  - `addToRunTimeSelectionTable(...)` when the module is runtime-selectable
- Preserve dictionary-driven behaviour and runtime configurability.
- Keep field naming and dimension handling consistent with existing solver code.
- Avoid introducing dependencies that break `wmake` or `Allwmake`.
- When adding source files, update the relevant `Make/files` in the same
  module.

Runtime-selectable surfaces in this repository include:

- `beamModel` implementations under `src/wireBunchingModels/beamModels`
- `crossSectionModel` implementations under
  `src/wireBunchingModels/beamModels/beamModel/crossSections`
- `beamMomentumContribution` implementations under
  `src/wireBunchingModels/beamMomentumContribution`
- `morisonForce` implementations under `src/morisonForce`

Rules when adding a new runtime-selectable class:

- Add `TypeName("...")` in the header.
- Register the type in the matching source file.
- Add the source file to the relevant `Make/files`.
- Ensure the dictionary `type` string matches the registered class name.

## 5) Solver and Utility Patterns

The main solver creates a `beamModel` from `constant/beamProperties`, then
loops over:

1. `beam().evolve()`
2. `beam().updateTotalFields()`
3. `beam().writeFields()`

When editing solver-side behaviour:

- Preserve the solver/library split: `beamFoam.C` should stay thin, with beam
  mechanics inside the library.
- Keep `beamProperties` as the primary user-facing configuration surface.
- Preserve the primary unknowns and output fields: beam displacement `W`,
  rotation `Theta`, and derived fields such as `pointW`.
- Be careful with time-step control logic in
  `applications/solvers/beamFoam/readBeamFoamControls.H` and
  `applications/solvers/beamFoam/setDeltaT.H`.

When editing beam setup utilities:

- `createBeamMesh` is the canonical mesh generator for this repo; the solver
  does not rely on `blockMesh`.
- `setInitialPositionBeam` is used when a case starts from a rotated or
  translated beam reference configuration.
- Many tutorials expect `createBeamMesh` first and optionally
  `setInitialPositionBeam` before running `beamFoam`.

When editing tutorials:

- Keep `0/`, `constant/`, and `system/` inputs consistent with solver
  expectations.
- Do not change case dictionaries or plotting scripts unless needed for the
  requested fix.
- Preserve tutorial helper workflow unless the user asked for interface changes.

## 6) Beam-Specific Patterns

- The default beam model currently built in the repo is
  `coupledTotalLagNewtonRaphsonBeam`.
- Cross-section construction is driven by the `beams` entries in
  `constant/beamProperties`, which are consumed by `createBeamMesh` and the
  beam library.
- Existing cross-section models include shapes such as `circle` and
  `rectangle`; follow those implementations when adding another section type.
- External loading and environmental effects are layered through
  `beamMomentumContribution` classes such as `morisonDrag` and `groundContact`.
- Custom patch fields under `src/wireBunchingModels/beamModels/fvPatchFields`
  are part of the solver interface and should remain dictionary-driven.

## 7) Build and Test Workflow

Preferred workflow:

1. Read nearby code and match local patterns.
2. Make the minimum patch necessary.
3. Update `Make/files` or `Make/options` if required.
4. Build the narrowest relevant target first when possible.
5. Run the most relevant tutorial, utility, or compile check available.

Typical build entry points:

- Root build: `./Allwmake`
- Library build: `src/wireBunchingModels/Allwmake`
- Solver build: `wmake` in `applications/solvers/beamFoam`
- Utility build: `wmake all` in `applications/utilities`

Relevant validation surfaces usually include tutorials under `tutorials/`,
especially:

- `tutorials/cantilever`
- `tutorials/3dDynamicCantilever`
- `tutorials/cantileverFollowerForce`
- `tutorials/freeFlexibleBeam`
- `tutorials/beamAndSpringInSeries`

## 8) Minimise Changes

- Only modify files necessary for the requested task.
- Do not reformat unrelated code.
- Do not rename files or symbols unless required.
- Do not change solver behaviour outside requested scope.
- Prefer targeted fixes over opportunistic cleanup.

Before finalizing, verify:

- Build impact is localized.
- Only relevant files changed.
- No registration or `Make/files` updates were missed.
- Dictionary interfaces remain compatible unless the user asked to change them.

## 9) What to Avoid

- Large-scale refactors without explicit request.
- Moving core beam mechanics out of the library layer without a clear need.
- Introducing new coding conventions inconsistent with the repo.
- Silent changes to tutorial or dictionary interfaces.
- Replacing the custom beam-mesh workflow with generic OpenFOAM mesh tools.

## 10) Reference Files

- Main solver: `applications/solvers/beamFoam/beamFoam.C`
- Solver controls:
  `applications/solvers/beamFoam/readBeamFoamControls.H`
  `applications/solvers/beamFoam/setDeltaT.H`
- Beam model base:
  `src/wireBunchingModels/beamModels/beamModel/beamModel.H`
  `src/wireBunchingModels/beamModels/beamModel/beamModel.C`
  `src/wireBunchingModels/beamModels/beamModel/newBeamModel.C`
- Main beam model:
  `src/wireBunchingModels/beamModels/coupledTotalLagNewtonRaphsonBeam/coupledTotalLagNewtonRaphsonBeam.H`
  `src/wireBunchingModels/beamModels/coupledTotalLagNewtonRaphsonBeam/coupledTotalLagNewtonRaphsonBeam.C`
- Cross-section factory:
  `src/wireBunchingModels/beamModels/beamModel/crossSections/crossSection/crossSectionModel.H`
  `src/wireBunchingModels/beamModels/beamModel/crossSections/crossSection/newCrossSectionModel.C`
- Beam momentum contributions:
  `src/wireBunchingModels/beamMomentumContribution/beamMomentumContribution.H`
  `src/wireBunchingModels/beamMomentumContribution/beamMomentumContribution.C`
- Library build list: `src/wireBunchingModels/Make/files`
- Mesh utility: `applications/utilities/createBeamMesh/createBeamMesh.C`
- Initial-position utility:
  `applications/utilities/setInitialPositionBeam/setInitialPositionBeam.C`
