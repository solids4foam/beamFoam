# Agent instructions for beamFoam

Read `.agents/skills/beamfoam/SKILL.md` before writing or reviewing code in
this repository. It holds the full rules; the most important are:

- Make the smallest correct change. Do not reformat, rename or clean up code
  the task does not need.
- Reuse before you write: search `src/wireBunchingModels/numerics` (see its
  `README.md`), `beamModels/beamModel` and OpenFOAM's `src` before adding a
  helper.
- New sources go in the module's `Make/files`; runtime-selectable classes
  need `TypeName`, `defineTypeNameAndDebug` and
  `addToRunTimeSelectionTable`.
- Follow the OpenFOAM coding style in added and changed lines, and run
  `bin/beamFoamUpdateHeaders` on new files.
- Build with `./Allwmake -j`, run the affected tutorials, and state what was
  and was not verified. Say whether tutorial results changed.
- Pull requests target `main` and use `.github/pull_request_template.md`.
  See `CONTRIBUTING.md`.
