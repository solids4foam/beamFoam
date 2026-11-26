# beamFoam pull request

## Problem, solution, and verification

<!-- Use short sentences. Explain terms that other developers may not know.
Link related issues. Keep this summary useful without the technical details. -->

**Problem:**

**Solution:**

**Verification:** State what you tested and the result. Say what is untested.

Closes #<!-- issue number, or delete this line -->

## Before review

<!-- Reviewers check these. Tick them or say why they do not apply. Keep the
PR in draft until they are done. -->

- [ ] New or renamed source files are listed in the module's `Make/files`,
  and new runtime-selectable classes are registered.
- [ ] Existing helpers are reused (see
  `src/wireBunchingModels/numerics/README.md`); any new helper is explained
  below and, if in `numerics`, listed in that README.
- [ ] Added and changed lines follow the
  [OpenFOAM coding style](https://openfoam.org/dev/coding-style-guide/).
- [ ] Tutorial results: unchanged (say how you checked), or the changed
  cases and the size of the change are listed below.
- [ ] Dictionary or boundary condition changes are documented, with defaults
  that keep existing cases working.
- [ ] If an AI agent wrote or reviewed part of this change, a person has read
  every changed line and checked the results.

## Technical details

<!-- Keep only the details that apply. Use “Not applicable” where needed. -->

- **Implementation:** Explain the main changes and why you chose this approach.
- **Numerical considerations:** State assumptions and effects on accuracy,
  stability, conservation, and convergence. Link the method or reference.
- **Evidence:** Give case names, commands, OpenFOAM versions, reference
  results, and error tolerances. Link logs or plots where useful.
- **Compatibility:** Describe changes to dictionaries, defaults, dependencies,
  existing cases, coupled solvers (moorFV, solids4foam), or solver behaviour.
- **Review focus:** Identify uncertainties and decisions that need human review.

See the [contribution guide][contributing] for how to work on a change.

[contributing]: https://github.com/solids4foam/beamFoam/blob/main/CONTRIBUTING.md
