# rigidBodyEnd test cases

Fast beamFoam-only tests for a rigid body attached to the end of a beam
(Phase 0 of `codexLogs/monolithicCouplingPlan.pdf`). No fluid; each case runs
in a few seconds.

```sh
./Allrun      # runs all three cases, then checkRigidBodyEnd.py
./Allclean
```

`checkRigidBodyEnd.py` reads the parameters from each case and compares the
body history (`postProcessing/rigidBodyEnd/0/rigidBodyEnd.dat`) with analytic
results. It exits non-zero if a check fails.

## Cases

All use a rod along +x, 2 m long, pinned at `left`, with a 10 kg point mass
on `right`. Gravity acts on the body and the rod.

| Case | Setup | Checked against |
|---|---|---|
| `hangingMass` | Rod hanging along +x; damped Newmark on the body (β = 0.49, γ = 0.9) so it settles | End force = −mg; stretch = (mg + ½ rod weight)·L/EA |
| `axialOscillation` | As above, undamped (β = 0.25, γ = 0.5), released unstretched | Period of the exact rod-plus-end-mass mode (z tan z = rod mass/end mass), with the trapezoidal-rule period elongation; mean and peak displacement |
| `pendulum` | Thin rod clamped at the top, gravity tilted 5° | Pendulum period with finite-amplitude correction; effective length includes the stretch and the bending boundary layer √(EI/T) at the clamp |

The pendulum rod is clamped at the top because a rod pinned at both ends can
twist freely about its own axis, which leaves the beam solve almost singular.
It is thin enough that the clamp's bending stiffness hardly changes the period.

## Coupling settings

The body is coupled to the beam with `rigidBodyCoupling partitioned`: each
time step the body is moved with the beam force, the beam is solved with the
body position as its end condition, and this repeats up to
`nCouplingIterations` times until the displacement change is below
`couplingTolerance`.

**One pass per step is not stable for these cases.** With
`nCouplingIterations 1` the body always uses the beam force from the previous
solve; in `axialOscillation` the amplitude then grows by about 15% per period.
The cases therefore iterate to convergence (2–6 iterations per step), which
gives the reference that the monolithic solve (Phase 1) must reproduce.
