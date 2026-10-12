# rigidBodyEndTwist

Robustness sweeps comparing partitioned and monolithic coupling of a rod and
a rigid body **welded** to its end, so the body twists the rod. The
beamFoam-only sweeps take a few seconds per run; `moorFV/` runs the same kind
of case through moorFV's `FvBeamNewmark` (partitioned) and `beamFoamCoupled`
(monolithic) solvers.

```sh
./Allrun [moment|stiffness]   # generate, run (NJOBS at once, default 6), summarise
./summarise.py                # summary only; also written to summary.txt
./Allclean
moorFV/Allrun                 # moorFV cases, then moorFV/compare.py
```

## The welded joint

`rigidBodyEnd` takes `joint pin;` (the default, unchanged) or `joint welded;`.
A welded body needs a moment of inertia and a `fixedValue` Theta condition on
its patch. The beam's end face then turns with the body (by the body's step
rotation, from the end face's rotation at the start of the step), and the
beam's end moment acts back on the body:

    J alpha + omega x (J omega) + c omega = externalMoment + r x beamForce + beamMoment

with beamMoment = −M on the end face.

- **Partitioned:** each coupling iteration sets the end Theta so the next beam
  solve turns the end face with the body, and passes −M back with −Q. A body
  that only turns, such as one with its centre of mass at the tip, is
  measured against its radius of gyration in the coupling test.
- **Monolithic:** the end-face rotation enters the BlockEigen system through
  the body's rotation columns. The beam rows get the fixedValue-Theta
  coefficients times the rotation Jacobian, and the body moment row gets M
  and its derivatives. After each Newton solve the end face is turned to
  exactly the body's rotation. `jacobianCheck` covers the new terms; they
  match finite differences to round-off.

In moorFV, the `finiteVolumeBeam` restraint takes the same `joint` keyword.
`beamFoamCoupled` passes it to beamFoam. With `FvBeamNewmark` the restraint
turns the end face with the body (`beamModel::setEndStepRotation`, using the
body's rotation over the step) and returns −M as the restraint moment. A
welded joint needs `beamAttachmentKinematics prescribedByRigidBody` (the
default).

## Setup (`base`)

A 2 m rod along +x (radius 0.05 m, E = 5.4326e7 Pa, G = 2.0895e7 Pa, rho = 1,
40 segments) is clamped at `left`. A 1 kg body with J = 0.05 kg m² is welded
to `right`, with its centre of mass at the tip. There is no gravity. A step
moment M about the rod axis acts against `angularDamping` c = 1 N m s. The
oscillation decays at c/(2J) = 10 s⁻¹ whatever the stiffness, so each run
settles on the static twist M L/(G Ip). For a circular rod that answer is
exact at any twist.

| Base value | |
|---|---|
| Torsional stiffness kt = G Ip/L | 102.57 N m/rad |
| Torsional mode sqrt(kt/J) | 45.3 rad/s |
| M | 102.57 N m (static twist 1 rad) |
| Time | dt = 1e-3 s, 1.5 s, trapezoidal Newmark (β = 0.25, γ = 0.5) |
| Partitioned | `nCouplingIterations 50`, `couplingTolerance 1e-10` |

The undamped step response checks the dynamics as well: the twist follows
θs(1 − cos ωt) to 3.4e-4 rad (the trapezoidal rule's period error), the end
moment equals −kt θ, and the two couplings agree to 1e-10 rad.

## Sweeps (`makeCases`)

| Sweep | Varies | Values |
|---|---|---|
| `moment` | M, as the static twist it gives | 0.1, 0.5, 1, 2, 3, 4, 6 rad |
| `stiffness` | G (so kt) relative to base, M as base | 0.25, 1, 4, 16, 64, 256, 1024, 4096, 16384 |

The stiffness sweep raises the partitioned stability number β dt² kt/J from
1.3e-4 to 8.4. Above roughly 4.1 rad of static twist (Greenhill's torque
π EI/L, the same for any rod since E/G = 2.6), a clamped rod under an end
torque does not stay straight.

## Results

Run 12 Oct 2026, 6 cases at once (`summary.txt` has the full table; run
times vary with the load on the machine).

**Moment sweep (kt = 102.57 N m/rad).** Up to 4 rad, both couplings finish:

- The settled twist matches M L/(G Ip) to 3e-7, and the end moment matches M
  to 3e-6.
- Partitioned takes 2.3–2.4 coupling iterations and about 4.6 beam Newton
  iterations per step; monolithic takes 2.
- The two twist histories agree to 5e-7 rad.

At 6 rad (1.47 × Greenhill) the straight rod is unstable. A sideways tip
displacement grows tenfold every 0.1 s from round-off. The partitioned run
buckles and crashes at t = 1.33 s. The monolithic run is still straight at
1.5 s only because its round-off seed is smaller. Run on, it buckles at about
t = 1.8 s and crashes too: once the rod buckles, the constant end torque
spins the body up without limit. This case therefore does not separate the
two couplings.

**Stiffness sweep (M = 102.57 N m).**

| G factor | Stability number | Partitioned coupling iterations / step (max) | Partitioned Newton / step | Monolithic Newton / step | Runtime part. / mono. [s] |
|---|---|---|---|---|---|
| 0.25  | 1.3e-4 | 2.15 (3)  | 4.23  | 2.00 | 2.1 / 1.1 |
| 1     | 5.1e-4 | 2.30 (4)  | 4.58  | 2.00 | 2.3 / 1.2 |
| 4     | 2.1e-3 | 3.04 (4)  | 6.08  | 2.00 | 3.0 / 1.2 |
| 16    | 8.2e-3 | 3.57 (5)  | 7.13  | 2.00 | 3.5 / 1.2 |
| 64    | 0.033  | 5.14 (7)  | 10.28 | 2.00 | 4.7 / 1.2 |
| 256   | 0.13   | 8.53 (12) | 17.14 | 2.00 | 7.8 / 1.2 |
| 1024  | 0.53   | 28.15 (36) | 56.66 | 2.00 | 24.4 / 1.2 |
| 4096  | 2.1    | diverges in step 2 | – | 2.00 | – / 1.2 |
| 16384 | 8.4    | diverges in step 3 | – | 2.00 | – / 1.2 |

- Partitioned needs more coupling iterations as the rod stiffens and fails
  between stability numbers 0.53 and 2.1.
- Monolithic takes 2 Newton iterations per step at every stiffness.
- Where both finish, the settled twist matches theory to 4e-5 and the
  couplings agree to 7e-7 rad.
- At G × 4096 and × 16384 the monolithic runs are not yet settled at 1.5 s
  (twist 0.6% and 14.5% high). With ω dt = 2.9 and 5.8, the trapezoidal rule
  damps that mode much more slowly than c/(2J); the leftover amplitude matches
  that slower decay.

## moorFV (`moorFV/`)

The same rod, 2.6565 m long to reach the body centre (kt = 77.22 N m/rad,
60 segments), is welded to moorFV's floating box (1 kg, J = 0.05 kg m², centre
of mass at the rod end). The box is in interFoam with a coarse mesh (24k
cells), all air, `test true` (fluid forces ignored) and no gravity. The body
starts spinning about the rod axis at Omega0 = 7.86 rad/s, so the twist is a
free torsional oscillation (Omega0/omega) sin(omega t) with amplitude 0.2 rad
and omega = 39.3 rad/s. dt = 0.002 s, 0.5 s, 3 PIMPLE outer correctors with
`moveMeshOuterCorrectors`. `Allrun` meshes once (`mesh/`), runs
`partitioned` (`FvBeamNewmark`, `accelerationRelaxation 0.4`, as in
floatRigidBody_partioned) and `monolithic` (`beamFoamCoupled`) side by side,
then `compare.py` writes `comparison.txt`.

| Coupling | Peak twist / theory, per half period | Largest deviation from theory | Run time |
|---|---|---|---|
| monolithic, 3 outer correctors | 0.9998 every half period | 1.0% of the amplitude (period error) | 170 s |
| partitioned, 3 outer correctors | 1.028, 1.091, 1.157, 1.227, 1.301, 1.380 | 38% | 175 s |
| partitioned, 8 outer correctors (0.25 s) | 1.0016, 1.0052, 1.0089 | | 205 s |
| partitioned, 20 outer correctors (0.25 s) | 0.9998, 0.9998, 0.9998 | | 452 s |

The partitioned coupling iterates only through the PIMPLE outer correctors,
and the body is driven by the beam state of the previous corrector. Not
iterated to convergence, it feeds energy into the torsional mode: the
amplitude grows by 6–7% per half period with 3 correctors and by 0.4% with
8. With 20 it matches the monolithic result, so the welded coupling is
consistent; the monolithic solve gets there with 3 correctors.

Two fixes were needed for these runs, and they apply to pinned lines too:

- `rigidBodyEnd` now starts its first time step from its initial state.
  moorFV creates the body during the first time step, which used to lose an
  initial velocity or spin.
- BlockEigen no longer counts unused rigid-body rows in the beam residual.
  A body supplied by moorFV without BlockEigen coupling left dt times its
  angular momentum there, and a beam that starts undeformed then could not
  converge.
