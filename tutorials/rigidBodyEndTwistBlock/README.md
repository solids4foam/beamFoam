# rigidBodyEndTwistBlock

Robustness sweeps comparing partitioned and monolithic coupling of a rod and
a rigid body **welded** to its end, so the body twists the rod. Each run
takes a few seconds and needs only beamFoam.

```sh
./Allrun [moment|stiffness]   # generate, run (NJOBS at once, default 6), summarise
./summarise.py                # summary only; also written to summary.txt
./Allclean
```

## Viewing in ParaView

beamFoam writes the rigid body only to
`postProcessing/rigidBodyEnd/0/rigidBodyEnd.dat`, so `case.foam` shows the
beam alone. `Allrun` runs `./writeRigidBodyVTK <case> --size 0.3` after each
case. The body's centre of mass is on the rod end, so it is drawn as a 0.3 m
cube (three rod diameters) centred there and turned with the body, at every
time directory:

- open `case.foam`, then Warp By Vector with `pointW` (the beam);
- open `VTK/rigidBody.vtk.series` (the body; point fields `U`, `omega`).

The animation cases are not in the repository. They were made by hand from
`base`, not by `makeCases`, all with
`writeInterval 4` (every 0.004 s) so they animate smoothly:

- `animation/{monolithic,partitioned}`: base with M = 205.14 N m (2 rad).
  Both finish.
- `animation/Gx4096/{monolithic,partitioned}`: E and G × 4096 (E/G kept, so
  the Greenhill twist stays 4.1 rad), M = 84025.344 N m (0.2 rad); the
  partitioned stability number is 2.1. Monolithic finishes (twist 0.2012 rad
  at 1.5 s); partitioned diverges in its first step and writes nothing.
  At 2 rad monolithic diverges in step 3 too (the body turns 2.7 rad in one
  step); at 0.5 rad a sideways whirl grows from round-off and it diverges at
  0.98 s. The whirl also grows at 0.2 rad, but only to 1.6e-7 m by 1.5 s.
- `animation/twist6/{monolithic,partitioned}`: `moment/twist6` (M = 615.42 N m,
  1.47 × Greenhill). Partitioned buckles and crashes at 1.328 s, monolithic is
  still straight at 1.5 s. This is the physical instability, not the coupling:
  the sideways growth rate is the same in both and partitioned only starts from
  a larger round-off seed.

## The welded joint

`rigidBodyEnd` takes `joint pin;` (the default, unchanged) or `joint welded;`.
A welded body needs a moment of inertia and a `fixedValue` Theta condition on
its patch. The beam's end face then turns with the body (by the body's step
rotation, from the end face's rotation at the start of the step), and the
beam's end moment acts back on the body:

```text
J alpha + omega x (J omega) + c omega
    = externalMoment + r x beamForce + beamMoment
```

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

## Setup (`base`)

A 2 m rod along +x (radius 0.05 m, E = 5.4326e7 Pa, G = 2.0895e7 Pa, rho = 1,
40 segments) is clamped at `left`. A 1 kg body with J = 0.05 kg m² is welded
to `right`, with its centre of mass at the tip. There is no gravity. A step
moment M about the rod axis acts against `angularDamping` c = 1 N m s. The
oscillation decays at c/(2J) = 10 s⁻¹ whatever the stiffness, so each run
settles on the static twist M L/(G Ip). For a circular rod that answer is
exact at any twist.

| Quantity | Base value |
| --- | --- |
| Torsional stiffness kt = G Ip/L | 102.57 N m/rad |
| Torsional mode sqrt(kt/J) | 45.3 rad/s |
| M | 102.57 N m (static twist 1 rad) |
| Time | dt = 1e-3 s, 1.5 s, trapezoidal Newmark (β = 0.25, γ = 0.5) |
| Partitioned | `nCouplingIterations 50`, `couplingTolerance 1e-10` |

The undamped step response checks the dynamics as well: the twist follows
θs(1 − cos ωt) to 3.4e-4 rad (the trapezoidal rule's period error), the end
moment equals −kt θ, and the two couplings agree to 1e-10 rad.

## Sweeps (`makeCases`)

- `moment`: M, as the static twist it gives. Values: 0.1, 0.5, 1, 2, 3, 4, 6 rad
- `stiffness`: G (so kt) relative to base, M as base. Values: 0.25, 1, 4, 16,
  64, 256, 1024, 4096, 16384

The stiffness sweep raises the partitioned stability number β dt² kt/J from
1.3e-4 to 8.4. Above roughly 4.1 rad of static twist (Greenhill's torque
π EI/L, the same for any rod since E/G = 2.6), a clamped rod under an end
torque does not stay straight.

## Results

Run 10 Nov 2026, 5 cases at once (`summary.txt` has the full table; run
times vary with the load on the machine).

The Newton counts are the `N: Converged` numbers in `log.beamFoam`. Since the
merge of `main` (2 Nov 2026) the Newton loop counter is post-incremented, so
each beam solve reports one fewer than in the 12 Oct runs (monolithic 1
instead of 2). The solves and results are unchanged.

**Moment sweep (kt = 102.57 N m/rad).** Up to 4 rad, both couplings finish:

- The settled twist matches M L/(G Ip) to 3e-7, and the end moment matches M
  to 3e-6.
- Partitioned takes 2.3–2.4 coupling iterations and about 2.3 beam Newton
  iterations per step; monolithic takes 1.
- The two twist histories agree to 5e-7 rad.

At 6 rad (1.47 × Greenhill) the straight rod is unstable. A sideways tip
displacement grows tenfold every 0.1 s from round-off. The partitioned run
buckles and crashes at t = 1.33 s. The monolithic run is still straight at
1.5 s only because its round-off seed is smaller. Run on, it buckles at about
t = 1.8 s and crashes too: once the rod buckles, the constant end torque
spins the body up without limit. This case therefore does not separate the
two couplings.

**Stiffness sweep (M = 102.57 N m).**

| G     | Stab. no. | Coupl. it. (max)   | Newton P | Newton M | Run P / M [s] |
| ----- | --------- | ------------------ | -------- | -------- | ------------- |
| 0.25  | 1.3e-4    | 2.15 (3)           | 2.08     | 1.00     | 2.1 / 1.1     |
| 1     | 5.1e-4    | 2.30 (4)           | 2.29     | 1.00     | 2.2 / 1.3     |
| 4     | 2.1e-3    | 3.04 (4)           | 3.04     | 1.00     | 2.6 / 1.1     |
| 16    | 8.2e-3    | 3.57 (5)           | 3.56     | 1.00     | 3.0 / 1.1     |
| 64    | 0.033     | 5.14 (7)           | 5.14     | 1.00     | 4.0 / 1.1     |
| 256   | 0.13      | 8.53 (12)          | 8.61     | 1.00     | 6.6 / 1.1     |
| 1024  | 0.53      | 28.15 (36)         | 28.50    | 1.00     | 20.7 / 1.2    |
| 4096  | 2.1       | diverges in step 3 | –        | 1.00     | – / 1.1       |
| 16384 | 8.4       | diverges in step 2 | –        | 1.00     | – / 1.1       |

G: G relative to base. Stab. no.: stability number. Coupl. it.: partitioned
coupling iterations per step (max). Newton P / M: Newton iterations per step,
partitioned / monolithic. Run P / M: runtime, partitioned / monolithic.

- Partitioned needs more coupling iterations as the rod stiffens and fails
  between stability numbers 0.53 and 2.1.
- Monolithic takes 1 Newton iteration per step at every stiffness.
- Where both finish, the settled twist matches theory to 4e-5 and the
  couplings agree to 7e-7 rad.
- At G × 4096 and × 16384 the monolithic runs are not yet settled at 1.5 s
  (twist 0.6% and 14.5% high). With ω dt = 2.9 and 5.8, the trapezoidal rule
  damps that mode much more slowly than c/(2J); the leftover amplitude matches
  that slower decay.
