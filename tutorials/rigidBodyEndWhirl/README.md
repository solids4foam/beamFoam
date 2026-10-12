# rigidBodyEndWhirl

beamFoam-only robustness sweeps comparing partitioned and monolithic
`rigidBodyEnd` coupling under an applied rotation. No fluid; each run takes
under a minute.

```sh
./Allrun [moment|stiffness]   # generate, run (NJOBS at once, default 6), summarise
./summarise.py                # summary only; also written to summary.txt
./Allclean
```

## Setup (`base`)

A 2 m rod along +x is clamped at `left`. A 1 kg body is pinned to `right`,
with its centre of mass 0.1 m from the tip, perpendicular to the rod
(J = 0.05, 0.01, 0.01 kg m², so x is the largest-inertia axis about the pin).
There is no gravity. A constant `externalMoment` M about the rod axis acts
against `angularDamping` c = 0.03 N m s. The body spins up towards M/c and
whirls the rod tip round the x axis.

`rigidBodyEnd` is a pin: the beam passes the body only a force. The applied
moment reaches the rod through the arm from the tip to the centre of mass,
as the centrifugal pull of the whirling body.

| Base value | |
|---|---|
| Rod | radius 0.05 m, E = 5.4326e7 Pa, G = E/2.6, rho = 1 (rod mass 0.016 kg), 40 segments |
| Bending stiffness EI | 266.67 N m² (tip stiffness 3EI/L³ = 100 N/m) |
| Critical speed ω_c | sqrt(3EI/L³ / m) = 10 rad/s |
| M | 0.225 N m (terminal spin 7.5 rad/s = 0.75 ω_c) |
| Time | dt = 1e-3 s, 15 s, trapezoidal Newmark (β = 0.25, γ = 0.5) |
| Partitioned | `nCouplingIterations 50`, `couplingTolerance 1e-10` |

The rod is thick so its axial stiffness, EA/L = 2.13e5 N/m, stays low. The
partitioned coupling iterations converge only while roughly
β dt² EA/(L m_pin) < 1, where m_pin = 1/(1/m + a²/J) = 0.5 kg is the mass
the pin sees (the body is pinned off its centre of mass, so a force at the
pin also turns it). The number is 0.107 at base. A first try with radius
0.01 m at the same EI (EA/L = 5.3e6 N/m, number 10.6) diverged in the
first partitioned time step.

## Sweeps (`makeCases`)

Each case is written to `<sweep>/<label>/{partitioned,monolithic}`.

| Sweep | Varies | Values | Terminal spin M/(c ω_c) |
|---|---|---|---|
| `moment` | M [N m] | 0.075, 0.15, 0.225, 0.27, 0.3, 0.45, 0.9 | 0.25 – 3 |
| `stiffness` | EI relative to base (E and G scaled together) | 0.25, 0.5, 0.75, 1, 2, 4, 16, 64 | 1.5 – 0.094 |

The stiffness sweep also raises the axial stiffness, and with it the
partitioned stability number β dt² EA/(L m_pin): from 0.027 to 6.8.

## Checks (`summarise.py`)

- Linear whirl theory: at spin rate ω the tip runs on a circle of radius
  a r²/|1 − r²|, r = ω/ω_c (a = 0.1 m; rod mass lumped as 33/140 of it at
  the tip). This is checked at the final time and skipped near resonance
  (|1 − r²| < 0.3).
- Partitioned against monolithic: the largest tip separation over the run,
  relative to the largest tip radius.
- Cost and robustness: status, time reached, beam Newton iterations per step,
  Newton failures, partitioned coupling iterations (mean, max, and steps that
  hit the cap), and run time.

In the base case the tip radius follows the linear theory to within 3%
through the spin-up, and the two couplings agree to 6e-7 m.

## Results

Run 12 Oct 2026, 6 cases at once (`summary.txt` has the full table).

**Moment sweep (EI = 266.67 N m², ω_c = 10 rad/s).** Every case finishes
with both couplings, and the two tip paths agree to 1.2e-4 of the tip
radius or better. Partitioned takes 3.4–3.5 coupling iterations and 6.8–9.6
beam Newton iterations per step at every moment; monolithic takes 3.0–3.9
Newton iterations and runs 2.2–2.4 times faster (14–19 s against 32–46 s).
Where the check applies, the final tip radius is within 5.3% of linear
theory (−0.8% at M = 0.225 N m).

The spin does not reach M/c above about 0.8 ω_c. The applied moment goes
into the whirl instead and the spin stalls just below the critical speed
(the Sommerfeld effect). At M = 0.9 N m (M/c = 3 ω_c) the spin is 10.3 rad/s
at 15 s, and the tip radius has grown to 0.82 m on the 2 m rod. Neither
coupling has trouble with this, so the moment sweep tests large rotations
and deflections more than it tests the coupling.

**Stiffness sweep (M = 0.225 N m).**

| EI factor | Stability number | Partitioned coupling iterations / step (max) | Partitioned Newton / step | Monolithic Newton / step | Runtime part. / mono. [s] |
|---|---|---|---|---|---|
| 0.25 | 0.027 | 2.21 (4)  | 6.90  | 3.00 | 32 / 15 |
| 0.5  | 0.053 | 2.46 (5)  | 7.17  | 3.00 | 33 / 14 |
| 0.75 | 0.080 | 3.31 (6)  | 8.41  | 3.86 | 39 / 19 |
| 1    | 0.107 | 3.46 (6)  | 8.74  | 3.90 | 40 / 19 |
| 2    | 0.213 | 3.79 (8)  | 9.41  | 3.92 | 42 / 19 |
| 4    | 0.427 | 6.50 (14) | 14.83 | 3.94 | 64 / 19 |
| 16   | 1.707 | diverges in step 2 | – | 3.96 | – / 19 |
| 64   | 6.827 | diverges in step 1 | – | 3.97 | – / 18 |

Partitioned coupling slows as the stability number rises and fails between
0.43 and 1.7: the coupling iterations grow instead of shrinking, and the
beam Newton solve then diverges ("Residual grew excessively"). Monolithic
cost does not change with stiffness: 3.0–4.0 Newton iterations per step,
and the tip radius stays within 1.5% of linear theory. Where both finish,
the tip paths agree to 4.7e-4 of the tip radius or better.
