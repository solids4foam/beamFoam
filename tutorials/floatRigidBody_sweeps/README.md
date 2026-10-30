# floatRigidBody_sweeps

Parameter sweeps of `floatRigidBody`, comparing moorFV's partitioned
(`FvBeamNewmark`) and monolithic (`beamFoamCoupled`) coupling of the box and
its mooring line. Like `rigidBodyEndTwist`, each case is run with both
couplings.

```sh
./Allrun [stiffness|period|animation]   # generate, run (NJOBS at once, default 10), summarise
./summarise.py                          # summary only; also written to summary.txt
./Allclean
../plot_floatRigidBody_sweeps.py        # plots (Spyder, or SAVE_FIGURES = True)
```

Each case runs serially on one core and takes about 75 minutes with ten
running at once.

## Setup (`base`)

`base` is `../floatRigidBody/floatRigidBody_monolithic`: a
0.1 × 0.05 × 0.05 m box floats half submerged in a 1 m tank of water
0.15 m deep. It is free to surge, heave and pitch. One line runs from the tank
floor at (0.25, −0.15) to the box at (0.45, −0.025):
- length 0.236 m
- radius 1 mm
- E = 5e6 Pa, so EA/L = 66.6 N/m
- no pretension

StokesII waves are 0.02 m high with a 1 s period, ramped up over 2 s. The run
lasts 10 s, with 3 PIMPLE outer correctors and the time step set by Courant
number (at most 2 ms).

`makeCases` changes only:
- **Solver type:** the `solver` type in `constant/dynamicMeshDict`.
- **Line stiffness:** E and G = E/2.8 in `beamInit/makeBeams`. `beamInit`
  rebuilds the line for each case.
- **Wave period:** `wavePeriod` in `constant/waveProperties`.
- **Writes:** the sweep cases write fields at 5 and 10 s only, because all
  results come from the function objects in `postProcessing`. The animation
  cases write compressed fields every 0.05 s.

The box mass in `constant/dynamicMeshDict` balances the base line's weight on
the box (0.00447 N). It is not changed for the other stiffnesses, because the
line's weight on the box changes very little with E.

## Sweeps

| Sweep | Varies | Values |
|---|---|---|
| `stiffness` | line E (and G) relative to base, 1 s waves | 0.25, 1, 4, 16, 64, 256, 1024 |
| `period` | wave period (s), base line | 0.5, 0.6, 0.8, 1.0, 1.25, 1.5 |
| `animation` | none: base, for ParaView (`case.foam`) | |

- **Stiffness sweep:** the line's axial mode on the box, sqrt(EA/(L m)), goes
  from 12 to 740 rad/s. At the top of the range that is 8.5 ms per period, or
  4 to 8 time steps.
- **Period sweep:** the wavelength goes from 0.39 to 1.7 m, starting near the
  box's heave period (roughly 0.4 s). The longest period stops at an Ursell
  number of about 18, inside StokesII's range.

## Output

`summary.txt` (from `summarise.py`) gives, for each case and coupling:
- status and the time reached
- time steps
- beam Newton iterations per step (mean and max, over all outer correctors)
  and in total
- run time
- surge, heave and pitch amplitudes for t ≥ 5 s
- the largest line tension
- the largest partitioned−monolithic heave difference for t ≥ 5 s, relative to
  the heave amplitude

The Newton counts come from the `beamSolverIterations` function object.

## Results (30 Oct 2026)

`summary.txt` has the full table. 19 of 28 runs reached 10 s:

| Case | Partitioned | Monolithic |
|---|---|---|
| Ex0.25 | finished, 12.1 Newton/step | finished, 12.8 |
| Ex1 (= T1) | crashed 2.61 s | finished, 12.6 |
| Ex4 | finished, 11.9 | crashed 6.75 s |
| Ex16 | crashed 6.29 s | crashed 3.28 s |
| Ex64 | crashed 1.97 s | crashed 6.09 s |
| Ex256 | finished, 11.5 | finished, 13.6 |
| Ex1024 | finished, 11.2 | finished, 13.5 |
| T0.5 | crashed 6.61 s | crashed 3.74 s |
| T0.6 | finished, 12.2 | finished, 14.0 |
| T0.8 | finished, 12.1 | finished, 13.5 |
| T1.25 | finished, 11.9 | finished, 12.5 |
| T1.5 | finished, 11.9 | finished, 12.1 |
| animation (base) | finished, 12.0 | finished, 12.6 |

Every crash is the beam's Newton solve failing, either hitting 200
iterations or diverging. 8 of the 9 happen while the line is touching the
tank floor:
- **Contact chatter:** the number of line cells in contact flips between
  two values (for example 3, 4, 3) from one Newton iteration to the next.
- **Line on the floor:** most of the line is lying on the floor. This
  happens with the stiffer lines (Ex16, Ex64 monolithic) and in the short
  waves (T0.5). The line has no pretension, so it buckles when the box
  moves back towards the anchor.

Only Ex64 partitioned fails with no floor contact.

Whether a borderline case survives depends on the exact time steps:
- Ex1 partitioned crashes here but finishes in `animation/partitioned` and
  in `../floatRigidBody`. The only difference is the write interval, and
  `adjustableRunTime` shortens steps to land on write times.
- Ex1 and T1 are the same case and crash identically, so the runs are
  deterministic.

In the finished runs, monolithic takes 2–20% more Newton iterations per step
than partitioned. The difference is largest for the stiffest lines and the
0.6 s waves. Run times are within 7% of each other because the
CFD dominates. Where both couplings finish, the heave differs by 0.2–1.5% of
its amplitude, except at T0.6 (17%, the largest pitch response).
