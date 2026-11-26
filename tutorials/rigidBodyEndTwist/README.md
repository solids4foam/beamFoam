# rigidBodyEndTwist: a rod twisted by a welded rigid body

## Overview

A rod is clamped at one end, and a rigid body (`rigidBodyEnd`) is **welded**
to the other: the rod's end face turns with the body, and the rod's end
moment acts back on the body. A step moment about the rod axis, applied to
the body, twists the rod. Angular damping on the body lets the motion settle
on the static twist, which for a circular rod is exact at any twist:

```text
theta = M L/(G Ip)
```

The body and the beam are solved in turn (partitioned coupling), with
coupling iterations in each time step until the body motion changes by less
than `couplingTolerance`.

```sh
./Allrun
./Allclean
```

---

## Setup

| Quantity | Value |
| --- | --- |
| Rod | length 2 m along +x, radius 0.05 m, 40 cells |
| Rod material | E = 5.4326e7 Pa, G = 2.0895e7 Pa, rho = 1 kg/m³ |
| Torsional stiffness kt = G Ip/L | 102.57 N m/rad |
| Body | mass 1 kg, J = 0.05 kg m², centre of mass at the rod end |
| Applied moment M | (102.57 0 0) N m: static twist 1 rad |
| Angular damping c | 1 N m s: decays at c/(2J) = 10 1/s |
| Time | dt = 1e-3 s to 1.5 s, trapezoidal Newmark (beta 0.25, gamma 0.5) |
| Coupling | `nCouplingIterations 20`, `couplingTolerance 1e-10` |

The left end is clamped (`fixedValue` W and Theta). The right end has
`fixedValue` W and Theta conditions, which the body sets in every coupling
iteration. The partitioned stability number beta dt² kt/J is 5.1e-4, so each
time step converges in 2 to 4 coupling iterations (see `rigidBodyEnd.H`).

The body is advanced with the beam's Newmark parameters, `newmarkBeta(W)` and
`newmarkGamma(W)` in `d2dt2Schemes` in `system/fvSchemes`.

---

## Output

`postProcessing/rigidBodyEnd/0/rigidBodyEnd.dat` has one line per time step:
the body's displacement, velocity and acceleration, the force from the beam,
the coupling iterations used and whether they converged, the orientation
`Q` (9 entries, body to global axes), `rotationSum` (the sum of the step
rotation vectors, which for this rotation about x is the twist angle, and
does not wrap at pi), the angular velocity, and the torque and moment from
the beam on the body.

The body state is also written to `<time>/uniform/rigidBodyEnd` at every
write time, and read back when the run is restarted from that time.

`postProcessing/0/beamConvergenceData.dat` has the solver work per time step:
`nBeamSolves` is the number of coupling iterations (beam solves), and
`nNewton` the Newton iterations summed over them.

---

## Results

At t = 1.5 s:

- the twist (`rotationSumX`) is 1.00002 rad, against 1 rad static;
- the moment from the beam on the body (`beamMomentX`) is -102.570 N m,
  balancing the applied moment;
- every time step converged, in 2.3 coupling iterations on average (4 at
  most).

Restarted at t = 0.76 s and run on to 1.5 s, the run agrees with the
uninterrupted one to the precision of the written files (1e-7 rad in the
twist).
