# rigidBodyEndPendulum: a rigid body pinned to a rod end

## Overview

A compound pendulum: a rigid body (`rigidBodyEnd`) is **pinned** to the end of
a short, stiff rod and swings under gravity. The pin passes a force but no
moment, so the body turns freely about the rod end, and the force acts on the
body through the arm from its centre of mass to the pin. For small angles the
period is

```text
T = 2 pi sqrt((J + m d^2)/(m g d))
```

with d the distance from the pin to the centre of mass. The case checks the
lever-arm terms of the coupling (offset centre of mass), gravity on the body,
and the pin joint.

```sh
./Allrun
./Allclean
```

---

## Setup

x points down. The rod hangs along +x from its clamp, and the body hangs from
the rod end, with its centre of mass further down. It swings in the x-y plane.

| Quantity | Value |
| --- | --- |
| Rod | length 0.5 m along +x, radius 0.05 m, 20 cells |
| Rod material | E = 5e8 Pa, G = 1.923e8 Pa, rho = 1 kg/m³ |
| Body | m = 10 kg, J = 10 kg m², centre of mass (1.5 0 0): d = 1 m |
| Gravity | g = (9.81 0 0) m/s² |
| Initial motion | angular velocity 0.1 rad/s about z, amplitude 0.045 rad |
| Small-angle period T | 2.8370 s |
| Time | dt = 1e-3 s to 9 s, trapezoidal Newmark (beta 0.25, gamma 0.5) |
| Coupling | `nCouplingIterations 20`, `couplingTolerance 1e-10` |

The rod end needs a zero moment Theta condition (`momentBeamRotationNR` with
`moment uniform (0 0 0)`), so that it is free to rotate. A pinned body with a
condition that fixes the rotation stops with an error. The initial velocity
of the centre of mass, `angularVelocity ^ (1 0 0)`, keeps the pin at rest.

The partitioned stability number beta dt² k/m at the pin is 0.2 along the rod
(k = E A/L) and 0.003 across it (k = 3 E I/L³, m = 1/(1/m + d²/J)), so the
coupling converges in about 10 iterations per time step.

---

## Results

`rotationSumZ` in `postProcessing/rigidBodyEnd/0/rigidBodyEnd.dat` is the
swing angle. Its zero crossings over 3 periods give T = 2.83796 s, 3.4e-4
above the small-angle value. The amplitude (1.3e-4) and the rod's small
flexibility (the pin moves by 4e-5 m) account for the difference.

Restarted at t = 4.5 s and run on to 9 s, the body motion agrees with the
uninterrupted run to the precision of the written files (1e-9 rad).
