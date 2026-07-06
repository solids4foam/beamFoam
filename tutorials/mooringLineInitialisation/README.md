# Mooring Line Initialisation — Catenary Settling with Seabed Contact

## Overview
This tutorial demonstrates the **static initialisation of a catenary mooring line**: a finite-volume beam is placed on the straight chord between the anchor and the fairlead, and then dynamically settles under gravity, buoyancy, seabed contact and hydrodynamic damping into its static catenary shape.

The resulting equilibrium configuration (hanging catenary + grounded segment resting on the seabed) is the standard starting point for moored floating body simulations with the `moorFV` coupling library, where the initialised beam is used as a mooring-line restraint on a 6-DoF rigid body [Taran et al. (2025)].

The case exercises two runtime-selectable **beam momentum contributions**:
- `groundContact` — penalty spring/damper seabed contact,
- `morisonDrag` — Morison drag through still water (pure damping here).

---

## Initialisation Procedure
The `Allrun` script performs the three stages sketched below:

1. **`createBeamMesh`** — creates a straight beam mesh of length `L = 1.45 m` along the x-axis (60 control volumes, from `constant/beamProperties`).
2. **`setInitialPositionBeam`** — rigidly rotates the beam from the x-axis onto the anchor–fairlead chord direction and translates its first end onto the anchor point on the seabed.
3. **`beamFoam`** — transient settling run from `t = 0` to `t = 2 s`: the line sags under its effective (submerged) weight, the touchdown section is caught by the seabed contact model, Morison drag damps the transients, and the fairlead end is simultaneously ramped from the chord end to its target fairlead position by a prescribed displacement series.

```
 z = 0   ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~   free surface
                                                    o  <-- (2) chord end
                                                   *o  <-- (3) fairlead end ramped
                                               _.-* |          to target over 0..2 s
                        (2) straight        _.-*    |
                            chord        _.-*      /
                            after     _.-*        |
              setInitialPositionBeam *         (3) line settles into a catenary:
                                _.-*            |   gravity + buoyancy, damped by
                             _.-*              /        morisonDrag contribution
                          _.-*               _/
                       _.-*               __/
                    _.-*             ___-*
                 _.-*          ____-*
 z = -0.5   ___o=========----**_______________________________  seabed
            anchor   grounded segment                      (groundContact)
           (clamped)  (~28 of 60 CVs in contact at t = 2 s)

           (1) createBeamMesh first builds the beam straight along x:
               o==========================o
               x = 0                  x = L = 1.45 m
```

- **Anchor** (left patch): `(-1.385, -0.423, -0.5)` — clamped on the seabed.
- **Chord direction:** unit vector `(0.904, 0.276, 0.326)`.
- **Fairlead ramp** (right patch): total displacement `(-0.024, -0.0773, -0.051) m` relative to the chord end, applied linearly over `t ∈ [0, 2 s]` via `constant/timeVsDisplacement`.

---

## Geometry and Line Properties
The line represents a 1:80 model-scale studless mooring chain, modelled as an equivalent cylinder [Taran et al. (2025)]:

- **Unstretched length:** `L = 1.45 m`
- **Cross-section:** circular, radius `R = 1.828 mm` (equivalent diameter `d = 3.656 mm`)
- **Discretisation:** `60` beam CVs
- **Young's modulus:** `E = 1.826 MPa`
- **Shear modulus:** `G = 0.913 MPa`
- **Density:** `ρ = 5782.1 kg/m³` (mass per unit length `≈ 0.0607 kg/m`)
- **Water density:** `ρ_f = 1000 kg/m³` (buoyancy / effective weight)
- **Water depth:** `0.5 m` (seabed at `z = -0.5 m`)
- **Gravity:** `g = (0, 0, -9.81) m/s²`

---

## Boundary Conditions
- **`left` patch (anchor):** displacement `W` fixed to zero (`fixedValue`), i.e. the anchor is pinned to the seabed.
- **`right` patch (fairlead):** `fixedDisplacement` with a `displacementSeries` interpolated from `constant/timeVsDisplacement`.
- **Rotation `Theta`:** `momentBeamRotationNR` on both ends with a zero moment series (`constant/timeVsMoment`), i.e. moment-free ends.

---

## Seabed Contact and Damping
Both models are configured in `constant/beamMomentumContributionProperties` (not in `beamProperties`):

- **`groundContact`:** beam cells penetrating below `groundZ = -0.5 m` receive a normal penalty force with spring stiffness `kNormal = 1e3` and damping `cNormal = 1.0` (tangential stiffness and friction disabled here). The penalty is deliberately soft, giving a stable penetration of a few millimetres.
- **`morisonDrag`:** drag coefficients `Cdn = 1.6` / `Cdt = 0.05` acting on the beam's own velocity through still water. With no dissipation the line would bounce indefinitely instead of settling; the Morison drag provides the hydrodynamic damping that brings it to rest as a catenary by `t = 2 s`.

---

## Numerical Setup
- **Time integration:** Euler (first-order, implicit); `Δt = 0.002 s`, `t_end = 2 s`
- **Linear solver:** Eigen (direct)
- **Newton–Raphson tolerances:** `absoluteTol = 1e-7`, `residualTol = 1e-5`

---

## Running the Case
```
./Allclean
./Allrun
```

---

## Post-processing
- Open the case in **ParaView** (`case.foam`).
- Use the **WarpByVector** filter with the `pointW` field to visualise the deformed line (the mesh itself stays straight along the chord; `pointW` carries the displacement).
- The catenary shape at `t = 2 s` can be compared with the quasi-static lumped-mass solution (e.g. MoorDyn) and the experimental static tensions reported in Taran et al. (2025).

---

## Expected Results
- The line settles into a **static catenary with a grounded segment**: about `28` of the `60` CVs are in seabed contact at `t = 2 s`, with a stable penetration of a few millimetres (soft penalty).
- Newton–Raphson converges in **~4 iterations per time step**; total run time is a few seconds on a single core.
- The settled state (`2/` time directory) can be copied into a `moorFV` CFD case as the mooring-line initial condition. When doing so, change the `right` patch of `W` from `fixedDisplacement` to a `fixedValue` condition (the fairlead is then driven by the coupled 6-DoF body motion instead of the prescribed series).

---

## References
- Taran, A., Bali, S., Tuković, Ž., Pakrashi, V., & Cardiff, P. (2025). *A finite volume Simo-Reissner beam method for moored floating body dynamics.* Applied Ocean Research, 165, 104845. https://doi.org/10.1016/j.apor.2025.104845
