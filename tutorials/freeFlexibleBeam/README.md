# Tutorial case 3: Flexible beam with free ends - Kayak-rowing motion 

## Overview
This benchmark case reproduces the well-known **kayak-rowing pattern** of a flexible beam, investigated in several works
[Simo et al. (1988)], [Jelenić & Crisfield (1998)], and [Boyer et al. (2004)].  

The kayak-rowing motion arises from the **strongly coupled 3-D response** of a beam where the free ends trace paddle-like
trajectories due to nonlinear interaction of bending in two planes with oscillatory forward motion.

## Geometry and Setup
- Beam length: `L = 10 m`  
- Initial configuration: inclined beam in the xy-plane (see Figure: undeformed view)  
- Discretisation: `10` beam CVs  
- Time step: `Δt = 0.01 s`  
- Time integration: Newmark-β scheme  

### Cross-section and Material
- Circular cross-section radius: `R = 0.447214 m`  
- Density: `ρ = 1.59155 kg/m³`  
- Mechanical properties (from Simo et al., 1988):  
  - `EA = GA = 1 × 10⁴ N`  
  - `EI = 500 Nm`  
  - `Aρ = 1 kg/m`  
  - `I₂ρ = I₃ρ = 10 kg·m`  
- Artificial rotary inertia scaling factor: `k̃ = 200.39427` (to match reference values).

## Boundary Conditions
- Both ends are free (Neumann conditions).  
- At endpoint A (right patch):  
  - Force: `q(t) = (8, 0, 0) N` for `t ∈ [0, 2.5 s]`  
  - Moments: `m(t) = (mₓ(t), mᵧ(t), 0)` with `|m| = 80 Nm`  
- After `t = 2.5 s`, all external loads are removed, and the beam enters **free vibration** until `t = 40 s`.

## Results

### Deformation Pattern
- The free end exhibits a **paddle-like trajectory** (kayak-rowing motion).  
- Endpoint A oscillates in `y` and `z` while moving forward in `x`.  
- Residuals converge in ~3 iterations per timestep.  

### Figures
- ![Beam deformation pattern](figures/fig_case_flexBeam.eps)  
- ![Displacement vs time](figures/dispFlexibleBeamPlot.eps)  

### Energy Conservation
- Energy is conserved in free vibration (t > 2.5 s).  
- Comparison:  
  - Euler scheme → **numerical damping** (energy decay).  
  - Newmark-β scheme → **energy conservation** for Δt = 0.01 s and Δt = 0.1 s.  

- ![Energy stability](figures/case3_energy.png)

### Long-Time Simulation
- Simulated up to `t = 1000 s`.  
- Energy remains stable (small oscillations for Δt = 0.01 s).  

## How to Run
```bash
# 1. Create initial beam
setInitialPositionBeam

# 2. Run simulation
beamFoam
```

Energy can be extracted using the `beamEnergyData` functionObject in `controlDict`.

## References
- Simo, J.C., Vu-Quoc, L. (1988). *On the dynamics in space of rods undergoing large motions — A geometrically exact approach*. International Journal for Numerical Methods in Engineering.  
- Jelenić, G., Crisfield, M.A. (1998). *Geometrically exact 3D beam theory: implementation of a strain-invariant finite element for arbitrary nonlinear behaviour*. Computer Methods in Applied Mechanics and Engineering.  
- Boyer, F., et al. (2004). *Geometrically exact beam model for slender elastic rods: Numerical methods and applications*. Computers & Structures.
