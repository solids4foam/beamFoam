# Case 4: In-Plane Bending of a Cantilever Beam Subjected to Concentrated Moment

## Overview
This benchmark test investigates the **in-plane pure flexural bending** of a cantilever beam subjected to a concentrated moment at its free end.  
The problem has been widely studied in the literature [Simo & Vu-Quoc (1986), Ibrahimbegovic (1995)] and serves to validate beam formulations for large in-plane rotations.

---

## Geometry and Setup
- **Beam length:** `L = 10 m`  
- **Initial configuration:** straight cantilever, clamped at the left end  
- **Discretisation:** `10 CVs` (tested also with 5, 20, 40 CVs for convergence)  
---
## Material Properties
- Axial rigidity: `EA = 1 × 10⁴ N`  
- Shear rigidity: `GA₂ = GA₃ = 5000 N`  
- Flexural rigidity: `EI₂ = EI₃ = 100 Nm²`  
- Torsional rigidity: `GJ = 100 Nm²`  
- Cross-section radius: `r = 0.2 m`  
- Young’s modulus: `E = 7.95 × 10⁴ Pa`  
- Poisson’s ratio: `ν = 0`  

---

## Boundary Conditions
- **Left end:** clamped (fixed displacements and rotations).  
- **Right end (free end):** subjected to a **concentrated moment** about the z-axis:  
  - Example: `M_z = 20π Nm` applied in 4 increments.  
- Analytical solution (Euler formula):  
  \[
  ψ_z = \frac{M_z L}{EI}, \quad w_x = L - \frac{L}{ψ_z/2} \sin(ψ_z/2) \cos(ψ_z/2), \quad w_y = \frac{L}{ψ_z/2} \left( \sin(ψ_z/2) \right)^2
  \]

---

## Numerical Setup
- Discretisation: `10 CVs` (refined up to 40 CVs for convergence).  
- Newton–Raphson solver with average of 6 iterations per increment.  
- Execution time: `< 0.5 s` for 10 CVs.  

---

## Results

### Deformation Pattern
- For `M_z = 20π Nm (ψ_z = 2π)`, the beam bends into a full circle.

---

## How to Run
1. Execute:
```
./Allclean
./Allrun
```

Post-process in ParaView (`WarpByVector` using `pointW`).

---

## References
- Simo, J.C., & Vu-Quoc, L. (1986). *A three-dimensional finite-strain rod model. Part II: Computational aspects.* Computer Methods in Applied Mechanics and Engineering, 58, 79–116.  
- Ibrahimbegovic, A. (1995). *Computational aspects of vector-like parametrization of three-dimensional finite rotations.* International Journal for Numerical Methods in Engineering, 38(21), 3653–3673.  
