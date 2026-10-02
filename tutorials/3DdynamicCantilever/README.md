# Test case 1 - 3D Dynamic Cantilever Beam

## Case description  
This benchmark reproduces the large-deformation 3-D cantilever column studied in Cardiff et al. (2025) using the `solids4Foam` toolbox and ABAQUS.  
- A vertical cantilever beam with **square cross-section (0.2 m × 0.2 m)** and **length 2 m** is subjected to a sudden constant tip load.  
- The applied load is  
\[
\mathbf{q} = (2000, 2000, 0)\ \text{N},
\]  
acting on the top face from t=0 to t=1 s.  
- The case tests the ability of `beamFoam` to capture large dynamic deformations and finite strains.

**Reference:**  
- Cardiff et al., *A Jacobian-Free Newton–Krylov Hyperelastic Solid Solver in OpenFOAM* (2025).  
- ABAQUS simulations with C3D8 elements.

---

## Geometry and setup  
- Beam length: **2 m**  
- Cross-section: **Square, side = 0.2 m**  
- Material properties:  
  - Young’s modulus: E = 15.293 MPa  
  - Poisson’s ratio: ν = 0.3  
  - Density: ρ = 1000 kg/m^3  
- Boundary condition: **Clamped at the base, loaded at free end**  
- Mesh: tested with **5, 10, 20 segments**  

---

## Solver settings  
- Time integration:  
  - **Newmark-β scheme (second-order, implicit)**  
  - Optional: **Backward Euler (first-order, implicit)**  
- Default timestep: Δt = 0.001 s  
- Simulation time: t_end = 1.0 s  
- Convergence: ~4 Newton–Raphson iterations per step  

---

## Running a case  
1. To run and individual case, execute:
```
./Allclean
./Allrun
```
2. To perform a mesh discretisation study for a particular time scheme (Euler or Newmark), execute:
```
./runSweep.sh
```
3. To perform a time discretisation study, change timeSchemes in system/fvSchemes to either Euler or Newmark), and execute:
```
./runTimeSweep.sh
```

## Plots in gnuplot
# Run the command
gnuplot allPlots.gnuplot - generates diplacementPlot.pdf and energyPlot.pdf
---

## Post-processing  
- Open case in **ParaView** (`case.foam`)  
- Use **WarpByVector** filter with `pointW` field to visualise deformation.  
- Displacement of the top patch can be compared with reference ABAQUS and `solids4Foam` results.  
- Time-step sensitivity can be investigated by modifying `deltaT` in `controlDict`.  

---

## Expected results  
- **Displacement vs. time** at the top face shows good agreement with ABAQUS and `solids4Foam`.  
- Even coarse meshes (5 cells) capture the response well; with 20 segments, results converge to the ABAQUS solution.  
- **Performance:** `beamFoam` runs **orders of magnitude faster** (1.26 s vs. ~12,000 s for `solids4Foam`).  
- Time-integration comparison:  
  - Euler → higher dissipation, accurate only with very small Δt.  
  - Newmark → captures oscillatory motion better, stable for larger Δt.  

---

## 7. References  
- Cardiff, P., et al. (2025). *A jacobian-free newton-krylov method for cell-centred finite volume solid mechanics, arXiv preprint arXiv:2502.17217*  
- ABAQUS Documentation, C3D8 solid elements.  
