# Case Name: In-plane cantilever subjected to end follower load

## Case Description
This benchmark is an in-plane quasi-static problem originally studied by Argyris et al. (1981) and later by Simo and Vu-Quoc (1986).
A straight cantilever beam of length **L = 100 m**, clamped at the left end, is subjected at the free end to a **circulatory follower load** (load that follows instananeous deformation) of magnitude:

\[
\| \mathbf{q} \| = 134 \ \text{kN}
\]

## Beam and Material Properties
- Flexural rigidity: \( EI = 3.5 \times 10^{7} \ \text{N/m}^2 \)
- Shear rigidity: \( GA = 1.61538 \times 10^{8} \ \text{N/m}^2 \)
- Cross-section: circular, radius \( R \approx 0.658 \ \text{m} \)

The load of 134 kN is applied in the **z-direction** in 1000 increments (pseudo time step \( \Delta t = 0.001 \ \text{s} \)).

## Boundary Conditions
- **Displacement field `W`:**
  - Left end (clamped): `fixedValue uniform (0 0 0)`
  - Right end: `followerForceBeamDisplacementNR`
- **Rotation field `Theta`:**
  - Left end: clamped (`fixedValue uniform (0 0 0)`)
  - Right end: `momentBeamRotationNR` with time-varying zero `momentSeries`

Example boundary condition setup is shown in `0/W` and `0/Theta`.

## Solver Settings
- Newton–Raphson convergence tolerance for both `W` and `Theta` fields set to:
  \[
  1 \times 10^{-10}
  \]
  (default code tolerance is \( 1 \times 10^{-6} \))
- Typical convergence is achieved in **4–6 iterations per load step**.

## Running the Case
1. Ensure the beam mesh and case files are set up.
2. Run the solver with:
   ```bash
   beamFoam
   ```
3. Post-process results using ParaView:
   - Use the field `pointW` with the **Warp By Vector** filter to visualise the deformed configurations.
   - Displacement fields (`W`) can also be visualised with colour maps.

## Results
- **Displacement response:**
  Figure `cant-fol-force-disp` shows \( w_x \) and \( w_y \) displacements at the free end of the beam for meshes with 5, 10, and 20 cells. Results are in good agreement with Simo & Vu-Quoc (1986).
- **Mesh convergence:**
  A quadratic rate of convergence is observed when refining from 10 → 80 CVs (Figure `cant-fol-force-mesh-conv`).
- Residual convergence:
  \[
  \| \phi \| = \max (\| \Delta \mathbf{w} \|, \| \Delta \mathbf{R}_V \|) < 1 \times 10^{-10}
  \]

---

## References
- Argyris, J. et al. (1981). *Nonlinear finite element analysis of elastic systems under nonconservative loading-natural formulation. part i. quasistatic problems, Computer Methods in Applied Mechanics and  Engineering*.
- Simo, J. C., & Vu-Quoc, L. (1986). *A three-dimensional finite-strain rod model. Part II: Computational aspects.* Computer Methods in Applied Mechanics and Engineering, 58, 79–116.
