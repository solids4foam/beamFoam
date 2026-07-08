# 3D Dynamic Cantilever Multibeam Density Test

This tutorial is a multibeam variant of `tutorials/3dDynamicCantilever`.

It creates three independent cantilever beams with identical geometry, stiffness,
time integration, and tip loading, but different densities:

- `beam_0`: `rho = 500`
- `beam_1`: `rho = 1000`
- `beam_2`: `rho = 2000`

The purpose of the case is to expose multibeam inertia handling. With the same
applied end load on all three beams, the lighter beam should respond with larger
accelerations and a faster transient response, while the heavier beam should be
more inertially damped.

## Run

```bash
./Allclean
./Allrun
```

## Notes

- `createBeamMesh` generates three separate beam cell-zones from
  `constant/beamProperties`.
- The right-end load from `constant/timeVsForce` is applied to all three beam
  tips.
- The three beams are independent but co-located by default; compare them using
  the patch-specific history outputs for `right_0`, `right_1`, and `right_2`.

## Post-processing

Open the case in ParaView and inspect `W`, or compare the patch history data
written by the three `beamDisplacements` function objects.
