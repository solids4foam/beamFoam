# Beam Attachment Kinematics Ownership

Date: 2026-05-07

This note records the current ownership decision for the `rigidBodyAndBeam`
coupling work.

## Current Working Mode

The current stable mode is:

```cpp
blockEigenForceCoupling true;
blockEigenKinematicCoupling false;
blockEigenAuthoritativeMotion true;
```

In this mode, ownership is split as follows:

```text
finiteVolumeBeam restraint
    owns the beam attachment boundary displacement before beam.evolve()

BlockEigenSolverOF
    owns the rigid-body update after beam.evolve()

sixDoFRigidBodyMotionFvBeam
    owns storage of the authoritative rigid-body state
```

The beam attachment point is still imposed on the beam boundary from the current
rigid-body transform before the beam solve. The beam reaction is then included in
the BlockEigen rigid-body RHS using `blockEigenForceCoupling`. After the solve,
the BlockEigen rigid-body state is applied to `sixDoFRigidBodyMotionFvBeam`.

When this authoritative path is active, `finiteVolumeBeam` returns zero force and
zero moment to the old restraint accumulator. This avoids double-counting the
beam force.

## Decision For Future Kinematic Coupling

When `blockEigenKinematicCoupling` is developed further, the beam attachment
kinematics should be owned by the block system, not by both the boundary
condition and the matrix coupling.

The intended future ownership is:

```text
BlockEigen matrix
    owns compatibility between beam-end DOFs and rigid-body DOFs

finiteVolumeBeam restraint
    supplies rigid-body state data and beam reaction data
    does not independently prescribe the same beam-end motion

sixDoFRigidBodyMotionFvBeam
    remains the storage owner for the accepted rigid-body state
```

The key rule is:

```text
Do not prescribe the beam-end displacement through W.boundaryFieldRef()
and simultaneously enforce the same attachment kinematics through
blockEigenKinematicCoupling.
```

That double ownership is the likely source of the earlier instability when
kinematic coupling was enabled.

## Proposed Implementation Direction

Add an explicit attachment kinematics mode, for example:

```cpp
beamAttachmentKinematics prescribedByRigidBody;
beamAttachmentKinematics solvedByBlockEigen;
```

For the current mode:

```text
prescribedByRigidBody
```

`finiteVolumeBeam` keeps setting the beam attachment boundary from the rigid-body
transform.

For the future fully coupled mode:

```text
solvedByBlockEigen
```

`finiteVolumeBeam` should stop imposing the attachment displacement directly, and
`blockEigenKinematicCoupling` should provide the compatibility coefficients that
connect the beam-end rows to the rigid-body translation and rotation DOFs.

## Next Work Item

The next code experiment should be a small, switch-gated
`solvedByBlockEigen` branch in `finiteVolumeBeam`. It should disable the direct
beam-end displacement prescription only when `blockEigenKinematicCoupling` is
active and the attachment compatibility rows are populated consistently.
