#!/usr/bin/env python3
"""
Check the beamFoamCoupled case after ./Allrun:

- the body state moorFV stores (sixDoF history) equals the state beamFoam
  solved for (rigidBodyEnd history), at every write time;
- for information, how far the trajectory is from rigidBodyAndBeam_loop.

Uses only the Python standard library. Exit status is non-zero if the state
check fails.
"""

import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LOOP = os.path.join(os.path.dirname(HERE), "rigidBodyAndBeam_loop")


def sixdof_history(case):
    """time -> (centre of rotation, velocity, omega) from the moorFV history"""
    path = os.path.join(case, "postProcessing", "sixDoF_History", "0", "sixDoFRigidBodyStateFvBeam.dat")
    rows = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = float(line.split()[0])
            vectors = [[float(x) for x in v.split()] for v in re.findall(r"\(([^)]*)\)", line)]
            rows[t] = (vectors[0], vectors[3], vectors[4])
    return rows


def rigid_body_end_history(case):
    """time -> row of the beamFoam rigidBodyEnd history (last row per time)"""
    path = sorted(glob.glob(os.path.join(case, "postProcessing", "rigidBodyEnd", "*", "rigidBodyEnd.dat")))[0]
    rows = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            values = [float(v) for v in line.split()]
            rows[values[0]] = values
    return rows


def initial_centre_of_mass(case):
    text = open(os.path.join(case, "constant", "dynamicMeshDict")).read()
    return [float(v) for v in re.search(r"centreOfMass\s+\(([^)]*)\)", text).group(1).split()]


def norm(v):
    return math.sqrt(sum(x*x for x in v))


def closest(rows, t):
    key = min(rows, key=lambda k: abs(k - t))
    return rows[key] if abs(key - t) < 1e-9 else None


moorfv = sixdof_history(HERE)
beamfoam = rigid_body_end_history(HERE)
c0 = initial_centre_of_mass(HERE)

worst = {"centre of rotation": 0.0, "velocity": 0.0, "angular velocity": 0.0}
scale = {"centre of rotation": 0.0, "velocity": 0.0, "angular velocity": 0.0}
compared = 0

for t, (cor, vel, omega) in sorted(moorfv.items()):
    row = closest(beamfoam, t)
    if row is None:
        continue
    compared += 1
    expected = {
        "centre of rotation": [c + d for c, d in zip(c0, row[1:4])],
        "velocity": row[4:7],
        "angular velocity": row[20:23],
    }
    actual = {"centre of rotation": cor, "velocity": vel, "angular velocity": omega}
    for key in worst:
        worst[key] = max(worst[key], norm([a - e for a, e in zip(actual[key], expected[key])]))
        scale[key] = max(scale[key], norm(expected[key]) if key != "centre of rotation"
                         else norm([e - c for e, c in zip(expected[key], c0)]))

print(f"moorFV body state against beamFoam's solved state, {compared} write times")
ok = compared > 0
for key in worst:
    rel = worst[key]/max(scale[key], 1e-30)
    passed = rel <= 1e-6
    ok = ok and passed
    print(f"  {'PASS' if passed else 'FAIL'}  {key:<20} max difference {worst[key]:.3e} "
          f"({rel:.2e} of peak motion, tol 1e-06)")

if os.path.isdir(LOOP):
    loop = sixdof_history(LOOP)
    diffs, peak = [], 0.0
    for t, (cor, _, _) in moorfv.items():
        if t in loop:
            diffs.append(norm([a - b for a, b in zip(cor, loop[t][0])]))
            peak = max(peak, norm([b - c for b, c in zip(loop[t][0], c0)]))
    if diffs:
        print(f"  info  against rigidBodyAndBeam_loop: max centre-of-rotation difference "
              f"{max(diffs):.3e} m ({max(diffs)/peak:.2e} of peak displacement)")

sys.exit(0 if ok else 1)
