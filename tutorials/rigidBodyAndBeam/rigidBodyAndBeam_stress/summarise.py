#!/usr/bin/env python3
"""
Summarise the rigidBodyAndBeam stress tests (run ./Allrun first).

For each variant and coupling method: whether the run finished, the simulated
time reached, the largest difference from the converged reference
(monolithic_nOuter8) in body displacement and beam force, the run time, the
number of time steps and the number of beam Newton iterations.

Uses only the Python standard library.

Usage: ./summarise.py [variant ...]
"""

import glob
import math
import os
import re
import bisect
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
METHODS = ["loop_nOuter1", "loop_nOuter8", "monolithic_nOuter1", "monolithic_nOuter8"]
REFERENCES = ["monolithic_nOuter8", "loop_nOuter8"]

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?|nan|inf"
VECTOR = rf"\(\s*({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*\)"


def read_log(case):
    info = {"status": "not run", "time": 0.0, "steps": 0, "runtime": math.nan,
            "newton": 0, "failed": 0}
    path = os.path.join(case, "log.interFoam")
    if not os.path.isfile(path):
        return info

    info["status"] = "running"
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("Time = "):
                info["time"] = float(line.split()[2])
                info["steps"] += 1
            elif line.startswith("ExecutionTime"):
                info["runtime"] = float(line.split()[2])
            elif re.match(r"^\s*\d+: Converged", line):
                info["newton"] += int(line.split(":")[0])
            elif re.match(r"^\s*\d+: Failed", line):
                info["newton"] += int(line.split(":")[0])
                info["failed"] += 1
            elif line.startswith("End"):
                info["status"] = "finished"
            elif "FOAM FATAL" in line or line.startswith("#0 ") or "sigFpe : Floating" in line:
                info["status"] = "crashed"
    return info


def read_motion(case):
    path = os.path.join(case, "postProcessing", "sixDoF_History", "0",
                        "sixDoFRigidBodyStateFvBeam.dat")
    if not os.path.isfile(path):
        return None
    pattern = re.compile(rf"^\s*({NUMBER})\s+{VECTOR}")
    rows = {}
    with open(path) as f:
        for line in f:
            match = pattern.search(line)
            if match:
                values = [float(v) for v in match.groups()]
                rows[values[0]] = values[1:4]
    if not rows:
        return None
    time = sorted(rows)
    text = open(os.path.join(case, "constant", "dynamicMeshDict")).read()
    c0 = [float(v) for v in re.search(rf"centreOfMass\s+{VECTOR}", text).groups()]
    return time, [[a - b for a, b in zip(rows[t], c0)] for t in time]


def read_force(case):
    path = os.path.join(case, "postProcessing", "0", "forcebeam.dat")
    if not os.path.isfile(path):
        return None
    rows = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            values = [float(v) for v in line.split()[:4]]
            rows[values[0]] = values[1:4]
    if not rows:
        return None
    time = sorted(rows)
    return time, [rows[t] for t in time]


def norm(v):
    return math.sqrt(sum(x*x for x in v))


def interpolate(time, values, t):
    i = min(max(bisect.bisect_left(time, t), 1), len(time) - 1)
    w = (t - time[i - 1])/(time[i] - time[i - 1])
    return [(1 - w)*a + w*b for a, b in zip(values[i - 1], values[i])]


def max_difference(series, reference):
    """Largest |case - reference| on the common time range, and the reference peak"""
    if series is None or reference is None:
        return math.nan, math.nan
    time, values = series
    ref_time, ref_values = reference
    peak = max(norm(v) for v in ref_values)
    if not all(math.isfinite(x) for v in values for x in v):
        return math.inf, peak
    diffs = [
        norm([a - b for a, b in zip(v, interpolate(ref_time, ref_values, t))])
        for t, v in zip(time, values) if ref_time[0] <= t <= ref_time[-1]
    ]
    return (max(diffs) if diffs else math.nan), peak


variants = sys.argv[1:] or sorted(
    d for d in os.listdir(HERE) if os.path.isdir(os.path.join(HERE, d, METHODS[0]))
)

header = (f"{'method':<20} {'status':<9} {'t end':>7} {'steps':>6} {'run (s)':>8} "
          f"{'Newton':>7} {'fail':>4} {'disp diff':>10} {'force diff':>10}")

for variant in variants:
    cases = {m: os.path.join(HERE, variant, m) for m in METHODS}
    reference = next((m for m in REFERENCES if read_log(cases[m])["status"] == "finished"), None)
    ref_motion = read_motion(cases[reference]) if reference else None
    ref_force = read_force(cases[reference]) if reference else None

    print("")
    print(f"{variant}: differences relative to the peak of {reference or 'no finished reference'}")
    print(header)
    print("-"*len(header))
    for method, case in cases.items():
        info = read_log(case)
        if method == reference:
            disp = force = "reference"
        else:
            d, dp = max_difference(read_motion(case), ref_motion)
            f, fp = max_difference(read_force(case), ref_force)
            disp = f"{d/dp:10.2e}" if math.isfinite(d/dp) else f"{d/dp:>10}"
            force = f"{f/fp:10.2e}" if math.isfinite(f/fp) else f"{f/fp:>10}"
        print(f"{method:<20} {info['status']:<9} {info['time']:7.4f} {info['steps']:6d} "
              f"{info['runtime']:8.0f} {info['newton']:7d} {info['failed']:4d} "
              f"{disp:>10} {force:>10}")
