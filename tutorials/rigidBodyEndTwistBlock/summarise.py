#!/usr/bin/env python3
"""
Summarise the rigidBodyEndTwist sweeps (run ./Allrun first).

For each case and coupling: whether the run finished, the time reached, the
beam Newton iterations, the partitioned coupling iterations, the run time,
the final twist against the static twist M L/(G Ip) (exact for a circular
rod at any twist), the final beam end moment against M, the largest
sideways tip displacement (non-zero only if the rod buckles: Greenhill's
limit is about pi EI/L of torque) and the largest difference in twist
between the partitioned and monolithic runs.

Uses only the Python standard library. Writes summary.txt as well.

Usage: ./summarise.py [sweep ...]
"""

import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SWEEPS = ["moment", "stiffness"]
COUPLINGS = ["partitioned", "monolithic"]


def read_text(path):
    with open(path) as f:
        return re.sub(r"//.*", "", f.read())


def scalar_entry(text, key):
    match = re.search(r"^\s*" + key + r"\s+([-+0-9.eE]+)\s*;", text, flags=re.M)
    if not match:
        raise KeyError(f"{key} not found")
    return float(match.group(1))


def vector_entry(text, key):
    match = re.search(r"^\s*" + key + r"\s+\(([^)]*)\)\s*;", text, flags=re.M)
    if not match:
        raise KeyError(f"{key} not found")
    return [float(v) for v in match.group(1).split()]


def read_log(case):
    info = {"status": "not run", "time": 0.0, "steps": 0, "runtime": math.nan,
            "newton": 0, "failed": 0}
    path = os.path.join(case, "log.beamFoam")
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
            elif "FOAM FATAL" in line or line.startswith("#0 ") or "sigFpe" in line:
                info["status"] = "crashed"
    return info


def parameters(case):
    beam = read_text(os.path.join(case, "constant", "beamProperties"))
    body = beam[beam.index("rigidBodyEnd\n"):]
    L = scalar_entry(beam, "length")
    R = scalar_entry(beam, "radius")
    E = scalar_entry(beam, "E")
    G = scalar_entry(beam, "G")
    Ip = math.pi*R**4/2
    J = vector_entry(body, "momentOfInertia")[0]
    dt = scalar_entry(read_text(os.path.join(case, "system", "controlDict")), "deltaT")
    beta = scalar_entry(body, "newmarkBeta")
    kt = G*Ip/L
    return {
        "M": vector_entry(body, "externalMoment")[0],
        "kt": kt,
        "omega": math.sqrt(kt/J),
        "greenhill": math.pi*E*(Ip/2)/L,
        # Partitioned coupling iterations converge only while this is < 1
        "stability": beta*dt**2*kt/J,
        "nCoupling": int(scalar_entry(body, "nCouplingIterations")),
    }


def read_body(case):
    """List of (time, twist, tip y, tip z, beam moment x, coupling iterations),
    the twist unwrapped from the rotation vector"""
    path = os.path.join(case, "postProcessing", "rigidBodyEnd", "0", "rigidBodyEnd.dat")
    rows = []
    if not os.path.isfile(path):
        return rows
    twist = 0.0
    previous = 0.0
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            v = line.split()
            if len(v) < 29:
                continue
            try:
                angle = float(v[17])
                twist += math.remainder(angle - previous, 2*math.pi)
                previous = angle
                rows.append((float(v[0]), twist, float(v[15]), float(v[16]),
                             float(v[26]), int(v[13])))
            except ValueError:
                pass
    return rows


def summarise_case(case):
    p = parameters(case)
    log = read_log(case)
    body = read_body(case)
    s = {"p": p, "log": log, "body": body}

    iterations = [row[5] for row in body]
    s["coupling_mean"] = sum(iterations)/len(iterations) if iterations else math.nan
    s["coupling_max"] = max(iterations, default=0)
    s["coupling_capped"] = sum(1 for n in iterations if n >= p["nCoupling"])
    s["sideways"] = max((math.hypot(r[2], r[3]) for r in body), default=math.nan)

    s["twist"] = s["twist_error"] = s["moment_error"] = math.nan
    if body and log["status"] == "finished":
        last = body[-1]
        static = p["M"]/p["kt"]
        s["twist"] = last[1]
        s["twist_error"] = last[1]/static - 1
        # The beam acts on the body with -M at rest
        s["moment_error"] = -last[4]/p["M"] - 1
    return s


def twist_difference(a, b):
    """Largest twist difference over the common times"""
    tb = {r[0]: r[1] for r in b}
    diffs = [abs(r[1] - tb[r[0]]) for r in a if r[0] in tb]
    return max(diffs) if diffs else math.nan


def fmt(value, spec):
    return "-" if isinstance(value, float) and math.isnan(value) else format(value, spec)


def summarise(sweep):
    root = os.path.join(HERE, sweep)
    if not os.path.isdir(root):
        return [f"{sweep}: no cases (run ./Allrun {sweep})"]

    labels = sorted(os.listdir(root), key=lambda s: float(re.sub(r"^[A-Za-z]+", "", s)))
    lines = [
        f"Sweep: {sweep}",
        f"{'case':<10} {'M':>8} {'kt':>9} {'omega':>7} {'M/Mgh':>6} {'stab':>7} "
        f"{'coupling':<12} {'status':<9} {'t':>6} {'Newton/step':>11} {'failed':>6} "
        f"{'coupl mean/max/capped':>21} {'runtime':>8} {'twist':>8} {'twist err':>9} "
        f"{'moment err':>10} {'sideways':>9} {'part-mono':>9}",
    ]
    for label in labels:
        results = {}
        for coupling in COUPLINGS:
            case = os.path.join(root, label, coupling)
            if os.path.isdir(case):
                results[coupling] = summarise_case(case)
        difference = math.nan
        if all(
            c in results and results[c]["log"]["status"] == "finished"
            for c in COUPLINGS
        ):
            difference = twist_difference(
                results["partitioned"]["body"], results["monolithic"]["body"]
            )

        for coupling, s in results.items():
            p, log = s["p"], s["log"]
            per_step = log["newton"]/log["steps"] if log["steps"] else math.nan
            coupling_text = (
                f"{s['coupling_mean']:.2f}/{s['coupling_max']}/{s['coupling_capped']}"
                if coupling == "partitioned" else ""
            )
            lines.append(
                f"{label:<10} {p['M']:8.4g} {p['kt']:9.4g} {p['omega']:7.1f} "
                f"{p['M']/p['greenhill']:6.2f} {p['stability']:7.2g} {coupling:<12} "
                f"{log['status']:<9} {log['time']:6.3f} {fmt(per_step, '11.2f')} "
                f"{log['failed']:6d} {coupling_text:>21} {fmt(log['runtime'], '8.2f')} "
                f"{fmt(s['twist'], '8.4f')} {fmt(s['twist_error'], '+9.2e')} "
                f"{fmt(s['moment_error'], '+10.2e')} {fmt(s['sideways'], '9.2e')} "
                f"{fmt(difference, '9.1e') if coupling == 'monolithic' else '':>9}"
            )
    return lines


def main():
    sweeps = sys.argv[1:] or SWEEPS
    output = []
    for sweep in sweeps:
        output += summarise(sweep) + [""]

    output += [
        "kt: torsional stiffness G Ip/L [N m/rad]; omega: torsional mode",
        "sqrt(kt/J) [rad/s]; M/Mgh: moment over the Greenhill buckling torque",
        "pi EI/L; stab: beta dt^2 kt/J, partitioned coupling iterations converge",
        "only below about 1; Newton/step: beam Newton iterations per time step;",
        "failed: steps where Newton did not converge; coupl mean/max/capped:",
        "partitioned coupling iterations per step, and steps that hit",
        "nCouplingIterations unconverged; twist: final twist [rad]; twist err:",
        "against M L/(G Ip); moment err: final beam end moment against M;",
        "sideways: largest sideways tip displacement [m] (buckling); part-mono:",
        "largest partitioned-monolithic twist difference [rad] (both finished only).",
    ]
    text = "\n".join(output)
    print(text)
    with open(os.path.join(HERE, "summary.txt"), "w") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
