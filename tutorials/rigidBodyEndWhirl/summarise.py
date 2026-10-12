#!/usr/bin/env python3
"""
Summarise the rigidBodyEndWhirl sweeps (run ./Allrun first).

For each case and coupling: whether the run finished, the time reached, the
beam Newton iterations, the partitioned coupling iterations, the run time,
the final spin rate relative to the critical speed, the error in the whirl
radius against linear theory, and the largest difference between the
partitioned and monolithic tip paths.

Whirl theory (small deflection, rod mass lumped at the tip): a body of mass m
whose centre of mass is a from the tip, spinning at omega about the rod axis,
pulls the tip out to radius

    a r^2/|1 - r^2|,  r = omega/omega_c,  omega_c^2 = (3EI/L^3)/(m + 33/140 m_rod)

The check is made at the final time and skipped near resonance
(|1 - r^2| < 0.3), where the deflection is not small.

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


def parameters(case):
    beam = read_text(os.path.join(case, "constant", "beamProperties"))
    body = beam[beam.index("rigidBodyEnd\n"):]
    L = scalar_entry(beam, "length")
    R = scalar_entry(beam, "radius")
    E = scalar_entry(beam, "E")
    I = math.pi*R**4/4
    A = math.pi*R**2
    m = scalar_entry(body, "mass")
    rod_mass = scalar_entry(beam, "rho")*A*L
    com = vector_entry(body, "centreOfMass")
    # Mass the pin sees: the body is pinned off its centre of mass, so a
    # force at the pin also turns it (smallest moment of inertia: worst case)
    J = min(vector_entry(body, "momentOfInertia"))
    a = math.hypot(com[1], com[2])
    pin_mass = 1/(1/m + a**2/J)
    dt = scalar_entry(read_text(os.path.join(case, "system", "controlDict")), "deltaT")
    beta = scalar_entry(body, "newmarkBeta")
    return {
        "EI": E*I,
        "M": vector_entry(body, "externalMoment")[0],
        "c": scalar_entry(body, "angularDamping"),
        "a": a,
        "omega_c": math.sqrt(3*E*I/L**3/(m + 33.0/140.0*rod_mass)),
        # Partitioned coupling iterations converge only while this is < 1
        "stability": beta*dt**2*E*A/L/pin_mass,
        "nCoupling": int(scalar_entry(body, "nCouplingIterations")),
    }


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


def read_body(case):
    """Time -> (tip y, tip z, omega x, coupling iterations)"""
    path = os.path.join(case, "postProcessing", "rigidBodyEnd", "0", "rigidBodyEnd.dat")
    rows = {}
    if not os.path.isfile(path):
        return rows
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            v = line.split()
            if len(v) < 26:
                continue
            try:
                rows[v[0]] = (float(v[15]), float(v[16]), float(v[20]), int(v[13]))
            except ValueError:
                pass
    return rows


def summarise_case(case):
    p = parameters(case)
    log = read_log(case)
    body = read_body(case)
    s = {"p": p, "log": log, "body": body}

    iterations = [row[3] for row in body.values()]
    s["coupling_mean"] = sum(iterations)/len(iterations) if iterations else math.nan
    s["coupling_max"] = max(iterations, default=0)
    s["coupling_capped"] = sum(1 for n in iterations if n >= p["nCoupling"])

    s["ratio"] = s["whirl_error"] = math.nan
    if body and log["status"] == "finished":
        last = body[max(body, key=float)]
        r = last[2]/p["omega_c"]
        s["ratio"] = r
        if abs(1 - r**2) >= 0.3:
            predicted = p["a"]*r**2/abs(1 - r**2)
            s["whirl_error"] = math.hypot(last[0], last[1])/predicted - 1
    return s


def path_difference(a, b):
    """Largest tip separation over the common times, relative to the
    largest tip radius"""
    times = set(a) & set(b)
    if not times:
        return math.nan
    diff = max(math.hypot(a[t][0] - b[t][0], a[t][1] - b[t][1]) for t in times)
    radius = max(math.hypot(a[t][0], a[t][1]) for t in times)
    return diff/radius if radius > 0 else math.nan


def fmt(value, spec):
    return "-" if isinstance(value, float) and math.isnan(value) else format(value, spec)


def summarise(sweep):
    root = os.path.join(HERE, sweep)
    if not os.path.isdir(root):
        return [f"{sweep}: no cases (run ./Allrun {sweep})"]

    labels = sorted(os.listdir(root), key=lambda s: float(re.sub(r"^[A-Za-z]+x?", "", s)))
    lines = [
        f"Sweep: {sweep}",
        f"{'case':<10} {'M':>6} {'EI':>8} {'wc':>6} {'stab':>6} "
        f"{'coupling':<12} {'status':<9} {'t':>6} {'Newton/step':>11} {'failed':>6} "
        f"{'coupl mean/max/capped':>21} {'runtime':>8} {'w/wc':>6} {'whirl err':>9} "
        f"{'part-mono':>9}",
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
            difference = path_difference(
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
                f"{label:<10} {p['M']:6.3g} {p['EI']:8.4g} {p['omega_c']:6.2f} "
                f"{p['stability']:6.3f} {coupling:<12} {log['status']:<9} "
                f"{log['time']:6.3f} {fmt(per_step, '11.2f')} {log['failed']:6d} "
                f"{coupling_text:>21} {fmt(log['runtime'], '8.1f')} "
                f"{fmt(s['ratio'], '6.3f')} {fmt(s['whirl_error'], '+9.2%')} "
                f"{fmt(difference, '9.1e') if coupling == 'monolithic' else '':>9}"
            )
    return lines


sweeps = sys.argv[1:] or SWEEPS
output = []
for sweep in sweeps:
    output += summarise(sweep) + [""]

output += [
    "wc: critical speed [rad/s]; stab: beta dt^2 EA/(L m_pin), m_pin =",
    "1/(1/m + a^2/J): partitioned coupling iterations converge only below about",
    "1; Newton/step: beam Newton iterations per time step; failed: steps",
    "where Newton did not converge;",
    "coupl mean/max/capped: partitioned coupling iterations per step, and steps",
    "that hit nCouplingIterations unconverged; w/wc: final spin rate over wc;",
    "whirl err: final tip radius against linear theory (- near resonance or",
    "not finished);",
    "part-mono: largest partitioned-monolithic tip separation over the largest",
    "tip radius (both finished only).",
]
text = "\n".join(output)
print(text)
with open(os.path.join(HERE, "summary.txt"), "w") as f:
    f.write(text + "\n")
