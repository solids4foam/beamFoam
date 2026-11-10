#!/usr/bin/env python3
"""
Compare the moorFV welded twist runs (run ./Allrun first).

The body is welded to the rod end and starts spinning about the rod axis at
Omega0, with no gravity and the fluid forces ignored, so the twist is a free
torsional oscillation:

    theta(t) = (Omega0/omega) sin(omega t),  omega^2 = (G Ip/L)/(J + rho Ip L/3)

For each coupling: whether the run finished, the time reached, the run time,
the beam Newton iterations per time step, the largest difference from theory
and the amplitude, relative to the theoretical amplitude; then the largest
difference between the two couplings.

Uses only the Python standard library. Writes comparison.txt as well.
"""

import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
COUPLINGS = ["partitioned", "monolithic"]

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?|nan|inf"
VECTOR = rf"\(\s*({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*\)"


def read_text(path):
    with open(path) as f:
        return re.sub(r"//.*", "", f.read())


def scalar_entry(text, key):
    return float(re.search(rf"^\s*{key}\s+({NUMBER})\s*;", text, re.M).group(1))


def vector_entry(text, key):
    return [float(v) for v in re.search(rf"^\s*{key}\s+{VECTOR}\s*;", text, re.M).groups()]


def theory():
    beam = read_text(os.path.join(HERE, "base", "constant", "beam", "beamProperties"))
    body = read_text(os.path.join(HERE, "base", "constant", "dynamicMeshDict"))
    L = scalar_entry(beam, "length")
    R = scalar_entry(beam, "radius")
    G = scalar_entry(beam, "G")
    rho = scalar_entry(beam, "rho")
    Ip = math.pi*R**4/2
    J = vector_entry(body, "momentOfInertia")[0]
    pi0 = vector_entry(body, "angularMomentum")[0]
    omega = math.sqrt(G*Ip/L/(J + rho*Ip*L/3))
    return pi0/J, omega


def read_log(case):
    info = {"status": "not run", "time": 0.0, "steps": 0, "runtime": math.nan,
            "newton": 0}
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
            elif re.match(r"^\s*\d+: (Converged|Failed)", line):
                info["newton"] += int(line.split(":")[0])
            elif line.startswith("End"):
                info["status"] = "finished"
            elif "FOAM FATAL" in line or line.startswith("#0 ") or "sigFpe" in line:
                info["status"] = "crashed"
    return info


def read_twist(case):
    """Time -> twist about the rod axis, unwrapped"""
    path = os.path.join(case, "postProcessing", "sixDoF_History", "0",
                        "sixDoFRigidBodyStateFvBeam.dat")
    rows = {}
    if not os.path.isfile(path):
        return rows
    pattern = re.compile(rf"^\s*({NUMBER})\s+{VECTOR}\s+{VECTOR}\s+{VECTOR}")
    twist = previous = 0.0
    with open(path) as f:
        for line in f:
            match = pattern.search(line)
            if not match:
                continue
            values = [float(v) for v in match.groups()]
            angle = values[7]
            twist += math.remainder(angle - previous, 2*math.pi)
            previous = angle
            # The last line written at a time is the accepted state
            rows[round(values[0], 9)] = twist
    return rows


def main():
    Omega0, omega = theory()
    amplitude = Omega0/omega

    lines = [
        f"Free torsional oscillation: Omega0 = {Omega0:.4g} rad/s, "
        f"omega = {omega:.4g} rad/s, amplitude {amplitude:.4g} rad",
        "",
        f"{'coupling':<12} {'status':<9} {'t':>6} {'steps':>6} {'runtime':>8} "
        f"{'Newton/step':>11} {'amplitude':>9} {'max |theta - theory|':>21}",
    ]
    twists = {}
    for coupling in COUPLINGS:
        case = os.path.join(HERE, coupling)
        log = read_log(case)
        twist = read_twist(case)
        twists[coupling] = twist
        peak = max((abs(v) for v in twist.values()), default=math.nan)
        error = max(
            (abs(v - amplitude*math.sin(omega*t)) for t, v in twist.items()),
            default=math.nan,
        )
        per_step = log["newton"]/log["steps"] if log["steps"] else math.nan
        lines.append(
            f"{coupling:<12} {log['status']:<9} {log['time']:6.3f} {log['steps']:6d} "
            f"{log['runtime']:8.1f} {per_step:11.2f} {peak/amplitude:9.4f} "
            f"{error/amplitude:21.2e}"
        )

    common = set(twists["partitioned"]) & set(twists["monolithic"])
    if common:
        difference = max(abs(twists["partitioned"][t] - twists["monolithic"][t]) for t in common)
        lines += ["", f"Largest partitioned - monolithic twist difference: "
                  f"{difference/amplitude:.2e} of the amplitude"]

    lines += [
        "",
        "amplitude and max |theta - theory| are relative to the theoretical",
        "amplitude Omega0/omega. The trapezoidal rule lengthens the period by",
        "about (omega dt)^2/12, which accounts for a phase error that grows",
        "with time.",
    ]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(HERE, "comparison.txt"), "w") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
