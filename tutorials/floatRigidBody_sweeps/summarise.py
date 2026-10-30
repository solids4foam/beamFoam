#!/usr/bin/env python3
"""
Summarise the floatRigidBody sweeps (run ./Allrun first).

For each case and coupling: whether the run finished, the time reached, the
time steps, the beam Newton iterations per step (mean and max) and in all,
the run time, the box's surge, heave and pitch amplitudes once the waves
are up (half the range over t >= T_SETTLED), the largest line tension, and
the largest partitioned-monolithic heave difference over the same window,
relative to the monolithic heave amplitude.

Uses only the Python standard library. Writes summary.txt as well.

Usage: ./summarise.py [sweep ...]
"""

import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SWEEPS = ["stiffness", "period", "animation"]
COUPLINGS = ["partitioned", "monolithic"]

DEPTH = 0.15     # Water depth [m]
G_ACC = 9.81
T_SETTLED = 5.0  # Waves ramp up over 2 s and take ~2 s to cross to the box

NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
VECTOR = r"\(\s*(" + NUMBER + r")\s+(" + NUMBER + r")\s+(" + NUMBER + r")\s*\)"


def read_text(path):
    with open(path) as f:
        return re.sub(r"//.*", "", f.read())


def scalar_entry(text, key):
    """key value; or key key [dimensions] value;"""
    match = re.search(
        r"^\s*" + key + r"\s+(?:" + key + r"\s+\[[^\]]*\]\s+)?(" + NUMBER + r")\s*;",
        text, flags=re.M,
    )
    if not match:
        raise KeyError(f"{key} not found")
    return float(match.group(1))


def wavelength(period, depth=DEPTH):
    """Linear dispersion omega^2 = g k tanh(k h), by Newton's method"""
    omega = 2*math.pi/period
    k = omega**2/G_ACC
    for _ in range(50):
        f = G_ACC*k*math.tanh(k*depth) - omega**2
        df = G_ACC*(math.tanh(k*depth) + k*depth/math.cosh(k*depth)**2)
        k -= f/df
    return 2*math.pi/k


def parameters(case):
    beam = read_text(os.path.join(case, "constant", "beamone", "beamProperties"))
    E = scalar_entry(beam, "E")
    R = scalar_entry(beam, "R")
    L = scalar_entry(beam, "L")
    mass = scalar_entry(read_text(os.path.join(case, "constant", "dynamicMeshDict")), "mass")
    waves = read_text(os.path.join(case, "constant", "waveProperties"))
    period = scalar_entry(waves, "wavePeriod")
    k_axial = E*math.pi*R**2/L
    return {
        "E": E,
        "kAxial": k_axial,
        # Line's axial mode with the box on its end (no added mass)
        "omega": math.sqrt(k_axial/mass),
        "period": period,
        "wavelength": wavelength(period),
        "height": scalar_entry(waves, "waveHeight"),
    }


def read_log(case):
    info = {"status": "not run", "time": 0.0, "steps": 0, "runtime": math.nan}
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
            elif line.startswith("End"):
                info["status"] = "finished"
            elif "FOAM FATAL" in line or line.startswith("#0 ") or "sigFpe" in line:
                info["status"] = "crashed"
    return info


def read_rows(path, pattern):
    """Rows matching pattern, the last one per time"""
    rows = {}
    if not os.path.isfile(path):
        return []
    with open(path, errors="replace") as f:
        for line in f:
            match = pattern.match(line)
            if match:
                values = [float(v) for v in match.groups()]
                rows[values[0]] = values
    return [rows[t] for t in sorted(rows)]


def read_motion(case):
    """(time, surge x, heave z, pitch) of the box"""
    path = os.path.join(
        case, "postProcessing", "sixDoF_History", "0", "sixDoFRigidBodyStateFvBeam.dat"
    )
    pattern = re.compile(r"^\s*(" + NUMBER + r")\s+" + r"\s+".join([VECTOR]*5))
    # Columns: time, centre of rotation, centre of mass, rotation, ...
    return [(r[0], r[4], r[6], r[8]) for r in read_rows(path, pattern)]


def read_columns(path, ncols):
    pattern = re.compile(r"^\s*" + r"\s+".join([f"({NUMBER})"]*(ncols + 1)))
    return read_rows(path, pattern)


def amplitude(rows, column):
    values = [r[column] for r in rows if r[0] >= T_SETTLED]
    return (max(values) - min(values))/2 if values else math.nan


def summarise_case(case):
    p = parameters(case)
    log = read_log(case)
    motion = read_motion(case)
    post = os.path.join(case, "postProcessing", "0")
    s = {"p": p, "log": log, "motion": motion}

    # Time nPimple nBeamSolves nNewton totalNewton
    work = read_columns(os.path.join(post, "beamSolverIterations.dat"), 4)
    newton = [r[3] for r in work]
    s["newton_mean"] = sum(newton)/len(newton) if newton else math.nan
    s["newton_max"] = max(newton, default=math.nan)
    s["newton_total"] = work[-1][4] if work else math.nan

    for i, key in enumerate(["surge", "heave", "pitch"]):
        s[key] = amplitude(motion, i + 1)

    # Time anchorAxial anchorShear attachmentAxial attachmentShear
    forces = read_columns(os.path.join(post, "axialForcebeamone.dat"), 4)
    s["tension"] = max((max(r[1], r[3]) for r in forces), default=math.nan)
    return s


def heave_difference(a, b):
    """Largest heave difference over t >= T_SETTLED, matching the 0.01 s
    samples by time"""
    tb = {round(r[0], 2): r[2] for r in b}
    diffs = [
        abs(r[2] - tb[round(r[0], 2)])
        for r in a if r[0] >= T_SETTLED and round(r[0], 2) in tb
    ]
    return max(diffs) if diffs else math.nan


def fmt(value, spec):
    return "-" if isinstance(value, float) and math.isnan(value) else format(value, spec)


def label_value(label):
    try:
        return float(re.sub(r"^[A-Za-z]+", "", label))
    except ValueError:
        return 0.0


def summarise(sweep):
    root = os.path.join(HERE, sweep)
    if not os.path.isdir(root):
        return [f"{sweep}: no cases (run ./Allrun {sweep})"]

    if sweep == "animation":
        groups = {"base": root}
    else:
        groups = {
            label: os.path.join(root, label)
            for label in sorted(os.listdir(root), key=label_value)
        }

    lines = [
        f"Sweep: {sweep}",
        f"{'case':<8} {'E':>8} {'EA/L':>8} {'omega':>6} {'T':>5} {'lambda':>6} "
        f"{'coupling':<12} {'status':<9} {'t':>6} {'steps':>6} "
        f"{'Newton mean/max':>15} {'total':>7} {'runtime':>8} "
        f"{'surge':>8} {'heave':>8} {'pitch':>8} {'tension':>8} {'part-mono':>9}",
    ]
    for label, directory in groups.items():
        results = {}
        for coupling in COUPLINGS:
            case = os.path.join(directory, coupling)
            try:
                results[coupling] = summarise_case(case)
            except (FileNotFoundError, KeyError):
                pass
        difference = math.nan
        if all(c in results for c in COUPLINGS):
            mono = results["monolithic"]
            difference = heave_difference(
                results["partitioned"]["motion"], mono["motion"]
            )/mono["heave"]

        for coupling, s in results.items():
            p, log = s["p"], s["log"]
            lines.append(
                f"{label:<8} {p['E']:8.3g} {p['kAxial']:8.3g} {p['omega']:6.1f} "
                f"{p['period']:5.2f} {p['wavelength']:6.3f} {coupling:<12} "
                f"{log['status']:<9} {log['time']:6.3f} {log['steps']:6d} "
                f"{fmt(s['newton_mean'], '9.2f')}/{fmt(s['newton_max'], '<5.0f')} "
                f"{fmt(s['newton_total'], '7.0f')} {fmt(log['runtime'], '8.0f')} "
                f"{fmt(s['surge'], '8.2e')} {fmt(s['heave'], '8.2e')} "
                f"{fmt(s['pitch'], '8.2e')} {fmt(s['tension'], '8.2e')} "
                f"{fmt(difference, '9.1e') if coupling == 'monolithic' else '':>9}"
            )
    return lines


def main():
    sweeps = sys.argv[1:] or SWEEPS
    output = []
    for sweep in sweeps:
        output += summarise(sweep) + [""]

    output += [
        "E: line Young's modulus [Pa]; EA/L: line axial stiffness [N/m]; omega:",
        "the line's axial mode sqrt(EA/(L m)) on the box [rad/s]; T: wave period",
        f"[s]; lambda: wavelength in {DEPTH} m of water [m]; t: time reached [s];",
        "steps: time steps; Newton mean/max: beam Newton iterations per time step",
        "(all PIMPLE outer correctors); total: Newton iterations in the run;",
        "runtime: ExecutionTime [s]; surge, heave [m], pitch [rad]: half the",
        f"range for t >= {T_SETTLED:g} s; tension: largest axial line force [N];",
        "part-mono: largest partitioned-monolithic heave difference for",
        f"t >= {T_SETTLED:g} s over the monolithic heave amplitude.",
    ]
    text = "\n".join(output)
    print(text)
    with open(os.path.join(HERE, "summary.txt"), "w") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
