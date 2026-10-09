#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare the moored floating box in waves (floatRigidBody) between three
mooring models:

- partitioned: floatRigidBody/floatRigidBody_partioned (moorFV solver
  FvBeamNewmark; the body and beamFoam take turns)
- monolithic: floatRigidBody/floatRigidBody_monolithic (moorFV solver
  beamFoamCoupled; beamFoam solves the line and the body together)
- MoorDyn: floatRigidBody/floatRigidBody_moorDyn (foamMooring with MoorDyn v2
  lumped-mass line; copied from run/foamMooring/floatRigidBody_moorDyn)

All three use 3 PIMPLE outer correctors, the same mesh, waves, box and
constraints (surge, heave and pitch only), with the still water at z = 0.
All three have seabed contact for the line: the beamFoam cases at groundZ
-0.15 m with kNormal 1e4 (constant/beamMomentumContributionProperties),
MoorDyn with kBot 3e6.

Three groups of plots:
- motion: surge, heave and pitch, their difference from REFERENCE_CASE, the
  slow drift and the wave-frequency oscillation
- mooring tension: line tension at the box (fairlead) and at the anchor,
  and, for the beamFoam cases, the magnitude of the end force |Q|
- efficiency: run time, time-step size, linear solver iterations and, for the
  beamFoam cases, beam Newton iterations, all from log.interFoam

Motion is taken relative to each case's initial centreOfMass.

Line tension:
- MoorDyn: FairTen1 and AnchTen1 in Mooring/lines.out (tension in the end
  segments; never negative, a slack line has zero tension)
- beamFoam: the axial force, i.e. the end-face force Q along the line
  tangent, positive in tension and negative in compression (the beam has
  bending stiffness, so a slack line can push). Read from
  postProcessing/0/axialForcebeamone.dat, written every time step by the
  finiteVolumeBeam restraint in moorFV (since 7 Oct 2026). Runs from before
  that have no such file: the axial force is then worked out from the beam
  fields in the time directories (Q, W and refW on the end faces), so only
  at the write interval (0.05 s) and the snap load is under-resolved.

The end-face force magnitude |Q| in postProcessing/0/{attachment,anchor}
Forcebeamone.dat (what the AssessmentOfCoupledFVFramework scripts plot as
tension) is also plotted. It equals the tension only while the line is taut
and straight; once the line is slack and bends, most of it is shear.

Run from Spyder or from the beamFoam/tutorials directory. Cases that have not
been run are skipped, and a case that is still running is plotted up to its
latest output.
"""

import os
import re

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MaxNLocator


# =============================================================
# CASES
# =============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# (case directory, legend label, mooring model ("beam" or "moorDyn"), colour,
# line style). Directories are relative to this script.
CASES = [
    ("floatRigidBody/floatRigidBody_partioned", "Partitioned (beamFoam)", "beam", "b", "-"),
    ("floatRigidBody/floatRigidBody_monolithic", "Monolithic (beamFoam)", "beam", "r", "--"),
    ("floatRigidBody/floatRigidBody_moorDyn", "MoorDyn", "moorDyn", "k", ":"),
]

# Differences are taken against this case (its label)
REFERENCE_CASE = "MoorDyn"

# Beam region of the mooring line in the beamFoam cases
BEAM_NAME = "beamone"


# =============================================================
# USER SETTINGS
# =============================================================

T_START = None
T_END = None

# Wave period (constant/waveProperties): the moving average over one period
# splits the motion into a slow drift and the wave-frequency oscillation
WAVE_PERIOD = 1.0

# Statistics are taken from this time on, after the 2 s wave ramp
STEADY_START = 6.0

# Window for the zoomed wave-frequency plots
ZOOM_WINDOW = (8.0, 10.0)

# The snap load at t = 2.3 s (about 0.1 N) hides the tension after it (about
# 0.001 N): the tension is also plotted from this time on
TENSION_ZOOM_START = 3.0

SAVE_FIGURES = False
FIGURE_DIR = os.path.join(SCRIPT_DIR, "comparison_plots", "floatRigidBody_moorDyn")


# =============================================================
# READ DATA
# =============================================================

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
VECTOR = rf"\(\s*({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*\)"

TIME_PATTERN = re.compile(rf"^\s*Time\s*=\s*({NUMBER})\s*$")
DELTA_T_PATTERN = re.compile(rf"^\s*deltaT\s*=\s*({NUMBER})")
EXECUTION_PATTERN = re.compile(
    rf"^\s*ExecutionTime\s*=\s*({NUMBER})\s*s\s+ClockTime\s*=\s*({NUMBER})\s*s"
)
NEWTON_PATTERN = re.compile(r"^\s*(\d+)\s*:\s*(?:Converged|Failed)\b")
LINEAR_SOLVER_PATTERN = re.compile(r"Solving for\s+([^,]+),.*No Iterations\s+(\d+)")


def read_initial_centre_of_mass(case_dir):
    """Initial centre of mass from constant/dynamicMeshDict"""
    with open(os.path.join(case_dir, "constant", "dynamicMeshDict")) as handle:
        text = handle.read()

    match = re.search(rf"centreOfMass\s+{VECTOR}", text)
    if not match:
        raise ValueError(f"No centreOfMass in {case_dir}/constant/dynamicMeshDict")

    return np.array([float(value) for value in match.groups()])


def last_row_per_time(data):
    """A time can appear more than once (outer correctors, restarts); keep
    the last row"""
    _, last = np.unique(data[::-1, 0], return_index=True)
    return data[::-1][last]


def read_sixdof_history(filename):
    """Time, centre of rotation, rotation, velocity and angular velocity"""
    line_pattern = re.compile(
        rf"^\s*({NUMBER})\s+" + r"\s+".join([VECTOR]*5)
    )

    rows = []
    with open(filename, "r", encoding="utf-8") as handle:
        for line in handle:
            match = line_pattern.search(line)
            if match:
                rows.append([float(value) for value in match.groups()])

    if not rows:
        raise ValueError(f"No motion rows found in {filename}")

    data = last_row_per_time(np.asarray(rows))

    return {
        "time": data[:, 0],
        "centre_of_rotation": data[:, 1:4],
        "rotation": data[:, 7:10],
    }


def read_columns(filename, ncols, skip_header=0):
    """Time and ncols values per row, last row per time"""
    data = np.genfromtxt(
        filename, comments="#", skip_header=skip_header, invalid_raise=False
    )

    if data.size == 0:
        raise ValueError(f"No rows found in {filename}")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    data = last_row_per_time(data[:, :ncols + 1])
    return data[:, 0], data[:, 1:]


def read_log(filename):
    """Per time step: time, deltaT, cumulative ExecutionTime and ClockTime,
    linear solver iterations and beam Newton iterations (summed over all
    beam solves in the step; NaN when there are none, e.g. MoorDyn)"""
    rows = []

    time = None
    delta_t = np.nan
    linear = 0
    newton = np.nan

    with open(filename, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = DELTA_T_PATTERN.search(line)
            if match:
                delta_t = float(match.group(1))
                continue

            match = TIME_PATTERN.search(line)
            if match:
                time = float(match.group(1))
                linear = 0
                newton = np.nan
                continue

            if time is None:
                continue

            match = NEWTON_PATTERN.search(line)
            if match:
                newton = np.nansum([newton, float(match.group(1))])
                continue

            match = LINEAR_SOLVER_PATTERN.search(line)
            if match:
                linear += int(match.group(2))
                continue

            match = EXECUTION_PATTERN.search(line)
            if match:
                rows.append([
                    time, delta_t,
                    float(match.group(1)), float(match.group(2)),
                    linear, newton,
                ])
                time = None

    if not rows:
        raise ValueError(f"No time steps found in {filename}")

    data = np.asarray(rows)
    return {
        "time": data[:, 0],
        "delta_t": data[:, 1],
        "execution_time": data[:, 2],
        "clock_time": data[:, 3],
        "step_execution_time": np.diff(np.concatenate(([0.0], data[:, 2]))),
        "linear_iterations": data[:, 4],
        "newton_iterations": data[:, 5],
    }


def apply_time_window(time, *arrays):
    mask = np.ones_like(time, dtype=bool)

    if T_START is not None:
        mask &= time >= T_START

    if T_END is not None:
        mask &= time <= T_END

    return (time[mask],) + tuple(array[mask] for array in arrays)


def read_field_values(filename, patches):
    """Internal vector values and the uniform value on each named patch of
    an OpenFOAM field file"""
    with open(filename, "r", encoding="utf-8") as handle:
        text = handle.read()

    start = text.index("internalField")
    boundary = text.index("boundaryField", start)
    internal = np.array(
        [[float(value) for value in match]
         for match in re.findall(VECTOR, text[start:boundary])]
    )

    # The first value after the patch name (a patch can hold nested
    # dictionaries, e.g. displacementSeries)
    values = {}
    for patch in patches:
        start = re.search(rf"\n\s*{patch}\s*\n\s*\{{", text[boundary:])
        match = start and re.compile(rf"value\s+uniform\s+{VECTOR}").search(
            text, boundary + start.end()
        )
        if not match:
            raise ValueError(f"No uniform value on patch {patch} in {filename}")
        values[patch] = np.array([float(value) for value in match.groups()])

    return internal, values


def axial_force_from_fields(path):
    """Axial force at the anchor (patch left) and attachment (patch right)
    from the beam fields in the time directories, for runs without
    axialForce<beam>.dat. Q on an end face along the outward line tangent,
    positive in tension, as in the finiteVolumeBeam restraint. The tangent
    is taken from the end-cell centre to the end face.

    The beam mesh is straight along x from 0 to L0 in its own coordinates
    (createCircularBeamMesh); refW moves it to the reference line and W is
    the displacement from there."""
    points_file = os.path.join(path, "constant", BEAM_NAME, "polyMesh", "points")
    ref_file = os.path.join(path, "0", BEAM_NAME, "refW")
    if not (os.path.isfile(points_file) and os.path.isfile(ref_file)):
        return None

    with open(points_file, "r", encoding="utf-8") as handle:
        points = np.array(
            [[float(value) for value in match]
             for match in re.findall(VECTOR, handle.read())]
        )
    length = points[:, 0].max()
    if points[:, 0].min() < -1e-9*length or np.abs(points[:, 1:]).max() > 0.1*length:
        raise ValueError(f"{points_file}: beam mesh is not straight along x")

    ref_cells, ref_faces = read_field_values(ref_file, ("left", "right"))
    n = len(ref_cells)
    cell_x = (np.arange(n) + 0.5)*length/n
    local_cells = np.column_stack([cell_x, np.zeros(n), np.zeros(n)])
    local_faces = {"left": np.zeros(3), "right": np.array([length, 0, 0])}
    end_cells = {"left": 0, "right": n - 1}

    times = []
    for entry in os.listdir(path):
        try:
            time = float(entry)
        except ValueError:
            continue
        if os.path.isfile(os.path.join(path, entry, BEAM_NAME, "Q")):
            times.append((time, entry))

    rows = []
    for time, entry in sorted(times):
        W_cells, W_faces = read_field_values(
            os.path.join(path, entry, BEAM_NAME, "W"), ("left", "right")
        )
        _, Q_faces = read_field_values(
            os.path.join(path, entry, BEAM_NAME, "Q"), ("left", "right")
        )
        row = [time]
        for patch in ("left", "right"):
            i = end_cells[patch]
            face = local_faces[patch] + ref_faces[patch] + W_faces[patch]
            cell = local_cells[i] + ref_cells[i] + W_cells[i]
            tangent = (face - cell)/np.linalg.norm(face - cell)
            axial = Q_faces[patch] @ tangent
            row += [axial, np.linalg.norm(Q_faces[patch] - axial*tangent)]
        rows.append(row)

    if not rows:
        return None

    data = np.asarray(rows)
    # Same columns as axialForce<beam>.dat
    return data[:, 0], data[:, 1:]


def tension_files(path, model):
    """Fairlead and anchor tension (time and a single-column array each) and
    where they come from; for the beamFoam cases also the end-force
    magnitude |Q| at the box and at the anchor"""
    if model == "moorDyn":
        # Time, FairTen1, AnchTen1 after a name row and a unit row
        lines_file = os.path.join(path, "Mooring", "lines.out")
        if not os.path.isfile(lines_file):
            return None
        time, tension = read_columns(lines_file, 2, skip_header=2)
        return {
            "fairlead": (time, tension[:, [0]]),
            "anchor": (time, tension[:, [1]]),
            "source": "Mooring/lines.out",
        }

    post = os.path.join(path, "postProcessing", "0")
    attachment_file = os.path.join(post, f"attachmentForce{BEAM_NAME}.dat")
    anchor_file = os.path.join(post, f"anchorForce{BEAM_NAME}.dat")
    axial_file = os.path.join(post, f"axialForce{BEAM_NAME}.dat")
    if not os.path.isfile(attachment_file):
        return None

    result = {}
    for name, filename in (("fairlead_Q", attachment_file), ("anchor_Q", anchor_file)):
        time, force = read_columns(filename, 3)
        result[name] = (time, np.linalg.norm(force, axis=1)[:, None])

    # Columns: anchor axial, anchor shear, attachment axial, attachment shear
    if os.path.isfile(axial_file):
        axial = read_columns(axial_file, 4)
        result["source"] = os.path.relpath(axial_file, path)
    else:
        axial = axial_force_from_fields(path)
        result["source"] = "beam fields in the time directories"
        if axial is None:
            print(f"{path}: no axial force output or beam fields")
            return None

    time, values = axial
    result["anchor"] = (time, values[:, [0]])
    result["fairlead"] = (time, values[:, [2]])
    return result


def load_case(case_dir, label, model, colour, style):
    path = os.path.join(SCRIPT_DIR, case_dir)
    post = os.path.join(path, "postProcessing")
    motion_name = (
        "sixDoFRigidBodyState.dat" if model == "moorDyn"
        else "sixDoFRigidBodyStateFvBeam.dat"
    )
    motion_file = os.path.join(post, "sixDoF_History", "0", motion_name)
    height_file = os.path.join(post, "interfaceHeight1", "0", "height.dat")
    log_file = os.path.join(path, "log.interFoam")

    if not os.path.isfile(motion_file):
        print(f"Skipping {label}: {motion_file} not found")
        return None

    tensions = tension_files(path, model)
    if tensions is None:
        print(f"Skipping {label}: no mooring line force output")
        return None

    motion = read_sixdof_history(motion_file)
    time, displacement, rotation = apply_time_window(
        motion["time"],
        motion["centre_of_rotation"] - read_initial_centre_of_mass(path),
        motion["rotation"],
    )

    case = {
        "label": label,
        "model": model,
        "colour": colour,
        "style": style,
        "time": time,
        "displacement": displacement,
        "rotation": rotation,
    }

    case["tension_source"] = tensions["source"]
    for end in ("fairlead", "anchor"):
        case[f"{end}_time"], case[f"{end}_tension"] = apply_time_window(*tensions[end])
        if f"{end}_Q" in tensions:
            case[f"{end}_Q_time"], case[f"{end}_Q"] = apply_time_window(
                *tensions[f"{end}_Q"]
            )

    # Free-surface probes at x = 0.25 and 0.75 m: columns are height above
    # the bottom and above the probe location for each probe
    if os.path.isfile(height_file):
        height_time, height = read_columns(height_file, 4)
        case["height_time"], case["height"] = apply_time_window(
            height_time, height[:, [0, 2]]
        )

    if os.path.isfile(log_file):
        log = read_log(log_file)
        window = apply_time_window(log["time"], *(log[k] for k in log if k != "time"))
        case["log"] = dict(zip(log.keys(), window))

    return case


cases = [case for case in (load_case(*entry) for entry in CASES) if case is not None]

if not cases:
    raise RuntimeError("None of the cases in CASES have results")

reference = next((case for case in cases if case["label"] == REFERENCE_CASE), None)
if reference is None:
    print(f"Reference case {REFERENCE_CASE} has no results: difference plots skipped")


def difference_from_reference(case, time_key, value_key):
    """Case values minus the reference, interpolated onto the case times"""
    ref_time = reference[time_key]
    ref_values = reference[value_key]
    time = case[time_key]

    inside = (time >= ref_time[0]) & (time <= ref_time[-1])
    interpolated = np.column_stack(
        [np.interp(time[inside], ref_time, ref_values[:, i])
         for i in range(ref_values.shape[1])]
    )
    return time[inside], case[value_key][inside] - interpolated


# =============================================================
# PRINT SUMMARY: MOTION AND TENSION
# =============================================================

print("")
print("Cases")
for case in cases:
    final = case["displacement"][-1]
    print(
        f"  {case['label']:<24} {len(case['time']):5d} motion rows to "
        f"t = {case['time'][-1]:.3f} s, final displacement "
        f"({final[0]: .4e} {final[1]: .4e} {final[2]: .4e}) m"
    )
    print(f"  {'':<24} tension from {case['tension_source']}")

if reference is not None:
    comparisons = [
        ("time", "displacement", "displacement", "m"),
        ("time", "rotation", "rotation", "rad"),
        ("fairlead_time", "fairlead_tension", "fairlead tension", "N"),
        ("anchor_time", "anchor_tension", "anchor tension", "N"),
    ]
    if "height" in reference:
        comparisons.append(("height_time", "height", "surface height", "m"))

    print("")
    print(f"Maximum difference from {reference['label']}")
    for time_key, value_key, name, unit in comparisons:
        peak = np.max(np.linalg.norm(reference[value_key], axis=1))
        print(f"  {name} (reference peak {peak:.4e} {unit})")
        for case in cases:
            if case is reference or value_key not in case:
                continue
            _, diff = difference_from_reference(case, time_key, value_key)
            if diff.size == 0:
                continue
            max_diff = np.max(np.linalg.norm(diff, axis=1))
            print(
                f"    {case['label']:<22} {max_diff:.3e} {unit} "
                f"({max_diff/peak:.2e} of peak)"
            )


def split_drift(time, values, dt=0.005):
    """Uniform time, moving average over one wave period (drift) and the
    remainder (wave-frequency oscillation). Half a period at each end is
    dropped, where the moving average is incomplete"""
    uniform = np.arange(time[0], time[-1], dt)
    n = max(int(round(WAVE_PERIOD/dt)), 1)

    # Not yet a full wave period (e.g. a run that has just started)
    if len(uniform) <= n:
        return None

    resampled = np.interp(uniform, time, values)
    drift = np.convolve(resampled, np.ones(n)/n, mode="same")

    keep = slice(n//2, len(uniform) - n//2)
    return uniform[keep], drift[keep], (resampled - drift)[keep]


# (name, unit, time key, values from a case)
SPLIT_QUANTITIES = [
    ("surge", "m", "time", lambda c: c["displacement"][:, 0]),
    ("heave", "m", "time", lambda c: c["displacement"][:, 2]),
    ("pitch", "rad", "time", lambda c: c["rotation"][:, 1]),
    ("fairlead tension", "N", "fairlead_time", lambda c: c["fairlead_tension"][:, 0]),
    ("anchor tension", "N", "anchor_time", lambda c: c["anchor_tension"][:, 0]),
]

for case in cases:
    case["split"] = {}
    for name, unit, time_key, values in SPLIT_QUANTITIES:
        split = split_drift(case[time_key], values(case))
        if split is not None:
            case["split"][name] = split

print("")
print(
    f"Drift (mean) and wave-frequency amplitude (half peak-to-peak) "
    f"from t = {STEADY_START} s"
)
for name, unit, _, _ in SPLIT_QUANTITIES:
    print(f"  {name} ({unit})")
    for case in cases:
        if name not in case["split"]:
            continue
        time, drift, oscillation = case["split"][name]
        steady = time >= STEADY_START
        if not np.any(steady):
            print(f"    {case['label']:<24} no data after t = {STEADY_START} s")
            continue
        print(
            f"    {case['label']:<24} mean {np.mean(drift[steady]): .4e}, "
            f"amplitude {0.5*np.ptp(oscillation[steady]):.4e}"
        )

print("")
print("Peak and minimum line tension (negative: compression)")
for case in cases:
    for name, time_key, value_key in (
        ("fairlead", "fairlead_time", "fairlead_tension"),
        ("anchor", "anchor_time", "anchor_tension"),
    ):
        tension = case[value_key][:, 0]
        i = np.argmax(tension)
        j = np.argmin(tension)
        print(
            f"  {case['label']:<24} {name:<9} peak {tension[i]: .4e} N at t = "
            f"{case[time_key][i]:.3f} s, minimum {tension[j]: .4e} N at t = "
            f"{case[time_key][j]:.3f} s"
        )


# =============================================================
# PRINT SUMMARY: EFFICIENCY
# =============================================================

def efficiency_metrics(log):
    simulated = log["time"][-1] - log["time"][0]
    newton = log["newton_iterations"]
    finite_newton = newton[np.isfinite(newton)]

    return {
        "time steps": float(log["time"].size),
        "last simulated time (s)": float(log["time"][-1]),
        "mean deltaT (s)": float(np.nanmean(log["delta_t"])),
        "ExecutionTime (s)": float(log["execution_time"][-1]),
        "ClockTime (s)": float(log["clock_time"][-1]),
        "ExecutionTime / simulated s": (
            float(log["execution_time"][-1])/simulated if simulated > 0 else np.nan
        ),
        "mean step ExecutionTime (s)": float(np.mean(log["step_execution_time"])),
        "max step ExecutionTime (s)": float(np.max(log["step_execution_time"])),
        "mean linear solver iterations": float(np.mean(log["linear_iterations"])),
        "total beam Newton iterations": (
            float(np.sum(finite_newton)) if finite_newton.size else np.nan
        ),
        "mean beam Newton iterations": (
            float(np.mean(finite_newton)) if finite_newton.size else np.nan
        ),
    }


logged = [case for case in cases if "log" in case]
for case in logged:
    case["efficiency"] = efficiency_metrics(case["log"])

if logged:
    ref_log = reference if reference is not None and "log" in reference else None
    width = max(16, *(len(case["label"]) + 2 for case in logged))

    print("")
    print("Computational efficiency (log.interFoam)")
    header = f"{'metric':32s}" + "".join(f"{case['label']:>{width}s}" for case in logged)
    if ref_log is not None:
        header += f"   ExecutionTime / {ref_log['label']}"
    print(header)
    for metric in logged[0]["efficiency"]:
        row = f"{metric:32s}"
        for case in logged:
            value = case["efficiency"][metric]
            row += f"{'-':>{width}s}" if np.isnan(value) else f"{value:{width}.6g}"
        print(row)

    if ref_log is not None:
        ref_time = ref_log["efficiency"]["ExecutionTime (s)"]
        print("")
        print(f"ExecutionTime relative to {ref_log['label']}")
        for case in logged:
            print(
                f"  {case['label']:<24} "
                f"{case['efficiency']['ExecutionTime (s)']/ref_time:.4f}"
            )


# =============================================================
# PLOT DATA
# =============================================================

def finish_plot(ax, ylabel, xlabel="Time (s)", bottom=None):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if bottom is not None:
        ax.set_ylim(bottom=bottom)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(which="major", ls="--", lw=0.6, alpha=0.75)
    ax.grid(which="minor", ls=":", lw=0.4, alpha=0.45)
    ax.legend()
    plt.tight_layout()


def save_or_show(filename):
    if SAVE_FIGURES:
        os.makedirs(FIGURE_DIR, exist_ok=True)
        path = os.path.join(FIGURE_DIR, filename)
        plt.savefig(path, dpi=300)
        print("Saved:", path)

    plt.show()


def plot_lines(series, ylabel, filename, title=None, bottom=None):
    """series: (case, time, values) for each line"""
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()

    for case, time, values in series:
        ax.plot(
            time, values, color=case["colour"], linestyle=case["style"],
            linewidth=1.4, label=case["label"],
        )

    if title:
        ax.set_title(title)
    finish_plot(ax, ylabel, bottom=bottom)
    save_or_show(filename)


def plot_component(time_key, value_key, component, ylabel, filename, title=None):
    plot_lines(
        [(c, c[time_key], c[value_key][:, component])
         for c in cases if value_key in c],
        ylabel, filename, title,
    )


def plot_difference(time_key, value_key, component, ylabel, filename):
    if reference is None:
        return

    series = []
    for case in cases:
        if case is reference or value_key not in case:
            continue
        time, difference = difference_from_reference(case, time_key, value_key)
        series.append((case, time, difference[:, component]))

    plot_lines(series, ylabel, filename, f"Difference from {reference['label']}")


# ---- Motion ----

for label, component in (("Surge", 0), ("Heave", 2)):
    plot_component(
        "time", "displacement", component,
        f"Box {label.lower()} displacement (m)",
        f"box_{label.lower()}.png",
    )
    plot_difference(
        "time", "displacement", component,
        f"Box {label.lower()} difference (m)",
        f"box_{label.lower()}_difference.png",
    )

plot_component("time", "rotation", 1, "Box pitch (rad)", "box_pitch.png")
plot_difference("time", "rotation", 1, "Box pitch difference (rad)", "box_pitch_difference.png")

# Free surface: checks the waves reaching the box are the same in all runs
for label, component in (("x = 0.25 m", 0), ("x = 0.75 m", 1)):
    tag = label.replace(" ", "").replace("=", "").replace(".", "p")
    plot_component(
        "height_time", "height", component,
        f"Free-surface height at {label} (m)",
        f"surface_height_{tag}.png",
    )

# ---- Mooring tension ----

plot_component(
    "fairlead_time", "fairlead_tension", 0,
    "Line tension at the box (N)", "fairlead_tension.png",
)
plot_difference(
    "fairlead_time", "fairlead_tension", 0,
    "Line tension at the box, difference (N)", "fairlead_tension_difference.png",
)
plot_component(
    "anchor_time", "anchor_tension", 0,
    "Line tension at the anchor (N)", "anchor_tension.png",
)
plot_difference(
    "anchor_time", "anchor_tension", 0,
    "Line tension at the anchor, difference (N)", "anchor_tension_difference.png",
)

# After the snap load
for end, ylabel in (("fairlead", "Line tension at the box (N)"),
                    ("anchor", "Line tension at the anchor (N)")):
    series = []
    for case in cases:
        time = case[f"{end}_time"]
        after = time >= TENSION_ZOOM_START
        series.append((case, time[after], case[f"{end}_tension"][after, 0]))
    plot_lines(
        series, ylabel, f"{end}_tension_after_snap.png",
        f"From t = {TENSION_ZOOM_START} s (negative: compression)",
    )

# End-force magnitude |Q| of the beamFoam cases (the tension only while the
# line is taut), with the MoorDyn tension for reference
for end, ylabel in (("fairlead", "End force at the box (N)"),
                    ("anchor", "End force at the anchor (N)")):
    series = []
    for case in cases:
        key = f"{end}_Q" if f"{end}_Q" in case else f"{end}_tension"
        series.append((case, case[f"{key}_time" if key.endswith("_Q") else f"{end}_time"], case[key][:, 0]))
    plot_lines(
        series, ylabel, f"{end}_force_magnitude.png",
        "beamFoam: |Q| (axial + shear); MoorDyn: tension",
    )

# Slow drift (moving average over one wave period)
for name, ylabel in (
    ("surge", "Box surge drift (m)"),
    ("heave", "Box heave drift (m)"),
    ("fairlead tension", "Line tension at the box, drift (N)"),
    ("anchor tension", "Line tension at the anchor, drift (N)"),
):
    plot_lines(
        [(c, c["split"][name][0], c["split"][name][1])
         for c in cases if name in c["split"]],
        ylabel, f"{name.replace(' ', '_')}_drift.png",
        f"Moving average over {WAVE_PERIOD} s",
    )

# Wave-frequency oscillation (drift removed) in ZOOM_WINDOW: amplitude and
# phase
for name, ylabel in (
    ("surge", "Box surge oscillation (m)"),
    ("heave", "Box heave oscillation (m)"),
    ("pitch", "Box pitch oscillation (rad)"),
    ("fairlead tension", "Line tension at the box, oscillation (N)"),
):
    series = []
    for case in cases:
        if name not in case["split"]:
            continue
        time, _, oscillation = case["split"][name]
        window = (time >= ZOOM_WINDOW[0]) & (time <= ZOOM_WINDOW[1])
        series.append((case, time[window], oscillation[window]))
    plot_lines(
        series, ylabel, f"{name.replace(' ', '_')}_oscillation_zoom.png",
        "Drift removed",
    )

# ---- Efficiency ----

if logged:
    plot_lines(
        [(c, c["log"]["time"], c["log"]["execution_time"]) for c in logged],
        "Cumulative ExecutionTime (s)", "cumulative_execution_time.png",
        bottom=0,
    )

    # Moving average over 50 steps: single steps are noisy
    def smoothed(values, n=50):
        if values.size <= n:
            return values
        return np.convolve(values, np.ones(n)/n, mode="same")

    plot_lines(
        [(c, c["log"]["time"], smoothed(c["log"]["step_execution_time"]))
         for c in logged],
        "Step ExecutionTime, 50-step average (s)", "step_execution_time.png",
        bottom=0,
    )
    plot_lines(
        [(c, c["log"]["time"], c["log"]["delta_t"]) for c in logged],
        "Time-step size (s)", "time_step_size.png", bottom=0,
    )
    plot_lines(
        [(c, c["log"]["time"], smoothed(c["log"]["linear_iterations"]))
         for c in logged],
        "Linear solver iterations per step, 50-step average",
        "linear_solver_iterations.png", bottom=0,
    )

    beam_logged = [
        c for c in logged if np.any(np.isfinite(c["log"]["newton_iterations"]))
    ]
    if beam_logged:
        plot_lines(
            [(c, c["log"]["time"], np.nancumsum(c["log"]["newton_iterations"]))
             for c in beam_logged],
            "Cumulative beam Newton iterations", "cumulative_newton_iterations.png",
            title="beamFoam cases (MoorDyn has no Newton iterations)", bottom=0,
        )

    # Total run time
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()
    x = np.arange(len(logged))
    width = 0.36
    ax.bar(
        x - width/2, [c["efficiency"]["ExecutionTime (s)"] for c in logged],
        width, label="ExecutionTime",
    )
    ax.bar(
        x + width/2, [c["efficiency"]["ClockTime (s)"] for c in logged],
        width, label="ClockTime",
    )
    ax.set_xticks(x)
    ax.set_xticklabels([c["label"] for c in logged])
    finish_plot(ax, "Total run time (s)", xlabel="", bottom=0)
    save_or_show("total_run_time.png")
