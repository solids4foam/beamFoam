#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare the floatRigidBody cases between coupling frameworks:

- partitioned: floatRigidBody/floatRigidBody_partioned (moorFV solver
  FvBeamNewmark; the body and beamFoam take turns)
- monolithic: floatRigidBody/floatRigidBody_monolithic (moorFV solver
  beamFoamCoupled; beamFoam solves the line and the body together)

each with 3 PIMPLE outer correctors and with 8 (the _nOuter8 cases), to see
how far the coupling has converged within each time step.

Both use the same plane/axis constraints, so the box only surges, heaves and
pitches. Sway, roll and yaw are still plotted as a check that the
constraints hold (they should be exactly zero).

Differences are plotted against REFERENCE_CASE.

Run from Spyder or from the beamFoam/tutorials directory after the cases have
been run with ./Allrun. Cases that have not been run are skipped, and a case
that is still running is plotted up to its latest output.
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

# (case directory, legend label, colour, line style)
CASES = [
    ("floatRigidBody/floatRigidBody_partioned_nOuter8", "Partitioned, 8 outer correctors", "k", "-"),
    ("floatRigidBody/floatRigidBody_partioned", "Partitioned, 3 outer correctors", "b", "-"),
    ("floatRigidBody/floatRigidBody_monolithic", "Monolithic, 3 outer correctors", "r", "--"),
    ("floatRigidBody/floatRigidBody_monolithic_nOuter8", "Monolithic, 8 outer correctors", "m", ":"),
]

# Differences are taken against this case (the best-converged partitioned run)
REFERENCE_CASE = "floatRigidBody/floatRigidBody_partioned_nOuter8"

# Beam region of the mooring line
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

SAVE_FIGURES = False
FIGURE_DIR = os.path.join(SCRIPT_DIR, "comparison_plots", "floatRigidBody")


# =============================================================
# READ DATA
# =============================================================

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
VECTOR = rf"\(\s*({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*\)"


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
        "velocity": data[:, 10:13],
        "omega": data[:, 13:16],
    }


def read_columns(filename, ncols):
    """Time and ncols values per row, last row per time"""
    data = np.genfromtxt(filename, comments="#", invalid_raise=False)

    if data.size == 0:
        raise ValueError(f"No rows found in {filename}")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    data = last_row_per_time(data[:, :ncols + 1])
    return data[:, 0], data[:, 1:]


def apply_time_window(time, *arrays):
    mask = np.ones_like(time, dtype=bool)

    if T_START is not None:
        mask &= time >= T_START

    if T_END is not None:
        mask &= time <= T_END

    return (time[mask],) + tuple(array[mask] for array in arrays)


def load_case(case_dir, label, colour, style):
    path = os.path.join(SCRIPT_DIR, case_dir)
    post = os.path.join(path, "postProcessing")
    motion_file = os.path.join(
        post, "sixDoF_History", "0", "sixDoFRigidBodyStateFvBeam.dat"
    )
    attachment_file = os.path.join(post, "0", f"attachmentForce{BEAM_NAME}.dat")
    anchor_file = os.path.join(post, "0", f"anchorForce{BEAM_NAME}.dat")
    height_file = os.path.join(post, "interfaceHeight1", "0", "height.dat")

    for filename in (motion_file, attachment_file):
        if not os.path.isfile(filename):
            print(f"Skipping {case_dir}: {os.path.relpath(filename, SCRIPT_DIR)} not found")
            return None

    motion = read_sixdof_history(motion_file)
    centre_of_mass = read_initial_centre_of_mass(path)

    time, displacement, rotation, velocity, omega = apply_time_window(
        motion["time"],
        motion["centre_of_rotation"] - centre_of_mass,
        motion["rotation"],
        motion["velocity"],
        motion["omega"],
    )

    case = {
        "name": case_dir,
        "label": label,
        "colour": colour,
        "style": style,
        "time": time,
        "displacement": displacement,
        "rotation": rotation,
        "velocity": velocity,
        "omega": omega,
    }

    # Line tension at the box (attachment) and at the seabed (anchor)
    case["attachment_time"], case["attachment_force"] = apply_time_window(
        *read_columns(attachment_file, 3)
    )
    case["attachment_tension"] = np.linalg.norm(case["attachment_force"], axis=1)[:, None]

    if os.path.isfile(anchor_file):
        case["anchor_time"], case["anchor_force"] = apply_time_window(
            *read_columns(anchor_file, 3)
        )
        case["anchor_tension"] = np.linalg.norm(case["anchor_force"], axis=1)[:, None]

    # Free-surface probes at x = 0.25 and 0.75 m: columns are height above
    # the bottom and above the probe location for each probe
    if os.path.isfile(height_file):
        height_time, height = read_columns(height_file, 4)
        case["height_time"], case["height"] = apply_time_window(
            height_time, height[:, [0, 2]]
        )

    return case


cases = [case for case in (load_case(*entry) for entry in CASES) if case is not None]

if not cases:
    raise RuntimeError("None of the cases in CASES have results")

reference = next((case for case in cases if case["name"] == REFERENCE_CASE), None)
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
# PRINT SUMMARY
# =============================================================

print("")
print("Cases")
for case in cases:
    final = case["displacement"][-1]
    print(
        f"  {case['label']:<32} {len(case['time']):5d} motion rows to "
        f"t = {case['time'][-1]:.3f} s, final displacement "
        f"({final[0]: .4e} {final[1]: .4e} {final[2]: .4e}) m"
    )

print("")
print("Out-of-plane motion (constrained to zero)")
for case in cases:
    print(
        f"  {case['label']:<32} max |y| {np.max(np.abs(case['displacement'][:, 1])):.3e} m, "
        f"max |roll| {np.max(np.abs(case['rotation'][:, 0])):.3e} rad, "
        f"max |yaw| {np.max(np.abs(case['rotation'][:, 2])):.3e} rad"
    )

if reference is not None:
    comparisons = [
        ("time", "displacement", "displacement", "m"),
        ("time", "rotation", "rotation", "rad"),
        ("attachment_time", "attachment_force", "attachment force", "N"),
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
                f"    {case['label']:<28} {max_diff:.3e} {unit} "
                f"({max_diff/peak:.2e} of peak)"
            )


# =============================================================
# DRIFT AND WAVE-FREQUENCY MOTION
# =============================================================

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
    ("tension", "N", "attachment_time", lambda c: c["attachment_tension"][:, 0]),
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
            print(f"    {case['label']:<32} no data after t = {STEADY_START} s")
            continue
        print(
            f"    {case['label']:<32} mean {np.mean(drift[steady]): .4e}, "
            f"amplitude {0.5*np.ptp(oscillation[steady]):.4e}"
        )

print("")
print("Peak line tension at the box")
for case in cases:
    tension = case["attachment_tension"][:, 0]
    i = np.argmax(tension)
    print(
        f"  {case['label']:<32} {tension[i]:.4e} N at t = "
        f"{case['attachment_time'][i]:.3f} s"
    )


# =============================================================
# PLOT DATA
# =============================================================

def finish_plot(ax, ylabel):
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(ylabel)
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


def plot_component(time_key, value_key, component, ylabel, filename, title=None):
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()

    for case in cases:
        if value_key not in case:
            continue
        ax.plot(
            case[time_key],
            case[value_key][:, component],
            color=case["colour"],
            linestyle=case["style"],
            linewidth=1.4,
            label=case["label"],
        )

    if title:
        ax.set_title(title)
    finish_plot(ax, ylabel)
    save_or_show(filename)


def plot_difference(time_key, value_key, component, ylabel, filename):
    if reference is None:
        return

    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()

    for case in cases:
        if case is reference or value_key not in case:
            continue
        time, difference = difference_from_reference(case, time_key, value_key)
        ax.plot(
            time,
            difference[:, component],
            color=case["colour"],
            linestyle=case["style"],
            linewidth=1.4,
            label=case["label"],
        )

    ax.set_title(f"Difference from {reference['label']}")
    finish_plot(ax, ylabel)
    save_or_show(filename)


# Surge, heave and pitch: the in-plane motion both solvers resolve
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

# Sway, roll and yaw: zero when the constraints hold
plot_component(
    "time", "displacement", 1, "Box sway displacement (m)", "box_sway.png",
    title="Out-of-plane motion (zero when constrained)",
)
for label, component in (("Roll", 0), ("Yaw", 2)):
    plot_component(
        "time", "rotation", component, f"Box {label.lower()} (rad)",
        f"box_{label.lower()}.png",
        title="Out-of-plane motion (zero when constrained)",
    )

# Line force on the box
for label, component in (("X", 0), ("Z", 2)):
    plot_component(
        "attachment_time", "attachment_force", component,
        f"Line attachment force, {label} (N)",
        f"attachment_force_{label.lower()}.png",
    )
    plot_difference(
        "attachment_time", "attachment_force", component,
        f"Line attachment force difference, {label} (N)",
        f"attachment_force_{label.lower()}_difference.png",
    )

plot_component(
    "attachment_time", "attachment_tension", 0,
    "Line tension at the box (N)", "attachment_tension.png",
)
plot_component(
    "anchor_time", "anchor_tension", 0,
    "Line tension at the anchor (N)", "anchor_tension.png",
)

# Free surface: checks the waves reaching the box are the same in both runs
for label, component in (("x = 0.25 m", 0), ("x = 0.75 m", 1)):
    tag = label.replace(" ", "").replace("=", "").replace(".", "p")
    plot_component(
        "height_time", "height", component,
        f"Free-surface height at {label} (m)",
        f"surface_height_{tag}.png",
    )


# Slow drift (moving average over one wave period)
for name, ylabel in (
    ("surge", "Box surge drift (m)"),
    ("heave", "Box heave drift (m)"),
    ("tension", "Line tension drift (N)"),
):
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()
    for case in cases:
        if name not in case["split"]:
            continue
        time, drift, _ = case["split"][name]
        ax.plot(
            time, drift, color=case["colour"], linestyle=case["style"],
            linewidth=1.4, label=case["label"],
        )
    ax.set_title(f"Moving average over {WAVE_PERIOD} s")
    finish_plot(ax, ylabel)
    save_or_show(f"{name}_drift.png")

# Wave-frequency oscillation (drift removed) in ZOOM_WINDOW: amplitude and
# phase
for name, ylabel in (
    ("surge", "Box surge oscillation (m)"),
    ("heave", "Box heave oscillation (m)"),
    ("pitch", "Box pitch oscillation (rad)"),
):
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()
    for case in cases:
        if name not in case["split"]:
            continue
        time, _, oscillation = case["split"][name]
        window = (time >= ZOOM_WINDOW[0]) & (time <= ZOOM_WINDOW[1])
        ax.plot(
            time[window], oscillation[window], color=case["colour"],
            linestyle=case["style"], linewidth=1.4, label=case["label"],
        )
    ax.set_title("Drift removed")
    finish_plot(ax, ylabel)
    save_or_show(f"{name}_oscillation_zoom.png")
