#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare rigidBodyAndBeam results between coupling frameworks:

- monolithic: rigidBodyAndBeam/rigidBodyAndBeam_beamFoamCoupled (moorFV
  solver beamFoamCoupled; beamFoam solves the beam and the body together)
- loop: rigidBodyAndBeam/rigidBodyAndBeam_loop (moorFV solver FvBeamNewmark;
  the body and beamFoam take turns)

each with 1 PIMPLE outer corrector and with 8 (converged coupling within each
time step). Differences are plotted against REFERENCE_CASE. Set CASE_SET to
plot one of the stress-test variants in rigidBodyAndBeam/rigidBodyAndBeam_stress
instead, or "floatRigidBody" for the moored floating box in waves
(floatRigidBody/floatRigidBody_partioned_gravity against
floatRigidBody/floatRigidBody_monolithic_gravity; for more detail on that
case see plot_monolithic_vs_partitioned_floatRigidBody.py).

Run from Spyder or from the beamFoam/tutorials directory after the cases have
been run with ./Allrun. Cases that have not been run are skipped.
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

# "base" for the rigidBodyAndBeam cases, one of the Phase 4 stress-test
# variants in rigidBodyAndBeam/rigidBodyAndBeam_stress: "stiffLine",
# "lightBody", "stiffLight", "stiffLargeDt", or "floatRigidBody"
CASE_SET = "floatRigidBody"

# Force from the beam on the body, in postProcessing/0 (named after the beam
# region)
FORCE_FILE = "forcebeam.dat"

if CASE_SET == "base":
    # (case directory, legend label, colour, line style)
    CASES = [
        ("rigidBodyAndBeam/rigidBodyAndBeam_loop_nOuter8", "Loop, 8 outer correctors", "k", "-"),
        ("rigidBodyAndBeam/rigidBodyAndBeam_loop", "Loop, 1 outer corrector", "b", "-"),
        ("rigidBodyAndBeam/rigidBodyAndBeam_beamFoamCoupled", "Monolithic, 1 outer corrector", "r", "--"),
        ("rigidBodyAndBeam/rigidBodyAndBeam_beamFoamCoupled_nOuter8", "Monolithic, 8 outer correctors", "m", ":"),
    ]

    # Differences are taken against this case (the converged loop)
    REFERENCE_CASE = "rigidBodyAndBeam/rigidBodyAndBeam_loop_nOuter8"
elif CASE_SET == "floatRigidBody":
    CASES = [
        ("floatRigidBody/floatRigidBody_partioned_gravity", "Partitioned, 3 outer correctors", "b", "-"),
        ("floatRigidBody/floatRigidBody_monolithic_gravity", "Monolithic, 3 outer correctors", "r", "--"),
    ]

    # Differences are taken against the partitioned run
    REFERENCE_CASE = "floatRigidBody/floatRigidBody_partioned_gravity"

    # The mooring line region is beamone
    FORCE_FILE = "forcebeamone.dat"
else:
    stress = os.path.join("rigidBodyAndBeam/rigidBodyAndBeam_stress", CASE_SET)
    CASES = [
        (os.path.join(stress, "monolithic_nOuter8"), "Monolithic, 8 outer correctors", "k", "-"),
        (os.path.join(stress, "loop_nOuter8"), "Loop, 8 outer correctors", "m", ":"),
        (os.path.join(stress, "loop_nOuter1"), "Loop, 1 outer corrector", "b", "-"),
        (os.path.join(stress, "monolithic_nOuter1"), "Monolithic, 1 outer corrector", "r", "--"),
    ]

    # The stress tests are compared against the converged monolithic run
    REFERENCE_CASE = os.path.join(stress, "monolithic_nOuter8")


# =============================================================
# USER SETTINGS
# =============================================================

T_START = None
T_END = None

SAVE_FIGURES = False
FIGURE_DIR = os.path.join(SCRIPT_DIR, "comparison_plots", CASE_SET)


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

    data = np.asarray(rows)

    # A time can appear more than once (restarts); keep the last row
    _, last = np.unique(data[::-1, 0], return_index=True)
    data = data[::-1][last]

    return {
        "time": data[:, 0],
        "centre_of_rotation": data[:, 1:4],
        "rotation": data[:, 7:10],
        "velocity": data[:, 10:13],
        "omega": data[:, 13:16],
    }


def read_beam_force(filename):
    """Force from the beam on the body (FORCE_FILE)"""
    data = np.genfromtxt(filename, comments="#")

    if data.size == 0:
        raise ValueError(f"No force rows found in {filename}")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    # The restraint writes once per update; keep the last row per time
    _, last = np.unique(data[::-1, 0], return_index=True)
    data = data[::-1][last]

    return data[:, 0], data[:, 1:4]


def apply_time_window(time, *arrays):
    mask = np.ones_like(time, dtype=bool)

    if T_START is not None:
        mask &= time >= T_START

    if T_END is not None:
        mask &= time <= T_END

    return (time[mask],) + tuple(array[mask] for array in arrays)


def load_case(case_dir, label, colour, style):
    path = os.path.join(SCRIPT_DIR, case_dir)
    motion_file = os.path.join(
        path, "postProcessing", "sixDoF_History", "0", "sixDoFRigidBodyStateFvBeam.dat"
    )
    force_file = os.path.join(path, "postProcessing", "0", FORCE_FILE)

    for filename in (motion_file, force_file):
        if not os.path.isfile(filename):
            print(f"Skipping {case_dir}: {os.path.relpath(filename, SCRIPT_DIR)} not found")
            return None

    motion = read_sixdof_history(motion_file)
    force_time, force = read_beam_force(force_file)
    centre_of_mass = read_initial_centre_of_mass(path)

    time, displacement, rotation = apply_time_window(
        motion["time"],
        motion["centre_of_rotation"] - centre_of_mass,
        motion["rotation"],
    )
    force_time, force = apply_time_window(force_time, force)

    return {
        "name": case_dir,
        "label": label,
        "colour": colour,
        "style": style,
        "time": time,
        "displacement": displacement,
        "rotation": rotation,
        "force_time": force_time,
        "force": force,
    }


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
        [np.interp(time[inside], ref_time, ref_values[:, i]) for i in range(3)]
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
        f"  {case['label']:<34} {len(case['time']):4d} motion rows, "
        f"final displacement ({final[0]: .6e} {final[1]: .6e} {final[2]: .6e}) m"
    )

if reference is not None:
    peak_displacement = np.max(np.linalg.norm(reference["displacement"], axis=1))
    peak_force = np.max(np.linalg.norm(reference["force"], axis=1))

    print("")
    print(f"Maximum difference from {reference['label']}")
    print(f"  (peak displacement {peak_displacement:.4e} m, peak force {peak_force:.4e} N)")
    for case in cases:
        if case is reference:
            continue
        _, d_disp = difference_from_reference(case, "time", "displacement")
        _, d_force = difference_from_reference(case, "force_time", "force")
        max_disp = np.max(np.linalg.norm(d_disp, axis=1))
        max_force = np.max(np.linalg.norm(d_force, axis=1))
        print(
            f"  {case['label']:<34} displacement {max_disp:.3e} m "
            f"({max_disp/peak_displacement:.2e} of peak), "
            f"force {max_force:.3e} N ({max_force/peak_force:.2e} of peak)"
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


def plot_component(time_key, value_key, component, ylabel, filename, scale=1.0):
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()

    for case in cases:
        ax.plot(
            case[time_key],
            scale*case[value_key][:, component],
            color=case["colour"],
            linestyle=case["style"],
            linewidth=1.6,
            label=case["label"],
        )

    finish_plot(ax, ylabel)
    save_or_show(filename)


def plot_difference(time_key, value_key, component, ylabel, filename):
    if reference is None:
        return

    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()

    for case in cases:
        if case is reference:
            continue
        time, difference = difference_from_reference(case, time_key, value_key)
        ax.plot(
            time,
            difference[:, component],
            color=case["colour"],
            linestyle=case["style"],
            linewidth=1.6,
            label=case["label"],
        )

    ax.set_title(f"Difference from {reference['label']}")
    finish_plot(ax, ylabel)
    save_or_show(filename)


for label, component in (("X", 0), ("Height", 2)):
    plot_component(
        "time", "displacement", component,
        f"Body {label.lower()} displacement (m)",
        f"rigid_body_{label.lower()}_displacement.png",
    )
    plot_difference(
        "time", "displacement", component,
        f"Body {label.lower()} displacement difference (m)",
        f"rigid_body_{label.lower()}_displacement_difference.png",
    )

plot_component(
    "time", "rotation", 1,
    "Body rotation about y (rad)",
    "rigid_body_rotation_y.png",
)

for label, component in (("X", 0), ("Z", 2)):
    plot_component(
        "force_time", "force", component,
        f"Beam force on body, {label} (N)",
        f"beam_force_{label.lower()}.png",
    )
    plot_difference(
        "force_time", "force", component,
        f"Beam force difference, {label} (N)",
        f"beam_force_{label.lower()}_difference.png",
    )
