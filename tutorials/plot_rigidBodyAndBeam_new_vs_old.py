#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare rigidBodyAndBeam results from the new BlockEigen-coupled path against
the old coupled-solver case.

Run from Spyder or from the beamFoam/tutorials directory.
"""

import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MaxNLocator


# =============================================================
# INPUT FILES
# =============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

NEW_CASE = os.path.join(SCRIPT_DIR, "rigidBodyAndBeam")
OLD_CASE = os.path.join(SCRIPT_DIR, "rigidBodyAndBeam_coupledSolver")

NEW_MOTION_FILE = sorted(
    glob.glob(
        os.path.join(
            NEW_CASE,
            "postProcessing",
            "BlockEigenSolve_rigidBodyMotion",
            "*",
            "BlockEigenSolve_rigidBodyMotion.dat",
        )
    )
)[0]

NEW_FORCE_COUPLING_FILE = sorted(
    glob.glob(
        os.path.join(
            NEW_CASE,
            "postProcessing",
            "BlockEigenSolve_rigidBodyForceCoupling",
            "*",
            "BlockEigenSolve_rigidBodyForceCoupling.dat",
        )
    )
)[0]

NEW_BEAM_FORCE_FILE = os.path.join(
    NEW_CASE,
    "postProcessing",
    "0",
    "forcebeam.dat",
)

OLD_MOTION_FILE = sorted(
    glob.glob(
        os.path.join(
            OLD_CASE,
            "postProcessing",
            "sixDoF*",
            "*",
            "sixDoFRigidBodyStateFvBeam.dat",
        )
    )
)[0]

OLD_BEAM_FORCE_FILE = os.path.join(
    OLD_CASE,
    "postProcessing",
    "0",
    "forcebeam.dat",
)

print("Input files")
print("New BlockEigen motion:       ", NEW_MOTION_FILE)
print("New BlockEigen force coupling:", NEW_FORCE_COUPLING_FILE)
print("New beam force:              ", NEW_BEAM_FORCE_FILE)
print("Old sixDoF motion:           ", OLD_MOTION_FILE)
print("Old beam force:              ", OLD_BEAM_FORCE_FILE)


# =============================================================
# USER SETTINGS
# =============================================================

REFERENCE_CENTRE_OF_ROTATION = np.array([2.6565, 0.0, 0.0])

# BlockEigen writes one row for every nonlinear solve. The last row at each
# time is the converged value.
BLOCKEIGEN_TIME_REDUCTION = "last"  # "last" or "all"

# Use "beam" to compare forcebeam.dat from both cases.
# Use "rigid_body_reaction" to compare the new force applied to the rigid body
# against the old restraint force, which is the negative of forcebeam.dat.
FORCE_COMPARISON = "beam"  # "beam" or "rigid_body_reaction"

T_START = None
T_END = None

SAVE_FIGURES = False
FIGURE_DIR = os.path.join(SCRIPT_DIR, "comparison_plots")


# =============================================================
# READ DATA
# =============================================================

def reduce_to_last_row_per_time(time, values):
    keep = []

    for t in np.unique(time):
        matching_rows = np.where(time == t)[0]
        keep.append(matching_rows[-1])

    keep = np.asarray(keep, dtype=int)
    return time[keep], values[keep]


def read_blockeigen_position(filename):
    data = np.genfromtxt(filename, names=True)

    if data.ndim == 0:
        data = np.array([data], dtype=data.dtype)

    time = np.asarray(data["Time"], dtype=float)
    displacement = np.column_stack(
        (
            np.asarray(data["disp_x"], dtype=float),
            np.asarray(data["disp_y"], dtype=float),
            np.asarray(data["disp_z"], dtype=float),
        )
    )

    position = displacement + REFERENCE_CENTRE_OF_ROTATION

    if BLOCKEIGEN_TIME_REDUCTION == "all":
        return time, position

    if BLOCKEIGEN_TIME_REDUCTION != "last":
        raise ValueError("BLOCKEIGEN_TIME_REDUCTION must be 'last' or 'all'")

    return reduce_to_last_row_per_time(time, position)


def read_sixdof_centre_of_rotation(filename):
    number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
    vector = rf"\(\s*({number})\s+({number})\s+({number})\s*\)"

    line_pattern = re.compile(
        rf"^\s*({number})\s+"
        rf"{vector}\s+"
        rf"{vector}\s+"
        rf"{vector}\s+"
        rf"{vector}\s+"
        rf"{vector}"
    )

    time = []
    centre_of_rotation = []

    with open(filename, "r", encoding="utf-8") as handle:
        for line in handle:
            match = line_pattern.search(line)

            if not match:
                continue

            values = [float(value) for value in match.groups()]
            time.append(values[0])
            centre_of_rotation.append(values[1:4])

    if not time:
        raise ValueError(f"No motion rows found in {filename}")

    return np.asarray(time), np.asarray(centre_of_rotation)


def read_beam_force(filename):
    data = np.genfromtxt(filename, comments="#")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    time = data[:, 0]
    force = data[:, 1:4]
    return time, force


def read_blockeigen_force_coupling(filename):
    data = np.genfromtxt(filename, names=True)

    if data.ndim == 0:
        data = np.array([data], dtype=data.dtype)

    time = np.asarray(data["Time"], dtype=float)
    force = np.column_stack(
        (
            np.asarray(data["force_x"], dtype=float),
            np.asarray(data["force_y"], dtype=float),
            np.asarray(data["force_z"], dtype=float),
        )
    )

    return reduce_to_last_row_per_time(time, force)


def apply_time_window(time, values):
    mask = np.ones_like(time, dtype=bool)

    if T_START is not None:
        mask &= time >= T_START

    if T_END is not None:
        mask &= time <= T_END

    return time[mask], values[mask]


new_motion_time, new_position = read_blockeigen_position(NEW_MOTION_FILE)
old_motion_time, old_position = read_sixdof_centre_of_rotation(OLD_MOTION_FILE)

if FORCE_COMPARISON == "beam":
    new_force_time, new_force = read_beam_force(NEW_BEAM_FORCE_FILE)
    old_force_time, old_force = read_beam_force(OLD_BEAM_FORCE_FILE)
    force_ylabel_prefix = "beam attachment force"
elif FORCE_COMPARISON == "rigid_body_reaction":
    new_force_time, new_force = read_blockeigen_force_coupling(
        NEW_FORCE_COUPLING_FILE
    )
    old_force_time, old_beam_force = read_beam_force(OLD_BEAM_FORCE_FILE)
    old_force = -old_beam_force
    force_ylabel_prefix = "rigid-body reaction force"
else:
    raise ValueError("FORCE_COMPARISON must be 'beam' or 'rigid_body_reaction'")

new_motion_time, new_position = apply_time_window(
    new_motion_time,
    new_position,
)
old_motion_time, old_position = apply_time_window(
    old_motion_time,
    old_position,
)
new_force_time, new_force = apply_time_window(new_force_time, new_force)
old_force_time, old_force = apply_time_window(old_force_time, old_force)


# =============================================================
# PRINT CHECKS
# =============================================================

def print_position_summary(name, time, position):
    displacement = position - REFERENCE_CENTRE_OF_ROTATION
    rows = [0]

    if len(time) > 1:
        rows.append(len(time) - 1)

    print("")
    print(name)
    print(
        "row  time"
        "  x_disp  y_disp  z_disp"
        "  x_position  y_position  z_position"
    )

    for row in rows:
        print(
            f"{row:3d}"
            f"  {time[row]:.10g}"
            f"  {displacement[row, 0]: .10e}"
            f"  {displacement[row, 1]: .10e}"
            f"  {displacement[row, 2]: .10e}"
            f"  {position[row, 0]: .10e}"
            f"  {position[row, 1]: .10e}"
            f"  {position[row, 2]: .10e}"
        )


def print_force_summary(name, time, force):
    rows = [0]

    if len(time) > 1:
        rows.append(len(time) - 1)

    print("")
    print(name)
    print("row  time  force_x  force_y  force_z")

    for row in rows:
        print(
            f"{row:3d}"
            f"  {time[row]:.10g}"
            f"  {force[row, 0]: .10e}"
            f"  {force[row, 1]: .10e}"
            f"  {force[row, 2]: .10e}"
        )


print("")
print("Rows plotted")
print("New motion:", len(new_motion_time))
print("Old motion:", len(old_motion_time))
print("New force: ", len(new_force_time))
print("Old force: ", len(old_force_time))

print_position_summary("New BlockEigen centre-of-rotation", new_motion_time, new_position)
print_position_summary("Old coupled-solver centre-of-rotation", old_motion_time, old_position)
print_force_summary("New force", new_force_time, new_force)
print_force_summary("Old force", old_force_time, old_force)


# =============================================================
# PLOT DATA
# =============================================================

def set_informative_y_ticks(ax, values):
    finite_values = values[np.isfinite(values)]

    if finite_values.size == 0:
        return

    y_min = np.min(finite_values)
    y_max = np.max(finite_values)

    if np.isclose(y_min, y_max):
        padding = max(abs(y_min)*0.1, 1.0e-6)
    else:
        padding = 0.08*(y_max - y_min)

    ax.set_ylim(y_min - padding, y_max + padding)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def finish_plot(ax, values):
    set_informative_y_ticks(ax, values)
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


def plot_displacement_component(label, component):
    new_displacement = new_position - REFERENCE_CENTRE_OF_ROTATION
    old_displacement = old_position - REFERENCE_CENTRE_OF_ROTATION

    plt.figure(figsize=(8, 4))
    ax = plt.gca()

    ax.plot(
        old_motion_time,
        old_displacement[:, component],
        "b-",
        linewidth=1.8,
        label="old coupled solver",
    )
    ax.plot(
        new_motion_time,
        new_displacement[:, component],
        "r--",
        linewidth=1.5,
        label="new BlockEigen force coupling",
    )

    ax.set_xlabel("Time (s)")
    ax.set_ylabel(f"{label} displacement (m)")
    ax.set_title(f"Centre-of-rotation {label} displacement")

    combined_values = np.concatenate(
        (old_displacement[:, component], new_displacement[:, component])
    )
    finish_plot(ax, combined_values)
    save_or_show(f"rigid_body_{label.lower()}_displacement_new_vs_old.png")


def plot_force_component(label, component):
    plt.figure(figsize=(8, 4))
    ax = plt.gca()

    ax.plot(
        old_force_time,
        old_force[:, component],
        "b-",
        linewidth=1.8,
        label="old coupled solver",
    )
    ax.plot(
        new_force_time,
        new_force[:, component],
        "r--",
        linewidth=1.5,
        label="new BlockEigen force coupling",
    )

    ax.set_xlabel("Time (s)")
    ax.set_ylabel(f"{label} {force_ylabel_prefix} (N)")
    ax.set_title(f"{force_ylabel_prefix}: {label} component")

    combined_values = np.concatenate((old_force[:, component], new_force[:, component]))
    finish_plot(ax, combined_values)
    save_or_show(f"{force_ylabel_prefix.replace(' ', '_')}_{label.lower()}_new_vs_old.png")


for component_label, component_index in (("X", 0), ("Y", 1), ("Z", 2)):
    plot_displacement_component(component_label, component_index)

for component_label, component_index in (("X", 0), ("Y", 1), ("Z", 2)):
    plot_force_component(component_label, component_index)
