#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare BlockEigen and sixDoF centre-of-rotation positions.

Run from Spyder or from this case directory.
"""

import os
import re
import glob

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator


# =============================================================
# INPUT FILES
# =============================================================

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
POST_DIR = os.path.join(CASE_DIR, "postProcessing")

BLOCKEIGEN_FILE = sorted(
    glob.glob(
        os.path.join(
            POST_DIR,
            "BlockEigenSolve_rigidBodyMotion",
            "*",
            "BlockEigenSolve_rigidBodyMotion.dat",
        )
    )
)[0]

SIXDOF_FILE = sorted(
    glob.glob(
        os.path.join(
            POST_DIR,
            "sixDoF*",
            "*",
            "sixDoFRigidBodyStateFvBeam.dat",
        )
    )
)[0]

print("Input files")
print("BlockEigen:", BLOCKEIGEN_FILE)
print("sixDoF:    ", SIXDOF_FILE)


# =============================================================
# USER SETTINGS
# =============================================================

# BlockEigen writes displacement. Add this reference point to convert to
# centre-of-rotation position.
REFERENCE_CENTRE_OF_ROTATION = np.array([2.6565, 0.0, 0.0])

# BlockEigen writes one row for every nonlinear solve. Usually the last row at
# each time is the useful converged value.
BLOCKEIGEN_TIME_REDUCTION = "last"  # "last" or "all"

# Optional time window. Set both to None for the full run.
T_START = None
T_END = None

SAVE_FIGURES = False
FIGURE_DIR = "plots"


# =============================================================
# READ DATA
# =============================================================

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

    keep = []
    for t in np.unique(time):
        matching_rows = np.where(time == t)[0]
        keep.append(matching_rows[-1])

    keep = np.asarray(keep, dtype=int)
    return time[keep], position[keep]


def read_sixdof_centre_of_rotation(filename):
    number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
    vector = rf"\(\s*({number})\s+({number})\s+({number})\s*\)"

    # Columns:
    # Time centreOfRotation centreOfMass rotation velocity omega
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

    with open(filename, "r") as handle:
        for line in handle:
            match = line_pattern.search(line)
            if not match:
                continue

            values = [float(value) for value in match.groups()]
            time.append(values[0])

            # centreOfRotation is the first vector after time.
            centre_of_rotation.append(values[1:4])

    if not time:
        raise ValueError(f"No motion rows found in {filename}")

    return np.asarray(time), np.asarray(centre_of_rotation)


def apply_time_window(time, values):
    mask = np.ones_like(time, dtype=bool)

    if T_START is not None:
        mask &= time >= T_START

    if T_END is not None:
        mask &= time <= T_END

    return time[mask], values[mask]


block_time, block_position = read_blockeigen_position(BLOCKEIGEN_FILE)
sixdof_time, sixdof_position = read_sixdof_centre_of_rotation(SIXDOF_FILE)

block_time, block_position = apply_time_window(block_time, block_position)
sixdof_time, sixdof_position = apply_time_window(sixdof_time, sixdof_position)

print("")
print("Rows plotted")
print("BlockEigen:", len(block_time))
print("sixDoF:    ", len(sixdof_time))


def print_xyz_table(title, time, position):
    displacement = position - REFERENCE_CENTRE_OF_ROTATION

    print("")
    print(title)
    print(
        "row  time"
        "  x_position  y_position  z_position"
        "  x_displacement  y_displacement  z_displacement"
    )

    rows_to_print = [0]

    if len(time) > 1:
        rows_to_print.append(len(time) - 1)

    for row in rows_to_print:
        print(
            f"{row:3d}"
            f"  {time[row]:.10g}"
            f"  {position[row, 0]: .10e}"
            f"  {position[row, 1]: .10e}"
            f"  {position[row, 2]: .10e}"
            f"  {displacement[row, 0]: .10e}"
            f"  {displacement[row, 1]: .10e}"
            f"  {displacement[row, 2]: .10e}"
        )


print_xyz_table("BlockEigen centre-of-rotation position/displacement", block_time, block_position)
print_xyz_table("sixDoF centre-of-rotation position/displacement", sixdof_time, sixdof_position)


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


def plot_position_component(label, component):
    plt.figure(figsize=(8, 4))
    ax = plt.gca()

    plt.plot(
        sixdof_time,
        sixdof_position[:, component],
        "b-",
        linewidth=1.8,
        label="sixDoF centreOfRotation",
    )

    plt.plot(
        block_time,
        block_position[:, component],
        "r--",
        linewidth=1.5,
        label="BlockEigen centreOfRotation",
    )

    plt.xlabel("Time (s)")
    plt.ylabel(f"{label} position (m)")
    plt.title(f"Centre of rotation: {label} position")
    plt.legend()

    combined_values = np.concatenate(
        (
            sixdof_position[:, component],
            block_position[:, component],
        )
    )
    set_informative_y_ticks(ax, combined_values)

    ax.grid(which="major", ls="--", lw=0.6, alpha=0.75)
    ax.grid(which="minor", ls=":", lw=0.4, alpha=0.45)

    plt.tight_layout()

    if SAVE_FIGURES:
        plot_dir = os.path.join(CASE_DIR, FIGURE_DIR)
        os.makedirs(plot_dir, exist_ok=True)

        figure_path = os.path.join(
            plot_dir,
            f"centre_of_rotation_{label.lower()}_position_comparison.png",
        )
        plt.savefig(figure_path, dpi=300)
        print("Saved:", figure_path)

    plt.show()


plot_position_component("X", 0)
plot_position_component("Y", 1)
plot_position_component("Z", 2)
