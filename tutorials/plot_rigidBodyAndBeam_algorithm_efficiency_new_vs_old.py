#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare computational efficiency of the new BlockEigen-coupled framework
against the original coupled-solver framework.

Run from Spyder or from the beamFoam/tutorials directory after both cases have
completed and contain log.interFoam files.
"""

import glob
import os
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MaxNLocator


# =============================================================
# INPUT FILES
# =============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

NEW_CASE = os.path.join(SCRIPT_DIR, "rigidBodyAndBeam_fullTime")
OLD_CASE = os.path.join(SCRIPT_DIR, "rigidBodyAndBeam_compareSixDoF_Eigen")

NEW_LABEL = "New Solver"
OLD_LABEL = "Original Solver"


def find_one(pattern, description):
    matches = sorted(glob.glob(pattern))

    if not matches:
        raise FileNotFoundError(f"No {description} file found for {pattern}")

    return matches[0]


NEW_LOG_FILE = find_one(os.path.join(NEW_CASE, "log.interFoam"), "new log")
OLD_LOG_FILE = find_one(os.path.join(OLD_CASE, "log.interFoam"), "old log")

print("Input files")
print("New framework log:     ", NEW_LOG_FILE)
print("Original framework log:", OLD_LOG_FILE)


# =============================================================
# USER SETTINGS
# =============================================================

T_START = None
T_END = None

SAVE_FIGURES = False
FIGURE_DIR = os.path.join(SCRIPT_DIR, "algorithm_efficiency_plots")


# =============================================================
# READ DATA
# =============================================================

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
TIME_PATTERN = re.compile(rf"^\s*Time\s*=\s*({NUMBER})\s*$")
EXECUTION_PATTERN = re.compile(
    rf"^\s*ExecutionTime\s*=\s*({NUMBER})\s*s\s+ClockTime\s*=\s*({NUMBER})\s*s"
)
CONVERGENCE_PATTERN = re.compile(r"^\s*(\d+)\s*:\s*(?:Converged|Failed)\b")
LINEAR_SOLVER_PATTERN = re.compile(
    r"Solving for\s+([^,]+),.*No Iterations\s+(\d+)"
)
PIMPLE_PATTERN = re.compile(r"^\s*PIMPLE:\s*iteration\s+(\d+)")


def read_optional_beam_convergence(case_dir):
    pattern = os.path.join(
        case_dir,
        "postProcessing",
        "*",
        "beamConvergenceData.dat",
    )
    matches = sorted(glob.glob(pattern))

    if not matches:
        return None

    filename = matches[0]
    data = np.genfromtxt(filename, names=True)

    if data.size == 0:
        return None

    if data.ndim == 0:
        data = np.array([data], dtype=data.dtype)

    names = data.dtype.names
    if names is None or len(names) < 2:
        raise ValueError(f"Could not read convergence columns from {filename}")

    return {
        "filename": filename,
        "time": np.asarray(data[names[0]], dtype=float),
        "newton_iterations": np.asarray(data[names[1]], dtype=float),
    }


def match_values_by_time(target_time, source_time, source_values):
    matched = np.full_like(target_time, np.nan, dtype=float)

    for i, time_value in enumerate(target_time):
        nearest = np.argmin(np.abs(source_time - time_value))

        if np.isclose(source_time[nearest], time_value, rtol=1.0e-8, atol=1.0e-10):
            matched[i] = source_values[nearest]

    return matched


def parse_solver_log(filename):
    records = []
    all_solver_fields = set()

    current_time = None
    current_newton_iterations = np.nan
    current_linear_iterations = 0
    current_linear_by_field = defaultdict(int)
    current_pimple_iterations = 0

    with open(filename, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            time_match = TIME_PATTERN.search(line)
            if time_match:
                current_time = float(time_match.group(1))
                current_newton_iterations = np.nan
                current_linear_iterations = 0
                current_linear_by_field = defaultdict(int)
                current_pimple_iterations = 0
                continue

            if current_time is None:
                continue

            convergence_match = CONVERGENCE_PATTERN.search(line)
            if convergence_match:
                current_newton_iterations = float(convergence_match.group(1))
                continue

            pimple_match = PIMPLE_PATTERN.search(line)
            if pimple_match:
                current_pimple_iterations = max(
                    current_pimple_iterations,
                    int(pimple_match.group(1)),
                )
                continue

            linear_solver_match = LINEAR_SOLVER_PATTERN.search(line)
            if linear_solver_match:
                field = linear_solver_match.group(1).strip()
                n_iterations = int(linear_solver_match.group(2))
                current_linear_iterations += n_iterations
                current_linear_by_field[field] += n_iterations
                all_solver_fields.add(field)
                continue

            execution_match = EXECUTION_PATTERN.search(line)
            if execution_match:
                records.append(
                    {
                        "time": current_time,
                        "newton_iterations": current_newton_iterations,
                        "linear_iterations": current_linear_iterations,
                        "linear_by_field": dict(current_linear_by_field),
                        "pimple_iterations": current_pimple_iterations,
                        "execution_time": float(execution_match.group(1)),
                        "clock_time": float(execution_match.group(2)),
                    }
                )
                current_time = None

    if not records:
        raise ValueError(f"No time-step records found in {filename}")

    time = np.asarray([record["time"] for record in records], dtype=float)
    execution_time = np.asarray(
        [record["execution_time"] for record in records],
        dtype=float,
    )
    clock_time = np.asarray(
        [record["clock_time"] for record in records],
        dtype=float,
    )
    newton_iterations = np.asarray(
        [record["newton_iterations"] for record in records],
        dtype=float,
    )
    linear_iterations = np.asarray(
        [record["linear_iterations"] for record in records],
        dtype=float,
    )
    pimple_iterations = np.asarray(
        [record["pimple_iterations"] for record in records],
        dtype=float,
    )

    step_execution_time = np.diff(np.concatenate(([0.0], execution_time)))
    step_clock_time = np.diff(np.concatenate(([0.0], clock_time)))

    linear_by_field = {}
    for field in sorted(all_solver_fields):
        linear_by_field[field] = np.asarray(
            [
                record["linear_by_field"].get(field, 0)
                for record in records
            ],
            dtype=float,
        )

    return {
        "filename": filename,
        "time": time,
        "newton_iterations": newton_iterations,
        "linear_iterations": linear_iterations,
        "linear_by_field": linear_by_field,
        "pimple_iterations": pimple_iterations,
        "execution_time": execution_time,
        "clock_time": clock_time,
        "step_execution_time": step_execution_time,
        "step_clock_time": step_clock_time,
    }


def apply_convergence_file_if_present(case_dir, case_data):
    convergence_data = read_optional_beam_convergence(case_dir)

    if convergence_data is None:
        return case_data

    matched = match_values_by_time(
        case_data["time"],
        convergence_data["time"],
        convergence_data["newton_iterations"],
    )
    replace = np.isfinite(matched)

    if np.any(replace):
        case_data = dict(case_data)
        case_data["newton_iterations"] = np.asarray(
            case_data["newton_iterations"],
            dtype=float,
        )
        case_data["newton_iterations"][replace] = matched[replace]
        case_data["beam_convergence_file"] = convergence_data["filename"]

    return case_data


def apply_time_window(case_data):
    mask = np.ones_like(case_data["time"], dtype=bool)

    if T_START is not None:
        mask &= case_data["time"] >= T_START

    if T_END is not None:
        mask &= case_data["time"] <= T_END

    windowed = {}

    for key, value in case_data.items():
        if isinstance(value, np.ndarray) and value.shape == mask.shape:
            windowed[key] = value[mask]
        elif key == "linear_by_field":
            windowed[key] = {
                field: values[mask]
                for field, values in value.items()
            }
        else:
            windowed[key] = value

    if windowed["time"].size == 0:
        raise ValueError("The selected time window contains no time steps")

    return windowed


def add_derived_metrics(case_data):
    case_data = dict(case_data)
    newton_iterations = case_data["newton_iterations"]
    step_execution_time = case_data["step_execution_time"]
    step_clock_time = case_data["step_clock_time"]

    with np.errstate(divide="ignore", invalid="ignore"):
        case_data["execution_time_per_newton"] = np.where(
            newton_iterations > 0,
            step_execution_time / newton_iterations,
            np.nan,
        )

    case_data["cumulative_newton_iterations"] = np.cumsum(
        np.nan_to_num(newton_iterations, nan=0.0)
    )
    case_data["cumulative_step_execution_time"] = np.cumsum(
        step_execution_time
    )
    case_data["cumulative_step_clock_time"] = np.cumsum(
        case_data["step_clock_time"]
    )
    return case_data


new_data = parse_solver_log(NEW_LOG_FILE)
old_data = parse_solver_log(OLD_LOG_FILE)

new_data = apply_convergence_file_if_present(NEW_CASE, new_data)
old_data = apply_convergence_file_if_present(OLD_CASE, old_data)

new_data = add_derived_metrics(apply_time_window(new_data))
old_data = add_derived_metrics(apply_time_window(old_data))

if np.all(~np.isfinite(new_data["newton_iterations"])):
    raise ValueError(f"No Newton iteration counts found in {NEW_LOG_FILE}")

if np.all(~np.isfinite(old_data["newton_iterations"])):
    raise ValueError(f"No Newton iteration counts found in {OLD_LOG_FILE}")


# =============================================================
# PRINT CHECKS
# =============================================================

def finite_sum(values):
    return float(np.nansum(values))


def finite_mean(values):
    finite_values = values[np.isfinite(values)]

    if finite_values.size == 0:
        return np.nan

    return float(np.mean(finite_values))


def finite_max(values):
    finite_values = values[np.isfinite(values)]

    if finite_values.size == 0:
        return np.nan

    return float(np.max(finite_values))


def total_simulation_time(case_data):
    return float(case_data["time"][-1] - case_data["time"][0])


def total_runtime(case_data, runtime_key):
    return float(case_data[runtime_key][-1])


def print_time_taken():
    print("")
    print("Time taken")
    print(
        f"{OLD_LABEL}: ExecutionTime = "
        f"{total_runtime(old_data, 'execution_time'):.6g} s, "
        f"ClockTime = {total_runtime(old_data, 'clock_time'):.6g} s"
    )
    print(
        f"{NEW_LABEL}:      ExecutionTime = "
        f"{total_runtime(new_data, 'execution_time'):.6g} s, "
        f"ClockTime = {total_runtime(new_data, 'clock_time'):.6g} s"
    )

    old_execution_time = total_runtime(old_data, "execution_time")
    new_execution_time = total_runtime(new_data, "execution_time")

    if new_execution_time > 0:
        print(
            "ExecutionTime speed-up "
            f"({OLD_LABEL}/{NEW_LABEL}) = "
            f"{old_execution_time/new_execution_time:.6g}x"
        )


def summary_metrics(case_data):
    execution_time = total_runtime(case_data, "execution_time")
    simulation_time = float(case_data["time"][-1])

    return {
        "time steps": float(case_data["time"].size),
        "last simulated time (s)": simulation_time,
        "ExecutionTime (s)": execution_time,
        "ClockTime (s)": total_runtime(case_data, "clock_time"),
        "mean step ExecutionTime (s)": finite_mean(
            case_data["step_execution_time"]
        ),
        "max step ExecutionTime (s)": finite_max(
            case_data["step_execution_time"]
        ),
        "total Newton iterations": finite_sum(
            case_data["newton_iterations"]
        ),
        "mean Newton iterations": finite_mean(
            case_data["newton_iterations"]
        ),
        "max Newton iterations": finite_max(
            case_data["newton_iterations"]
        ),
        "total linear solver iterations": finite_sum(
            case_data["linear_iterations"]
        ),
        "mean linear solver iterations": finite_mean(
            case_data["linear_iterations"]
        ),
        "ExecutionTime / simulated s": (
            execution_time / simulation_time
            if simulation_time > 0
            else np.nan
        ),
        "ExecutionTime / Newton iteration": (
            execution_time / finite_sum(case_data["newton_iterations"])
            if finite_sum(case_data["newton_iterations"]) > 0
            else np.nan
        ),
    }


def print_summary_table():
    old_metrics = summary_metrics(old_data)
    new_metrics = summary_metrics(new_data)

    print("")
    print("Computational efficiency summary")
    print(
        f"{'metric':36s}"
        f"{OLD_LABEL:>20s}"
        f"{NEW_LABEL:>20s}"
        f"{'new/original':>18s}"
    )

    for metric in old_metrics:
        old_value = old_metrics[metric]
        new_value = new_metrics[metric]

        if np.isfinite(old_value) and old_value != 0.0:
            ratio = new_value / old_value
        else:
            ratio = np.nan

        print(
            f"{metric:36s}"
            f"{old_value:20.6g}"
            f"{new_value:20.6g}"
            f"{ratio:18.6g}"
        )


print("")
print("Rows plotted")
print(f"{NEW_LABEL}: {new_data['time'].size}")
print(f"{OLD_LABEL}: {old_data['time'].size}")

if "beam_convergence_file" in new_data:
    print("New convergence data:     ", new_data["beam_convergence_file"])

if "beam_convergence_file" in old_data:
    print("Original convergence data:", old_data["beam_convergence_file"])

print_time_taken()
print_summary_table()


# =============================================================
# PLOT DATA
# =============================================================

def set_informative_y_ticks(ax, values, integer=False, bottom=None):
    finite_values = values[np.isfinite(values)]

    if finite_values.size == 0:
        return

    y_min = np.min(finite_values)
    y_max = np.max(finite_values)

    if np.isclose(y_min, y_max):
        padding = max(abs(y_min)*0.1, 1.0e-6)
    else:
        padding = 0.08*(y_max - y_min)

    if bottom is not None:
        ax.set_ylim(bottom=bottom)
    else:
        ax.set_ylim(y_min - padding, y_max + padding)

    ax.yaxis.set_major_locator(MaxNLocator(nbins=7, integer=integer))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def finish_plot(ax, values, integer_y=False, bottom=None):
    set_informative_y_ticks(ax, values, integer=integer_y, bottom=bottom)
    ax.grid(which="major", ls="--", lw=0.6, alpha=0.75)
    ax.grid(which="minor", ls=":", lw=0.4, alpha=0.45)
    ax.legend()


def save_current_figure(filename):
    if SAVE_FIGURES:
        os.makedirs(FIGURE_DIR, exist_ok=True)
        path = os.path.join(FIGURE_DIR, filename)
        plt.savefig(path, dpi=300)
        print("Saved:", path)


def combined_values(key):
    return np.concatenate((old_data[key], new_data[key]))


def plot_line_comparison(
    key,
    ylabel,
    filename,
    title=None,
    integer_y=False,
    bottom=None,
):
    plt.figure(figsize=(8, 4))
    ax = plt.gca()

    ax.plot(
        old_data["time"],
        old_data[key],
        "b-",
        linewidth=1.8,
        label=OLD_LABEL,
    )
    ax.plot(
        new_data["time"],
        new_data[key],
        "r--",
        linewidth=1.5,
        label=NEW_LABEL,
    )

    ax.set_xlabel("Simulation time (s)")
    ax.set_ylabel(ylabel)

    if title is not None:
        ax.set_title(title)

    finish_plot(
        ax,
        combined_values(key),
        integer_y=integer_y,
        bottom=bottom,
    )
    plt.tight_layout()
    save_current_figure(filename)


def plot_newton_iterations():
    plt.figure(figsize=(8, 4))
    ax = plt.gca()

    ax.plot(
        old_data["time"],
        old_data["newton_iterations"],
        "bo-",
        linewidth=1.3,
        markersize=3.0,
        label=OLD_LABEL,
    )
    ax.plot(
        new_data["time"],
        new_data["newton_iterations"],
        "rs--",
        linewidth=1.2,
        markersize=3.0,
        label=NEW_LABEL,
    )

    ax.set_xlabel("Simulation time (s)")
    ax.set_ylabel("Newton iterations per time step")
    finish_plot(
        ax,
        combined_values("newton_iterations"),
        integer_y=True,
        bottom=0,
    )
    plt.tight_layout()
    save_current_figure("newton_iterations_per_time_step_new_vs_old.png")


def plot_total_time_bars():
    labels = [OLD_LABEL, NEW_LABEL]
    execution_values = [
        total_runtime(old_data, "execution_time"),
        total_runtime(new_data, "execution_time"),
    ]
    clock_values = [
        total_runtime(old_data, "clock_time"),
        total_runtime(new_data, "clock_time"),
    ]

    x = np.arange(len(labels))
    width = 0.36

    plt.figure(figsize=(7, 4))
    ax = plt.gca()
    ax.bar(x - width/2, execution_values, width, label="ExecutionTime")
    ax.bar(x + width/2, clock_values, width, label="ClockTime")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Elapsed time (s)")
    finish_plot(ax, np.asarray(execution_values + clock_values), bottom=0)
    plt.tight_layout()
    save_current_figure("total_elapsed_time_new_vs_old.png")


def plot_runtime_subplots():
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    for ax, key, ylabel in (
        (axes[0], "step_execution_time", "Step ExecutionTime (s)"),
        (axes[1], "step_clock_time", "Step ClockTime (s)"),
    ):
        ax.plot(
            old_data["time"],
            old_data[key],
            "b-",
            linewidth=1.8,
            label=OLD_LABEL,
        )
        ax.plot(
            new_data["time"],
            new_data[key],
            "r--",
            linewidth=1.5,
            label=NEW_LABEL,
        )
        ax.set_ylabel(ylabel)
        finish_plot(ax, combined_values(key), bottom=0)

    axes[1].set_xlabel("Simulation time (s)")
    plt.tight_layout()
    save_current_figure("per_step_runtime_new_vs_old.png")


def plot_cumulative_runtime_subplots():
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    for ax, key, ylabel in (
        (axes[0], "execution_time", "Cumulative ExecutionTime (s)"),
        (axes[1], "clock_time", "Cumulative ClockTime (s)"),
    ):
        ax.plot(
            old_data["time"],
            old_data[key],
            "b-",
            linewidth=1.8,
            label=OLD_LABEL,
        )
        ax.plot(
            new_data["time"],
            new_data[key],
            "r--",
            linewidth=1.5,
            label=NEW_LABEL,
        )
        ax.set_ylabel(ylabel)
        finish_plot(ax, combined_values(key), bottom=0)

    axes[1].set_xlabel("Simulation time (s)")
    plt.tight_layout()
    save_current_figure("cumulative_runtime_new_vs_old.png")


def plot_linear_solver_field_breakdown(case_data, label, filename):
    fields = sorted(case_data["linear_by_field"])

    if not fields:
        return

    plt.figure(figsize=(8, 4))
    ax = plt.gca()

    bottom_values = np.zeros_like(case_data["time"], dtype=float)
    for field in fields:
        values = case_data["linear_by_field"][field]
        ax.fill_between(
            case_data["time"],
            bottom_values,
            bottom_values + values,
            step="mid",
            alpha=0.65,
            label=field,
        )
        bottom_values += values

    ax.set_xlabel("Simulation time (s)")
    ax.set_ylabel("Linear solver iterations per time step")
    ax.set_title(label)
    finish_plot(ax, bottom_values, integer_y=True, bottom=0)
    plt.tight_layout()
    save_current_figure(filename)


plot_newton_iterations()
plot_line_comparison(
    "cumulative_newton_iterations",
    "Cumulative Newton iterations",
    "cumulative_newton_iterations_new_vs_old.png",
    integer_y=True,
    bottom=0,
)
plot_runtime_subplots()
plot_cumulative_runtime_subplots()
plot_line_comparison(
    "execution_time_per_newton",
    "Step ExecutionTime per Newton iteration (s)",
    "execution_time_per_newton_iteration_new_vs_old.png",
    bottom=0,
)
plot_line_comparison(
    "linear_iterations",
    "Linear solver iterations per time step",
    "linear_solver_iterations_per_time_step_new_vs_old.png",
    integer_y=True,
    bottom=0,
)
plot_total_time_bars()
plot_linear_solver_field_breakdown(
    old_data,
    OLD_LABEL,
    "linear_solver_field_breakdown_original_solver.png",
)
plot_linear_solver_field_breakdown(
    new_data,
    NEW_LABEL,
    "linear_solver_field_breakdown_new_solver.png",
)

if SAVE_FIGURES:
    print("Saved figures to:", FIGURE_DIR)

plt.show()
