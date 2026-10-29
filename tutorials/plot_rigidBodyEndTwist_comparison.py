#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare the partitioned and monolithic rigidBodyEnd couplings over the
rigidBodyEndTwist sweeps: a rod twisted by a rigid body welded to its end.

beamFoam-only sweeps (rigidBodyEndTwist/moment and rigidBodyEndTwist/stiffness,
run with rigidBodyEndTwist/Allrun):
    - twist histories, partitioned against monolithic, one panel per case
    - iterations per time step (beam Newton; partitioned coupling)
    - settled twist error against M L/(G Ip), and the partitioned-monolithic
      difference
    - run time
    - grouped bars of total and average beam Newton iterations against
      applied moment M (N m) and torsional stiffness kt (N m/rad), including
      partial runs (hatched bars)
The moment sweep is plotted against the static twist, the stiffness sweep
against the partitioned stability number beta dt^2 kt/J. A cross marks the
point a run crashed at.

moorFV cases (rigidBodyEndTwist/moorFV, run with its Allrun): the free
torsional oscillation against theory for FvBeamNewmark (partitioned, with
3, 8 and 20 outer correctors where run) and beamFoamCoupled (monolithic), and
the peak twist in each half period.

Beam Newton iterations (nNewton, BlockEigenSolverOF::solve calls per time
step, summed over outer correctors and coupling iterations) come from
postProcessing/0/beamSolverIterations.dat (beamSolverIterations function
object) where a case has it, else from the log. Plotted for both:
    - nNewton per time step against time, one panel per sweep case
    - average nNewton per time step against the sweep parameter
    - moorFV: nNewton per time step and its running average against time

Run from Spyder or from the beamFoam/tutorials directory after the cases have
completed. Uses the readers in rigidBodyEndTwist/summarise.py and
rigidBodyEndTwist/moorFV/compare.py.
"""

import math
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MaxNLocator


# =============================================================
# INPUT FILES
# =============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TWIST_DIR = os.path.join(SCRIPT_DIR, "rigidBodyEndTwist")
MOORFV_DIR = os.path.join(TWIST_DIR, "moorFV")

sys.path.insert(0, SCRIPT_DIR)
sys.path.insert(0, TWIST_DIR)
sys.path.insert(0, MOORFV_DIR)
import summarise  # noqa: E402  (rigidBodyEndTwist/summarise.py)
import compare    # noqa: E402  (rigidBodyEndTwist/moorFV/compare.py)
from read_beam_solver_iterations import read_beam_solver_iterations  # noqa: E402

COUPLINGS = ["partitioned", "monolithic"]

# moorFV cases: (directory in rigidBodyEndTwist/moorFV, legend label, coupling)
# Cases that have not been run are skipped
MOORFV_CASES = [
    ("monolithic", "Monolithic, 3 outer correctors", "monolithic"),
    ("partitioned", "Partitioned, 3 outer correctors", "partitioned"),
    ("partitioned_nOuter8", "Partitioned, 8 outer correctors", "partitioned"),
    ("partitioned_nOuter20", "Partitioned, 20 outer correctors", "partitioned"),
]


# =============================================================
# USER SETTINGS
# =============================================================

SAVE_FIGURES = False
FIGURE_DIR = os.path.join(SCRIPT_DIR, "comparison_plots", "rigidBodyEndTwist")

# One colour per coupling, the same in every figure. Partitioned runs with
# more outer correctors use the same colour with a different line style
COLOURS = {"monolithic": "#2a78d6", "partitioned": "#eb6834"}
THEORY_COLOUR = "0.35"
PARTITIONED_STYLES = ["-", "--", ":"]


# =============================================================
# READ DATA
# =============================================================

def sweep_labels(sweep):
    root = os.path.join(TWIST_DIR, sweep)
    if not os.path.isdir(root):
        return []
    return sorted(
        os.listdir(root),
        key=lambda s: float(re.sub(r"^[A-Za-z]+", "", s)),
    )


def newton_counts(case, log):
    """beamSolverIterations output (or None), total and average nNewton per
    time step: from the function object if the case has it, else the log"""
    solver = read_beam_solver_iterations(case)
    if solver is not None and solver["time"].size:
        return solver, solver["n_newton"].sum(), solver["n_newton"].mean()
    steps = log["steps"]
    return None, log["newton"], log["newton"]/steps if steps else math.nan


def load_sweep(sweep):
    """List of {label, partitioned, monolithic, difference} for a sweep"""
    rows = []
    for label in sweep_labels(sweep):
        row = {"label": label}
        for coupling in COUPLINGS:
            case = os.path.join(TWIST_DIR, sweep, label, coupling)
            if os.path.isdir(case):
                s = summarise.summarise_case(case)
                body = s["body"]
                s["time"] = np.array([r[0] for r in body])
                s["twist"] = np.array([r[1] for r in body])
                s["coupling_iterations"] = np.array([r[5] for r in body])
                s["solver"], s["newton_total"], s["newton_per_step"] = (
                    newton_counts(case, s["log"])
                )
                s["finished"] = s["log"]["status"] == "finished"
                row[coupling] = s
        if all(c in row for c in COUPLINGS):
            row["difference"] = (
                summarise.twist_difference(
                    row["partitioned"]["body"], row["monolithic"]["body"]
                )
                if row["partitioned"]["finished"] and row["monolithic"]["finished"]
                else math.nan
            )
            rows.append(row)
    return rows


def static_twist(s):
    return s["p"]["M"]/s["p"]["kt"]


def load_moorfv():
    cases = []
    for name, label, coupling in MOORFV_CASES:
        case = os.path.join(MOORFV_DIR, name)
        if not os.path.isdir(case):
            continue
        twist = compare.read_twist(case)
        if not twist:
            continue
        time = np.array(sorted(twist))
        log = compare.read_log(case)
        solver, newton_total, newton_per_step = newton_counts(case, log)
        cases.append({
            "label": label,
            "coupling": coupling,
            "time": time,
            "twist": np.array([twist[t] for t in time]),
            "log": log,
            "solver": solver,
            "newton_total": newton_total,
            "newton_per_step": newton_per_step,
        })
    return cases


moment = load_sweep("moment")
stiffness = load_sweep("stiffness")
moorfv = load_moorfv()

print("Cases read")
print(f"moment sweep: {len(moment)} cases")
print(f"stiffness sweep: {len(stiffness)} cases")
print(f"moorFV: {', '.join(c['label'] for c in moorfv) or 'none'}")


# =============================================================
# PLOTTING HELPERS
# =============================================================

def finish_plot(ax, ylabel, xlabel, legend=True, log_y=False):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if not log_y:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(which="major", ls="--", lw=0.6, alpha=0.75)
    ax.grid(which="minor", ls=":", lw=0.4, alpha=0.45)
    if legend:
        ax.legend()
    plt.tight_layout()


def save_or_show(filename):
    if SAVE_FIGURES:
        os.makedirs(FIGURE_DIR, exist_ok=True)
        path = os.path.join(FIGURE_DIR, filename)
        plt.savefig(path, dpi=300)
        print("Saved:", path)

    plt.show()


def coupling_label(coupling):
    return coupling.capitalize()


def plot_finished(ax, x, y, rows, coupling, **kwargs):
    """Line through the runs that finished; crashed runs are left out of the
    line (see mark_crashes)"""
    keep = [i for i, row in enumerate(rows) if row[coupling]["finished"]]
    ax.plot([x[i] for i in keep], [y[i] for i in keep], **kwargs)


def mark_crashes(ax, x, y, rows, coupling):
    """Cross at the x position of each crashed run, on the curve's last value"""
    for xi, yi, row in zip(x, y, rows):
        if not row[coupling]["finished"]:
            ax.plot(
                xi, yi, marker="x", markersize=10, markeredgewidth=2,
                color=COLOURS[coupling], linestyle="none",
            )


# =============================================================
# SWEEP PLOTS
# =============================================================

def plot_histories(rows, sweep, title_of, filename):
    """Twist over the static twist against time, one panel per case"""
    n = len(rows)
    if n == 0:
        return
    ncols = 3
    nrows = math.ceil(n/ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(11, 2.8*nrows), sharex=True, squeeze=False
    )

    for ax, row in zip(axes.flat, rows):
        reference = static_twist(row["monolithic"])

        # y range from the monolithic run, so a partitioned run that blows up
        # does not flatten the panel
        mono = row["monolithic"]["twist"]/reference
        low, high = min(mono.min(), 0.0), max(mono.max(), 1.0)
        margin = 0.1*(high - low)
        low, high = low - margin, high + margin

        # Monolithic underneath, partitioned dashed on top: where they agree
        # both stay visible
        for coupling in ["monolithic", "partitioned"]:
            s = row[coupling]
            if len(s["time"]) == 0:
                continue
            ax.plot(
                s["time"], s["twist"]/reference, color=COLOURS[coupling],
                linewidth=2.2 if coupling == "monolithic" else 1.1,
                linestyle="-" if coupling == "monolithic" else (0, (4, 3)),
                label=coupling_label(coupling),
            )
            if not s["finished"]:
                t_crash = s["time"][-1]
                y_crash = min(max(s["twist"][-1]/reference, low), high)
                ax.plot(
                    t_crash, y_crash, marker="x", markersize=9,
                    markeredgewidth=2, color=COLOURS[coupling],
                    clip_on=False,
                )
                ax.text(
                    0.97, 0.95, f"x: {coupling} crashed at t = {t_crash:.3g} s",
                    transform=ax.transAxes, fontsize=8, ha="right",
                    va="top", color=COLOURS[coupling],
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
                )
        ax.set_ylim(low, high)
        ax.axhline(1, color=THEORY_COLOUR, ls=":", lw=1.0, label="Static twist")
        ax.set_title(title_of(row), fontsize=10)
        ax.grid(which="major", ls="--", lw=0.6, alpha=0.75)

    for ax in axes.flat[n:]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Time (s)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Twist / static twist")

    axes.flat[0].legend(fontsize=8)
    fig.suptitle(f"{sweep.capitalize()} sweep: twist histories")
    fig.tight_layout()
    save_or_show(filename)


def plot_iterations(rows, x, xlabel, log_x, sweep, filename, stability_line=False):
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()

    for coupling in COUPLINGS:
        y = [row[coupling]["newton_per_step"] for row in rows]
        plot_finished(
            ax, x, y, rows, coupling, color=COLOURS[coupling], marker="o",
            markersize=6, linewidth=1.6,
            label=f"{coupling_label(coupling)}: beam Newton iterations",
        )
        mark_crashes(ax, x, y, rows, coupling)

    y = [row["partitioned"]["coupling_mean"] for row in rows]
    plot_finished(
        ax, x, y, rows, "partitioned", color=COLOURS["partitioned"], marker="s",
        markersize=6, linewidth=1.2, linestyle="--",
        label="Partitioned: coupling iterations",
    )
    mark_crashes(ax, x, y, rows, "partitioned")

    if stability_line:
        ax.axvline(1, color=THEORY_COLOUR, ls=":", lw=1.0)
        ax.text(
            1.1, 0.3, "stability\nnumber 1", transform=ax.get_xaxis_transform(),
            color=THEORY_COLOUR, fontsize=8, va="center", ha="left",
        )
    if log_x:
        ax.set_xscale("log")
    # Log scale only where the iterations span decades (the stiffness sweep)
    if log_x:
        ax.set_yscale("log")
    else:
        ax.set_ylim(bottom=0)
    ax.set_title(
        f"{sweep.capitalize()} sweep: iterations per time step"
        "\n(x: crashed; its iterations up to the crash)", fontsize=10,
    )
    finish_plot(ax, "Iterations per time step", xlabel, log_y=log_x)
    save_or_show(filename)


def plot_accuracy(rows, x, xlabel, log_x, sweep, filename):
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()

    for coupling in ["monolithic", "partitioned"]:
        y = [abs(row[coupling]["twist_error"]) for row in rows]
        ax.plot(
            x, y, color=COLOURS[coupling],
            marker="o", markersize=9 if coupling == "monolithic" else 6,
            markerfacecolor=COLOURS[coupling] if coupling == "monolithic" else "white",
            linewidth=2.0 if coupling == "monolithic" else 1.0,
            linestyle="-" if coupling == "monolithic" else (0, (4, 3)),
            label=f"{coupling_label(coupling)}: |settled twist error|",
        )

    # With omega dt well above 1 the trapezoidal rule damps the torsional
    # mode slowly: those runs have not settled by the end time
    for xi, row in zip(x, rows):
        s = row["monolithic"]
        omega_dt = math.sqrt(s["p"]["stability"]/0.25)
        if s["finished"] and omega_dt > 2:
            ax.annotate(
                f"not settled\n(omega dt = {omega_dt:.2g})",
                (xi, abs(s["twist_error"])), xytext=(-8, 0),
                textcoords="offset points", fontsize=8, ha="right", va="center",
            )

    y = [
        row["difference"]/static_twist(row["monolithic"])
        for row in rows
    ]
    ax.plot(
        x, y, color=THEORY_COLOUR, marker="^", markersize=6, linewidth=1.0,
        label="|Partitioned - monolithic| (max over time)",
    )

    if log_x:
        ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(
        f"{sweep.capitalize()} sweep: error relative to the static twist"
        "\n(missing points: a run crashed)", fontsize=10,
    )
    finish_plot(ax, "Relative to the static twist", xlabel, log_y=True)
    save_or_show(filename)


def plot_runtime(rows, x, xlabel, log_x, sweep, filename):
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()

    for coupling in COUPLINGS:
        y = [row[coupling]["log"]["runtime"] for row in rows]
        plot_finished(
            ax, x, y, rows, coupling, color=COLOURS[coupling], marker="o",
            markersize=6, linewidth=1.6, label=coupling_label(coupling),
        )
        mark_crashes(ax, x, y, rows, coupling)

    if log_x:
        ax.set_xscale("log")
    ax.set_ylim(bottom=0)
    ax.set_title(
        f"{sweep.capitalize()} sweep: run time\n(x: crashed; time to the crash)",
        fontsize=10,
    )
    finish_plot(ax, "Execution time (s)", xlabel)
    save_or_show(filename)


def plot_newton_bars(rows, sweep, x_of, xlabel, average=False):
    """Two bars per sweep case, using the same logged counts as the summaries.

    Categories keep every case readable despite the uneven parameter spacing.
    Hatched bars show iterations up to an early stop, not a full-length run.
    """
    if not rows:
        return

    fig, ax = plt.subplots(figsize=(max(8, 0.9*len(rows)), 4.8))
    positions = np.arange(len(rows))
    width = 0.36
    metric = "average" if average else "total"
    ylabel = (
        "Average beam Newton iterations per time step"
        if average else "Total beam Newton iterations"
    )
    partial = False
    for offset, coupling in zip((-width/2, width/2), ("monolithic", "partitioned")):
        values = [
            row[coupling]["newton_per_step"] if average
            else row[coupling]["newton_total"]
            for row in rows
        ]
        bars = ax.bar(
            positions + offset, values, width, color=COLOURS[coupling],
            label=coupling_label(coupling), zorder=3,
        )
        for bar, row in zip(bars, rows):
            if not row[coupling]["finished"]:
                partial = True
                bar.set_hatch("///")
                bar.set_edgecolor("0.25")

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{x_of(row):.6g}" for row in rows], rotation=30, ha="right")
    ax.set_ylim(bottom=0)
    ax.set_title(f"{sweep.capitalize()} sweep: {metric} beam Newton iterations")
    finish_plot(ax, ylabel, xlabel)
    # Vertical grid lines have no quantitative meaning on a categorical axis.
    ax.xaxis.grid(False, which="both")
    if partial:
        fig.text(
            0.5, 0.01,
            "Hatched bars: run stopped early; counts cover only the logged run.",
            ha="center", fontsize=8, color="0.35",
        )
        fig.tight_layout(rect=(0, 0.06, 1, 1))
    save_or_show(f"{sweep}_newton_{metric}_bars.png")


def plot_newton_histories(rows, sweep, title_of, filename):
    """nNewton per time step against time, one panel per case (cases without
    beamSolverIterations output are left out)"""
    rows = [
        row for row in rows
        if all(row[c]["solver"] is not None for c in COUPLINGS)
    ]
    if not rows:
        return
    ncols = 3
    nrows = math.ceil(len(rows)/ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(11, 2.8*nrows), sharex=True, squeeze=False
    )

    for ax, row in zip(axes.flat, rows):
        for coupling in ["partitioned", "monolithic"]:
            solver = row[coupling]["solver"]
            ax.step(
                solver["time"], solver["n_newton"], where="post",
                color=COLOURS[coupling], linewidth=1.0,
                label=coupling_label(coupling),
            )
            if not row[coupling]["finished"]:
                ax.plot(
                    solver["time"][-1], solver["n_newton"][-1], marker="x",
                    markersize=9, markeredgewidth=2, color=COLOURS[coupling],
                    clip_on=False,
                )
        ax.set_ylim(bottom=0)
        ax.set_title(title_of(row), fontsize=10)
        ax.grid(which="major", ls="--", lw=0.6, alpha=0.75)

    for ax in axes.flat[len(rows):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Time (s)")
    for ax in axes[:, 0]:
        ax.set_ylabel("nNewton per time step")

    axes.flat[0].legend(fontsize=8)
    fig.suptitle(
        f"{sweep.capitalize()} sweep: beam Newton iterations per time step "
        "(x: crashed)"
    )
    fig.tight_layout()
    save_or_show(filename)


def plot_average_newton(rows, x, xlabel, log_x, sweep, filename, stability_line=False):
    """Average nNewton per time step against the sweep parameter"""
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()

    for coupling in COUPLINGS:
        y = [row[coupling]["newton_per_step"] for row in rows]
        plot_finished(
            ax, x, y, rows, coupling, color=COLOURS[coupling], marker="o",
            markersize=6, linewidth=1.6, label=coupling_label(coupling),
        )
        mark_crashes(ax, x, y, rows, coupling)
        for xi, yi, row in zip(x, y, rows):
            if row[coupling]["finished"]:
                ax.annotate(
                    f"{yi:.3g}", (xi, yi), xytext=(0, 6),
                    textcoords="offset points", fontsize=7, ha="center",
                    color=COLOURS[coupling],
                )

    if stability_line:
        ax.axvline(1, color=THEORY_COLOUR, ls=":", lw=1.0)
    if log_x:
        ax.set_xscale("log")
        ax.set_yscale("log")
    else:
        ax.set_ylim(bottom=0)
    ax.set_title(
        f"{sweep.capitalize()} sweep: average beam Newton iterations per time step"
        "\n(x: crashed; its average up to the crash)", fontsize=10,
    )
    finish_plot(ax, "Average nNewton per time step", xlabel, log_y=log_x)
    save_or_show(filename)


def plot_sweep(rows, sweep, x_of, xlabel, log_x, title_of, stability_line=False):
    if not rows:
        print(f"No {sweep} sweep results: run rigidBodyEndTwist/Allrun {sweep}")
        return
    x = [x_of(row) for row in rows]
    plot_histories(rows, sweep, title_of, f"{sweep}_histories.png")
    plot_iterations(
        rows, x, xlabel, log_x, sweep, f"{sweep}_iterations.png", stability_line
    )
    plot_accuracy(rows, x, xlabel, log_x, sweep, f"{sweep}_accuracy.png")
    plot_runtime(rows, x, xlabel, log_x, sweep, f"{sweep}_runtime.png")
    plot_newton_histories(rows, sweep, title_of, f"{sweep}_newton_per_step.png")
    plot_average_newton(
        rows, x, xlabel, log_x, sweep, f"{sweep}_newton_average.png",
        stability_line,
    )


plot_sweep(
    moment, "moment",
    x_of=lambda row: static_twist(row["monolithic"]),
    xlabel="Static twist M L/(G Ip) (rad)",
    log_x=False,
    title_of=lambda row: f"static twist {static_twist(row['monolithic']):.3g} rad",
)

plot_sweep(
    stiffness, "stiffness",
    x_of=lambda row: row["monolithic"]["p"]["stability"],
    xlabel="Partitioned stability number beta dt^2 kt/J",
    log_x=True,
    title_of=lambda row: (
        f"G x {row['label'][2:]}, stability {row['monolithic']['p']['stability']:.2g}"
    ),
    stability_line=True,
)

for average in (False, True):
    plot_newton_bars(
        moment, "moment", x_of=lambda row: row["monolithic"]["p"]["M"],
        xlabel="Applied moment M (N m)", average=average,
    )

for average in (False, True):
    plot_newton_bars(
        stiffness, "stiffness", x_of=lambda row: row["monolithic"]["p"]["kt"],
        xlabel="Torsional stiffness kt = G Ip/L (N m/rad)", average=average,
    )


# =============================================================
# MOORFV PLOTS
# =============================================================

def moorfv_styles(cases):
    styles = {}
    partitioned = 0
    for case in cases:
        if case["coupling"] == "monolithic":
            styles[case["label"]] = ("-", 1.6)
        else:
            styles[case["label"]] = (
                PARTITIONED_STYLES[partitioned % len(PARTITIONED_STYLES)], 1.2
            )
            partitioned += 1
    return styles


if moorfv:
    Omega0, omega = compare.theory()
    amplitude = Omega0/omega
    half_period = math.pi/omega
    styles = moorfv_styles(moorfv)

    # Twist against theory
    plt.figure(figsize=(8, 4.5))
    ax = plt.gca()
    t_end = max(c["time"][-1] for c in moorfv)
    t = np.linspace(0, t_end, 2000)
    ax.plot(
        t, np.sin(omega*t), color=THEORY_COLOUR, ls=":", lw=1.2, zorder=5,
        label="Theory (Omega0/omega) sin(omega t)",
    )
    for case in moorfv:
        ls, lw = styles[case["label"]]
        ax.plot(
            case["time"], case["twist"]/amplitude, color=COLOURS[case["coupling"]],
            ls=ls, lw=lw, label=case["label"],
        )
    ax.set_title("moorFV: free torsional oscillation of the welded body", fontsize=10)
    finish_plot(ax, "Twist / theoretical amplitude", "Time (s)", legend=False)
    ax.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=9,
        frameon=False,
    )
    plt.tight_layout()
    save_or_show("moorFV_twist.png")

    # Peak twist in each half period
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()
    for case in moorfv:
        ls, lw = styles[case["label"]]
        n = int(case["time"][-1]/half_period + 1e-9)
        peaks = []
        for k in range(n):
            mask = (case["time"] >= k*half_period) & (case["time"] < (k + 1)*half_period)
            if mask.any():
                peaks.append(np.max(np.abs(case["twist"][mask]))/amplitude)
        ax.plot(
            np.arange(1, len(peaks) + 1), peaks, color=COLOURS[case["coupling"]],
            ls=ls, lw=lw, marker="o", markersize=6, label=case["label"],
        )
    ax.axhline(1, color=THEORY_COLOUR, ls=":", lw=1.0, label="Theory")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title(
        "moorFV: peak twist in each half period\n"
        "(lines that match the monolithic run lie under it)", fontsize=10,
    )
    finish_plot(ax, "Peak twist / theoretical amplitude", "Half period")
    save_or_show("moorFV_peaks.png")

    # Cost: beam Newton iterations per time step and run time per simulated second
    labels = [c["label"] for c in moorfv]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    newton = [c["newton_per_step"] for c in moorfv]
    runtime = [
        c["log"]["runtime"]/c["log"]["time"] if c["log"]["time"] else math.nan
        for c in moorfv
    ]
    colours = [COLOURS[c["coupling"]] for c in moorfv]
    y = np.arange(len(moorfv))
    for ax, values, xlabel in (
        (axes[0], newton, "Beam Newton iterations per time step"),
        (axes[1], runtime, "Execution time per simulated second (s)"),
    ):
        ax.barh(y, values, color=colours, height=0.6)
        for yi, v in zip(y, values):
            ax.text(
                v, yi, f" {v:.3g}" if v < 1000 else f" {v:.0f}",
                va="center", fontsize=9,
            )
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.invert_yaxis()
        ax.set_xlim(0, max(v for v in values if not math.isnan(v))*1.2)
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", ls="--", lw=0.6, alpha=0.75)
    axes[1].set_yticklabels([])
    fig.suptitle("moorFV: cost", fontsize=11)
    fig.tight_layout()
    save_or_show("moorFV_cost.png")

    # nNewton per time step and its running average
    solver_cases = [c for c in moorfv if c["solver"] is not None]
    for key, ylabel, title, filename in (
        (
            "n_newton", "nNewton per time step",
            "moorFV: beam Newton iterations per time step (all outer correctors)",
            "moorFV_newton_per_step.png",
        ),
        (
            "running_mean_newton", "Average nNewton per time step",
            "moorFV: running average of the beam Newton iterations per time step",
            "moorFV_newton_running_average.png",
        ),
    ):
        if not solver_cases:
            break
        plt.figure(figsize=(8, 4.5))
        ax = plt.gca()
        for case in solver_cases:
            ls, lw = styles[case["label"]]
            ax.plot(
                case["solver"]["time"], case["solver"][key],
                color=COLOURS[case["coupling"]], ls=ls, lw=lw,
                label=case["label"],
            )
        ax.set_ylim(bottom=0)
        ax.set_title(title, fontsize=10)
        finish_plot(ax, ylabel, "Time (s)")
        save_or_show(filename)
else:
    print("No moorFV results: run rigidBodyEndTwist/moorFV/Allrun")
