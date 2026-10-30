#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare the partitioned (FvBeamNewmark) and monolithic (beamFoamCoupled)
couplings over the floatRigidBody sweeps: a floating box moored by one line
in waves (floatRigidBody_sweeps, run with floatRigidBody_sweeps/Allrun).

For each sweep (stiffness: the line's E; period: the wave period):
    - surge, heave and pitch histories, partitioned against monolithic, one
      panel per case
    - surge, heave and pitch amplitudes and the largest line tension against
      the sweep parameter (heave also as a ratio to the wave amplitude H/2
      for the period sweep)
    - the largest partitioned-monolithic heave difference, relative to the
      monolithic heave amplitude
    - run time and total beam Newton iterations
    - bar charts, one pair of bars (monolithic, partitioned) per case, of
      the total, average and largest nNewton per time step and the
      execution time per simulated second (hatched: crashed run)
Beam Newton iterations (nNewton, BlockEigenSolverOF::solve calls per time
step over all PIMPLE outer correctors) from postProcessing/0/
beamSolverIterations.dat:
    - nNewton per time step and its running average against time, one panel
      per case
    - average nNewton per time step against the sweep parameter
Animation cases (base, written every 0.05 s for ParaView): motion and
nNewton histories, and bars of the same four costs.

A cross marks a run that crashed, at the last value it reached.

Run from Spyder or from the beamFoam/tutorials directory after the cases have
completed. Uses the readers in floatRigidBody_sweeps/summarise.py.
"""

import math
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import AutoMinorLocator, MaxNLocator
import numpy as np


# =============================================================
# INPUT FILES
# =============================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SWEEP_DIR = os.path.join(SCRIPT_DIR, "floatRigidBody_sweeps")

sys.path.insert(0, SCRIPT_DIR)
sys.path.insert(0, SWEEP_DIR)
import summarise  # noqa: E402  (floatRigidBody_sweeps/summarise.py)
from read_beam_solver_iterations import read_beam_solver_iterations  # noqa: E402

COUPLINGS = ["partitioned", "monolithic"]

# Initial centre of mass of the box (base/constant/dynamicMeshDict)
CENTRE_OF_MASS = (0.5, 0.0)


# =============================================================
# USER SETTINGS
# =============================================================

SAVE_FIGURES = False
FIGURE_DIR = os.path.join(SCRIPT_DIR, "comparison_plots", "floatRigidBody_sweeps")

# One colour per coupling, the same in every figure
COLOURS = {"monolithic": "#2a78d6", "partitioned": "#eb6834"}
REFERENCE_COLOUR = "0.35"

# (key in the motion arrays, axis label) for the history panels
MOTION = [
    ("surge", "Surge (m)"),
    ("heave", "Heave (m)"),
    ("pitch", "Pitch (rad)"),
]


# =============================================================
# READ DATA
# =============================================================

def load_case(case):
    s = summarise.summarise_case(case)
    motion = np.array(s["motion"]).reshape(-1, 4)
    s["time"] = motion[:, 0]
    s["surge_history"] = motion[:, 1] - CENTRE_OF_MASS[0]
    s["heave_history"] = motion[:, 2] - CENTRE_OF_MASS[1]
    s["pitch_history"] = motion[:, 3]
    s["solver"] = read_beam_solver_iterations(case)
    s["finished"] = s["log"]["status"] == "finished"
    return s


def load_group(directory, label):
    """{label, partitioned, monolithic, difference}, or None if either
    coupling is missing"""
    row = {"label": label}
    for coupling in COUPLINGS:
        case = os.path.join(directory, coupling)
        try:
            row[coupling] = load_case(case)
        except (FileNotFoundError, KeyError):
            return None
    mono = row["monolithic"]
    row["difference"] = summarise.heave_difference(
        row["partitioned"]["motion"], mono["motion"]
    )/mono["heave"]
    return row


def load_sweep(sweep):
    root = os.path.join(SWEEP_DIR, sweep)
    if not os.path.isdir(root):
        return []
    rows = []
    for label in sorted(os.listdir(root), key=summarise.label_value):
        row = load_group(os.path.join(root, label), label)
        if row is not None:
            rows.append(row)
    return rows


stiffness = load_sweep("stiffness")
period = load_sweep("period")
animation = load_group(os.path.join(SWEEP_DIR, "animation"), "base")


# =============================================================
# PLOT HELPERS
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


def plot_over_cases(ax, x, rows, value_of, coupling, annotate=None):
    """Line through the runs that finished, a cross at each crashed run"""
    y = [value_of(row[coupling]) for row in rows]
    keep = [i for i, row in enumerate(rows) if row[coupling]["finished"]]
    ax.plot(
        [x[i] for i in keep], [y[i] for i in keep], color=COLOURS[coupling],
        marker="o", markersize=5, linewidth=1.6, label=coupling_label(coupling),
    )
    for xi, yi, row in zip(x, y, rows):
        if row[coupling]["log"]["status"] == "crashed" and not math.isnan(yi):
            ax.plot(
                xi, yi, marker="x", markersize=10, markeredgewidth=2,
                color=COLOURS[coupling], linestyle="none",
            )
        elif annotate and row[coupling]["finished"]:
            ax.annotate(
                format(yi, annotate), (xi, yi), xytext=(0, 6),
                textcoords="offset points", fontsize=7, ha="center",
                color=COLOURS[coupling],
            )


def case_panels(rows, title_of, suptitle, ylabel, draw, filename):
    """One panel per case; draw(ax, case, coupling) plots one coupling"""
    if not rows:
        return
    ncols = 3
    nrows = math.ceil(len(rows)/ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(11, 2.8*nrows), sharex=True, squeeze=False
    )

    for ax, row in zip(axes.flat, rows):
        # Monolithic underneath, partitioned dashed on top: where they agree
        # both stay visible
        for coupling in ["monolithic", "partitioned"]:
            case = row[coupling]
            draw(ax, case, coupling)
            if case["log"]["status"] == "crashed":
                ax.text(
                    0.97, 0.95,
                    f"x: {coupling} crashed at t = {case['log']['time']:.3g} s",
                    transform=ax.transAxes, fontsize=8, ha="right", va="top",
                    color=COLOURS[coupling],
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
                )
        ax.set_title(title_of(row), fontsize=10)
        ax.grid(which="major", ls="--", lw=0.6, alpha=0.75)

    for ax in axes.flat[len(rows):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Time (s)")
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel)

    axes.flat[0].legend(fontsize=8)
    fig.suptitle(suptitle)
    fig.tight_layout()
    save_or_show(filename)


def history_drawer(key):
    def draw(ax, case, coupling):
        t, y = case["time"], case[f"{key}_history"]
        if t.size == 0:
            return
        ax.plot(
            t, y, color=COLOURS[coupling],
            linewidth=1.6 if coupling == "monolithic" else 0.9,
            linestyle="-" if coupling == "monolithic" else (0, (4, 3)),
            label=coupling_label(coupling),
        )
        if case["log"]["status"] == "crashed":
            ax.plot(
                t[-1], y[-1], marker="x", markersize=9, markeredgewidth=2,
                color=COLOURS[coupling], clip_on=False,
            )
    return draw


def solver_drawer(key):
    def draw(ax, case, coupling):
        solver = case["solver"]
        if solver is None or solver["time"].size == 0:
            return
        ax.plot(
            solver["time"], solver[key], color=COLOURS[coupling],
            linewidth=1.0, label=coupling_label(coupling),
        )
        if case["log"]["status"] == "crashed":
            ax.plot(
                solver["time"][-1], solver[key][-1], marker="x", markersize=9,
                markeredgewidth=2, color=COLOURS[coupling], clip_on=False,
            )
        ax.set_ylim(bottom=0)
    return draw


# =============================================================
# SWEEP PLOTS
# =============================================================

def plot_response(rows, x, xlabel, log_x, sweep, filename):
    """Surge, heave, pitch amplitudes and largest tension against the sweep
    parameter"""
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    panels = [
        ("surge", "Surge amplitude (m)"),
        ("heave", "Heave amplitude (m)"),
        ("pitch", "Pitch amplitude (rad)"),
        ("tension", "Largest line tension (N)"),
    ]
    for ax, (key, ylabel) in zip(axes.flat, panels):
        for coupling in COUPLINGS:
            plot_over_cases(ax, x, rows, lambda s, k=key: s[k], coupling)
        if key == "heave" and sweep == "period":
            # Wave amplitude H/2 for reference
            half_height = rows[0]["monolithic"]["p"]["height"]/2
            ax.axhline(half_height, color=REFERENCE_COLOUR, ls=":", lw=1.0,
                       label="Wave amplitude H/2")
        if log_x:
            ax.set_xscale("log")
        finish_plot(ax, ylabel, xlabel if ax in axes[-1] else "",
                    legend=ax is axes.flat[0] or key == "heave")
    fig.suptitle(
        f"{sweep.capitalize()} sweep: box response for t >= "
        f"{summarise.T_SETTLED:g} s (half the range; x: crashed)", fontsize=11,
    )
    fig.tight_layout()
    save_or_show(filename)


def plot_difference(rows, x, xlabel, log_x, sweep, filename):
    y = [row["difference"] for row in rows]
    if not any(v > 0 for v in y):
        return
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()
    ax.plot(x, y, color=REFERENCE_COLOUR, marker="o", markersize=5, linewidth=1.6)
    if log_x:
        ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(
        f"{sweep.capitalize()} sweep: largest partitioned-monolithic heave "
        f"difference\nfor t >= {summarise.T_SETTLED:g} s, over the monolithic "
        "heave amplitude", fontsize=10,
    )
    finish_plot(ax, "Heave difference / heave amplitude", xlabel,
                legend=False, log_y=True)
    save_or_show(filename)


def plot_cost(rows, x, xlabel, log_x, sweep, filename):
    """Run time and total nNewton against the sweep parameter"""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, key, ylabel in (
        (axes[0], lambda s: s["log"]["runtime"], "Run time (s)"),
        (axes[1], lambda s: s["newton_total"], "Total nNewton"),
    ):
        for coupling in COUPLINGS:
            plot_over_cases(ax, x, rows, key, coupling)
        if log_x:
            ax.set_xscale("log")
        ax.set_ylim(bottom=0)
        finish_plot(ax, ylabel, xlabel)
    fig.suptitle(f"{sweep.capitalize()} sweep: cost (x: crashed, up to the crash)",
                 fontsize=11)
    fig.tight_layout()
    save_or_show(filename)


def plot_average_newton(rows, x, xlabel, log_x, sweep, filename):
    plt.figure(figsize=(7, 4.5))
    ax = plt.gca()
    for coupling in COUPLINGS:
        plot_over_cases(
            ax, x, rows, lambda s: s["newton_mean"], coupling, annotate=".3g"
        )
    if log_x:
        ax.set_xscale("log")
    ax.set_ylim(bottom=0)
    ax.set_title(
        f"{sweep.capitalize()} sweep: average beam Newton iterations per time "
        "step\n(all outer correctors; x: crashed, its average up to the crash)",
        fontsize=10,
    )
    finish_plot(ax, "Average nNewton per time step", xlabel)
    save_or_show(filename)


# (key, value of one run, y label, title, file suffix) for the bar charts
BAR_METRICS = [
    (
        "total", lambda s: s["newton_total"], "Total beam Newton iterations",
        "total beam Newton iterations", "newton_total_bars",
    ),
    (
        "average", lambda s: s["newton_mean"],
        "Average beam Newton iterations per time step",
        "average beam Newton iterations per time step", "newton_average_bars",
    ),
    (
        "max", lambda s: s["newton_max"],
        "Most beam Newton iterations in a time step",
        "most beam Newton iterations in a time step", "newton_max_bars",
    ),
    (
        "runtime",
        lambda s: s["log"]["runtime"]/s["log"]["time"] if s["log"]["time"] else math.nan,
        "Execution time per simulated second (s)",
        "execution time per simulated second", "runtime_bars",
    ),
]


def plot_bars(rows, name, tick_of, xlabel, value_of, ylabel, title, filename):
    """Two bars per case, monolithic and partitioned, labelled with their
    values. Categories keep every case readable despite the uneven parameter
    spacing; hatched bars are runs that crashed, counted up to the crash"""
    if not rows:
        return

    fig, ax = plt.subplots(figsize=(max(8, 1.1*len(rows)), 4.8))
    positions = np.arange(len(rows))
    width = 0.36
    crashed = False
    for offset, coupling in zip((-width/2, width/2), ("monolithic", "partitioned")):
        values = [value_of(row[coupling]) for row in rows]
        bars = ax.bar(
            positions + offset, values, width, color=COLOURS[coupling],
            label=coupling_label(coupling), zorder=3,
        )
        for bar, row, value in zip(bars, rows, values):
            if row[coupling]["log"]["status"] == "crashed":
                crashed = True
                bar.set_hatch("///")
                bar.set_edgecolor("0.25")
            if not math.isnan(value):
                ax.annotate(
                    f"{value:.3g}" if value < 1000 else f"{value:.0f}",
                    (bar.get_x() + bar.get_width()/2, value), xytext=(0, 3),
                    textcoords="offset points", ha="center", fontsize=7,
                    rotation=90 if len(rows) > 4 else 0,
                )

    ax.set_xticks(positions)
    ax.set_xticklabels([tick_of(row) for row in rows])
    ax.set_ylim(0, ax.get_ylim()[1]*1.12)
    ax.set_title(f"{name}: {title}", fontsize=11)
    finish_plot(ax, ylabel, xlabel, legend=False)
    # Plain patches: the bars' own handles carry the first bar's hatching
    ax.legend(handles=[
        Patch(color=COLOURS[c], label=coupling_label(c))
        for c in ("monolithic", "partitioned")
    ])
    # Vertical grid lines have no quantitative meaning on a categorical axis
    ax.xaxis.grid(False, which="both")
    if crashed:
        fig.text(
            0.5, 0.01, "Hatched bars: run crashed; values cover the run up to the crash.",
            ha="center", fontsize=8, color="0.35",
        )
        fig.tight_layout(rect=(0, 0.05, 1, 1))
    save_or_show(filename)


def plot_sweep_bars(rows, sweep, tick_of, xlabel):
    for _, value_of, ylabel, title, suffix in BAR_METRICS:
        plot_bars(
            rows, f"{sweep.capitalize()} sweep", tick_of, xlabel, value_of,
            ylabel, title, f"{sweep}_{suffix}.png",
        )


def plot_sweep(rows, sweep, x_of, xlabel, log_x, title_of):
    if not rows:
        print(f"No {sweep} sweep results: run floatRigidBody_sweeps/Allrun {sweep}")
        return
    x = [x_of(row) for row in rows]

    for key, ylabel in MOTION:
        case_panels(
            rows, title_of, f"{sweep.capitalize()} sweep: {key} histories",
            ylabel, history_drawer(key), f"{sweep}_{key}_histories.png",
        )
    plot_response(rows, x, xlabel, log_x, sweep, f"{sweep}_response.png")
    plot_difference(rows, x, xlabel, log_x, sweep, f"{sweep}_difference.png")
    plot_cost(rows, x, xlabel, log_x, sweep, f"{sweep}_cost.png")

    case_panels(
        rows, title_of,
        f"{sweep.capitalize()} sweep: beam Newton iterations per time step",
        "nNewton per time step", solver_drawer("n_newton"),
        f"{sweep}_newton_per_step.png",
    )
    case_panels(
        rows, title_of,
        f"{sweep.capitalize()} sweep: running average of nNewton per time step",
        "Average nNewton per time step", solver_drawer("running_mean_newton"),
        f"{sweep}_newton_running_average.png",
    )
    plot_average_newton(rows, x, xlabel, log_x, sweep, f"{sweep}_newton_average.png")


plot_sweep(
    stiffness, "stiffness",
    x_of=lambda row: row["monolithic"]["p"]["kAxial"],
    xlabel="Line axial stiffness EA/L (N/m)",
    log_x=True,
    title_of=lambda row: (
        f"E x {row['label'][2:]}, axial mode "
        f"{row['monolithic']['p']['omega']:.3g} rad/s"
    ),
)

plot_sweep(
    period, "period",
    x_of=lambda row: row["monolithic"]["p"]["period"],
    xlabel="Wave period (s)",
    log_x=False,
    title_of=lambda row: (
        f"T = {row['monolithic']['p']['period']:g} s, "
        f"wavelength {row['monolithic']['p']['wavelength']:.2f} m"
    ),
)

plot_sweep_bars(
    stiffness, "stiffness",
    tick_of=lambda row: f"x{row['label'][2:]}\n{row['monolithic']['p']['kAxial']:.3g}",
    xlabel="Line E relative to base, and axial stiffness EA/L (N/m)",
)

plot_sweep_bars(
    period, "period",
    tick_of=lambda row: f"{row['monolithic']['p']['period']:g}",
    xlabel="Wave period (s)",
)


# =============================================================
# ANIMATION CASES
# =============================================================

if animation is not None:
    fig, axes = plt.subplots(4, 1, figsize=(9, 10), sharex=True)
    for ax, (key, ylabel) in zip(axes, MOTION):
        for coupling in ["monolithic", "partitioned"]:
            history_drawer(key)(ax, animation[coupling], coupling)
        finish_plot(ax, ylabel, "", legend=key == "surge")
    for coupling in COUPLINGS:
        solver_drawer("n_newton")(axes[-1], animation[coupling], coupling)
    finish_plot(axes[-1], "nNewton per time step", "Time (s)", legend=False)
    fig.suptitle(
        "Animation cases (base): box motion and beam Newton iterations "
        f"(heave difference {animation['difference']:.2g} of the amplitude)",
        fontsize=11,
    )
    fig.tight_layout()
    save_or_show("animation_histories.png")

    # Cost of the two couplings for the base case: one panel per metric
    fig, axes = plt.subplots(1, len(BAR_METRICS), figsize=(13, 2.6), sharey=True)
    y = np.arange(len(COUPLINGS))
    for ax, (_, value_of, xlabel, _, _) in zip(axes, BAR_METRICS):
        values = [value_of(animation[c]) for c in COUPLINGS]
        ax.barh(y, values, color=[COLOURS[c] for c in COUPLINGS], height=0.6, zorder=3)
        for yi, v in zip(y, values):
            ax.text(v, yi, f" {v:.3g}" if v < 1000 else f" {v:.0f}",
                    va="center", fontsize=9)
        ax.set_xlim(0, max(values)*1.3)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.grid(axis="x", ls="--", lw=0.6, alpha=0.75, zorder=0)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([coupling_label(c) for c in COUPLINGS])
    axes[0].invert_yaxis()
    fig.suptitle("Animation cases (base): cost", fontsize=11)
    fig.tight_layout()
    save_or_show("animation_cost_bars.png")
else:
    print("No animation results: run floatRigidBody_sweeps/Allrun animation")
