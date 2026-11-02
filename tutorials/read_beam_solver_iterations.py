#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read the output of the beamSolverIterations function object
(libwireBunchingModels): postProcessing/<startTime>/<name>.dat with one row
per time step

    Time  nPimple  nBeamSolves  nNewton  totalNewton

nPimple      beam solves requested in the step: one per PIMPLE outer
             corrector in a coupled interFoam run, 1 in a beamFoam run
nBeamSolves  beam Newton loops (more than nPimple with beamFoam's
             partitioned rigidBodyEnd coupling iterations)
nNewton      Newton iterations (BlockEigenSolverOF::solve calls) in the step
totalNewton  nNewton summed from the start of the run

Used by the monolithic/partitioned comparison scripts in this directory.
"""

import glob
import os

import numpy as np


def read_beam_solver_iterations(case_dir, name="beamSolverIterations"):
    """Dict of arrays (time, n_pimple, n_beam_solves, n_newton,
    total_newton, running_mean_newton), or None if the case has no output.
    Files from restarts (several postProcessing/<startTime>) are joined; a
    time written twice keeps its last row"""
    files = glob.glob(os.path.join(case_dir, "postProcessing", "*", f"{name}.dat"))
    if not files:
        return None

    def start_time(path):
        try:
            return float(os.path.basename(os.path.dirname(path)))
        except ValueError:
            return float("inf")

    blocks = []
    for path in sorted(files, key=start_time):
        data = np.genfromtxt(path, comments="#", invalid_raise=False)
        if data.size:
            blocks.append(data.reshape(-1, 5))
    if not blocks:
        return None

    data = np.vstack(blocks)
    _, last = np.unique(data[::-1, 0], return_index=True)
    data = data[::-1][last]

    return window({
        "time": data[:, 0],
        "n_pimple": data[:, 1],
        "n_beam_solves": data[:, 2],
        "n_newton": data[:, 3],
        "total_newton": data[:, 4],
    })


def window(data, t_start=None, t_end=None):
    """The rows in [t_start, t_end] (all rows by default), with the running
    mean of nNewton taken from t_start"""
    if data is None:
        return None
    mask = np.ones_like(data["time"], dtype=bool)
    if t_start is not None:
        mask &= data["time"] >= t_start
    if t_end is not None:
        mask &= data["time"] <= t_end
    out = {key: value[mask] for key, value in data.items()}
    # Mean nNewton per time step up to each step
    out["running_mean_newton"] = (
        np.cumsum(out["n_newton"])/np.arange(1, len(out["n_newton"]) + 1)
    )
    return out
