#!/usr/bin/env python3
"""Regression check for the authoritative BlockEigen rigid-body path."""

from pathlib import Path
import math
import re
import sys


CASE_DIR = Path(__file__).resolve().parent

LOG_FILE = CASE_DIR / "log.interFoam"
BLOCK_EIGEN_FILE = (
    CASE_DIR
    / "postProcessing"
    / "BlockEigenSolve_rigidBodyMotion"
    / "0"
    / "BlockEigenSolve_rigidBodyMotion.dat"
)
SIX_DOF_FILE = (
    CASE_DIR
    / "postProcessing"
    / "sixDoF_History"
    / "0"
    / "sixDoFRigidBodyStateFvBeam.dat"
)
DYNAMIC_MESH_DICT = CASE_DIR / "constant" / "dynamicMeshDict"

POSITION_TOL = 1.0e-8
VELOCITY_TOL = 1.0e-6


def fail(message):
    print(f"FAIL: {message}")
    raise SystemExit(1)


def parse_vector(text):
    values = text.strip().strip("()").split()
    if len(values) != 3:
        raise ValueError(f"Expected a 3-vector, got: {text}")
    return tuple(float(value) for value in values)


def vector_error(a, b):
    return tuple(ai - bi for ai, bi in zip(a, b))


def max_abs(values):
    return max(abs(value) for value in values)


def read_initial_centre_of_rotation():
    text = DYNAMIC_MESH_DICT.read_text()

    match = re.search(r"initialCentreOfRotation\s+(\([^;]+\));", text)
    if not match:
        match = re.search(r"centreOfMass\s+(\([^;]+\));", text)

    if not match:
        fail(f"Could not read initial centre of rotation from {DYNAMIC_MESH_DICT}")

    return parse_vector(match.group(1))


def read_last_block_eigen_row():
    if not BLOCK_EIGEN_FILE.exists():
        fail(f"Missing BlockEigen output: {BLOCK_EIGEN_FILE}")

    data_lines = [
        line.strip()
        for line in BLOCK_EIGEN_FILE.read_text().splitlines()
        if line.strip() and not line.startswith("#") and not line.startswith("Time")
    ]

    if not data_lines:
        fail(f"No data rows found in {BLOCK_EIGEN_FILE}")

    fields = data_lines[-1].split()
    if len(fields) < 21:
        fail(f"Unexpected BlockEigen row format in {BLOCK_EIGEN_FILE}")

    return {
        "time": float(fields[0]),
        "displacement": tuple(float(value) for value in fields[3:6]),
        "velocity": tuple(float(value) for value in fields[9:12]),
    }


def read_last_six_dof_row():
    if not SIX_DOF_FILE.exists():
        fail(f"Missing sixDoF output: {SIX_DOF_FILE}")

    data_lines = [
        line.strip()
        for line in SIX_DOF_FILE.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]

    if not data_lines:
        fail(f"No data rows found in {SIX_DOF_FILE}")

    line = data_lines[-1]
    vectors = re.findall(r"\([^)]+\)", line)

    if len(vectors) < 4:
        fail(f"Unexpected sixDoF row format in {SIX_DOF_FILE}")

    return {
        "time": float(line.split()[0]),
        "centre_of_rotation": parse_vector(vectors[0]),
        "velocity": parse_vector(vectors[3]),
    }


def check_log():
    if not LOG_FILE.exists():
        fail(f"Missing solver log: {LOG_FILE}")

    text = LOG_FILE.read_text(errors="replace")
    lower_text = text.lower()

    if "end\n" not in lower_text:
        fail("log.interFoam does not end with a completed run marker")

    forbidden = ["foam fatal", "segmentation fault", "core dumped"]
    for pattern in forbidden:
        if pattern in lower_text:
            fail(f"log.interFoam contains failure marker: {pattern}")

    if re.search(r"(?<![a-z])nan(?![a-z])", lower_text):
        fail("log.interFoam contains NaN")

    authoritative_count = text.count(
        "finiteVolumeBeam using authoritative BlockEigen rigid-body solution"
    )

    if authoritative_count == 0:
        fail("authoritative BlockEigen motion path was not reported in log.interFoam")

    return authoritative_count


def main():
    print("Input files")
    print(f"Log:        {LOG_FILE}")
    print(f"BlockEigen: {BLOCK_EIGEN_FILE}")
    print(f"sixDoF:     {SIX_DOF_FILE}")
    print()

    authoritative_count = check_log()
    initial_cor = read_initial_centre_of_rotation()
    block_eigen = read_last_block_eigen_row()
    six_dof = read_last_six_dof_row()

    if not math.isclose(block_eigen["time"], six_dof["time"], abs_tol=1.0e-12):
        fail(
            "Final output times differ: "
            f"BlockEigen={block_eigen['time']}, sixDoF={six_dof['time']}"
        )

    block_eigen_cor = tuple(
        initial + displacement
        for initial, displacement in zip(initial_cor, block_eigen["displacement"])
    )

    position_difference = vector_error(block_eigen_cor, six_dof["centre_of_rotation"])
    velocity_difference = vector_error(block_eigen["velocity"], six_dof["velocity"])

    print("Final comparison")
    print(f"Time:                       {block_eigen['time']}")
    print(f"Initial centre of rotation: {initial_cor}")
    print(f"BlockEigen displacement:    {block_eigen['displacement']}")
    print(f"BlockEigen CoR:             {block_eigen_cor}")
    print(f"sixDoF CoR:                 {six_dof['centre_of_rotation']}")
    print(f"CoR difference:             {position_difference}")
    print(f"BlockEigen velocity:        {block_eigen['velocity']}")
    print(f"sixDoF velocity:            {six_dof['velocity']}")
    print(f"Velocity difference:        {velocity_difference}")
    print(f"Authoritative log count:    {authoritative_count}")
    print()

    if max_abs(position_difference) > POSITION_TOL:
        fail(
            "Final centre-of-rotation mismatch exceeds tolerance "
            f"{POSITION_TOL}"
        )

    if max_abs(velocity_difference) > VELOCITY_TOL:
        fail(f"Final velocity mismatch exceeds tolerance {VELOCITY_TOL}")

    print("PASS: authoritative BlockEigen rigid-body regression check")


if __name__ == "__main__":
    main()
