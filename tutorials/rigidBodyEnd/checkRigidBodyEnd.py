#!/usr/bin/env python3
"""
Compare the rigidBodyEnd test cases with analytic results, for partitioned
coupling (the case directories) and monolithic coupling (monolithic/<case>),
and compare the two couplings with each other.

Run after the cases have been run with ./Allrun:

    python3 checkRigidBodyEnd.py

Parameters (length, radius, E, rho, mass, g) are read from each case's own
files. Uses only the Python standard library. Exit status is non-zero if any
check fails.
"""

import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


# -----------------------------------------------------------------------------
# Reading
# -----------------------------------------------------------------------------

def read_text(case, *path):
    with open(os.path.join(HERE, case, *path)) as f:
        return re.sub(r"//.*", "", f.read())


def scalar_entry(text, key):
    match = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)\s*;", text)
    if not match:
        raise KeyError(f"{key} not found")
    return float(match.group(1))


def case_parameters(case):
    beam = read_text(case, "constant", "beamProperties")
    body = beam[beam.index("rigidBodyEnd"):]
    g_text = read_text(case, "constant", "g")
    g = [float(v) for v in re.search(r"value\s+\(([^)]*)\)", g_text).group(1).split()]

    radius = scalar_entry(beam, "radius")
    area = math.pi*radius**2
    return {
        "L": scalar_entry(beam, "length"),
        "A": area,
        "I": math.pi*radius**4/4.0,
        "E": scalar_entry(beam, "E"),
        "rho": scalar_entry(beam, "rho"),
        "m": scalar_entry(body, "mass"),
        "g": g,
        "dt": scalar_entry(read_text(case, "system", "controlDict"), "deltaT"),
    }


def read_history(case):
    path = os.path.join(HERE, case, "postProcessing", "rigidBodyEnd", "0", "rigidBodyEnd.dat")
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            rows.append([float(v) for v in line.split()])
    return rows


# -----------------------------------------------------------------------------
# Signal helpers
# -----------------------------------------------------------------------------

def crossings(times, values, level):
    """Times at which the signal crosses level, in either direction."""
    found = []
    for i in range(1, len(values)):
        a, b = values[i - 1] - level, values[i] - level
        if a*b < 0 or (a != 0 and b == 0):
            found.append(times[i - 1] + (times[i] - times[i - 1])*(-a)/(b - a))
    return found


def period_about(times, values, level):
    """Mean period from crossings of the equilibrium level."""
    found = crossings(times, values, level)
    if len(found) < 3:
        raise ValueError("fewer than three crossings")
    half_periods = len(found) - 1
    return 2.0*(found[-1] - found[0])/half_periods, half_periods/2.0


def mean_over_whole_periods(times, values, level):
    """Time average between the first and last full-period crossings."""
    found = crossings(times, values, level)
    n = len(found) - 1
    start, end = found[0], found[n - n % 2]
    total = 0.0
    for i in range(1, len(times)):
        t0, t1 = max(times[i - 1], start), min(times[i], end)
        if t1 > t0:
            total += 0.5*(values[i - 1] + values[i])*(t1 - t0)
    return total/(end - start)


# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

results = []


def check(name, value, expected, rel_tol):
    error = abs(value - expected)/abs(expected)
    ok = error <= rel_tol
    results.append(ok)
    print(f"  {'PASS' if ok else 'FAIL'}  {name:<34} {value: .8g}  expected {expected: .8g}"
          f"  (rel. error {error:.2e}, tol {rel_tol:.0e})")


def hanging_mass(case="hangingMass"):
    p = case_parameters(case)
    rows = read_history(case)
    g = p["g"][0]
    k = p["E"]*p["A"]/p["L"]
    weight = p["m"]*g
    rod_weight = p["rho"]*p["A"]*p["L"]*g

    # Static stretch: body weight plus half the rod weight on average
    stretch = (weight + 0.5*rod_weight)/k

    print(f"{case}: static equilibrium after damped settling (t = {rows[-1][0]:g} s)")
    check("end force on body = -mg", rows[-1][10], -weight, 1e-4)
    check("body displacement = static stretch", rows[-1][1], stretch, 2e-3)


def axial_oscillation(case="axialOscillation"):
    p = case_parameters(case)
    rows = read_history(case)
    g = p["g"][0]
    k = p["E"]*p["A"]/p["L"]
    c = math.sqrt(p["E"]/p["rho"])
    mu = p["rho"]*p["A"]*p["L"]/p["m"]

    # Exact first axial mode of a fixed rod with an end mass:
    #   z tan(z) = rod mass / end mass,  z = omega L / c
    z = math.sqrt(mu)
    for _ in range(50):
        z -= (z*math.tan(z) - mu)/(math.tan(z) + z/math.cos(z)**2)
    omega = z*c/p["L"]

    # Trapezoidal rule period elongation
    omega_h = 2.0/p["dt"]*math.atan(omega*p["dt"]/2.0)
    period = 2.0*math.pi/omega_h

    times = [r[0] for r in rows]
    x = [r[1] for r in rows]
    stretch = p["m"]*g/k
    measured, n = period_about(times, x, stretch)

    print(f"{case}: free axial oscillation, {n:g} periods measured")
    check("period", measured, period, 5e-3)
    check(
        "mean displacement = static stretch",
        mean_over_whole_periods(times, x, stretch), stretch, 2e-2
    )
    check("peak displacement = 2 x stretch", max(x), 2.0*stretch, 1e-2)


def pendulum(case="pendulum"):
    p = case_parameters(case)
    rows = read_history(case)
    gx, _, gz = p["g"]
    g = math.hypot(gx, gz)
    theta_eq = math.atan2(gz, gx)
    k = p["E"]*p["A"]/p["L"]

    # Simple pendulum with the finite-amplitude correction; rod mass is
    # negligible. Effective length: rod length, plus the mean stretch, minus
    # the bending boundary layer at the clamp, sqrt(EI/T): below it the
    # tensioned rod hangs like a string pinned about that far from the clamp
    tension = p["m"]*g
    boundary_layer = math.sqrt(p["E"]*p["I"]/tension)
    length = p["L"] + tension/k - boundary_layer
    amplitude = theta_eq
    period = (
        2.0*math.pi*math.sqrt(length/g)
       *(1.0 + amplitude**2/16.0 + 11.0*amplitude**4/3072.0)
    )

    times = [r[0] for r in rows]
    angle = [math.atan2(r[3], p["L"] + r[1]) for r in rows]
    measured, n = period_about(times, angle, theta_eq)

    print(f"{case}: small-angle swing, amplitude {math.degrees(amplitude):.1f} deg, "
          f"{n:g} periods measured")
    check("period", measured, period, 2e-3)

    # Angle is measured from the clamp, so it is smaller than the swing
    # angle about the effective pivot by (L - boundary layer)/L
    check(
        "peak angle (from clamp)",
        max(angle), 2.0*theta_eq*(p["L"] - boundary_layer)/p["L"], 5e-3
    )


def coupling_agreement(case, rel_tol=1e-4):
    """Monolithic and converged partitioned coupling solve the same equations"""
    partitioned = read_history(case)
    monolithic = read_history(os.path.join("monolithic", case))
    n = min(len(partitioned), len(monolithic))
    if len(partitioned) != len(monolithic):
        raise ValueError(f"{len(partitioned)} partitioned rows, {len(monolithic)} monolithic")

    print(f"{case}: monolithic against partitioned, {n} time steps")
    for col, name in ((1, "displacement x"), (3, "displacement z"), (10, "beam force x"), (12, "beam force z")):
        peak = max(abs(r[col]) for r in partitioned)
        if peak < 1e-8:
            continue
        diff = max(abs(monolithic[i][col] - partitioned[i][col]) for i in range(n))
        ok = diff <= rel_tol*peak
        results.append(ok)
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<34} max difference {diff/peak:.2e} of peak"
              f"  (tol {rel_tol:.0e})")

    def newton_per_step(path):
        total = 0
        with open(os.path.join(HERE, path, "log.beamFoam")) as f:
            for line in f:
                match = re.match(r"\s+(\d+): Converged", line)
                if match:
                    total += int(match.group(1))
        return total/n

    print(f"  info  beam Newton iterations per step: partitioned "
          f"{newton_per_step(case):.2f}, monolithic "
          f"{newton_per_step(os.path.join('monolithic', case)):.2f}")


def jacobian_check(case, tol=1e-5):
    """Finite-difference check of the monolithic body columns (jacobianCheck).

    One-sided differences with a small step: the physically significant
    columns agree to about 1e-9; the tolerance allows for round-off on
    columns that are nearly zero."""
    worst = 0.0
    count = 0
    with open(os.path.join(HERE, "monolithic", case, "log.beamFoam")) as f:
        for line in f:
            if line.strip().startswith("column"):
                values = re.findall(r"rows ([-+0-9.eE]+)", line)
                worst = max([worst] + [float(v) for v in values])
                count += 1
    if count == 0:
        raise ValueError("no Jacobian check output (set jacobianCheck)")
    ok = worst <= tol
    results.append(ok)
    print(f"  {'PASS' if ok else 'FAIL'}  {case + ' Jacobian columns':<34} worst relative difference "
          f"{worst:.2e} over {count} columns  (tol {tol:.0e})")


tests = [
    (hanging_mass, "hangingMass"),
    (axial_oscillation, "axialOscillation"),
    (pendulum, "pendulum"),
]

for coupling, prefix in (("partitioned", ""), ("monolithic", "monolithic")):
    print(f"=== {coupling} coupling ===")
    print()
    for test, case in tests:
        try:
            test(os.path.join(prefix, case) if prefix else case)
        except (OSError, KeyError, ValueError) as error:
            results.append(False)
            print(f"  FAIL  {test.__name__} ({coupling}): {error}")
        print()

print("=== monolithic Jacobian against finite differences ===")
print()
for case in ("axialOscillation", "pendulum"):
    try:
        jacobian_check(case)
    except (OSError, ValueError) as error:
        results.append(False)
        print(f"  FAIL  {case} Jacobian check: {error}")
print()

print("=== monolithic against partitioned ===")
print()
for _, case in tests:
    try:
        coupling_agreement(case)
    except (OSError, KeyError, ValueError) as error:
        results.append(False)
        print(f"  FAIL  {case}: {error}")
    print()

print(f"{sum(results)} of {len(results)} checks passed")
sys.exit(0 if all(results) else 1)
