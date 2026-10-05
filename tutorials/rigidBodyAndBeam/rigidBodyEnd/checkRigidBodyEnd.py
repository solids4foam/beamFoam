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


def vector_entry(text, key):
    match = re.search(r"\b" + key + r"\s+\(([^)]*)\)\s*;", text)
    if not match:
        raise KeyError(f"{key} not found")
    return [float(v) for v in match.group(1).split()]


def body_parameters(case):
    """Rotation parameters of the rigidBodyEnd dictionary"""
    beam = read_text(case, "constant", "beamProperties")
    body = beam[beam.index("rigidBodyEnd"):]
    p = case_parameters(case)
    attachment = [p["L"], 0.0, 0.0]
    com = vector_entry(body, "centreOfMass") if "centreOfMass" in body else attachment
    p["J"] = vector_entry(body, "momentOfInertia")
    p["arm0"] = [a - c for a, c in zip(attachment, com)]
    p["omega0"] = vector_entry(body, "angularVelocity") if "angularVelocity" in body else [0.0, 0.0, 0.0]
    return p


def rotation_matrix(v):
    """Rotation tensor of a rotation vector (Rodrigues)"""
    angle = math.sqrt(sum(x*x for x in v))
    if angle < 1e-14:
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    k = [x/angle for x in v]
    K = [[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]]
    K2 = [[sum(K[i][m]*K[m][j] for m in range(3)) for j in range(3)] for i in range(3)]
    return [[(1.0 if i == j else 0.0) + math.sin(angle)*K[i][j] + (1 - math.cos(angle))*K2[i][j]
             for j in range(3)] for i in range(3)]


def mat_vec(M, v):
    return [sum(M[i][j]*v[j] for j in range(3)) for i in range(3)]


def norm(v):
    return math.sqrt(sum(x*x for x in v))


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


def hanging_body(case="hangingBody"):
    """Body with an offset centre of mass swings 90 degrees and settles"""
    p = body_parameters(case)
    rows = read_history(case)
    last = rows[-1]
    g = p["g"]
    weight = p["m"]*norm(g)

    # Force from the beam balances the weight; the arm from the centre of
    # mass to the attachment points against gravity
    force = last[10:13]
    arm = mat_vec(rotation_matrix(last[17:20]), p["arm0"])
    cosine = -sum(a*gi for a, gi in zip(arm, g))/(norm(arm)*norm(g))
    angle = math.acos(max(-1.0, min(1.0, cosine)))

    print(f"{case}: equilibrium after swinging "
          f"{math.degrees(math.acos(-sum(a*gi for a, gi in zip(p['arm0'], g))/(norm(p['arm0'])*norm(g)))):.0f} deg "
          f"(t = {last[0]:g} s)")
    force_error = norm([f + p["m"]*gi for f, gi in zip(force, g)])/weight
    ok = force_error <= 1e-2
    results.append(ok)
    print(f"  {'PASS' if ok else 'FAIL'}  {'force from beam = -m g':<34} error {force_error:.2e} of mg  (tol 1e-02)")
    ok = angle <= 1e-3
    results.append(ok)
    print(f"  {'PASS' if ok else 'FAIL'}  {'centre of mass below attachment':<34} misalignment {angle:.2e} rad  (tol 1e-03)")


def free_rotation(case="freeRotation"):
    """Torque-free rotation: angular momentum and kinetic energy conserved"""
    p = body_parameters(case)
    rows = read_history(case)
    J = p["J"]

    def momentum_energy(row):
        R = rotation_matrix(row[17:20])
        omega = row[20:23]
        body = [sum(R[k][i]*omega[k] for k in range(3)) for i in range(3)]
        L = mat_vec(R, [J[i]*body[i] for i in range(3)])
        return L, 0.5*sum(omega[i]*L[i] for i in range(3))

    L0 = [J[i]*p["omega0"][i] for i in range(3)]
    E0 = 0.5*sum(J[i]*p["omega0"][i]**2 for i in range(3))
    dL = max(norm([a - b for a, b in zip(momentum_energy(r)[0], L0)]) for r in rows)/norm(L0)
    dE = max(abs(momentum_energy(r)[1] - E0) for r in rows)/E0
    turns = sum(norm(r[20:23]) for r in rows)*(rows[1][0] - rows[0][0])/(2*math.pi)

    print(f"{case}: torque-free rotation, about {turns:.1f} turns")
    check_abs = lambda name, value, tol: (results.append(value <= tol), print(
        f"  {'PASS' if value <= tol else 'FAIL'}  {name:<34} max change {value:.2e}  (tol {tol:.0e})"))
    check_abs("angular momentum (global)", dL, 1e-5)
    check_abs("kinetic energy", dE, 1e-6)


def compound_pendulum(case="compoundPendulum"):
    """Rod plus pinned body: the two linear pendulum mode frequencies"""
    p = body_parameters(case)
    rows = read_history(case)
    gx, _, gz = p["g"]
    g = math.hypot(gx, gz)
    m, L = p["m"], p["L"]
    d = norm(p["arm0"])
    Jc = p["J"][1]
    k = p["E"]*p["A"]/L
    tension = m*g
    length = L + tension/k - math.sqrt(p["E"]*p["I"]/tension)

    # Linearised double pendulum: rod angle and body angle
    M = [[m*length**2, m*length*d], [m*length*d, m*d*d + Jc]]
    K = [[m*g*length, 0.0], [0.0, m*g*d]]
    a = M[0][0]*M[1][1] - M[0][1]**2
    b = -(K[0][0]*M[1][1] + K[1][1]*M[0][0])
    c = K[0][0]*K[1][1]
    omega1, omega2 = (math.sqrt((-b - s*math.sqrt(b*b - 4*a*c))/(2*a)) for s in (1, -1))

    times = [r[0] for r in rows]
    rod = [math.atan2(r[16], L + r[14]) for r in rows]
    # Body angle about -y from its rotation vector (rotation stays in the xz plane)
    relative = [-r[18] - phi for r, phi in zip(rows, rod)]

    def residual(signal, omegas):
        n = len(omegas) + 1
        N = [[0.0]*n for _ in range(n)]
        y = [0.0]*n
        for t, s in zip(times, signal):
            f = [1.0] + [math.cos(w*t) for w in omegas]
            for i in range(n):
                y[i] += f[i]*s
                for j in range(n):
                    N[i][j] += f[i]*f[j]
        for i in range(n):
            for j in range(i + 1, n):
                factor = N[j][i]/N[i][i]
                for kk in range(n):
                    N[j][kk] -= factor*N[i][kk]
                y[j] -= factor*y[i]
        coeff = [0.0]*n
        for i in reversed(range(n)):
            coeff[i] = (y[i] - sum(N[i][kk]*coeff[kk] for kk in range(i + 1, n)))/N[i][i]
        return sum((s - coeff[0] - sum(ci*math.cos(w*t) for ci, w in zip(coeff[1:], omegas)))**2
                   for t, s in zip(times, signal))

    def best_frequency(f, lo, hi):
        ratio = (math.sqrt(5) - 1)/2
        x1, x2 = hi - ratio*(hi - lo), lo + ratio*(hi - lo)
        for _ in range(50):
            if f(x1) < f(x2):
                hi = x2
            else:
                lo = x1
            x1, x2 = hi - ratio*(hi - lo), lo + ratio*(hi - lo)
        return 0.5*(lo + hi)

    fit1 = best_frequency(lambda w: residual(rod, [w, omega2]), 0.95*omega1, 1.05*omega1)
    fit2 = best_frequency(lambda w: residual(relative, [omega1, w]), 0.95*omega2, 1.05*omega2)

    print(f"{case}: rod + body double pendulum, modes {omega1:.4f} and {omega2:.4f} rad/s")
    check("mode 1 frequency (rod angle)", fit1, omega1, 3e-3)
    check("mode 2 frequency (body - rod angle)", fit2, omega2, 1e-2)


def coupling_agreement(case, rel_tol=1e-4):
    """Monolithic and converged partitioned coupling solve the same equations"""
    partitioned = read_history(case)
    monolithic = read_history(os.path.join("monolithic", case))
    n = min(len(partitioned), len(monolithic))
    if len(partitioned) != len(monolithic):
        raise ValueError(f"{len(partitioned)} partitioned rows, {len(monolithic)} monolithic")

    print(f"{case}: monolithic against partitioned, {n} time steps")
    columns = [(1, "displacement x"), (3, "displacement z"), (10, "beam force x"), (12, "beam force z")]
    if len(partitioned[0]) > 17:
        columns += [(17, "rotation x"), (18, "rotation y"), (19, "rotation z"), (24, "beam torque y")]
    for col, name in columns:
        # Skip quantities that are round-off in both runs
        peak = max(abs(r[col]) for r in partitioned)
        if peak < 1e-6:
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

    Central differences with a small step: the physically significant
    entries agree to 1e-7 or better; the tolerance allows for round-off on
    nearly-zero entries."""
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
    (hanging_body, "hangingBody"),
    (free_rotation, "freeRotation"),
    (compound_pendulum, "compoundPendulum"),
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
for case in ("axialOscillation", "pendulum", "hangingBody", "freeRotation", "compoundPendulum"):
    try:
        jacobian_check(case)
    except (OSError, ValueError) as error:
        results.append(False)
        print(f"  FAIL  {case} Jacobian check: {error}")
print()

print("=== monolithic against partitioned ===")
print()
# The compound pendulum's undamped axial mode pumps energy into the body mode
# in beats, which amplifies small differences between the two couplings
agreement_tolerance = {"compoundPendulum": 5e-4}

for _, case in tests:
    try:
        coupling_agreement(case, agreement_tolerance.get(case, 1e-4))
    except (OSError, KeyError, ValueError) as error:
        results.append(False)
        print(f"  FAIL  {case}: {error}")
    print()

print(f"{sum(results)} of {len(results)} checks passed")
sys.exit(0 if all(results) else 1)
