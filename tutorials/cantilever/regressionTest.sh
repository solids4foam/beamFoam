#!/bin/bash
#------------------------------------------------------------------------------
# License
#     This file is part of beamFoam, licensed under the GNU General Public
#     License <http://www.gnu.org/licenses/>.
#
# Script
#     regressionTest.sh
#
# Description
#     Regression test for the cantilever tutorial: the end rotation is
#     increased to 2*pi with no end force, so the beam bends into a circle of
#     radius L/psi with a constant moment EI*psi/L.
#
#     Checks the tip path against the analytical solution, the root moment
#     against EI*psi/L, and the final results against references from
#     main (eae94a2, identical on OpenFOAM v2306/clang and v2312/gcc).
#
#     $ ./regressionTest.sh [--check-only]
#
#------------------------------------------------------------------------------

cd "${0%/*}" || exit 1
. ../RegressionFunctions

regression::run "$@"

disp=postProcessing/0/beamDisplacements_right.dat
forces=postProcessing/0/beamForcesMoments_left.dat

# Analytical tip displacement for an end rotation psi = 2*pi*t:
#     w_x = L*sin(psi)/psi - L,  w_z = L*(1 - cos(psi))/psi
# Largest distance between the computed and analytical tip, divided by L
tipError=$(awk -v L=10 '
    $1 ~ /^[0-9.]+([eE][-+]?[0-9]+)?$/ {
        psi = 2*3.14159265358979*$1
        ax = L*sin(psi)/psi - L
        az = L*(1 - cos(psi))/psi
        e = sqrt(($2 - ax)^2 + ($4 - az)^2)/L
        if (e > max) max = e
    }
    END { print max + 0 }
' "$disp")
regression::checkRange "Tip path error / L vs analytical" "$tipError" 0 1e-3

# Root moment EI*psi/L = 100*2*pi/10
regression::check "Root moment M_y vs EI*psi/L" \
    "$(regression::column "$forces" M_y final)" 62.8318531 1e-5

regression::check "Final tip w_x" \
    "$(regression::column "$disp" w_x final)" -10 1e-5
regression::check "Final tip w_z" \
    "$(regression::column "$disp" w_z final)" 0 0 1e-5
regression::check "Largest tip |w|" \
    "$(regression::column "$disp" 'mag(w)' max)" 12.5942 1e-5

regression::checkRange "Newton iterations (at most 10% more)" \
    "$(regression::column postProcessing/0/beamConvergenceData.dat \
        nNewton sum)" 0 330

regression::finish

#------------------------------------------------------------------------------
