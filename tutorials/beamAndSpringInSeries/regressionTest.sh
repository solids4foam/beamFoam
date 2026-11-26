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
#     Regression test for the beamAndSpringInSeries tutorial: a very stiff
#     beam, with a spring of stiffness 1000 N/m at its right end, is moved by
#     its left end to -1 m in x.
#
#     Checks that the beam moves with its left end, that the root force is
#     the spring force, 1000 N/m times the end displacement, and the final
#     results against references from main (eae94a2, identical on OpenFOAM
#     v2306/clang and v2312/gcc).
#
#     $ ./regressionTest.sh [--check-only]
#
#------------------------------------------------------------------------------

cd "${0%/*}" || exit 1
. ../RegressionFunctions

regression::run "$@"

disp=postProcessing/0/beamDisplacements_right.dat
forces=postProcessing/0/beamForcesMoments_left.dat

# The left end is moved by -t in x, and the beam hardly stretches
endError=$(awk '
    $1 ~ /^[0-9.]+([eE][-+]?[0-9]+)?$/ {
        e = $2 + $1
        if (e < 0) e = -e
        if (e > max) max = e
    }
    END { print max + 0 }
' "$disp")
regression::checkRange "Largest |w_x + t| at the right end" "$endError" 0 1e-6

regression::check "Root force Q_x at t = 0.5 vs spring" \
    "$(regression::column "$forces" Q_x at 0.5)" -500 1e-5
regression::check "Final root force Q_x vs spring" \
    "$(regression::column "$forces" Q_x final)" -1000 1e-5
regression::check "Final right end w_x" \
    "$(regression::column "$disp" w_x final)" -1 1e-5

regression::checkRange "Newton iterations (at most 10% more)" \
    "$(regression::column postProcessing/0/beamConvergenceData.dat \
        nNewton sum)" 0 22

regression::finish

#------------------------------------------------------------------------------
