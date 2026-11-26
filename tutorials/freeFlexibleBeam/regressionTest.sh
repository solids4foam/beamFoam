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
#     Regression test for the freeFlexibleBeam tutorial: a free beam is
#     forced at one end until 2.5 s, then moves freely until 40 s.
#
#     Checks that the total energy is conserved in the free motion (the
#     Newmark scheme conserves it), and the final results against references
#     from main (eae94a2, identical on OpenFOAM v2306/clang and v2312/gcc).
#
#     $ ./regressionTest.sh [--check-only]
#
#------------------------------------------------------------------------------

cd "${0%/*}" || exit 1
. ../RegressionFunctions

regression::run "$@"

disp=postProcessing/0/beamDisplacements_right.dat
energy=postProcessing/0/beamEnergyData.dat

# Spread of the total energy after the loads are removed, relative to it
eMax=$(regression::column "$energy" 'E_{tot}' max 2.5)
eMin=$(regression::column "$energy" 'E_{tot}' min 2.5)
regression::checkRange "Total energy drift after 2.5 s (relative)" \
    "$(awk -v a="$eMax" -v b="$eMin" 'BEGIN { print (a - b)/a }')" 0 1e-3

regression::check "Final total energy" \
    "$(regression::column "$energy" 'E_{tot}' final)" 130.968 1e-5
regression::check "Final tip w_x" \
    "$(regression::column "$disp" w_x final)" 74.787 1e-5
regression::check "Final tip w_y" \
    "$(regression::column "$disp" w_y final)" 8.42462 1e-5
regression::check "Final tip w_z" \
    "$(regression::column "$disp" w_z final)" -2.01391 1e-5
regression::check "Largest tip |w|" \
    "$(regression::column "$disp" 'mag(w)' max)" 75.4963 1e-5

regression::checkRange "Newton iterations (at most 10% more)" \
    "$(regression::column postProcessing/0/beamConvergenceData.dat \
        nNewton sum)" 0 13333

regression::finish

#------------------------------------------------------------------------------
