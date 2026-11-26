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
#     Regression test for the 3DdynamicCantileverMultiBeamDensity tutorial:
#     three independent cantilevers with densities 500, 1000 and 2000 under
#     the same sudden end force.
#
#     Without damping, a sudden constant load gives the same peak response
#     for any density, only later for a heavier beam. Checks that the three
#     peak tip displacements agree within 0.1%, and the final results against
#     references from main (eae94a2, identical on OpenFOAM v2306/clang and
#     v2312/gcc).
#
#     $ ./regressionTest.sh [--check-only]
#
#------------------------------------------------------------------------------

cd "${0%/*}" || exit 1
. ../RegressionFunctions

regression::run "$@"

pp=postProcessing/0

peak0=$(regression::column $pp/beamDisplacements_right_0.dat 'mag(w)' max)
peak1=$(regression::column $pp/beamDisplacements_right_1.dat 'mag(w)' max)
peak2=$(regression::column $pp/beamDisplacements_right_2.dat 'mag(w)' max)
regression::checkRange "Spread of the peak tip |w| (relative)" \
    "$(awk -v a="$peak0" -v b="$peak1" -v c="$peak2" 'BEGIN {
        hi = a; lo = a
        if (b > hi) hi = b; if (b < lo) lo = b
        if (c > hi) hi = c; if (c < lo) lo = c
        print (hi - lo)/hi
    }')" 0 1e-3

regression::check "Peak tip |w|, beam_1" "$peak1" 2.79738 1e-5
regression::check "Final tip |w|, beam_0 (rho 500)" \
    "$(regression::column $pp/beamDisplacements_right_0.dat 'mag(w)' final)" \
    1.1581 1e-5
regression::check "Final tip |w|, beam_1 (rho 1000)" \
    "$(regression::column $pp/beamDisplacements_right_1.dat 'mag(w)' final)" \
    2.69671 1e-5
regression::check "Final tip |w|, beam_2 (rho 2000)" \
    "$(regression::column $pp/beamDisplacements_right_2.dat 'mag(w)' final)" \
    0.464351 1e-5
regression::check "Final total energy" \
    "$(regression::column $pp/beamEnergyData.dat 'E_{tot}' final)" \
    9284.82 1e-5

regression::checkRange "Newton iterations (at most 10% more)" \
    "$(regression::column $pp/beamConvergenceData.dat nNewton sum)" 0 3583

regression::finish

#------------------------------------------------------------------------------
