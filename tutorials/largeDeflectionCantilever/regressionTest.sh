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
#     Regression test for the largeDeflectionCantilever tutorial: a thin steel
#     cantilever under its own weight and a 0.588 N end load
#     (constant/timeVsForce).
#
#     Checks the tip deflection against the experiment of Belendez et al.
#     (2003), 0.281 m, within 2%, and the final results against references
#     from main (eae94a2, identical on OpenFOAM v2306/clang and v2312/gcc).
#
#     $ ./regressionTest.sh [--check-only]
#
#------------------------------------------------------------------------------

cd "${0%/*}" || exit 1
. ../RegressionFunctions

regression::run "$@"

disp=postProcessing/0/beamDisplacements_right.dat

regression::checkRange "Tip w_z vs experiment, 0.281 m within 2%" \
    "$(regression::column "$disp" w_z final)" -0.28662 -0.27538

regression::check "Final tip w_x" \
    "$(regression::column "$disp" w_x final)" -0.143135 1e-5
regression::check "Final tip w_z" \
    "$(regression::column "$disp" w_z final)" -0.278404 1e-5
regression::check "Final strain energy" \
    "$(regression::column postProcessing/0/beamEnergyData.dat \
        'E_{int}' final)" 0.0670928 1e-5

regression::finish

#------------------------------------------------------------------------------
