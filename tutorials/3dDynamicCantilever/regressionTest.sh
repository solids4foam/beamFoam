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
#     Regression test for the 3dDynamicCantilever tutorial: a cantilever
#     column under a sudden constant end force, 0 to 1 s.
#
#     Checks the tip displacement history against the Abaqus solution
#     (cantilever_abaqus_50percent.txt) within 2% of its peak, and the final
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

# 2% of the Abaqus peak, 2.8007 m
regression::checkRange "Tip |w| vs Abaqus, largest difference (m)" \
    "$(python3 -I ../regressionCurveError "$disp" 'mag(w)' \
        cantilever_abaqus_50percent.txt 1 2)" 0 0.056

regression::check "Final tip w_x" \
    "$(regression::column "$disp" w_x final)" 1.20903 1e-5
regression::check "Final tip w_y" \
    "$(regression::column "$disp" w_y final)" 1.23644 1e-5
regression::check "Final tip w_z" \
    "$(regression::column "$disp" w_z final)" -2.06362 1e-5
regression::check "Largest tip |w|" \
    "$(regression::column "$disp" 'mag(w)' max)" 2.79644 1e-5
regression::check "Final root force Q_z" \
    "$(regression::column "$forces" Q_z final)" 6031.71 1e-5
regression::check "Final root moment M_x" \
    "$(regression::column "$forces" M_x final)" 4378.73 1e-5
regression::check "Final total energy" \
    "$(regression::column postProcessing/0/beamEnergyData.dat \
        'E_{tot}' final)" 4979.04 1e-5

regression::checkRange "Newton iterations (at most 10% more)" \
    "$(regression::column postProcessing/0/beamConvergenceData.dat \
        nNewton sum)" 0 3634

regression::finish

#------------------------------------------------------------------------------
