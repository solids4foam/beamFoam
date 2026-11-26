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
#     Regression test for the cantileverSpring tutorial: a static cantilever
#     with a linear spring at the free end
#     (linearSpringForceBeamDisplacementNR), solved in one step.
#
#     Checks the final results against references from main (eae94a2,
#     identical on OpenFOAM v2306/clang and v2312/gcc).
#
#     $ ./regressionTest.sh [--check-only]
#
#------------------------------------------------------------------------------

cd "${0%/*}" || exit 1
. ../RegressionFunctions

regression::run "$@"

disp=postProcessing/0/beamDisplacements_right.dat
forces=postProcessing/0/beamForcesMoments_left.dat

regression::check "Tip w_x" \
    "$(regression::column "$disp" w_x final)" -0.313862 1e-5
regression::check "Tip w_y" \
    "$(regression::column "$disp" w_y final)" -2.36527 1e-5
regression::check "Root force Q_y" \
    "$(regression::column "$forces" Q_y final)" 8107.79 1e-5
regression::check "Root moment M_z" \
    "$(regression::column "$forces" M_z final)" 30438.9 1e-5

regression::checkRange "Newton iterations (at most 10% more)" \
    "$(regression::column postProcessing/0/beamConvergenceData.dat \
        nNewton sum)" 0 642

regression::finish

#------------------------------------------------------------------------------
