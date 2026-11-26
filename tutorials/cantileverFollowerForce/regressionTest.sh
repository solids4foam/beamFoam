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
#     Regression test for the cantileverFollowerForce tutorial: a follower
#     force at the free end increases to 134 kN over pseudo-time 0 to 1.
#
#     Checks the tip displacements against Simo's results
#     (referenceResultsSimo.dat, digitised, so within 3% of L = 100 m), and
#     the final results against references from main (eae94a2, identical on
#     OpenFOAM v2306/clang and v2312/gcc).
#
#     $ ./regressionTest.sh [--check-only]
#
#------------------------------------------------------------------------------

cd "${0%/*}" || exit 1
. ../RegressionFunctions

regression::run "$@"

disp=postProcessing/0/beamDisplacements_right.dat
forces=postProcessing/0/beamForcesMoments_left.dat

# The load is 134 kN times the pseudo-time; the reference gives w_z against
# load in columns 1 and 2, and -w_x against load in columns 3 and 4
regression::checkRange "Tip w_z vs Simo, largest difference (m)" \
    "$(python3 -I ../regressionCurveError "$disp" w_z \
        referenceResultsSimo.dat 1 2 --x-scale 134)" 0 3
regression::checkRange "Tip -w_x vs Simo, largest difference (m)" \
    "$(python3 -I ../regressionCurveError "$disp" w_x \
        referenceResultsSimo.dat 3 4 --x-scale 134 --y-scale -1)" 0 3

regression::check "Final tip w_x" \
    "$(regression::column "$disp" w_x final)" -79.7138 1e-5
regression::check "Final tip w_z" \
    "$(regression::column "$disp" w_z final)" -49.3637 1e-5
regression::check "Largest tip |w|" \
    "$(regression::column "$disp" 'mag(w)' max)" 124.669 1e-5
regression::check "Final root force Q_x" \
    "$(regression::column "$forces" Q_x final)" 91471.7 1e-5
regression::check "Final root force Q_z" \
    "$(regression::column "$forces" Q_z final)" -97923.1 1e-5
regression::check "Final root moment M_y" \
    "$(regression::column "$forces" M_y final)" -2.52889e+06 1e-5

regression::checkRange "Newton iterations (at most 10% more)" \
    "$(regression::column postProcessing/0/beamConvergenceData.dat \
        nNewton sum)" 0 5040

regression::finish

#------------------------------------------------------------------------------
