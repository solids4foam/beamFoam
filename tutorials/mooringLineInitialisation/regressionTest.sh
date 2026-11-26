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
#     Regression test for the mooringLineInitialisation tutorial: a mooring
#     line settles onto the seabed (z = -0.5 m) under gravity and buoyancy,
#     with groundContact, morisonDrag and morisonAddedMass contributions.
#
#     Reads the final fields, at t = 2 s: the cell positions are refW + W.
#     Checks how many cells rest on the seabed (about 28 of 60, as in the
#     README), how far the penalty contact lets the line sink into it, and
#     the results against references from main (eae94a2, identical on
#     OpenFOAM v2306/clang and v2312/gcc).
#
#     $ ./regressionTest.sh [--check-only]
#
#------------------------------------------------------------------------------

cd "${0%/*}" || exit 1
. ../RegressionFunctions

regression::run "$@"

# Print one statistic of the final line: nContact, minZ, maxW or midW
lineStat()
{
    python3 -I - "$1" <<'EOF'
import re
import sys

groundZ = -0.5
radius = 1.828e-3


def cell_vectors(path):
    with open(path) as f:
        text = f.read()
    text = text[text.index("internalField"):]
    n = int(re.search(r"\n\s*(\d+)\s*\n\s*\(", text).group(1))
    number = r"([-+0-9.eE]+)"
    vectors = re.findall(r"\(" + r"\s+".join([number]*3) + r"\)", text)[:n]
    if len(vectors) != n:
        sys.exit("Could not read {} vectors from {}".format(n, path))
    return [tuple(float(x) for x in v) for v in vectors]


W = cell_vectors("2/W")
ref = cell_vectors("2/refW")
z = [r[2] + w[2] for r, w in zip(ref, W)]
mag = [sum(x*x for x in w)**0.5 for w in W]

stat = sys.argv[1]
if stat == "nContact":
    print(sum(1 for zi in z if zi - radius < groundZ))
elif stat == "minZ":
    print(min(z))
elif stat == "maxW":
    print(max(mag))
elif stat == "midW":
    print(mag[len(mag)//2])
EOF
}

nContact=$(lineStat nContact)
regression::checkRange "Cells resting on the seabed (README: about 28)" \
    "$nContact" 24 32
regression::check "Cells resting on the seabed" "$nContact" 28 0 0.5

# Penalty contact: a few mm into the seabed, not through it
regression::checkRange "Lowest cell z (m)" "$(lineStat minZ)" -0.51 -0.5
regression::check "Lowest cell z" "$(lineStat minZ)" -0.503255 1e-5

regression::check "Largest |W|" "$(lineStat maxW)" 0.246825 1e-5
regression::check "|W| of the middle cell" "$(lineStat midW)" 0.23935 1e-5

regression::finish

#------------------------------------------------------------------------------
