#!/bin/bash
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Authors
#     Philip Cardiff, UCD
#     Seevani Bali, UCD
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "| beamFoamCleanScripts start                                         |"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo

# 1. Removing files to set the initial configuration of beams
# The directory where the fields related to initial configuration of beams
# are located
DIR_NAME_1="0"

# Files
file_pattern1="ref*"
file_pattern2="pointW"

# If files from previously ran simulation are found, then remove them
if [[ -n $(find "${DIR_NAME_1}" -type f -name ${file_pattern1}) || -n $(find "${DIR_NAME_1}" -type f -name ${file_pattern2}) ]]
then
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "The following residual files to set the initial configuration of"; 
    echo "beams from the previous simulation are being removed from ${DIR_NAME_1}/"; echo
    find "${DIR_NAME_1}" -type f -name "${file_pattern1}" -print -exec rm {} \;
    find "${DIR_NAME_1}" -type f -name "${file_pattern2}" -print -exec rm {} \;
    echo
    echo "Successfully removed";
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
fi

# 2. Removing boundary field file from constant/polyMesh
DIR_NAME_2="constant/polyMesh"

if [[ -n $(find "${DIR_NAME_2}" -name boundary) ]]
then
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    find "${DIR_NAME_2}" -type f -name boundary -print -exec rm {} \;
    echo "Removing boundary field file from ${DIR_NAME_2}"; echo
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
fi

# 3. Removing the files stored in the history/ folder
DIR_NAME_3="history"
SUBDIR_1="0"

if [[ -d "${DIR_NAME_3}" ]]
then
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    find "${DIR_NAME_3}" -type d -name "${SUBDIR_1}" -exec echo "Removing directory: {}" \; -exec rm -rf {} +
    echo "Removing (if present) previous simulation data files from ${DIR_NAME_3}/${SUBDIR_1} folder";
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
fi

echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "| beamFoamCleanScripts end                                           |"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo