#!/bin/bash

#### This is a shell script to run the test case for varying mesh sizes

# List of mesh segment values to test
mesh_sizes=(5 10 20 40)

# File paths
beam_properties_file="constant/beamProperties"
log_file="log.beamFoam"
displacement_file="postProcessing/0/beamDisplacements_right.dat"

schemesfile="system/fvSchemes"

# Create a summary file for timing results
echo "nSegments ExecutionTime(s) ClockTime(s)" > timing_summary.txt

for n in "${mesh_sizes[@]}"
do
    echo "Running case with nSegments = $n"

    # Modify nSegments in beamProperties file
    sed -i "s/^\([ \t]*nSegments[ \t]*\)[0-9]\+;/\1$n;/" "$beam_properties_file"

    # Clean and run simulation
    ./Allclean

    # Mesh utility
    createBeamMesh

    # Secondary utility to set initial configuration
    setInitialPositionBeam -cellZone beam_0 -translate '(0 0 0)' -rotateAngle '((0 -1 0) 90)'

    # run beamFoam
    beamFoam

    # Extract execution and clock time from log
    exec_time=$(grep "ExecutionTime =" "$log_file" | tail -1 | awk '{print $3}')
    clock_time=$(grep "ClockTime =" "$log_file" | tail -1 | awk '{print $7}')

    # Save timing data
    echo "$n $exec_time $clock_time" >> timing_summary.txt

    # Create results folder if it doesn't exist
    mkdir -p dispResults

    # Save displacement data with mesh size in name
    if [ -f "$displacement_file" ]; then
        cp "$displacement_file" "dispResults/beamDisplacements_right_n$n.dat"
    else
        echo "Warning: Displacement file not found for nSegments = $n"
    fi

    echo "Finished nSegments = $n"
    echo "------------------------"
done

echo "All cases completed. See timing_summary.txt for performance."

# Plotting the displacement variation for various mesh sizes
gnuplot dispPlotMesh.gnuplot

echo "Finished plotting displacements vs time for varying mesh sizes."

open dispPlotVaryingMesh.pdf
