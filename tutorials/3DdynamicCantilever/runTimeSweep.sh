#!/bin/bash

# List of time dt values to test
dt_sizes=(0.1 0.05 0.01 0.001 0.0001)
write_interval=(1 1 10 100 1000)

# File paths
controlDict="system/controlDict"
log_file="log.beamFoam"
displacement_file="postProcessing/0/beamDisplacements_right.dat"

# Create a summary file for timing results
echo "dt ExecutionTime(s) ClockTime(s)" > dt_timing_summary.txt

for i in "${!dt_sizes[@]}"
do
    n=${dt_sizes[$i]}
    wInt=${write_interval[$i]}

    echo "Running case for dt = $n"

    # Modify deltaT in controlDict file
    sed -i "s/^\([ \t]*deltaT[ \t]*\)[0-9.eE+-]\+;/\1$n;/" "$controlDict"

    # Modify writeInterval in controlDict
    sed -i "s/^\([ \t]*writeInterval[ \t]*\)[0-9.eE+-]\+;/\1$wInt;/" "$controlDict"

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
    echo "$n $exec_time $clock_time" >> dt_timing_summary.txt

    # Create results folder if it doesn't exist
    mkdir -p dispResultsTime
    
    # Save displacement data with mesh size in name
    if [ -f "$displacement_file" ]; then
        cp "$displacement_file" "dispResultsTime/beamDisplacements_right_Dt_$n.dat"
    else
        echo "Warning: Displacement file not found for dt  = $n"
    fi

    echo "Finished for dt  = $n"
    echo "------------------------"
done

echo "All cases completed. See dt_timing_summary.txt for performance."

# Plotting displacement vs time for varying dt values
gnuplot dispPlotDT.gnuplot

# Displacement Plot name
open dispPlotVaryingDT.pdf
