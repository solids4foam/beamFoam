#!/bin/bash

# List of mesh segment values to test
mesh_sizes=(5 10 20 40)

# File paths
beam_properties_file="constant/beamProperties"
log_file="log.beamFoam"
displacement_file="postProcessing/0/beamDisplacements_right.dat"

# Create a summary file for timing results
echo "nSegments ExecutionTime(s) ClockTime(s)" > timing_summary.txt

for n in "${mesh_sizes[@]}"
do
    echo "Running case with nSegments = $n"

    # Modify nSegments in beamProperties file
    sed -i "s/^\([ \t]*nSegments[ \t]*\)[0-9]\+;/\1$n;/" "$beam_properties_file"

    # Clean and run simulation
    ./Allclean

    createBeamMesh

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

# Generate a Load.dat file for plotting results
# Input and output file
input="dispResults/beamDisplacements_right_n20.dat"
output="AppliedLoad.dat"
factor=134 # 134kN l=total load

# Extract only the load column into a new file
awk -v f=$factor '
NR==1 { print "load"; next }   # header
      { print $1 * f }         # load = time * factor
' "$input" > "$output"

# Plot the load vs displacement graph for different mesh sizes
gnuplot dispPlot.gnuplot

echo "Finished Plotting Load Vs Displacement Graph"

open dispComparison_Wx_MeshSizes.pdf
open dispComparison_Wz_MeshSizes.pdf
