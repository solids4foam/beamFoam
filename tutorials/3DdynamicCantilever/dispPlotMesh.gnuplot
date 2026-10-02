# === pdf OUTPUT ===
set term pdfcairo dashed enhanced
set output 'dispPlotVaryingMesh.pdf'

# --- Setting graph ticks and grid ---- #
set grid
set size ratio 0.75
set xrange [0:1]
# set yrange [0:3]
set xtics
set ytics

set xlabel 'Time (s)'
set ylabel 'Displacement (m)'

#set key inside top center
set key at 0.87, 2.93

plot \
    'dispResults/beamDisplacements_right_n5.dat' u 1:($5) every 10 w lp pt 6 ps 0.6 lc "blue" lw 2  t "nSeg = 5", \
    'dispResults/beamDisplacements_right_n10.dat' u 1:($5) every 10 w lp pt 6 ps 0.4 lc "red" lw 2  t "nSeg = 10", \
    'dispResults/beamDisplacements_right_n20.dat' u 1:($5) every 10 w lp pt 6 ps 0.2 lc "green" lw 2 t "nSeg = 20", \
   'solidPointDisplacement_pointDisp.dat' u 1:($5) every 10 w lp pt 6 ps 0.1 lc "orange" lw 2 t "solids4Foam-JFNK = 4.2 mm", \
   'cantilever_abaqus_50percent.txt' u 1:2 w l lc "magenta" lw 2 t "Abaqus (C3D8)" 
   
# === EPS OUTPUT ===
#set terminal postscript eps color enhanced font 'Arial,16'
#set output 'dispPlotVaryingMesh.eps'

#replot
