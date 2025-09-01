# === pdf OUTPUT ===
set term pdfcairo dashed enhanced
set output 'dispPlotVaryingDT.pdf'

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
    'dispResultsTime/beamDisplacements_right_Dt_0.1.dat' u 1:($5) every 1 w lp pt 6 ps 0.4 lc "red" lw 2  t "dt = 0.1", \
    'dispResultsTime/beamDisplacements_right_Dt_0.05.dat' u 1:($5) every 1 w lp pt 2 ps 0.4 lc "brown" lw 2 dt 2 t "dt = 0.05", \
    'dispResultsTime/beamDisplacements_right_Dt_0.01.dat' u 1:($5) every 1 w lp pt 6 ps 0.2 lc "dark-green" lw 2 dt 3 t "dt = 0.01", \
    'dispResultsTime/beamDisplacements_right_Dt_0.001.dat' u 1:($5) every 10 w lp pt 8 ps 0.3 lc "blue" lw 2 dt 2  t "dt = 0.001", \
    'dispResultsTime/beamDisplacements_right_Dt_0.0001.dat' u 1:($5) every 100 w l lc "violet" lw 2 t "dt = 0.0001"


# === EPS OUTPUT ===
#set terminal postscript eps color enhanced font 'Arial,16'
#set output 'dispPlotVaryingDT.eps'

#replot
