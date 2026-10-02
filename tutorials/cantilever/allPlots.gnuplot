#----------------------------------------
# Define file paths
dispDataFile   = 'postProcessing/0/beamDisplacements_right.dat'

#========================================
# 1. Plot Displacement Data — PNG
#========================================
set term pdfcairo dashed enhanced
set output 'displacementPlot.pdf'

set grid
#set title "Displacement components at tip vs Load"
set xlabel "Load (kN)"
set ylabel "Displacement (m)"
set key inside top right
plot \
    dispDataFile  u 1:($5) every 1 w lp pt 6 ps 0.4 lc "blue" lw 2  t '{w}_{mag}'

#========================================
# 2. Plot Displacement Data — EPS
#========================================
#set terminal postscript eps enhanced color font 'Arial,24' linewidth 2
#set output 'displacementPlot.eps'

#replot
