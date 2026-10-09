#----------------------------------------
# Define file paths
#dispDataFile   = 'postProcessing/0/beamDisplacements_right.dat'

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
    '< paste AppliedLoad.dat postProcessing/0/beamDisplacements_right.dat'  u 1:(-$3) every 10 w lp pt 6 ps 0.4 lc "blue" lw 2  t '{w}_{x}', \
    '< paste AppliedLoad.dat postProcessing/0/beamDisplacements_right.dat'  u 1:($5) every 10 w lp pt 2 ps 0.4 lc "red" lw 2  t '{w}_z', \
     'referenceResultsSimo.dat' u 3:4 w l lc "black" lw 2 t "ref w_x", \
     'referenceResultsSimo.dat' u 1:2 w l lc "black" lw 2 t "ref w_z"

#========================================
# 2. Plot Displacement Data — EPS
#========================================
#set terminal postscript eps enhanced color font 'Arial,24' linewidth 2
#set output 'displacementPlot.eps'

#replot
