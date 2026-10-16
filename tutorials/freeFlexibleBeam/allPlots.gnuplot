#----------------------------------------
# Define file paths
energyDataFile = 'postProcessing/0/beamEnergyData.dat'
dispDataFile   = 'postProcessing/0/beamDisplacements_right.dat'

#========================================
# 1. Plot Energy Data — PDF
#========================================
set term pdfcairo dashed enhanced
set output 'energyPlot.pdf'

#set title "Energy Components vs Time"
set xlabel "Time (s)"
set ylabel "Energy"
set key outside top right
plot \
    energyDataFile using 1:2 with lines dt 1 lc rgb "red"   lw 4 title 'E_{int}', \
    ''            using 1:3 with lines dt 1 lc rgb "blue"  lw 4 title 'E_{kin}', \
    ''            using 1:6 with lines dt 1 lc rgb "black" lw 6 title 'E_{tot}'

#========================================
# 2. Plot Energy Data — EPS
#========================================
#set terminal postscript eps enhanced color font 'Arial,24' linewidth 2
#set output 'energyPlot.eps'

#replot

#========================================
# 3. Plot Displacement Data — PNG
#========================================
set term pdfcairo dashed enhanced
set output 'displacementPlot.pdf'

set grid
#set title "Beam Displacement at Right End vs Time"
set xlabel "Time (s)"
set ylabel "Displacement at Point A (m)"
set key inside top left
plot \
    dispDataFile  u 1:($2) every 50 w lp pt 6 ps 0.4 lc "blue" lw 2  t '{w}_{x}', \
    dispDataFile  u 1:($3) every 50 w lp pt 4 ps 0.4 lc "dark-green" lw 2 dt 2 t '{w}_y', \
    dispDataFile  u 1:($4) every 50 w lp pt 2 ps 0.4 lc "red" lw 2 dt 3 t '{w}_z'

#========================================
# 4. Plot Displacement Data — EPS
#========================================
#set terminal postscript eps enhanced color font 'Arial,24' linewidth 2
#set output 'displacementPlot.eps'

#replot
