#----------------------------------------
# Define file paths
energyDataFile = 'postProcessing/0/beamEnergyData.dat'
dispDataFile   = 'postProcessing/0/beamDisplacements_right.dat'

#========================================
# 1. Plot Energy Data — PNG
#========================================
set term pdfcairo dashed enhanced
set output 'energyPlot.pdf'

set title "Energy Components vs Time"
set xlabel "Time (s)"
set ylabel "Energy"
set key outside top right
plot \
    energyDataFile using 1:2 with lines dt 2 lc rgb "red"   lw 4 title 'E_{int}', \
    ''            using 1:3 with lines dt 3 lc rgb "blue"  lw 4 title 'E_{kin}', \
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

set title "Beam Displacement at Right End vs Time"
set xlabel "Time (s)"
set ylabel "Displacement (m)"
set key outside top right
plot \
    dispDataFile using 1:5 with lines dt 1 lc rgb "red" lw 4 title 'mag(w)'

#========================================
# 4. Plot Displacement Data — EPS
#========================================
#set terminal postscript eps enhanced color font 'Arial,24' linewidth 2
#set output 'displacementPlot.eps'

#replot
