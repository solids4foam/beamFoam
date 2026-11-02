# Plot the beam tip streamwise displacement W_x versus time.
# Input : tipDisplacement.dat  (written by extractTipDisplacement.sh)
# Output: tipDisplacement.pdf

set terminal pdfcairo enhanced font 'Helvetica,12' size 6in,4in
set output 'tipDisplacement.pdf'

set grid
set xlabel "Time (s)"
set ylabel "Tip x-displacement W_x (m)"
set key off

plot 'tipDisplacement.dat' using 1:2 with linespoints \
     pointtype 7 pointsize 0.3 linecolor rgb '#1f77b4' linewidth 2
