# === pdf OUTPUT ===
set term pdfcairo dashed enhanced
set output 'dispComparison_Wz_MeshSizes.pdf'

# --- Setting graph ticks and grid ---- #
set grid
set size ratio 0.75
set xrange [0:140]
# set yrange [-60:120]
set xtics
set ytics

set xlabel 'Load (kN)'
set ylabel 'Displacement (m)'

set key inside top right
#set key at 120, 75

plot \
    '< paste AppliedLoad.dat dispResults/beamDisplacements_right_n5.dat' u 1:($5) every 20 w lp pt 6 ps 1 lc "blue" lw 2  t '{w}_{z}, nSeg = 5', \
    '< paste AppliedLoad.dat dispResults/beamDisplacements_right_n10.dat' u 1:($5) every 20 w lp pt 6 ps 0.7 lc "green" lw 2 t '{w}_z, nSeg = 10', \
    '< paste AppliedLoad.dat dispResults/beamDisplacements_right_n20.dat' u 1:($5) every 20 w lp pt 6 ps 0.4 lc "red" lw 2 t '{w}_z, nSeg = 20', \
    'referenceResultsSimo.dat' u 1:2 w l lc "black" lw 2 t "ref w_z"

# === pdf OUTPUT ===
set term pdfcairo dashed enhanced
set output 'dispComparison_Wx_MeshSizes.pdf'


plot \
    '< paste AppliedLoad.dat dispResults/beamDisplacements_right_n5.dat' u 1:(-$3) every 20 w lp pt 6 ps 1 lc "blue" lw 2  t '{w}_{x}, nSeg = 5', \
    '< paste AppliedLoad.dat dispResults/beamDisplacements_right_n10.dat' u 1:(-$3) every 20 w lp pt 6 ps 0.7 lc "green" lw 2 t '{w}_x, nSeg = 10', \
    '< paste AppliedLoad.dat dispResults/beamDisplacements_right_n20.dat' u 1:(-$3) every 20 w lp pt 6 ps 0.4 lc "red" lw 2 t '{w}_x, nSeg = 20', \
    'referenceResultsSimo.dat' u 3:4 w l lc "black" lw 2 t "ref w_x"

# === EPS OUTPUT ===
#set terminal postscript eps color enhanced font 'Arial,16'
#set output 'dispComparison_wz_MeshSizes.eps'

#replot
