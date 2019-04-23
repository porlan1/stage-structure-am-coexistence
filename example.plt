reset
set terminal svg
set multiplot layout 2,4
set tmargin at screen 0.85
#set rmargin 0.05
set label 'Unstructured population' at screen 0.25,0.95 center front
set label 'A' center at 75,6.3
set title 'Species 1 alone'
unset key
set border lw 1.2
#set tics no mirror
set yrange[0:6]
set xrange[0:80]
set size ratio 6/80
set xlabel 'Days'
set ylabel 'Biomass (mg)'
plot "CASE0_sp1_alone.dat" u 1:2 w l lw 2 lc rgb "black", \
"CASE0_sp1_alone.dat" u 1:3 w l lw 2 lc rgb "blue"

