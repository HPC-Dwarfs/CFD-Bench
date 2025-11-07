set terminal png size 1800,768 enhanced font ,12
set output 'residual.png'
set datafile separator whitespace
set xlabel "Timestep"
set ylabel "Residual"

set logscale y 2

plot 'residual.dat' using 1:2 title "Residual"