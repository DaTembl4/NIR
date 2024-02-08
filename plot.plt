set terminal png size 1920, 1080 font size 18; set output 'plot.png'
set xlabel 'time, h'; set ylabel 'magnitude'
plot 'source_prime.dat' with lines linetype 3 linewidth 2 title 'experimental data'
