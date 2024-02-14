set terminal png size 1920, 1080 font size 18; set output 'plot.png'
#set xlabel 'time, h'; set ylabel 'magnitude'
set xlabel 'frequency, 1/h'
#set xlabel 'time, h'
#plot 'DATA\t&m_interpolated.dat' with lines linetype 3 linewidth 1 title 'interpolated_data'
#plot 'DATA\periodogramm.dat' with lines linetype 3 linewidth 2 title 'periodogramm'
plot 'DATA\periodogramm_prime.dat' with lines linetype 3 linewidth 2 title 'periodogramm'
#plot 'DATA\correlogramm.dat' with lines linetype 3 linewidth 2 title 'correlogramm', 'DATA\correlogramm_prime.dat' with lines linetype 1 linewidth 2 title 'correlogramm_prime'





