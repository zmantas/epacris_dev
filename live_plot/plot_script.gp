set terminal png size 800,600
set output 'Results/55cnce/helios_benchmark/55cnce_comp5_f067//TP_profile_step_150.png'
set logscale y
set yrange [200:1e-8] reverse
set xrange [800:3200]
set xlabel 'Temperature (K)'
set ylabel 'Pressure (Bar)'
set title 'Temperature-Pressure Profile (Step 150)'
set grid
plot 'live_plot/temp_data.txt' using 1:3 with lines lw 3 lc rgb 'black' title 'T_{new}', \
     'live_plot/temp_data.txt' using 2:3 with lines lw 3 dt 2 lc rgb 'red' title 'T_{old}'
