set term pngcairo enhanced font "OpenSans" size 1500,1500
set output "Figure3.png"
# set multiplot layout 1, 2 title "Figure 3 - irregular oscillatory states"

unset border
set yrange[-1:21]
set ytics 5 nomirror
set xtics nomirror
set xlabel "Cell number"
plot 'spike_detector-22-0.gdf' using 2:1 with points pt 7 title "Raster plot";
