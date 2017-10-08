set term pngcairo enhanced font "OpenSans" size 1500,1500
f(x) = x > -40 ? 40 : x

set output outputfile."-rasters.png"
unset border
set yrange[-1:21]
set ytics 5 nomirror
set xtics nomirror
set xlabel "Cell number"
plot 'spike_detector-'.sd.'-0.gdf' using 2:1 with points pt 7 title "Raster plot";

set output outputfile."-voltage.png"
set multiplot layout 4, 1 title "Figure 1 - different cell spiking behaviours"
#set xrange [1000]
set ytics 40
set yrange [-100:60]
plot "voltmeter-".v1."-0.dat" using ($2):(f($3)) with lines lw 2 title "RE 1"
plot "voltmeter-".v2."-0.dat" using ($2):(f($3)) with lines lw 2 title "RE 2"
plot "voltmeter-".v3."-0.dat" using ($2):(f($3)) with lines lw 2 title "TC 1"
plot "voltmeter-".v4."-0.dat" using ($2):(f($3)) with lines lw 2 title "TC 2"
