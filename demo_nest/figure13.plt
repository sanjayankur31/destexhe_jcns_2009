set term pngcairo enhanced font "OpenSans" size 1500,1500
f(x) = x > -40 ? 40 : x

set output outputfile."-rasters.png"
unset border
set xrange[0:5000]
set yrange[-1:2500]
set ytics 5 nomirror
set xtics nomirror
set ytics 100
set ylabel "Cell number"
plot 'spike_detector-'.sd.'-0.gdf' using 2:1 with points pt 7 title "Raster plot";
