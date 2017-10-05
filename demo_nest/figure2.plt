set term pngcairo enhanced font "OpenSans" size 400,1200
set output "Figure2.png"
set multiplot layout 4, 1 title "Figure 2 - self sustained TC RE network"

unset border
set size square
set yrange [-120:0]

plot "voltmeter-05-0.dat" using 2:3 with lines lw 2 title "RE"

plot "voltmeter-06-0.dat" using 2:3 with lines lw 2 title "RE"

plot "voltmeter-07-0.dat" using 2:3 with lines lw 2 title "TC"

plot "voltmeter-08-0.dat" using 2:3 with lines lw 2 title "TC"
