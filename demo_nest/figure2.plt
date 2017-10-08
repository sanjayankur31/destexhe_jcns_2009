set term pngcairo enhanced font "OpenSans" size 1200,1200
set output outputfile
set multiplot layout 4, 1 title "Figure 2 - self sustained TC RE network"

f(x) = x > -40 ? 40 : x
unset border
set xtics 200
set ytics 40
set yrange [-100:60]

plot "voltmeter-05-0.dat" using ($2):(f($3)) with lines lw 2 title "RE"

plot "voltmeter-06-0.dat" using ($2):(f($3)) with lines lw 2 title "RE"

plot "voltmeter-07-0.dat" using ($2):(f($3)) with lines lw 2 title "TC"

plot "voltmeter-08-0.dat" using ($2):(f($3)) with lines lw 2 title "TC"
