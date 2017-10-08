set term pngcairo enhanced font "OpenSans" size 1500,1500
set output "Figure1.png"
set multiplot layout 3, 3 title "Figure 1 - different cell spiking behaviours"

f(x) = x > -40 ? 40 : x

unset border
set size square
#set xrange [1000]
set yrange [-100:60]
unset xtics
unset ytics

set title "(a) RS cell strong adaptation"
set arrow from 0,-80 to 200,-80 nohead lw 2
set arrow from 200,-70 to 600,-70 nohead lw 2
set arrow from 600,-80 to 1000,-80 nohead lw 2
set arrow from 200,-80 to 200,-70 nohead lw 2
set arrow from 600,-80 to 600,-70 nohead lw 2
set arrow from 750,-90 to 950,-90 nohead lw 2
set label "200ms" at 800,-95
plot "voltmeter-03-0.dat" using ($2):(f($3)) with lines lw 2 title "a = 0.001{/Symbol u}S, b = 0.04nA";

set title "(b) RS cell weak adaptation"
plot "voltmeter-04-0.dat" using ($2):(f($3)) with lines lw 2 title "a = 0.001{/Symbol u}S, b = 0.005nA";

set title "(c) FS cell"
plot "voltmeter-05-0.dat" using ($2):(f($3)) with lines lw 2 title "a = 0.001{/Symbol u}S, b = 0.nA";

set title "(d - top) LTS cell"
plot "voltmeter-06-0.dat" using ($2):(f($3)) with lines lw 2 title "a = 0.02{/Symbol u}S, b = 0.nA";

set title "(e - top) TC cell"
plot "voltmeter-07-0.dat" using ($2):(f($3)) with lines lw 2 title "a = 0.04{/Symbol u}S, b = 0nA";

set title "(f - top) RE cell"
plot "voltmeter-08-0.dat" using ($2):(f($3)) with lines lw 2 title "a = 0.08{/Symbol u}S, b = 0.03nA";
unset arrow
unset label

set title "(d - bottom) LTS cell"
set arrow from 0,-80 to 200,-80 nohead lw 2
set arrow from 200,-90 to 600,-90 nohead lw 2
set arrow from 600,-80 to 1000,-80 nohead lw 2
set arrow from 200,-90 to 200,-80 nohead lw 2
set arrow from 600,-90 to 600,-80 nohead lw 2
set arrow from 750,-90 to 950,-90 nohead lw 2
set label "200ms" at 800,-95
plot "voltmeter-09-0.dat" using ($2):(f($3)) with lines lw 2 title "a = 0.05{/Symbol u}S, b = 0.nA";

set title "(e - bottom) TC cell"
plot "voltmeter-10-0.dat" using ($2):(f($3)) with lines lw 2 title "a = 0.024{/Symbol u}S, b = 0nA";

set title "(f - bottom) RE cell"
plot "voltmeter-11-0.dat" using ($2):(f($3)) with lines lw 2 title "a = 0.04{/Symbol u}S, b = 0.03nA";
