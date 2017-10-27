reset
set grid
set log x
set log y
set xlabel "pressure [bar]"
set ylabel "Saturation limit x [ppm vol]"

set xrange [10:5000]

plot "iceExp1.dat" using ($2/1e5):3 title "Rabinovitch and Beketov (1995) I" pt 3
rep "iceExp2.dat" using ($2/1e5):3 title "Rabinovitch and Beketov (1995) II" pt 4

rep "icemap.dat" using ($1/1e5):($2*1e6) title "Extended BWR model" ls 1 w l
rep "icemap.dat" using ($1/1e5):($3*1e6) notitle ls 1 w l
rep "icemap.dat" using ($1/1e5):($4*1e6) notitle ls 1 w l
rep "icemap.dat" using ($1/1e5):($5*1e6) notitle ls 1 w l
rep "icemap.dat" using ($1/1e5):($6*1e6) notitle ls 1 w l
rep "icemap.dat" using ($1/1e5):($7*1e6) notitle ls 1 w l
 

set term fig lan mon big
set out "icemap.fig"
rep
set term x11
set out


