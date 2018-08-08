# Gnuplot script file for plotting data in file "force.dat"
# This file is called   force.p
set   autoscale                        # scale axes automatically
#unset log                              # remove any log-scaling
#unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
#set title "Force Deflection Data for a Beam and a Column"
#set xlabel "Deflection (meters)"
#set ylabel "Force (kN)"
#set key 0.01,100
#set label "Yield Point" at 0.003,260
#set arrow from 0.0028,250 to 0.003,280
#set xr [0.0:0.022]
#set yr [0:325]
#set datafile commentschars "()"
#plot    "Monitors/sd16_gas_flux_outlet.mon" using 1:1 title 'x' with linespoints , \
#      "Monitors/sd16_gas_flux_outlet.mon" using 1:1 title 'y' with points
plot    "Monitors/TMP/sd16_gas_flux_outlet.mon" every ::5
set term png
set out
set out "t.png"
rep


