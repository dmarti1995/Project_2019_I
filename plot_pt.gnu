set term postscript enhanced color eps 
set output "P_T.eps"

if (!exists("MP_LEFT"))   MP_LEFT = .1
if (!exists("MP_RIGHT"))  MP_RIGHT = .95
if (!exists("MP_BOTTOM")) MP_BOTTOM = .15
if (!exists("MP_TOP"))    MP_TOP = .95
if (!exists("MP_GAP"))    MP_GAP = 0.05
set key right bottom

set multiplot layout 2,1 columnsfirst  \
margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP

set xr [0:500]
set yr [0:600]

set ylabel "Temperature (K)"
unset xlabel
plot "data_EK_EP_T_P.dat" u 1:(($4)) w l lt 1 lw 2 t ""


set xr [0:500]
set yr [-1:4]
set ylabel "Pressure (Atm)"
set xlabel "Time (Picoseconds)"
plot "data_EK_EP_T_P.dat" u 1:($5) w l lt 1 lw 2 t ""




unset multiplot


