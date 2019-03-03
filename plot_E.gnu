set term postscript enhanced color eps 
set output "energies.eps"

if (!exists("MP_LEFT"))   MP_LEFT = .1
if (!exists("MP_RIGHT"))  MP_RIGHT = .95
if (!exists("MP_BOTTOM")) MP_BOTTOM = .15
if (!exists("MP_TOP"))    MP_TOP = .95
if (!exists("MP_GAP"))    MP_GAP = 0.05
set key right bottom

set multiplot layout 3,1 columnsfirst  \
margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP

set xr [0:500]
set yr [100:700]

set ylabel "E_{tot} (kJ/mol)"
unset xtics
plot "data_EK_EP_T_P.dat" u 1:(($2+$3)/1000) w l lt 1 lw 2 t ""


set xr [0:500]
set yr [100:700]
set ylabel "E_{kin} (kJ/mol)"
unset xtics
plot "data_EK_EP_T_P.dat" u 1:($2/1000) w l lt 1 lw 2 t ""


set xtics
set ylabel
set xr [0:500]
set yr [-2:20]
set ylabel "E_{pot} (kJ/mol)"
set xlabel "Time (picoseconds)"
plot "data_EK_EP_T_P.dat" u 1:($3/1000)  w l lt 3 lw 2 t ""
set xtics
unset multiplot


