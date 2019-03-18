set term postscript enhanced color eps 
set output "RDF_MSD.eps"

if (!exists("MP_LEFT"))   MP_LEFT = .2
if (!exists("MP_RIGHT"))  MP_RIGHT = .75
if (!exists("MP_BOTTOM")) MP_BOTTOM = .15
if (!exists("MP_TOP"))    MP_TOP = .95
if (!exists("MP_GAP"))    MP_GAP = 0.15
set key right bottom

set multiplot layout 2,1 columnsfirst  \
margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP

set xr [0:60]


set ylabel "RDF g(r)"
unset xtics
set xlabel "Distance (A)"
set xtics
plot "rad_dist_func.dat" u 1:(($2)) w l lt 1 lw 2 t ""


set xr [0:500]

set xtics
f(x)=a*x+b
fit f(x) "mean_square_disp.dat" u 1:($2) via a,b
set ylabel "Mean square displacement (A^2)"
set xlabel "Time (ps)"
plot "mean_square_disp.dat" u 1:($2) w l lt 1 lw 2 t "",f(x) t "Diff coef=0.89"




unset multiplot


