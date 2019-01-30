
set title 'Tube center positions'
set xlabel 'x [cm]'
set ylabel 'z [cm]'
set zlabel 'y [cm]'
set term qt 1
splot 'T1X.dat' u 1:3:2 t 'T1X', 'T1U.dat' u 1:3:2 t 'T1U', 'T2X.dat' u 1:3:2 t 'T2X', 'T2V.dat' u 1:3:2 t 'T2V'
set term qt 2
splot 'T3aX.dat' u 1:3:2 t 'T3aX', 'T3bX.dat' u 1:3:2 t 'T3bX', 'T3cX.dat' u 1:3:2 t 'T3cX', 'T3dX.dat' u 1:3:2 t 'T3dX', 'T4aX.dat' u 1:3:2 t 'T4aX', 'T4bX.dat' u 1:3:2 t 'T4bX', 'T4cX.dat' u 1:3:2 t 'T4cX', 'T4dX.dat' u 1:3:2 t 'T4dX'
set term qt 3
plot 'T3aX.dat' u 1:3 t 'T3aX', 'T3bX.dat' u 1:3 t 'T3bX', 'T3cX.dat' u 1:3 t 'T3cX', 'T3dX.dat' u 1:3 t 'T3dX', 'T4aX.dat' u 1:3 t 'T4aX', 'T4bX.dat' u 1:3 t 'T4bX', 'T4cX.dat' u 1:3 t 'T4cX', 'T4dX.dat' u 1:3 t 'T4dX'

