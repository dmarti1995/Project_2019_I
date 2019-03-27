#gfortran -c sub_fuerzas.f90
#mpif90 -c sub_periodic_paralel.f90
#mpif90 -c sub_integrate.f90
mpif90 -c sub_int_PBC_therm.f90
mpif90 -c sub_meansquare.f90
mpif90 -c sub_gr.f90
mpif90 -c sub_initialize.f90
mpif90 -c main.f90
mpif90 -c sub_fuerzas_paralelizadas.f90
mpif90 -c sub_kinetic_press.f90
#mpif90 sub_integrate.o sub_periodic_paralel.o sub_fuerzas.o  sub_initialize.o sub_gr.o sub_meansquare.o main.o
#mpif90 sub_int_PBC_therm.o sub_fuerzas.o sub_initialize.o sub_gr.o sub_meansquare.o main.o
mpif90 sub_int_PBC_therm.o sub_fuerzas_paralelizadas.o sub_initialize.o sub_gr.o sub_meansquare.o main.o sub_kinetic_press.o
