gfortran -c sub_distances.f90
gfortran -c sub_fuerzas.f90
gfortran -c sub_initialize_r.f90
#mpif90 -c sub_periodic_paralel.f90
#mpif90 -c sub_integrate.f90
mpif90 -c sub_int_PBC_therm.f90
mpif90 -c sub_meansquare.f90
mpif90 -c sub_gr.f90
mpif90 -c sub_initialize_v.f90
mpif90 -c main.f90

#mpif90 sub_integrate.o sub_periodic_paralel.o sub_fuerzas.o sub_initialize_r.o sub_initialize_v.o sub_gr.o sub_meansquare.o main.o
mpif90 sub_int_PBC_therm.o sub_fuerzas.o sub_initialize_r.o sub_initialize_v.o sub_gr.o sub_meansquare.o main.o
