
plot_files=./*eps
data_files=./*dat
mpif90=mpif90
mpirun=mpirun
F90=gfortran

ifeq ($(numproc),)
numproc :=2
endif


##
## to execute the program --> type:
##                                     >> make numproc=N
##                                     where N is the number of processors (the default value is 2)    
                             
execute: program.exe
	@echo 'Executing the main program with $(numproc) processors...'
	@$(mpirun) -np $(numproc) ./$^

##
## program.exe: obtain the program executable
program.exe: sub_int_PBC_therm.f90 sub_fuerzas.f90 sub_kinetic_press.f90 sub_initialize.f90 sub_gr.f90 sub_meansquare.f90 main.f90
	@$(mpif90) -o $@ $^

##
## binning: obtain statistics, the output is written in binning.dat

ifeq ($(wildcard $(data_files)),)
binning: 
	@echo 'There are not .dat files...'
	@echo 'Please, execute the main program to obtain them: type "make"'
else
binning:
	@rm -fr binning_dir
	@mkdir binning_dir
	@$(F90) binning.f90
	@./a.out
	@rm a.out
	@mv binning.dat binning_dir
endif

##
## backup: compress some results: .gnu, .dat and .eps files

ifeq ($(wildcard $(plot_files)),)
backup:
	@echo 'There are not .eps files'
	@echo 'Please, type "make plot" to get them'
else
backup:
	@rm -f RESULTS.tar
	@echo 'Compressing the following files'
	@tar -cvf RESULTS.tar *.eps *.dat *.gnu
endif

##
## plot: obtain the plots of the program
ifeq ($(wildcard $(data_files)),)
plot:
	@echo 'There are not .dat files...'
	@echo 'Please, execute the main program to obtain them: type "make"'
else
plot:
	@gnuplot *.gnu
endif

.PHONY: help
help:
	@sed -n 's/^##//p' Makefile

##        
## clean: remove auto-generated files
##
.PHONY: clean
clean:
	@rm -f *.o
	@rm -f *.dat
	@rm -f *.eps
	@rm -f *.exe
	@rm -f *.log
	@rm -f *.mod
