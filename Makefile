
plot_files=./*eps
F90=gfortran

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


## plot: obtain the plots of the program
plot: program.exe
	@gnuplot *.gnu

## program.exe: obtain the program executable
program.exe: sub_integrate.o sub_periodic.o sub_fuerzas.o sub_initialize.o sub_gr.o sub_meansquare.o main.o
	@$(F90) $^ -o $@
	@echo 'Executing the main program...'
	@./program.exe

%.o: %.f90
	@$(F90) -c -freal-4-real-8 $^

.PHONY: help
help:
	@sed -n 's/^##//p' Makefile
        
## clean: remove auto-generated files
.PHONY: clean
clean:
	@rm -f *.o
	@rm -f *.dat
	@rm -f *.eps
	@rm -f *.exe
	@rm -f *.log
