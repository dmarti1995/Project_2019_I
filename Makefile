
## backup: compress some results: .gnu, .dat and .pdf files
backup: plot
	@rm -f Results.tar
	@echo 'Compressing the following files'
	@tar -cvf Results.tar *.pdf *.dat *.gnu

## plot: obtain the plots of the program
plot: program.exe
	@./program.exe
	@gnuplot *.gnu

program.exe: sub_integrate.o periodic.o Fuerzas.o initialize.o main.o
	@gfortran $^ -o $@

%.o: %.f90
	@gfortran -c $^

.PHONY: help
help:
	@sed -n 's/^##//p' Makefile
        
## clean: remove auto-generated files
.PHONY: clean
clean:
	@rm -f *.o
	@rm -f *.dat
	@rm -f *.pdf
	@rm -f *.exe
