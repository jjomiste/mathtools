.SUFFIXES : .o .f90
.f90.o: ;  gfortran -c $<

all:: fft wigner dft windowdft dir

fft:: mtmodules.o fftcode.o
	gfortran  mtmodules.o fftcode.o -o $@

wigner:: mtmodules.o wignercode.o
	gfortran  mtmodules.o wignercode.o -o $@

dft:: mtmodules.o dftcode.o
	gfortran  mtmodules.o dftcode.o -o $@

windowdft:: mtmodules.o windowdftcode.o
	gfortran  mtmodules.o windowdftcode.o -o $@

dir:: 
	mkdir bin
	mv fft bin
	mv wigner bin
	mv dft bin
	mv windowdft bin

clean::
	rm -f *.o
	rm -f *.mod
	rm -rf bin

