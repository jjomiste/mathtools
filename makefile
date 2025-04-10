.SUFFIXES : .o .f90
.f90.o: ;  gfortran -c $<

fft:: mtmodules.o fftcode.o
	gfortran  mtmodules.o fftcode.o -o $@

fft:: mtmodules.o fftcode.o
	gfortran  mtmodules.o fftcode.o -o $@

wigner:: mtmodules.o wignercode.o
	gfortran  mtmodules.o wignercode.o -o $@

dft:: mtmodules.o dftcode.o
	gfortran  mtmodules.o dftcode.o -o $@

dir:: 
	mkdir bin
	mv fft bin
	mv wigner bin
	mv dft bin

clean::
	rm -f *.o
	rm -f *.mod
	rm -rf bin

all:: fft wigner dft dir
