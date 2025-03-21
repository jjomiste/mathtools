.SUFFIXES : .o .f90
.f90.o: ;  gfortran -c $<

fft:: mtmodules.o fftcode.o
	gfortran  mtmodules.o fftcode.o -o $@

wigner:: mtmodules.o wignercode.o
	gfortran  mtmodules.o wignercode.o -o $@

clean::
	rm -f *.o
	rm -f *.mod
	rm -f fft
	rm -f wigner

all:: fft wigner
