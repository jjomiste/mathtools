.SUFFIXES : .o .f90
.f90.o: ; gfortran -c $<

fft: mtmodules.o fftcode.o
	gfortran  mtmodules.o fftcode.o -o $@

clean::
	rm -f *.o
	rm -f *.mod
	rm -f fft

all::	fft
