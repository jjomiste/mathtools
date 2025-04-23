program windowdft
  use mtmodules
  implicit none

  character(len=100):: infile
  integer:: ifile,xcol,ycol,outfile, wdftfile
  real(8), allocatable :: x(:), y(:), omega(:), fre(:), fim(:)
  real(8), allocatable :: ywindow(:)
  integer:: ij, jk, ik,tstep
  real(8) :: windl!! duration of the window
  real(8) :: peak !! peak of the Hann window
  real(8) :: xx !! variable
  character(len=50) :: outfilechar
  
  call getarg(1,infile) !! get the input file name

  !! set the columns to analyze
  write(*,*) 'Which is the column for the x?'
  read(*,*) xcol
  
  write(*,*) 'Which is the column for the y?'
  read(*,*) ycol

  write(*,*) 'What is the duration of the window?'
  read(*,*) windl

  write(*,*) 'Time steps in the to be stored'
  read(*,*) tstep

  !!avoiding errors in the duration of the window
  if (windl.le.1d-3) then
     write(*,*)
     write(*,*) 'THE DURATION OF THE PULSE MUST BE LARGER THAN ZERO'
     write(*,*) 'Breaking the code...'
     write(*,*)
     stop
  end if  
  
  write(*,*) ''

  !!READING THE FILE
  call reading2list(infile,xcol,ycol,x,y)

  !!COMPUTING THE Windowed-DFT

  !! Open the new file
  write(outfilechar,'(A17,E9.3)') 'salida_wdft_file_',windl

  open(newunit=wdftfile, file=trim(outfilechar))

  write(wdftfile,*) '#Input file: ', trim(infile)
  write(wdftfile,*) ''
  write(wdftfile,*) '#Duration of the Hann Window: ', windl
  write(wdftfile,*) ''

  write(wdftfile,*) '# Note that the frequency units are the inverse of the input'
  write(wdftfile,'(A2,3X, A1,17X,A1,16X,A13,3X,A7,10X,A7,10X,A9)') ' #',  't','f', 'omega(2*pi*f)',  &
       'Re(FFT)',    'Im(FFT)',    '|WFFT|^2'
  write(wdftfile,*) '# '

  write(*,*) 'Computing the Windowed DFT'
  write(*,*)
  
  !! LOOP IN THE TIME POINTS
  Do ij=1,size(x),tstep

     peak=x(ij)
     !! Obtaining the time signal times the window
     call hannarray(x,y,peak,windl,ywindow)

     !! Compute the Fourier Transform of the tiem signal times the window
     call compute_wdft(x, ywindow,peak-windl*half,peak+windl*half, omega, fre, fim)

     do jk=1,size(omega)

        write(wdftfile,'(6E17.5)') x(ij), omega(jk)/(two*pi),omega(jk), fre(jk),fim(jk),fre(jk)**2.+fim(jk)**2.
     end do

     write(wdftfile,*) ''
  End Do
 
  
end program windowdft
