program windowdft
  use mtmodules
  implicit none

  character(len=100):: infile
  integer:: ifile,xcol,ycol,outfile, dftfile
  real(8), allocatable :: x(:), y(:), omega(:), fre(:), fim(:)
  integer:: ij, jk, ik
  real(8) :: windl
  
  call getarg(1,infile) !! get the input file name

  !! set the columns to analyze
  write(*,*) 'Which is the column for the x?'
  read(*,*) xcol
  
  write(*,*) 'Which is the column for the y?'
  read(*,*) ycol

  write(*,*) 'What is the duration of the window?'
  read(*,*) windl

  !!avoiding errors in the duration of the window
  if (windl.le.1d-3) then
     write(*,*)
     write(*,*) 'THE DURATION OF THE PULSE MUST BE LARGER THAN ZERO'
     write(*,*) 'Breaking the code...'
     write(*,*)
     stop
  end if  
  
  write(*,*) ''
  
  call reading2list(infile,xcol,ycol,x,y)

  call compute_dft(x, y, omega, fre, fim)

  open(newunit=dftfile, file='salida_dft_file')

  do ij=1,size(omega)
     write(dftfile,*) omega(ij)/(two*pi),omega(ij), fre(ij),fim(ij),fre(ij)**2.+fim(ij)**2.
  end do

  
  
end program windowdft
