program fftcode
    use mtmodules
!!$ This program is used to test the mathtools subroutines
  implicit none
  real(8), allocatable :: x(:),y(:), xy(:), y2(:)
  real(8), allocatable :: xnew(:), ynew(:)
  real(8), allocatable :: refft(:), imfft(:)
  real(8), allocatable :: xaux(:), yaux(:)
  real(8) :: xp, yp, delta
  real(8) :: tinicial !! Starting point for the FFT
  character(len=100) :: nombre
  character(len=2) :: columnach
  character(len=1) :: auxch
  integer :: nexp, stat, columna, ltotales, lcabecera
  integer :: outexp
  logical :: ex
  integer :: ij, jk, ik, nn

  !!FIRST CHECK THAT THE OUTPUT EXISTS
  
  inquire(file='salida.dat',exist=ex)
  open(newunit=outexp, file='salida.dat')
  
  if (ex) then
     write(*,*)
     write(*,*) '"salida.dat" already exists. Breaking the program...'
     write(*,*)
     stop
  else
     write(outexp,*) '# FFT of a signal given as a function of time in nanoseconds'
  end if
  
  call getarg(1,nombre) !! file with the
       write(outexp,*) '# File: ', trim(nombre)
!!$  call getarg(2,columnach) !! column where the quantity is
!!$  read(columnach,*,iostat=stat) columna !! set the ordenates column as 

  write(*,*) 
  write(*,*) 'Which column shall I use in the FFT?'
  read(*,*) columna
  write(*,*)
  write(outexp,*) '# Column transformed: ', columna
  write(*,*) 'Starting time for the FFT (ns)?'
  read(*,*) tinicial
  write(outexp,*) '# t0: ', tinicial
  write(*,*)

  open(newunit=nexp,file=trim(nombre)) !! open the file

  !! READ THE FILE

  !! initialize the counters for the heading
  lcabecera=-33
  ij=0
  
  Do     
     read(nexp,*,iostat=stat) auxch
     if (auxch.ne.'#'.and.lcabecera.lt.0) lcabecera=ij !! heading lines
     if (stat.ne.0) then !! total lines
        stat=0
        ltotales=ij
        exit
     end if
     ij=ij+1 !! set the lines
  End Do

  write(outexp,*) '# Lines in heading removed: ', lcabecera
  
  rewind(nexp) !! BEGINNING OF THE FILE
  
  !!SKIP THE HEADING

  Do ij=1,lcabecera
     read(nexp,*) auxch
  End Do

  
  !!READ THE SIGNAL FROM THE FILE
  allocate(x(1:ltotales-lcabecera))
  allocate(y(1:ltotales-lcabecera))
  allocate(xy(1:columna))
  x=0.0d0
  y=0.0d0
  xy=0.0d0
  ik=0

  Do ij=1,ltotales-lcabecera
     read(nexp,*) xy(1:columna)
     if (xy(1).lt.tinicial) cycle
     ik=ik+1 !! number of valid points
     x(ik)=xy(1)
     y(ik)=xy(columna)
  End Do
  close(nexp)
  deallocate(xy)

  write(outexp,*) '# Number of data points: ', ik
  !! Resize the data points
  allocate(xaux, source=x)
  allocate(yaux, source=y)

  deallocate(x)
  deallocate(y)
  allocate(x(1:ik), y(1:ik))
  x=0.0d0
  y=0.0d0

  do ij=1,ik
     x(ij)=xaux(ij)
     y(ij)=yaux(ij)
  end do
  
  deallocate(xaux)
  deallocate(yaux)
  !! END of resizing the data
  !! Now the signal is in x and y
  !! Now we construct the splines

  allocate(y2(1:size(y)))
  y2=0.d0

  nn=size(x)

  !! SET THE CUBIC SPLINE TO INTERPOLATE
  call spline(x,y, nn, 0.0d0, 0.0d0, y2)

  !! SET A NEW ARRAY WITH THE DIMENSION 2^N, LARGER THAN THE CURRENT ONE
  !! SET THE X ARRAY
  call lpower2(x,xnew)

  xnew(1)=x(1)
  xnew(size(xnew))=x(size(x))

  Do ij=2,size(xnew)-1
     xnew(ij)=xnew(ij-1)+(xnew(size(xnew))-xnew(1))/(size(xnew)-1)
  End Do
  
  !! SET THE Y ARRAY
  call lpower2(y,ynew)
  write(outexp,*) '# Number of points for FFT: ', size(ynew)
  write(outexp,*) '# '
  Do ij=1,size(xnew)
     call splint(x,y,y2,size(x),xnew(ij),yp)
     ynew(ij)=yp
  End Do

  !! NOW THE FOURIER TRANSFORM
  !! FIRST THE DISCRETE FOURIER TRANSFORM
  call realft(ynew,size(ynew),1)
  !! THEN THE FOURIER TRANSFORM
    delta=(xnew(2)-xnew(1))*1.d-9 !! in seconds
    ynew=ynew*delta
    Do ij=0,size(ynew)-1
       xnew(ij+1)=ij/(size(xnew)*delta)
    End Do

    write(outexp,*) '# f (Hz)     FFT'
    write(outexp,*) '# '

    call realft2realimag(ynew,refft,imfft)

    Do ij=1, size(refft)
       write(outexp,*) xnew(ij), refft(ij), imfft(ij), refft(ij)**2.0d0+imfft(ij)**2.0d0
    End Do


    
end program fftcode
