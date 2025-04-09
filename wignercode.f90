program wigner_transform
  use mtmodules
!!$ This program is used to test the mathtools subroutines
  implicit none
  real(8), allocatable :: x(:),y(:), xy(:), y2(:)
  real(8), allocatable :: xnew(:), ynew(:)
  real(8), allocatable :: refft(:), imfft(:)
  real(8), allocatable :: xaux(:), yaux(:)
  real(8) :: xp, yp, delta,dt
  real(8) :: tinicial !! Starting point for the FFT
  character(len=100) :: nombre
  character(len=2) :: columnach
  character(len=1) :: auxch
  integer :: nexp, stat, columna, ltotales, lcabecera
  integer :: outexp
  logical :: ex
  integer :: ij, jk, ik, nn

!!$ To perform the wigner transform
  integer:: N 
  real(8) :: x_min, x_max
  real(8) :: p_min , p_max 
  real(8) :: dx, dp
  real(8), allocatable :: x_vals(:), p_vals(:), W(:, :)
  integer :: i, j

  inquire(file='salida_wigner.dat',exist=ex)
  open(newunit=outexp, file='salida_wigner.dat')
  
  if (ex) then
     write(*,*)
     write(*,*) '"salida_wigner.dat" already exists. Breaking the program...'
     write(*,*)
     stop
  else
     write(outexp,*) '# Wigner transform of a signal given as a function of time'
  end if

    call getarg(1,nombre) !! file with the function to transform
       write(outexp,*) '# File: ', trim(nombre)

  write(*,*) 
  write(*,*) 'Which column shall I use to perform the Wigner transform?'
  read(*,*) columna
  write(*,*)
  write(outexp,*) '# Column to transform: ', columna

  open(newunit=nexp,file=trim(nombre)) !! open the file  

  write(*,*)


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
  
  rewind(nexp) !! BEGINNING OF THE FILE
  
  !!SKIP THE HEADING

  Do ij=1,lcabecera
     read(nexp,*) auxch
  End Do

  !!READ THE SIGNAL FROM THE FILE
  allocate(x(1:ltotales-lcabecera))  !! time points
  allocate(y(1:ltotales-lcabecera))  !! time-dependent function
  allocate(xy(1:columna)) !! single time entry to be stored
  x=0.0d0
  y=0.0d0
  xy=0.0d0
  ik=0 !! counter

  Do ij=1,ltotales-lcabecera
     read(nexp,*) xy(1:columna)
!!$     if (xy(1).lt.tinicial) cycle
     ik=ik+1 !! number of valid points
     x(ij)=xy(1)
     y(ij)=xy(columna)
  End Do
  close(nexp)
  deallocate(xy)

  write(outexp,*) '# Initial time: ', x(1)
  write(outexp,*) '# Final time: ', x(size(x))
  dt=x(2)-x(1)
  write(outexp,*) '# Time step: ', dt
  write(outexp,*) '# Number of time points: ', size(x)
  
  write(outexp,*) '# Lines in heading removed: ', lcabecera
  write(outexp,*) '# Number of data points: ', size(x)

  print*, x(1), x(size(x)),dt
  stop

  !!Arrange the limits for x and p
  
  x_min=x(1)
  x_max=x(size(x))
  dx=dt
  N=size(x)

  p_min=2*acos(-1.0)/(x_max-x_min)
  p_max=2*acos(-1.d0)/dx
  dp=(p_max-p_min)/dble(N)

  allocate(x_vals(1:N))
  allocate(p_vals(1:N))
  allocate(W(1:N,1:N))
  
  ! Initialize x and p values
  do i = 1, N
    x_vals(i) = x_min + (i - 1) * dx
    p_vals(i) = p_min + (i - 1) * dp
  end do

  ! Compute Wigner transform
  do i = 1, N
    do j = 1, N
      W(i, j) = compute_wigner(x_vals(i), p_vals(j))
    end do
  end do

  ! Output results to file
  open(unit=10, file="wigner_output.dat")
  do i = 1, N
    do j = 1, N
      write(10, *) x_vals(i), p_vals(j), W(i, j)
    end do
  end do
  close(10)

  print *, "Wigner transform computation completed. Data saved to wigner_output.dat"

contains

  function compute_wigner(x, p) result(W) !!This is only done for the gaussian case. DO IT IN GENERAL!!
    implicit none
    real(8), intent(in) :: x, p
    real(8) :: W, y, integral
    integer :: k, Ny
    real(8), parameter :: dy = 0.1d0, y_min = -5.0d0, y_max = 5.0d0
    Ny = int((y_max - y_min) / dy)
    integral = 0.0d0
    
    do k = 0, Ny
      y = y_min + k * dy
!!$      integral = integral + psi_gaussian(x + y) * psi_gaussian(x - y) * cos(-2.0d0 * p * y) * dy
    end do
    
    W = integral / acos(-1.0d0)
  end function compute_wigner

end program wigner_transform
