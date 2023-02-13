program test_subs
    use mtmodules
!!$ This program is used to test the mathtools subroutines
  implicit none
  real(8), allocatable :: x(:),y(:), xy(:), y2(:)
  real(8), allocatable :: xnew(:), ynew(:)
  real(8) :: xp, yp
  character(len=100) :: nombre
  character(len=2) :: columnach
  character(len=1) :: auxch
  integer :: nexp, stat, columna, ltotales, lcabecera
  integer :: outexp
  integer :: ij, jk, ik, nn

  call getarg(1,nombre) !! file with the
  call getarg(2,columnach) !! column where the quantity is

  read(columnach,*,iostat=stat) columna !! set the ordenates column as 
  
  open(newunit=nexp,file=trim(nombre)) !! open the file

  !! READ THE FILE

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
  allocate(x(lcabecera+1:ltotales))
  allocate(y(lcabecera+1:ltotales))
  allocate(xy(1:columna))
  x=0.0d0
  y=0.0d0
  xy=0.0d0

  Do ij=lcabecera+1,ltotales
     read(nexp,*) xy(1:columna)
     x(ij)=xy(1)
     y(ij)=xy(columna)
  End Do

  deallocate(xy)
  
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

  open(newunit=outexp, file='salida.dat')
  Do ij=1,size(xnew)
     call splint(x,y,y2,size(x),xnew(ij),yp)
     ynew(ij)=yp
     write(outexp,*) xnew(ij), ynew(ij)
  End Do

  
  
end program test_subs
