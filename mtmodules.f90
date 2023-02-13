module mtmodules
  implicit none
  real(8), parameter :: two=2.0d0, zero=0.0d0

  private :: two, zero
  public
  
  contains

    subroutine lpower2(vec)
      !! This subroutines expands the length of the vector up to the closest integer which is a power of 2.
      real(8), allocatable, intent(inout) :: vec(:)
      real(8), allocatable :: vecaux(:)
      integer :: ij, jk, newl,newlexp
      
      newlexp=int(log(real(size(vec)))/log(two))+1 !! exponent for the next 2^n

      call intpower(2,newlexp, newl)

      allocate(vecaux, source=vec)

      deallocate(vec)
      allocate(vec(newl))

      do ij=1, size(vecaux)
         vec(ij)=vecaux(ij)
      end do
      
      deallocate(vecaux)
      return
    end subroutine lpower2

    subroutine intpower(base, exponente, resultado)
      !! This subroutine computes base^exponente in term of integers
      integer, intent(in) :: base, exponente
      integer, intent(inout) :: resultado
      integer :: ij
      resultado =1

      Do ij=1,exponente
         resultado=resultado*base
      End Do
      
    end subroutine intpower

!!INTERPOLATION FUNCTIONS

    SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      !! TAKEN FROM NUMERICAL RECIPES
      !! NOTE IT IS IN SINGLE PRECISION

!!$Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
!!$x1 < x2 < ::: < xN, and given values yp1 and ypn for the rst derivative of the interpolating
!!$function at points 1 and n, respectively, this routine returns an array y2(1:n) of
!!$length n which contains the second derivatives of the interpolating function at the tabulated
!!$points xi. If yp1 and/or ypn are equal to 1  1030 or larger, the routine is signaled to set
!!$the corresponding boundary condition for a natural spline, with zero second derivative on
!!$that boundary.
!!$Parameter: NMAX is the largest anticipated value of n.
      
      INTEGER n,NMAX
      REAL(8) yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL(8) p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
       END !! END SUBROUTINE
!!$C  (C) Copr. 1986-92 Numerical Recipes Software ]k1">"@w.

       SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL(8) x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL(8) a,b,h

!!$Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
!!$xai 's in order), and given the array y2a(1:n), which is the output from spline above,
!!$and given a value of x, this routine returns a cubic-spline interpolated value y.
      
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
         write(*,*) 'bad xa input in splint'
         stop
      end if
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
    END SUBROUTINE splint
!!C  (C) Copr. 1986-92 Numerical Recipes Software ]k1">"@w.
    
end module mtmodules
