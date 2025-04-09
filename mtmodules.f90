module mtmodules
  implicit none
  real(8), parameter :: two=2.0d0, zero=0.0d0, pi=acos(-1.0d0)

  private :: two, zero,pi
  public
  
  contains

    subroutine lpower2(vec,vecout)
      !! This subroutines expands the length of the vector up to the closest integer which is a power of 2.
      real(8), allocatable, intent(in) :: vec(:)
      real(8), allocatable, intent(out) :: vecout(:)
      integer :: ij, jk, newl,newlexp
      
      newlexp=int(log(real(size(vec)))/log(two))+1 !! exponent for the next 2^n

      call intpower(2,newlexp, newl)

      allocate(vecout(newl))

      vecout=zero
      
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

!!$Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
!!$x1 < x2 < ::: < xN, and given values yp1 and ypn for the rst derivative of the interpolating
!!$function at points 1 and n, respectively, this routine returns an array y2(1:n) of
!!$length n which contains the second derivatives of the interpolating function at the tabulated
!!$points xi. If yp1 and/or ypn are equal to 1  1030 or larger, the routine is signaled to set
!!$the corresponding boundary condition for a natural spline, with zero second derivative on
!!$that boundary.
!!$Parameter: NMAX is the largest anticipated value of n.
      
      INTEGER n,NMAX
      REAL(8) :: yp1,ypn
      REAL(8), allocatable, intent(in) :: x(:),y(:)
      REAL(8), allocatable, intent(inout) :: y2(:)
      PARAMETER (NMAX=15000)
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
11      continue

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

    !!FOURIER TRANSFORM

      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      REAL(8) data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
        1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software ]k1">"@w.
    
    !! FOURIER TRANSFORM OF A REAL FUNCTION
    
          SUBROUTINE realft(data,n,isign)
      INTEGER isign,n
      REAL(8) data(n)
!!$    USES four1
      INTEGER i,i1,i2,i3,i4,n2p3
      REAL c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      if (isign.eq.1) then
        c2=-0.5
        call four1(data,n/2,+1)
      else
        c2=0.5
        theta=-theta
      endif
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
               call four1(data,n/2,-1)
      endif
      return
      END
!!  (C) Copr. 1986-92 Numerical Recipes Software ]k1">"@w.

      !! This subroutine gives the real and the imaginary part of the
      !! Fourier Transform coming from realft. FFT is normalized to unity.
      subroutine realft2realimag(realftinput,re,imag)
        implicit none
        real(8), allocatable, intent(in) :: realftinput(:)
        real(8), allocatable, intent(out) :: re(:), imag(:)
        integer :: dimft, ij, ik, jk
        real(8), allocatable :: fft2(:)
        real(8) :: maxfft2

        allocate(re(1:size(realftinput)/2))
        re=0.0d0
        allocate(imag,source=re)
        allocate(fft2, source=re)

        dimft=size(re)
        
        re(1)=realftinput(1)
        re(dimft)=realftinput(2)

        Do ij=2,dimft
           re(ij)=realftinput(2*ij-1)
           imag(ij)=realftinput(2*ij)
        End Do

        fft2=sqrt(re**2.0d0+imag**2.0d0)

        maxfft2=maxval(fft2)
        
        re=re/maxfft2
        imag=imag/maxfft2
        
        return
      end subroutine realft2realimag

!!! WIGNER TRANSFORM SUBROUTINES

      subroutine wigner(x,ff,p,w) !! compute the Wigner transformation of a real function ff
        !!INPUT
        !!ff is a function of x
        real(8), allocatable, intent(in) :: x(:), ff(:)
        !!OUTPUT
        !! p is the grid in p
        !! w is the wigner transform in the grid x,p
        real(8), allocatable, intent(inout) ::  p(:), w(:,:)
        !! AUXILIAR
        real(8) :: dx, dp,y,integral
        integer:: nn !! number of grid points
        integer:: ij, jk,ik

        allocate(w(1:nn,1:nn))
        w=zero
        
        !! setting the grid in x
        nn=size(x)
        dx=(x(nn)-x(1))/dble(nn)

        !! setting the grid in p
        allocate(p(1:nn))
        p=zero
        p(1)=2*pi/(x(nn)-x(1))
        p(nn)=2*pi/dx
        dp=(p(nn)-p(1))/dble(nn)

        do ij = 2, nn
           p(ij) = p(1) + (ij - 1) * dp
        end do

        !! computing the integral

        do ik=1,nn !! Running in p
           do ij=1, nn !! Running in x
              do jk = ij+1, nn-ij !! Integration in y
                 y = x(1) + (jk-1) * dx
                 integral = integral + ff(ij+jk) * ff(ij-jk) * cos(-two * p(ik) * y) * dx !!NOT DONE YET
              end do
              w(ij,ik)=integral
           end do
        end do
    
        
        return
      end subroutine wigner
      
end module mtmodules
