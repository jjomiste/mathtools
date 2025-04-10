module mtmodules
  implicit none
  real(8), parameter :: two=2.0d0, zero=0.0d0, pi=acos(-1.0d0), half=dble(0.5d0)

  public::   two, zero,pi


  
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

!!! DISCRETE FOURIER TRANSFORM

subroutine compute_dft(t, f, omega, F_re, F_im)
  implicit none
  ! Input
  real(8), allocatable, intent(in) :: t(:), f(:)

  ! Output
  real(8),allocatable, intent(out) :: F_re(:), F_im(:),omega(:)

  ! Internal variables
  integer :: k, j,n
  real(8) :: dt, arg
  

  ! Assume uniform spacing
  dt = t(2) - t(1)

  n=size(t)

  allocate(F_re(1:n))
  allocate(F_im(1:n))
  allocate(omega(1:n))
  F_re=zero
  F_im=zero
  omega=zero
  
  ! Loop over angular frequencies
  do k = 1, n
    omega(k) = 2*dble(k - 1) / (n * dt)

    ! Sum over time samples
    do j = 1, n
      arg = omega(k) * t(j)
      F_re(k) = F_re(k) + f(j) * cos(arg)
      F_im(k) = F_im(k) - f(j) * sin(arg)
    end do
  end do

  !! Include the differential and the 1/sqrt(2pi)
  
  F_re=F_re*dt/sqrt(two*pi)
  F_im=F_im*dt/sqrt(two*pi)
  
end subroutine compute_dft
      
      
!!! WIGNER TRANSFORM SUBROUTINES

      subroutine wignert(x,ff,psteps,outfile,p) !! compute the Wigner transformation of a real function ff
        !!INPUT
        !!ff is a function of x
        real(8), allocatable, intent(in) :: x(:), ff(:)
        integer :: psteps !! number of frequencies
        integer :: outfile !! output file to write on
        !!OUTPUT
        !! p is the grid in p
        !! w is the wigner transform in the grid x,p
        real(8), allocatable, intent(inout) ::  p(:)
        !! AUXILIAR
        real(8) :: dx, dp,y,integral
        integer:: nn !! number of grid points
        integer:: ij, jk,ik

        
        !! setting the grid in x
        nn=size(x)
        dx=x(2)-x(1)

        !! setting the grid in p
        allocate(p(1:psteps))
        p=zero
        p(1)=zero
        p(psteps)=two*pi/dx
        dp=(p(psteps)-p(1))/(psteps-1)

        do ij = 2, psteps
           p(ij) = p(1) + (ij - 1) * dp
        end do

        !! computing the integral

        do ik=1,psteps !! Running in p
           do ij=1, nn !! Running in x
              integral=zero
              do jk = 0,min(ij+1,nn-ij) !! Integration over y.
                 !! at each step y=jk*dx
                 !! ff(ij+jk) runs over x(ij)+jk*dx
                 !! ff(ij-jk) runs over x(ij)-jk*dx
                 y=jk*dx
                 integral = integral + ff(ij+jk) * ff(ij-jk) * cos(-two * p(ik) * y) * dx
              end do
              integral=integral/pi
              write(outfile,*) p(ik)*half/pi,x(ij),integral
              
           end do
           write(*,*) dble(ik)/dble(psteps)*100, '% done...'
           write(outfile,*) ''
           !! now we print out the result

        end do

        write(*,*) 'Wigner transform computed!'
        
        return
      end subroutine wignert


!! DISCRETE FOURIER TRANSFORM


      !! READING TWO LISTS OF NUMBERS

      subroutine reading2list(inpfile,colx,coly,x,y)
        !!INPUT
        !! inpfile: name of the input file
        !! colx: column for the x values
        !! coly: column for the y values
        !! OUTPUT:
        !! x: array with the x values starting in label 1
        !! y: array with the y values starting in label 1
        implicit none
        !!INPUT
        character(len=*), intent(in) :: inpfile
        integer, intent(in) :: colx, coly
        !!OUTPUT
        real(8), allocatable, intent(inout) :: x(:), y(:)
        !!AUXILIAR
        logical :: ex
        integer :: lcabecera, ij,stat,ltotales,inpnumber
        character(1) :: auxch
        integer :: ik
        real(8), allocatable :: xy(:)
                  
        !!Checking the file
        
        if (len(trim(inpfile)).gt.100) then
           write(*,*) 'The name of the input file is too long'
           write(*,*) 'Breaking the code...'
           write(*,*)
           stop
        end if

        inquire(file=trim(inpfile),exist=ex)

        if (.not.ex) then
           write(*,*)
           write(*,*) trim(inpfile), ' does not exists. Breaking the program...'
           write(*,*)
           stop
        else
           open(newunit=inpnumber,file=trim(inpfile))
        end if

        !!Starting the reading
        
        write(*,*) 'Reading file: ', trim(inpfile)
        write(*,*) ''
        write(*,*) 'Reading the columns ', colx, ' and ', coly
        write(*,*)

  !! READ THE FILE

        
        !! initialize the counters for the heading
        lcabecera=-33
        ij=0


        Do     
           read(inpnumber,*,iostat=stat) auxch
           if (auxch.ne.'#'.and.lcabecera.lt.0) lcabecera=ij !! heading lines     
           if (stat.ne.0) then !! total lines
              stat=0
              ltotales=ij
              exit
           end if
           ij=ij+1 !! set the lines
        End Do
  
        rewind(inpnumber) !! BEGINNING OF THE FILE

          Do ij=1,lcabecera
             read(inpnumber,*) auxch
          End Do

  !!READ THE SIGNAL FROM THE FILE
          allocate(x(1:ltotales-lcabecera))  !! time points
          allocate(y(1:ltotales-lcabecera))  !! time-dependent function
          allocate(xy(1:max(colx,coly))) !! single time entry to be stored
          x=0.0d0
          y=0.0d0
          xy=0.0d0
          ik=0 !! counter

          
          Do ij=1,ltotales-lcabecera
             read(inpnumber,*) xy(1:max(colx,coly))
             ik=ik+1 !! number of valid points
             x(ij)=xy(colx)
             y(ij)=xy(coly)
          End Do

        
        return
      end subroutine reading2list
      
end module mtmodules
