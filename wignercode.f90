program wigner_transform
  implicit none
  integer, parameter :: N = 100
  real(8), parameter :: x_min = -5.0d0, x_max = 5.0d0
  real(8), parameter :: p_min = -5.0d0, p_max = 5.0d0
  real(8), parameter :: dx = (x_max - x_min) / N, dp = (p_max - p_min) / N
  real(8) :: x_vals(N), p_vals(N), W(N, N)
  integer :: i, j

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

  function psi_gaussian(x) result(psi)
    implicit none
    real(8), intent(in) :: x
    real(8), parameter :: sigma = 1.0d0
    real(8) :: psi
    psi = (1.0d0 / (sigma * sqrt(3.141592653589793d0))) * exp(-x**2 / (2.0d0 * sigma**2))
  end function psi_gaussian

  function compute_wigner(x, p) result(W)
    implicit none
    real(8), intent(in) :: x, p
    real(8) :: W, y, integral
    integer :: k, Ny
    real(8), parameter :: dy = 0.1d0, y_min = -5.0d0, y_max = 5.0d0
    Ny = int((y_max - y_min) / dy)
    integral = 0.0d0
    
    do k = 0, Ny
      y = y_min + k * dy
      integral = integral + psi_gaussian(x + y) * psi_gaussian(x - y) * cos(-2.0d0 * p * y) * dy
    end do
    
    W = integral / 3.141592653589793d0
  end function compute_wigner

end program wigner_transform
