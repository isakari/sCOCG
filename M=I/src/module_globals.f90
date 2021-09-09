module module_globals
  implicit none

  ! parameters
  complex(8),parameter :: zero=cmplx(0.d0,0.d0,kind(1.d0)) !< 0
  complex(8),parameter :: ione=cmplx(0.d0,1.d0,kind(1.d0)) !< I
  complex(8),parameter :: one=cmplx(1.d0,0.d0,kind(1.d0)) !< 1
  real(8),parameter :: pi=acos(-1.d0) !< Pi
  real(8),parameter :: gma=0.57721566490153286d0 !< Euler's gamma
  real(8),parameter :: tiny=epsilon(1.d0) !< a tiny number

  !for module_pade
  real(8), parameter :: tol_order=1d-13
  real(8), parameter :: tol_GMRES=0
  integer, parameter :: DKA_max_itr = 100
  
  ! for GL quadrature
  integer,parameter :: maxs=10 !< GT公式の最大S
  integer,parameter :: maxng=100 !<  maximum number of integral points for GL quadrature
  
end module module_globals
