!-----------------------------------------------------------------------
!
!  Globals module.  Put all variables to be accessible everywhere
!    here. 
!
!-----------------------------------------------------------------------

module globals
! WP normally is 14
  integer, parameter     :: wp = selected_real_kind(p=14)   ! this sets the floating-point precision to double
  integer, parameter     :: int16 = selected_int_kind(16)   ! larger integers!
  complex(wp), parameter :: i = (0.0_wp, 1.0_wp)            ! sqrt (-1)
  real(wp), parameter    :: pi = 3.14159265358979_wp        ! pi! REAL: VERY WEIRD BUGS HAPPEN IF YOU SET PI=3

  ! for random number generator
  integer, parameter :: rand_pl_wp = wp
  integer, parameter :: rand_pl_mf = 100

end module globals
