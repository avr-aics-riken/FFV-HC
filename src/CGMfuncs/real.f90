function sizeof_real()
  implicit none
  integer                  :: sizeof_real
  real, parameter          :: x = 0.0
  sizeof_real = SIZEOF(x)
end function

