SUBROUTINE ElementUpdate(luh,lduh)
  USE typesDef
  IMPLICIT NONE 
  ! Argument list 
  REAL, INTENT(INOUT) :: lduh(nVar)      ! spatial degrees of freedom 
  REAL, INTENT(OUT)   :: luh(nVar)       ! nonlinear flux tensor in each space-time DOF 
  ! Local variables 
  INTEGER             :: i,l
  ! 
  ! Finally, sum the contribution to the spatial degrees of freedom 
  !
  luh = luh + dt*lduh 
  !
END SUBROUTINE ElementUpdate
    
    
