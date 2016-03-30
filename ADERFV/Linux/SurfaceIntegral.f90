SUBROUTINE ADERSurfaceIntegral(lduh,lFbnd)
  USE typesDef
  IMPLICIT NONE 
  ! Argument list 
  REAL, INTENT(IN)    :: lFbnd(nVar,2)            ! nonlinear flux tensor in each space-time DOF 
  REAL, INTENT(INOUT) :: lduh(nVar,nDOF(1))       ! spatial degrees of freedom 
  ! Local variables 
  INTEGER           :: i,l,iVar 
  REAL              :: aux(d) 
  !
  ! Now multiply the numerical fluxes on the surfaces with the test functions and compute the surface integrals 
  ! 
  ! Initialize the update
  lduh = 0.
  DO iVar = 1, nVar 
     lduh(iVar,:) = lduh(iVar,:) - 1.0/dx(1)*( lFbnd(iVar,2) - lFbnd(iVar,1) )      ! left flux minus right flux 
  ENDDO
  !
END SUBROUTINE ADERSurfaceIntegral
    
    
