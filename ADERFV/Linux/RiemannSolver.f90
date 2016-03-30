SUBROUTINE ADERRiemannSolver(lQbndL,lFbndL,lQbndR,lFbndR,nv)
  USE typesDef
  IMPLICIT NONE 
  ! Argument list 
  REAL, INTENT(IN)     :: lQbndL(nVar)              ! left space-time degrees of freedom 
  REAL, INTENT(IN)     :: lQbndR(nVar)              ! right space-time degrees of freedom 
  REAL, INTENT(IN)     :: nv(d)                                     ! normal vector 
  REAL, INTENT(INOUT)  :: lFbndL(nVar)              ! left flux 
  REAL, INTENT(INOUT)  :: lFbndR(nVar)              ! right flux 
  ! Local variables 
  INTEGER           :: i,j,k,l
  REAL              :: aux(d), QavL(nVar), QavR(nVar), smax
  REAL              :: LL(nVar), LR(nVar) 
  !
  ! Compute the average states from the left and the right, which we need to compute the numerical dissipation 
  QavL = 0. 
  QavR = 0. 
  QavL = QavL + lQbndL(:) 
  QavR = QavR + lQbndR(:) 
  !
  ! Here, we implement a very simple Rusanov scheme with scalar dissipation (smax*Id). 
  !
  CALL PDEEigenvalues(LL,QavL,nv) 
  CALL PDEEigenvalues(LR,QavR,nv) 
  smax = MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) ) 
  !
  ! We now compute the numerical flux. Note that the scheme is at the moment written in 
  ! CONSERVATION FORM => no fluctuations, but real fluxes. 
  ! Later, this will be converted into the left and right fluctuations. 
  !
  lFbndL(:) = 0.5*( lFbndR(:) + lFbndL(:) ) - 0.5*smax*( lQbndR(:) - lQbndL(:) ) 
  lFbndR(:) = lFbndL(:) 
  !
END SUBROUTINE ADERRiemannSolver
    
    
