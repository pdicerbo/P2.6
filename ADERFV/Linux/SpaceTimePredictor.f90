SUBROUTINE ADERSpaceTimePredictorNonlinear(lqhi,lFhi,lQbnd,lFbnd,lwh)
  USE typesDef
  IMPLICIT NONE
  ! Argument list
  REAL, INTENT(IN)  :: lwh(nVar,nDOF(1))              ! spatial degrees of freedom
  REAL, INTENT(OUT) :: lqhi(nVar,nDOF(1))             ! time-averaged space-time degrees of freedom
  REAL, INTENT(OUT) :: lFhi(nVar,d,nDOF(1))           ! time-averaged nonlinear flux tensor in each space-time DOF
  REAL, INTENT(OUT) :: lqbnd(nVar,2)                  ! time-averaged space-time degrees of freedom
  REAL, INTENT(OUT) :: lFbnd(nVar,2)                  ! time-averaged nonlinear flux tensor in each space-time DOF
  ! Local variables
  INTEGER :: i,l,iVar,iDim, iter
  REAL    :: rhs0(nVar,nDOF(1),nDOF(0))               ! contribution of the initial condition to the known right hand side
  REAL    :: rhs(nVar,nDOF(1),nDOF(0))                ! known right hand side
  REAL    :: lqh(nVar,nDOF(1),nDOF(0))                ! space-time degrees of freedom
  REAL    :: lFh(nVar,d,nDOF(1),nDOF(0))              ! nonlinear flux tensor in each space-time DOF
  REAL    :: aux(d), w                                                ! auxiliary variables
  REAL    :: lqhold(nVar,nDOF(1),nDOF(0))             ! old space-time degrees of freedom
  REAL    :: lqt(nVar,nDOF(1),nDOF(0))                ! time derivative qt of q
  REAL    :: res                                                      ! residual
  REAL, PARAMETER :: tol = 1e-7                                      ! tolerance
  !
  DO i = 1, N+1
     ! Trivial initial guess (can be significantly improved)
     DO iVar = 1, nVar
        lqh(iVar,i,:) = lwh(iVar,i)
     ENDDO
     ! Compute the contribution of the initial condition uh to the time update. 
     aux = wGPN(i)
     DO iVar = 1, nVar
        rhs0(iVar,i,:) = PRODUCT(aux(1:nDim))*F0(:)*lwh(iVar,i)
     ENDDO
     !
  ENDDO
  !
  ! Discrete Picard iterations. 
  DO iter = 1, N+1
     ! save old space-time DOF
     lqhold = lqh
     DO l = 1, nDOF(0) ! loop over DOF in time
        ! Compute the fluxes (once these fluxes are available, the subsequent operations are independent from each other)
        DO i = 1, nDOF(1)
           CALL PDEFlux(lFh(:,:,i,l),lqh(:,i,l))
        ENDDO
        ! Compute the "derivatives" (contributions of the stiffness matrix)
        ! x direction (independent from the y and z derivatives)
        rhs(:,:,l) = rhs0(:,:,l) - wGPN(l)*dt/dx(1)*MATMUL( lFh(:,1,:,l), Kxi )
        !
     ENDDO ! end loop over time DOF
     !
     ! Multiply with (K1)^(-1) to get the discrete time integral of the discrete Picard iteration
     !
     DO i = 1, nDOF(1)
        aux = wGPN(i)
        lqh(:,i,:) = 1./(PRODUCT(aux(1:nDim)))*MATMUL( rhs(:,i,:), TRANSPOSE(iK1) )
        !lqt(:,i,j,k,:) = 1.0/dt*MATMUL( lqh(:,i,j,k,:), TRANSPOSE(dudx) )         ! currently used only for debugging purposes, to check if derivatives are correctly computed
     ENDDO
     !
     ! We can stop the iterations if a certain tolerance has been reached. If you do not like this unpredictable feature (it depends on the solution of the PDE)
     ! simply comment the lines below, so each element will always do the same number of iterations in the predictor step, i.e. the same number of operations
     !
     res = SQRT(SUM((lqh-lqhold)**2))
     IF(res.LT.tol) THEN
        EXIT
     ENDIF
     !
  ENDDO
  !
  ! Immediately compute the time-averaged space-time polynomials
  !
  DO i = 1, nDOF(1)
     lqhi(:,i) = MATMUL( lqh(:,i,:), wGPN )
     DO iDim = 1, nDim
        lFhi(:,iDim,i) = MATMUL( lFh(:,iDim,i,:), wGPN )
     ENDDO
  ENDDO
  !
  ! Compute the bounday-extrapolated values for Q and F*n
  !
  lQbnd = 0.
  lFbnd = 0.
  ! x-direction: face 1 (left) and face 2 (right)
  lQbnd(:,1) = MATMUL( lqhi(:,:),   FLCoeff )   ! left
  lQbnd(:,2) = MATMUL( lqhi(:,:),   FRCoeff )   ! right
  lFbnd(:,1) = MATMUL( lFhi(:,1,:), FLCoeff )   ! left
  lFbnd(:,2) = MATMUL( lFhi(:,1,:), FRCoeff )   ! right
  !
  CONTINUE
  !
END SUBROUTINE ADERSpaceTimePredictorNonlinear








