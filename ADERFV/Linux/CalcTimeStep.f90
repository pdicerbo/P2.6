SUBROUTINE CalcTimeStep
  USE typesDef
  IMPLICIT NONE
  ! Local variables
  INTEGER :: iElem, iDim, i  
  REAL    :: lmax(d), Lambda(nVar), dtE 
  REAL    :: nv(d,d), denom  
  ! Normal vectors pointing into the three space dimensions 
  nv = 0. 
  DO i = 1, d
     nv(i,i) = 1. 
  ENDDO
  !
  ! This function computes the maximum admissible time step according to the CFL condition for a generic nonlinear PDE 
  !
  dt = 1e20 
  DO iElem = 1, nElem
     denom = 0.   
     DO iDim = 1, nDim
        CALL PDEEigenvalues(Lambda,uh(:,iElem),nv(:,iDim)) 
        denom = denom + MAXVAL(ABS(Lambda))/dx(iDim) 
     ENDDO
     dt = MIN(dt, CFL*PNPMTable(N)/denom )    
  ENDDO
  ! 
  IF(time>tend) THEN
     dt = tend - time
  ENDIF
  ! 
END SUBROUTINE CalcTimeStep
    
