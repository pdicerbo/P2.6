SUBROUTINE Reconstruction
  USE typesDef
  IMPLICIT NONE
  ! Local variables 
  INTEGER :: i, j, k, iw, iVar, jElem, iDim, iDOF, jDOF,ii, jj, kk, iii, jjj, kkk, neighp, neighm, counter, PositivityStatus 
  INTEGER :: sgn, ix(N+1)
  REAL    :: QR(nVar), QL(nVar), VR(nVar), VL(nVar), VE(nVar), umin, umax, vmin, vmax       
  REAL    :: InputVector( nVar, -ss:ss ) 
  REAL    :: wx(nVar,N+1), u0(nVar), aux(d)  
  REAL    :: luh(nVar,-ss:ss), lwh(nVar,N+1) 
  REAL    :: lx0(d), ldx(d), lxb(d), jx0(d), jdx(d), jxb(d)
  REAL    :: xGP(d),whTmp(nVar)
  REAL, PARAMETER    :: epsv = 1e-7, epsilon = 1e-1                
  !
  ! Pure WENO finite volume scheme
  ! ============================== 
  !
  DO i = 1, nElem

     DO ii = -ss, ss
        jElem = Neighbor(ii,i)  
        luh(:,ii) = uh(:,jElem) 
     ENDDO
     ! 
     lwh = 0. 
     ! Do WENO reconstruction in x-direction (1D)        
     DO j = -ss, ss 
        InputVector(:,j) = luh(:,j) 
     ENDDO

     CALL WENOReconstruction1D(wx,InputVector)
     lwh(:,1:N+1) = wx  

     wh(:,1:N+1,i) = lwh(:,1:N+1) 
     !
  ENDDO
  ! 

  ! 
END SUBROUTINE Reconstruction

SUBROUTINE WENOReconstruction1D(wh_WENO,InputVector) 
  USE typesDef
  IMPLICIT NONE
  ! Argument list  
  REAL    :: wh_WENO(nVar,N+1) 
  REAL    :: InputVector(nVar,-ss:ss) 
  ! Local variables 
  INTEGER :: i,j,k,l,iStencil, iDOF 
  REAL    :: wh_Stencil(nVar,N+1,WENO_nStencil), OI(nVar,WENO_nStencil), wh6(nVar,6)  
  REAL    :: omega(nVar,WENO_nStencil), omegasum(nVar) 
  REAL    :: b(N+1,nVar),centralweight     
  ! 
  wh6 = 0. 
  omegasum = 0.
  wh_Stencil = 0.  
  OI = 0. 
  DO iStencil = 1, WENO_nStencil
     !b = TRANSPOSE(InputVector(:,WENO_left(iStencil):WENO_right(iStencil)) )
     !wh_Stencil(:,:,iStencil) = TRANSPOSE(MATMUL(WENO_iMatrix(:,:,iStencil),b))   
     DO j = 1, N+1
        DO i = 1, N+1
           wh_Stencil(:,i,iStencil) = wh_Stencil(:,i,iStencil) + WENO_iMatrix(i,j,iStencil)*InputVector(:,WENO_left(iStencil)-1+j) 
        ENDDO
     ENDDO
     DO k = 1, N+1
        DO l = 1, N+1
           OI(:,iStencil) = OI(:,iStencil) + WENO_OIMatrix(k,l)*wh_Stencil(:,l,iStencil)*wh_Stencil(:,k,iStencil) 
        ENDDO
     ENDDO
     omega(:,iStencil) = WENO_lambda(iStencil)/(OI(:,iStencil)+WENO_epsilon)**WENO_r
     omegasum = omegasum + omega(:,iStencil) 
  ENDDO
  wh_WENO = 0. 
  DO iStencil = 1, WENO_nStencil
     omega(:,iStencil) = omega(:,iStencil)/omegasum
     DO iDOF = 1, N+1
        wh_WENO(:,iDOF) = wh_WENO(:,iDOF) + omega(:,iStencil)*wh_Stencil(:,iDOF,iStencil)
     ENDDO
  ENDDO
  ! 
END SUBROUTINE WENOReconstruction1D

