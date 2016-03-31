SUBROUTINE ADERFVInit
  USE typesDef
  IMPLICIT NONE
  ! Local variables 
  INTEGER :: i, j, k, ii, jj, kk, l, c, iGP, iElem, VMAX(d), count, cnt
  INTEGER :: iStencil, N2, counter
  REAL    :: A(N+1,N+1),iA(N+1,N+1),IntPhiPsi(1,N+1)
  REAL    :: phi(N+1), phi0(N+1), phi1(N+1), phi_xi(N+1)
  REAL    :: phi_i(N+1), phi_j(N+1), phi_k(N+1) 
  REAL    :: u0(nVar), xGP(d), x0(d), aux(d), nv(d), xGP1
  INTEGER, POINTER :: idxn(:), idxe(:)  
  INTEGER  :: tmpneighbor(-ss:ss)
  ! 
  ! ----------------- Some important preliminary stuff. Do not touch! ------------------- 
  dn(:) = 0 
  DO i = 1, nDim
     dn(i) = 1
  ENDDO
  ! According to the number of space dimensions, we set the number of degrees of freedom 
  ! i = 0 is the time dimension 
  nDOF(:) = 1 
  DO i = 0, nDim
     nDOF(i) = N+1
  ENDDO
  ! ------------------------------------------------------------------------------------- 
  !
  ! Some info about the PDE system 
  !
  EQN%gamma = 1.4                                                 ! ratio of specific heats for compressible Euler 
  !
  ! Here, you need to define the computational domain and the number of cells. 
  ! This is typically read from a parameter file 
  !
  xL = -0.5                                                       ! lower-left corner of the domain 
  xR =  0.5                                                       ! upper right corner of the domain 
  IMAX = 200 !36                                                       ! Number of elements in x,y,z direction 
  VMAX = IMAX                                                     ! Vector of the number of elements in each space dimension 
  dx = (xR-xL)/VMAX                                               ! Mesh spacing 
  NMAX = 100000                                                   ! Max. number of time steps 
  timestep = 0                                                    ! initial time step number 
  time = 0.                                                       ! initial time 
  tend = 0.1 ! 0.0 !25                                                     ! final time 
  Basefile = 'Test'                                               ! Base filename for writing results 
  !
  nElem = IMAX                                                    ! Number of elements 
  ALLOCATE(  uh(nVar, nElem) )               ! Allocate the cell averages 
  ALLOCATE( duh(nVar, nElem) )               ! Allocate the cell averages for the update 
  ALLOCATE( qhi(nVar, nDOF(1), nElem) )         ! Allocate the time-averaged space-time degrees of freedom 
  ALLOCATE( Fhi(nVar, d, nDOF(1), nElem) )      ! Allocate the time-averaged space-time degrees of freedom 
  ALLOCATE( qbnd(nVar, 2,  nElem) )              ! Allocate the time-averaged boundary-extrapolated values for Q 
  ALLOCATE( Fbnd(nVar, 2,  nElem) )              ! Allocate the time-averaged boundary-extrapolated values for the normal flux F * n 
  ALLOCATE(  wh(nVar,nDOF(1),nElem) )
  nNode = (IMAX+1)                                                ! number of nodes 
  ALLOCATE( x(d, nNode) )                                         ! Allocate the nodes 
  ALLOCATE( idxn(IMAX+dn(1))  )    
  ALLOCATE( BNDstatus(nElem))
  ! Define the node coordinates and the node numbers                                 
  count = 0 
  DO i = 1, IMAX+dn(1) 
     count = count + 1 
     x(:,count) = xL(:) + (i-1)*dx(:) 
     idxn(i) = count 
  ENDDO
  ! define the connectivity between the elements and the nodes. 
  ALLOCATE( tri(nVtx, nElem)     ) 
  ALLOCATE( neighbor(-ss:ss,nElem))
  ALLOCATE( idxe(0:IMAX+1) ) 
  BNDstatus       = 0
  count = 0  
  idxe = 0
  DO i = 1, IMAX 
     count = count + 1
     idxe(i) = count 
     tri(:,count) = (/ idxn(i), idxn(i+1) /) 
     !
     IF(i.LE.ss) THEN 
        BNDstatus(count) = 1
     ENDIF
     IF(i.GT.IMAX-ss) THEN
        BNDstatus(count) = 2
     ENDIF
  ENDDO
  idxe0 => idxe
  DEALLOCATE( idxn )
  ! define the connectivity between the faces and the elements 
  ! count how many faces we have 
  nFace = IMAX+1
  ALLOCATE( Face(nFace) ) 
  ! x faces qui
  count = 0 
  DO i = 1, IMAX+1  
     count = count + 1 
     IF(i.EQ.1) THEN
        Face(count)%Left  = 0 
        Face(count)%Right = idxe(i)
        ALLOCATE( Face(count)%qL(nVar) ) 
        ALLOCATE( Face(count)%FL(nVar) ) 
        Face(count)%qR => qBnd(:,1,Face(count)%Right) 
        Face(count)%FR => FBnd(:,1,Face(count)%Right) 
     ELSEIF(i.EQ.IMAX+1) THEN 
        Face(count)%Left  = idxe(i-1) 
        Face(count)%Right = 0 
        Face(count)%qL => qBnd(:,2,Face(count)%Left ) 
        Face(count)%FL => FBnd(:,2,Face(count)%Left ) 
        ALLOCATE( Face(count)%qR(nVar) ) 
        ALLOCATE( Face(count)%FR(nVar) ) 
     ELSE              
        Face(count)%Left  = idxe(i-1) 
        Face(count)%Right = idxe(i) 
        Face(count)%qL => qBnd(:,2,Face(count)%Left ) 
        Face(count)%qR => qBnd(:,1,Face(count)%Right) 
        Face(count)%FL => FBnd(:,2,Face(count)%Left ) 
        Face(count)%FR => FBnd(:,1,Face(count)%Right) 
     ENDIF
     Face(count)%nv = 1.0 ! set face normal vector 
  ENDDO
  !      
  !
  ! Initialize the neighbors
  DO i = 1, IMAX 
     neighbor(:,i) = i 
  ENDDO
  !
  DO iElem = 1, nElem
     CALL GetL0Neighbor(tmpneighbor,iElem)
     neighbor(:,iElem) = tmpneighbor(:)
  ENDDO
  DEALLOCATE( idxe ) 
  ! We now need to define our basis functions. This is done by choosing a set of distinct 1D nodes in the interval [0,1] 
  ! The basis functions will be the Lagrange interpolation polynomials running through these nodes 
  CALL gauleg(0.,1.,xiGPN,wGPN,N+1)
  xin = xiGPN  ! WE make the following choice: the basis functions run through the Gauss-Legendre nodes (=> really orthogonal basis) 
  !
  ! Now, let us compute some of the important matrices in the ADER-FV method 
  ! 
  MM   = 0.   ! Element mass matrix 
  Kxi  = 0.   ! Element stiffness matrix 
  dudx = 0.   ! discrete derivative operator, which projects the derivatives onto the basis  
  DO iGP = 1, N+1 
     CALL BaseFunc1D(phi,phi_xi,xiGPN(iGP))
     DO k = 1, N+1
        DO l = 1, N+1
           ! i) Mass-matrix 
           MM(k,l) = MM(k,l) + wGPN(iGP)*phi(k)*phi(l) 
           ! ii) Stiffness-matrix 
           Kxi(k,l) = Kxi(k,l) + wGPN(iGP)*phi_xi(k)*phi(l)  
        ENDDO
     ENDDO
  ENDDO
  CALL MatrixInverse(N+1,MM,iMM) 
  dudx = MATMUL( iMM, TRANSPOSE(Kxi) ) 
  CALL BaseFunc1D(phi0,phi_xi,0.0) ! Compute the basis functions on the left 
  CALL BaseFunc1D(phi1,phi_xi,1.0) ! Compute the basis function on the right 
  ! The flux matrices are all possible combinations of left and right 
  DO k = 1, N+1
     DO l = 1, N+1
        FLm(k,l) = phi0(k)*phi1(l)   ! Left contribution to the left flux matrix    (m = left  of the interface)  
        FLp(k,l) = phi0(k)*phi0(l)   ! Right contribution to the left flux matrix   (p = right of the interface) 
        FRm(k,l) = phi1(k)*phi1(l)   ! Left contribution to the right flux matrix   (m = left  of the interface) 
        FRp(k,l) = phi1(k)*phi0(l)   ! Right contribution to the right flux matrix  (p = right of the interface) 
     ENDDO
  ENDDO
  ! The time flux matrices for the ADER-FV predictor method are given by the principle of upwinding in time (causality principle) 
  F0 = phi0   ! upwinding in time = information comes from smaller times 
  F1 = FRm    ! upwinding in time = information comes from smaller times  
  K1 = F1 - Kxi 
  CALL MatrixInverse(N+1,K1,iK1)   
  FLcoeff = phi0  ! coefficients needed to extrapolate data onto the left  boundary 
  FRcoeff = phi1  ! coefficients needed to extrapolate data onto the right boundary 
  !
  ! Set the initial condition. Here, we assume a nodal basis. Otherwise, we would have to do formal L2 projection, 
  ! i.e. integration of the initial condition and multiplication with the inverse mass matrix. 
  !
  DO iElem = 1, nElem
     uh(:,iElem) = 0.0 
     x0 = x(:,tri(1,iElem)) ! get the coordinate of the lower left node 
     DO i = 1, N + 1
        xGP = x0 + xiGPN(i)*dx(:) 
        CALL InitialField(u0,xGP) 
        uh(:,iElem) = uh(:,iElem) + u0*wGPN(i) 
     ENDDO
  ENDDO
  !
  ! Set up the WENO matrices
  IF(MOD(N,2).EQ.0) THEN                     ! odd order scheme <=> even order polynomial
     WENO_left  = (/ -N/2, -N, 0, 0 /)       ! left  index for stencils
     WENO_right = (/  N/2,  0, N, 0 /)       ! right index for stencils
     WENO_nStencil = 3
     WENO_lambda = (/ WENO_lambdaC, 1., 1., 0. /) 
  ELSE 
     N2 = N/2 
     WENO_left  = (/ -N2-1, -N2,   -N, 0  /)  ! left  index for stencils
     WENO_right = (/  N2,    N2+1,  0, N  /)  ! right index for stencils
     WENO_nStencil = 4
     WENO_lambda = (/ WENO_lambdaC, WENO_lambdaC, 1., 1. /)     
  ENDIF
  IF(N.EQ.1) THEN
     WENO_nStencil = 2 
     WENO_left  = (/ -1, 0, 0, 0  /)  ! left  index for stencils
     WENO_right = (/  0, 1, 0, 0  /)  ! right index for stencils
     WENO_lambda = (/ 1., 1., 0., 0. /)             
  ENDIF
  ALLOCATE( WENO_iMatrix(N+1,N+1,WENO_nStencil) ) 
  DO iStencil = 1, WENO_nStencil
     counter = 1
     DO j = WENO_left(iStencil), WENO_right(iStencil)
        IntPhiPsi = 0. 
        DO iGP = 1, N+1
           xGP1 = REAL(j) + xiGPN(iGP) 
           CALL BaseFunc1D(phi,phi_xi,xGP1) 
           IntPhiPsi(1,:) = IntPhiPsi(1,:) + wGPN(iGP)*phi(:)   
        ENDDO
        A(counter,:) = IntPhiPsi(1,:)  
        counter = counter + 1
     ENDDO
     CALL MatrixInverse(N+1,A,iA)
     WENO_iMatrix(:,:,iStencil) = iA  
  ENDDO
  !
  CALL WriteDataGnuplot
  !
  CONTINUE 
  !
END SUBROUTINE ADERFVInit


SUBROUTINE InitialField(u0,xGP) 
  USE typesDef
  IMPLICIT NONE
  ! Argument list 
  REAL, INTENT(IN ) :: xGP(d)             ! spatial position vector 
  REAL, INTENT(OUT) :: u0(nVar)           ! initial data vector in terms of conserved variables 
  ! Local variables 
  REAL :: VBase(nVar), ampl(nVar), sigma(d) 
  REAL :: V0(nVar) 
  ! 
  
  ! Gaussian perturbation 
  ! sigma = (/ 0.05, 0.05, 0.05 /)       ! half-width
  ! VBase(:) = (/ 1., 0., 0., 0., 1. /)  ! base-state 
  ! ampl(:)  = 0.                        ! perturbation amplitude vector 
  ! ampl(5)   = 1e-3                     ! 
  ! V0(:) = VBase(:) + ampl(:)*EXP( -0.5*SUM(xGP(1:nDim)**2/sigma(1:nDim)**2) )    

  IF ( xGP(1) < 0.0) THEN
     
     V0(1) = 0.445
     V0(2) = 0.698
     V0(3) = 0.0
     V0(4) = 0.0
     V0(5) = 3.528
     
  ELSE

     V0(1) = 0.5
     V0(2) = 0.0
     V0(3) = 0.0
     V0(4) = 0.0
     V0(5) = 0.571

  ENDIF
  
  ! A simple debug check for the computation of derivatives 
  !u0 = 0. 
  !u0(1) = 0.123 !*xGP(1) 
  CALL PDEPrim2Cons(u0,V0) 
  !
END SUBROUTINE InitialField


SUBROUTINE GetL0Neighbor(tmpneighbor,iElem)
  USE TypesDef   
  IMPLICIT NONE 
  ! Argument list
  INTEGER :: iElem 
  INTEGER :: tmpneighbor(-ss:ss)
  ! Local Variables 
  INTEGER :: i,ii 
  INTEGER :: idx  
  !    
  i = iElem
  !
  DO ii = -ss, ss
     idx = i + ii
     IF(idx.LT.1) THEN 
        tmpneighbor(ii) = iElem
        CYCLE
     ENDIF
     IF((idx-IMAX).GT.0) THEN 
        tmpneighbor(ii) = iElem
        CYCLE 
     ENDIF
     IF(idxe0(idx).NE.0) THEN
        tmpneighbor(ii) = idxe0(idx) 
     ELSE
        tmpneighbor(ii) = iElem
     ENDIF
  ENDDO
  !    
END SUBROUTINE GetL0Neighbor
    
