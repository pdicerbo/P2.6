!----------------------------------------------------------------------!
!                                                                      !
!    This code has been in great part inspired by the routines         !
!    by E. F. Toro in his course on Riemann solvers and numerical      !
!    methods for fluid dynamics                                        !
!                                                                      !
!     Much of the theory can be found Toro, E. F., "Riemann Solvers    !
!                      and Numerical Methods for Fluid Dynamics"       !
!                      Springer-Verlag,  Second Edition, 1999          !
!                                                                      !
!     First-Order Godunov schemes for the one-dimensional              !
!                   Euler equations                                    !
!                                                                      !
!     Purpose: to solve the time-dependent one dimensional Euler       !
!              equations for an ideal gas by Godunov methods with      !
!              several approximate Rieman solvers of the               !
!              approximate flux type, namely                           !
!                                                                      !
!              INTFLX = 1: The HLL Riemann solver                      !
!              INTFLX = 2: The HLLC Riemann solver                     !
!                                                                      !
!     Input  file: file.input (initial data)                  !
!     Output files: *.dat (num. results for density, pressure, etc.)   !
!     Output files: *_exct.dat (exact results for density, etc.)       !
!                                                                      !
!----------------------------------------------------------------------!
!
!     Driver program
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      INTEGER INTFLX, CELLS, N, NFREQU, NTMAXI, i, IDIM
!
      REAL::    CFLCOE, PSCALE, TIME, TIMDIF, TIMEOU, TIMTOL
!
      COMMON /DRIVER/ CFLCOE, INTFLX, CELLS, NFREQU, NTMAXI, TIMEOU,    &
           &                PSCALE
!
      DATA TIME, TIMTOL /0.0, 1.0E-06/
!
 
      PARAMETER (IDIM = 3000)
!
      REAL FI
 
      DIMENSION FI(2, 3,-1:IDIM+2)
 
      COMMON /FLUXES/ FI
 
 
!     Parameters of problem are read in from file "finite_volume.input"
!
 
!     CHOOSE INITIAL CONDITIONS
!     TEST_NUMBER=1  (Modified Sod)
!     TEST_NUMBER=2  (123 problem)
!     TEST_NUMBER=3  (Left Woodward & Colella)
!     TEST_NUMBER=4  (Collision of 2 shocks)
!     TEST_NUMBER=5  (Stationary Contact ) 
!     TEST_NUMBER=6  (Woodward and Colella)
 
      
      CALL READER
 
!
!     Initial conditions are set up
!
      CALL INITIA(CELLS)
!
!     Output at time t=0 
      CALL OUTPUT(CELLS, PSCALE, TIME)
 
!     Time marching procedure
!
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)'   Time step N        TIME           TIMEOU'
      WRITE(6,*)'---------------------------------------------'
!
 
      DO 10 N = 1, NTMAXI
!
!        Boundary conditions are set
!
         CALL BCONDI(CELLS)
!
!        Courant-Friedrichs-Lewy (CFL) condition imposed
!
         CALL CFLCON(CFLCOE, CELLS, N, TIME, TIMEOU)
!
!        Intercell numerical fluxes are computed. Two
!        choices are available
!
         IF(INTFLX.EQ.1)CALL HLL(CELLS)
         IF(INTFLX.EQ.2)CALL HLLC(CELLS)
 
!
!        For comparison:
!        Intercell numerical fluxes are computed from Exact solution
!        of the Riemann problem         
 
         CALL RPGODU(CELLS, INTFLX)
 
!        Solution is updated according to conservative formula
!       
         CALL UPDATE(CELLS)
!
 
         IF(MOD(N,NFREQU).EQ.0)WRITE(6,20)N, TIME, TIMEOU
!
!        Check output time
!
         TIMDIF = ABS(TIME - TIMEOU)
!
         CALL OUTPUT(CELLS, PSCALE, TIME)
 
         IF(TIMDIF.LE.TIMTOL)THEN
 
         WRITE(6,*)'---------------------------------------------'
         WRITE(6,*)'   Number of time steps = ',N
         WRITE(6,*)'---------------------------------------------'
!
 
            GOTO 30
         ENDIF
!
 10   CONTINUE
 
!
 20   FORMAT(I12,6X,2(F12.7, 4X))
 21   FORMAT(' ')
 30   CONTINUE
!
      END
!----------------------------------------------------------------------*

!======================================================================*
      SUBROUTINE READER
!======================================================================*
!
!     Purpose: to read initial parameters of the problem
!
!     Input variables
!
!     DOMLEN    : Domain length
!     DIAPH1    : Position of diaphragm 1
!     CELLS     : Number of computing cells
!     GAMMA     : Ratio of specific heats
!     TIMEOU    : Output time
!     DLINIT    : Initial density  on left section of tube
!     ULINIT    : Initial velocity on left section of tube
!     PLINIT    : Initial pressure on left section of tube
!     DMINIT    : Initial density  on middle section of tube
!     UMINIT    : Initial velocity on middle section of tube
!     PMINIT    : Initial pressure on middle section of tube
!     DRINIT    : Initial density  on right section of tube
!     URINIT    : Initial velocity on right section of tube
!     PRINIT    : Initial pressure on right section of tube
!     DIAPH2    : Position of diaphragm 2
!     CFLCOE    : Courant number coefficient
!     IBCLEF    : Type of left boundary conditions
!     IBCRIG    : Type of right boundary conditions
!     NFREQU    : Output frequency to screen
!     NTMAXI    : Maximum number of time steps
!     PSCALE    : Pressure scaling factor
!     INTFLX    : Choice of intercell flux
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      INTEGER INTFLX , IBCLEF, IBCRIG, CELLS, NFREQU, NTMAXI,           &
           & TEST_NUMBER, TEST_ITER
      !
      REAL    CFLCOE, DOMLEN, DIAPH1, DIAPH2, PSCALE, TIMEOU,           &
           &  DLINIT, ULINIT, PLINIT, DMINIT, UMINIT, PMINIT, DRINIT,   &
           &  URINIT, PRINIT,                                           &
           &  GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      !
      COMMON /BOUNDA/ IBCLEF, IBCRIG
      COMMON /DOMAIN/ DOMLEN, DIAPH1, DIAPH2
      COMMON /DRIVER/ CFLCOE, INTFLX, CELLS, NFREQU, NTMAXI, TIMEOU,    &
           &                PSCALE
      COMMON /INISTA/ DLINIT, ULINIT, PLINIT, DMINIT, UMINIT, PMINIT,   &
           &                DRINIT, URINIT, PRINIT
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      !
      OPEN(UNIT = 1, FILE = 'inputfile.par', STATUS = 'UNKNOWN')
!
 
      READ(1,*)TEST_NUMBER
      READ(1,*)
      READ(1,*)
      READ(1,*)
      READ(1,*)
      READ(1,*)
      READ(1,*)
      READ(1,*)
      READ(1,*)
      READ(1,*)
 
      DO TEST_ITER=1,TEST_NUMBER
      READ(1,*)DOMLEN
      READ(1,*)DIAPH1
      READ(1,*)CELLS
      READ(1,*)GAMMA
      READ(1,*)TIMEOU
      READ(1,*)DLINIT
      READ(1,*)ULINIT
      READ(1,*)PLINIT
      READ(1,*)DMINIT
      READ(1,*)UMINIT
      READ(1,*)PMINIT
      READ(1,*)DRINIT
      READ(1,*)URINIT
      READ(1,*)PRINIT
      READ(1,*)DIAPH2
      READ(1,*)CFLCOE
      READ(1,*)IBCLEF
      READ(1,*)IBCRIG
      READ(1,*)NFREQU
      READ(1,*)NTMAXI
      READ(1,*)PSCALE
      READ(1,*)INTFLX
      READ(1,*)
      READ(1,*)
      READ(1,*)
      ENDDO
 
!
      CLOSE(1)
!
!     Input data is echoed to screen
!
      WRITE(6,*)
      WRITE(6,*)'Input data echoed to screen'
      WRITE(6,*)
      WRITE(6,*)'DOMLEN = ',DOMLEN
      WRITE(6,*)'DIAPH1 = ',DIAPH1
      WRITE(6,*)'CELLS  = ',CELLS
      WRITE(6,*)'GAMMA  = ',GAMMA
      WRITE(6,*)'TIMEOU = ',TIMEOU
      WRITE(6,*)'DLINIT = ',DLINIT
      WRITE(6,*)'ULINIT = ',ULINIT
      WRITE(6,*)'PLINIT = ',PLINIT
      WRITE(6,*)'DMINIT = ',DMINIT
      WRITE(6,*)'UMINIT = ',UMINIT
      WRITE(6,*)'PMINIT = ',PMINIT
      WRITE(6,*)'DRINIT = ',DRINIT
      WRITE(6,*)'URINIT = ',URINIT
      WRITE(6,*)'PRINIT = ',PRINIT
      WRITE(6,*)'DIAPH2 = ',DIAPH2
      WRITE(6,*)'CFLCOE = ',CFLCOE
      WRITE(6,*)'IBCLEF = ',IBCLEF
      WRITE(6,*)'IBCRIG = ',IBCRIG
      WRITE(6,*)'NFREQU = ',NFREQU
      WRITE(6,*)'NTMAXI = ',NTMAXI
      WRITE(6,*)'PSCALE = ',PSCALE
      WRITE(6,*)'INTFLX = ',INTFLX
!
      END
!----------------------------------------------------------------------*

!======================================================================*
      SUBROUTINE INITIA(CELLS)
!======================================================================*
!
!     Purpose: to set initial conditions
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      INTEGER I, CELLS, IDIM
      REAL    DOMLEN, DIAPH1, DIAPH2, DT, DX, D, U, P, CS,              &
           &   DLINIT, ULINIT, PLINIT, DMINIT, UMINIT, PMINIT, DRINIT,  &
           &   URINIT, PRINIT, XPOS,                                    &
           &   GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      !
      PARAMETER (IDIM = 3000)
!
      DIMENSION D(2,-1:IDIM+2), U(2,  -1:IDIM+2)
      DIMENSION P(2,-1:IDIM+2),CS(2,3,-1:IDIM+2)
!
      COMMON /DOMAIN/ DOMLEN, DIAPH1, DIAPH2
      COMMON /INISTA/ DLINIT, ULINIT, PLINIT, DMINIT, UMINIT, PMINIT,   &
           &                DRINIT, URINIT, PRINIT
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /PRIMIT/ D, U, P
      COMMON /CONSER/ CS
      COMMON /MESHPA/ DT, DX
!
!     Compute gamma related constants
!
      G1 = (GAMMA - 1.0)/(2.0*GAMMA)
      G2 = (GAMMA + 1.0)/(2.0*GAMMA)
      G3 = 2.0*GAMMA/(GAMMA - 1.0)
      G4 = 2.0/(GAMMA - 1.0)
      G5 = 2.0/(GAMMA + 1.0)
      G6 = (GAMMA - 1.0)/(GAMMA + 1.0)
      G7 = (GAMMA - 1.0)/2.0
      G8 = GAMMA - 1.0
!
!     Calculate mesh size DX
!
      DX = DOMLEN/REAL(CELLS)
!
!     Set initial data in tube of length DOMLEN, which is divided
!     into 3 sections by diaphragms at positions DIAPH1 and DIAPH2
!
      DO 10 I = 1, CELLS
!
         XPOS = (REAL(I) - 0.5)*DX
!
         IF(XPOS.LE.DIAPH1)THEN
!
!           Set initial values in left section of domaim
!
            D(:,I) = DLINIT
            U(:,I) = ULINIT
            P(:,I) = PLINIT
         ENDIF
!
         IF(XPOS.GT.DIAPH1.AND.XPOS.LE.DIAPH2)THEN
!
!           Set initial values in middle section of domaim
!
            D(:,I) = DMINIT
            U(:,I) = UMINIT
            P(:,I) = PMINIT
         ENDIF
!
         IF(XPOS.GT.DIAPH2)THEN
!
!           Set initial values in right section of domaim
!
            D(:,I) = DRINIT
            U(:,I) = URINIT
            P(:,I) = PRINIT
         ENDIF
!
!        Compute conserved variables
!
         CS(:,1,I) = D(:,I)
         CS(:,2,I) = D(:,I)*U(:,I)
         CS(:,3,I) = 0.5*CS(:,2,I)*U(:,I) + P(:,I)/G8
!
 10   CONTINUE
!
      END
!----------------------------------------------------------------------*
      
!======================================================================*
      SUBROUTINE BCONDI(CELLS)
!======================================================================*
!
!     Purpose: to set boundary conditions
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      INTEGER IBCLEF, IBCRIG, CELLS, IDIM
!
      REAL    D, U, P
!
      PARAMETER (IDIM = 3000)
!
      DIMENSION D(2,-1:IDIM+2), U(2,-1:IDIM+2), P(2,-1:IDIM+2)
!
      COMMON /PRIMIT/ D, U, P
      COMMON /BOUNDA/ IBCLEF, IBCRIG
!
      IF(IBCLEF.EQ.0)THEN
!
!        Transmissive boundary conditions on the left
!
         D(:,0)  =  D(:,1)
         U(:,0)  =  U(:,1)
         P(:,0)  =  P(:,1)
!
      ELSE
!
!        Reflective boundary conditions on the left
!
         D(:,0)  =  D(:,1)
         U(:,0)  = -U(:,1)
         P(:,0)  =  P(:,1)
!
      ENDIF
!
      IF(IBCRIG.EQ.0)THEN
!
!        Transmissive  boundary conditions on the right
!
         D(:,CELLS + 1) =  D(:,CELLS)
         U(:,CELLS + 1) =  U(:,CELLS)
         P(:,CELLS + 1) =  P(:,CELLS)
!
      ELSE
!
!        Reflective boundary conditions on the right
!
         D(:,CELLS + 1) =  D(:,CELLS)
         U(:,CELLS + 1) = -U(:,CELLS)
         P(:,CELLS + 1) =  P(:,CELLS)
!
      ENDIF
!
      END
!----------------------------------------------------------------------*

!======================================================================*
      SUBROUTINE CFLCON(CFLCOE, CELLS, N, TIME, TIMEOU)
!======================================================================*
!
!     Purpose: to apply the CFL condition to find a stable time
!              step size DT
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      INTEGER I, CELLS, IDIM, N
!
      REAL    C, CFLCOE, D, DT, DX, P, SMAX, SBEXTD, TIME,              &
           &        TIMEOU, U,                                          &
           &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      !
      PARAMETER (IDIM = 3000)
!
      DIMENSION D(2,-1:IDIM+2), U(2,-1:IDIM+2)
      DIMENSION P(2,-1:IDIM+2), C(2,-1:IDIM+2)
!
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /PRIMIT/ D, U, P
      COMMON /SOUNDS/ C
      COMMON /MESHPA/ DT, DX
!
      SMAX = 0.0
!
!     Find maximum velocity SMAX present in data
!
      DO 10 I = 0, CELLS + 1
!
!        Compute speed of sound
!
         C(:,I)   = SQRT(GAMMA*P(:,I)/D(:,I))
!
         SBEXTD  = ABS(U(1,I)) + C(1,I)
         IF(SBEXTD.GT.SMAX)SMAX = SBEXTD
 10   CONTINUE
!
!     Compute time step DT, for early times reduce its size
!
      DT = CFLCOE*DX/SMAX
!
!     For early times DT is reduced to compensate for approximate
!     calculation of SMAX
!
      IF(N.LE.5)DT = 0.2*DT
!
!     Check size of DT to avoid exceeding output time
!
      IF((TIME + DT).GT.TIMEOU)THEN
!
!        Recompute DT
!
         DT = TIMEOU - TIME
      ENDIF
!
!     Find current time
!
      TIME = TIME + DT
!
      END
!----------------------------------------------------------------------*

!======================================================================*
      SUBROUTINE OUTPUT(CELLS, PSCALE, TIME)
!======================================================================*
!
!     Purpose: to output the solution at a specified time TIMEOU
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      INTEGER I, CELLS, IDIM
!
      REAL    D, DT, DX, ENERGI, P, PSCALE, U, XPOS,                    &
           &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8, TIME
!
      PARAMETER (IDIM = 3000)
!
      DIMENSION D(2,-1:IDIM+2),U(2,-1:IDIM+2),P(2,-1:IDIM+2)
!
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /PRIMIT/ D, U, P
      COMMON /MESHPA/ DT, DX
!
      OPEN(UNIT = 1, FILE = 'out/rho.dat', STATUS = 'UNKNOWN',     &
           & POSITION = 'APPEND')
      !
      OPEN(UNIT = 2, FILE = 'out/vel.dat', STATUS = 'UNKNOWN',     &
           & POSITION = 'APPEND')
      !
      OPEN(UNIT = 3, FILE = 'out/prs.dat', STATUS = 'UNKNOWN',     &
           & POSITION = 'APPEND')
      !
      OPEN(UNIT = 4, FILE = 'out/ene.dat', STATUS = 'UNKNOWN',     &
           & POSITION = 'APPEND')
      
      
      OPEN(UNIT = 100, FILE = 'out/rho_exct.dat', STATUS = 'UNKNOWN',    &
           & POSITION = 'APPEND')
      !
      OPEN(UNIT = 200, FILE = 'out/vel_exct.dat', STATUS = 'UNKNOWN',    &
     & POSITION = 'APPEND')
!
      OPEN(UNIT = 300, FILE = 'out/prs_exct.dat', STATUS = 'UNKNOWN',    &
     & POSITION = 'APPEND')
      !
      OPEN(UNIT = 400, FILE = 'out/ene_exct.dat', STATUS = 'UNKNOWN',    &
           & POSITION = 'APPEND')
      
      !     TIME OF THE TIME-STEP
      WRITE(1,22) TIME
      WRITE(2,22) TIME
      WRITE(3,22) TIME
      WRITE(4,22) TIME
 
      WRITE(100,22) TIME
      WRITE(200,22) TIME
      WRITE(300,22) TIME
      WRITE(400,22) TIME
!
 
      DO I   = 1, CELLS
         XPOS   = (REAL(I) - 0.5)*DX
 
         ENERGI =  P(1,I)/D(1,I)/G8/PSCALE
         WRITE(1,20)XPOS, D(1,I)
         WRITE(2,20)XPOS, U(1,I)
         WRITE(3,20)XPOS, P(1,I)/PSCALE
         WRITE(4,20)XPOS, ENERGI
 
         ENERGI =  P(2,I)/D(2,I)/G8/PSCALE
         WRITE(100,20)XPOS, D(2,I)
         WRITE(200,20)XPOS, U(2,I)
         WRITE(300,20)XPOS, P(2,I)/PSCALE
         WRITE(400,20)XPOS, ENERGI
 
      ENDDO   
 
 
!     EMPTY LINE AFTER EACH TIME-STEP      
      WRITE(1,21)
      WRITE(2,21)
      WRITE(3,21)
      WRITE(4,21)
 
      WRITE(100,21)
      WRITE(200,21)
      WRITE(300,21)
      WRITE(400,21)
!
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
 
      CLOSE(100)
      CLOSE(200)
      CLOSE(300)
      CLOSE(400)
 
!
 20   FORMAT(2(F14.6,2X))
 21   FORMAT(' ')
 22   FORMAT('"Time = ',E13.6)
 
!
      END
!----------------------------------------------------------------------*

!======================================================================*
      SUBROUTINE UPDATE(CELLS)
!======================================================================*
!
!     Purpose: to update the solution according to the conservative
!              formula and compute physical variables
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      INTEGER I, K, CELLS, IDIM
!
      REAL    DT, DX, DTODX, D, U, P, CS, FI,                           &
           &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      !
      PARAMETER (IDIM = 3000)
!
      DIMENSION D(2,-1:IDIM+2), U(2,-1:IDIM+2), P(2,-1:IDIM+2),         &
           &          CS(2,3,-1:IDIM+2), FI(2,3,-1:IDIM+2)
      !
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /PRIMIT/ D, U, P
      COMMON /CONSER/ CS
      COMMON /FLUXES/ FI
      COMMON /MESHPA/ DT, DX
!
      DTODX = DT/DX
!
      DO I = 1, CELLS
!
         DO K = 1, 3
            CS(1,K,I) = CS(1,K,I) + DTODX*(FI(1,K,I-1) - FI(1,K,I))
         ENDDO
 
         DO K = 1, 3
            CS(2,K,I) = CS(2,K,I) + DTODX*(FI(2,K,I-1) - FI(2,K,I))
         ENDDO
!
      ENDDO
!
!     Compute physical variables
!
      DO I = 1, CELLS
         D(1,I) = CS(1,1,I)
         U(1,I) = CS(1,2,I)/D(1,I)
         P(1,I) = G8*(CS(1,3,I) - 0.5*CS(1,2,I)*U(1,I))
 
         D(2,I) = CS(2,1,I)
         U(2,I) = CS(2,2,I)/D(2,I)
         P(2,I) = G8*(CS(2,3,I) - 0.5*CS(2,2,I)*U(2,I))
      ENDDO
!
      END
!----------------------------------------------------------------------*

!======================================================================*
      SUBROUTINE HLLC(CELLS)
!======================================================================*
!
!     Purpose: to compute an intercell Godunov flux using
!              the HLLC approximate Riemann solver. See Chap. 10
!              Ref. 1
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      INTEGER  I, CELLS, IDIM, K
!
      REAL     C, CL, CR, CS, CSL, CSR, D, DL, DR, ENEL, ENER,          & 
           &         FD, FI, P, PL, PR, SL, SM, SR, U, UL, UR,          &
           &         GAMMA, G1, G2, G3, G4, G5, G6, G7, G8            
!
      PARAMETER (IDIM = 3000)
!
      DIMENSION D(2,-1:IDIM+2), U(2,-1:IDIM+2), P(2,-1:IDIM+2),         &
     &          C(2,-1:IDIM+2),                                         &
     &          CS(2,3,-1:IDIM+2), FD(3,-1:IDIM+2), FI(2,3,-1:IDIM+2),  &
     &          CSL(3), CSR(3)
!
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR
      COMMON /PRIMIT/ D, U, P
      COMMON /SOUNDS/ C
      COMMON /CONSER/ CS
      COMMON /FLUXES/ FI
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
!
!     Compute fluxes on data and conserved variables
!     in fictitious cells
!
      DO I = 0, CELLS + 1
!
         IF(I.LT.1.OR.I.GT.CELLS)THEN
            CS(1,1,I) = D(1,I)
            CS(1,2,I) = D(1,I)*U(1,I)
            CS(1,3,I) = 0.5* D(1,I)*U(1,I)*U(1,I) + P(1,I)/G8
         ENDIF
!
         FD(1,I) = CS(1,2,I)
         FD(2,I) = CS(1,2,I)*U(1,I)   + P(1,I)
         FD(3,I) = U(1,I)*(CS(1,3,I)  + P(1,I))
!
      ENDDO
!
!     Solve Riemann problem (i,i+1) and store quantities in I
!
      DO 20 I = 0, CELLS
!
         DL = D(1,I)
         UL = U(1,I)
         PL = P(1,I)
         CL = C(1,I)
!
         DR = D(1,I + 1)
         UR = U(1,I + 1)
         PR = P(1,I + 1)
         CR = C(1,I + 1)
!
!        Calculate estimates for wave speeds using adaptive
!        approximate-state Riemann solvers
!
         CALL ESTIME(SL, SM, SR)
!
         IF(SL.GE.0.0)THEN
!
!           Right-going supersonic flow
!
            DO K = 1, 3
               FI(1, K, I) = FD(K, I)
            ENDDO
!
         ENDIF
!
         IF(SL.LE.0.0.AND.SR.GE.0.0)THEN
!
!           Subsonic flow
!
            IF(SM.GE.0.0)THEN
!
!              Subsonic flow to the right
!
               ENEL   = CS(1, 3, I)/DL
               ENEL   = ENEL + (SM - UL)*(SM + PL/(DL*(SL - UL)))
               CSL(1) = DL*(SL - UL)/(SL - SM)
               CSL(2) = CSL(1)*SM
               CSL(3) = CSL(1)*ENEL
!
               DO K = 1, 3
                  FI(1, K, I) = FD(K, I) + SL*(CSL(K) - CS(1, K, I))
               ENDDO
!
            ELSE
!
!              Subsonic flow to the left
!
               ENER   = CS(1, 3,I+1)/DR
               ENER   = ENER + (SM - UR)*(SM + PR/(DR*(SR - UR)))
               CSR(1) = DR*(SR - UR)/(SR - SM)
               CSR(2) = CSR(1)*SM
               CSR(3) = CSR(1)*ENER
!
               DO K = 1, 3
                  FI(1, K, I) = FD(K, I + 1) + SR*(CSR(K) - CS(1, K, I + 1))
               ENDDO
            ENDIF
         ENDIF
!
         IF(SR.LE.0.0)THEN
!
!           Left-going supersonic flow
!
            DO K = 1, 3
               FI(1, K, I) = FD(K, I + 1)
            ENDDO
!
        ENDIF
 
 20   CONTINUE
!
      END
!----------------------------------------------------------------------*

!======================================================================*
      SUBROUTINE HLL(CELLS)
!======================================================================*
!
!     Purpose: to compute an intercell Godunov flux using
!              the HLL approximate Riemann solver. See Chap 10,
!              Ref. 1 and original references therein
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      INTEGER  I, CELLS, IDIM, K
!
      REAL     C, CL, CR, CS, D, DL, DR, FD, FI, HLLFLUX,               &
           &         P, PL, PR, SL, SM, SR, U, UL, UR,                  &
           &         GAMMA, G1, G2, G3, G4, G5, G6, G7, G8, XPOS, DX, DT              
!
      PARAMETER (IDIM = 3000)
!
      DIMENSION D(2, -1:IDIM+2), U(2, -1:IDIM+2), P(2, -1:IDIM+2),      &
           &    C(2, -1:IDIM+2),                                        &
           &    CS(2, 3,-1:IDIM+2), FD(3,-1:IDIM+2), FI(2, 3,-1:IDIM+2)  
      !
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR
      COMMON /PRIMIT/ D, U, P
      COMMON /SOUNDS/ C
      COMMON /CONSER/ CS
      COMMON /FLUXES/ FI
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /MESHPA/ DT, DX
!
!     Compute fluxes on data and conserved variables
!     in fictitious cells
!
      DO 10 I = 0, CELLS + 1
!
         IF(I.LT.1.OR.I.GT.CELLS)THEN
            CS(1, 1,I) = D(1, I)
            CS(1, 2,I) = D(1, I)*U(1, I)
            CS(1, 3,I) = 0.5* D(1, I)*U(1, I)*U(1, I) + P(1, I)/G8
         ENDIF
!
         FD(1,I) = CS(1, 2,I)
         FD(2,I) = CS(1, 2,I)*U(1, I)   + P(1, I)
         FD(3,I) = U(1, I)*(CS(1, 3,I)  + P(1, I))
!
 10   CONTINUE
!
!     Solve Riemann problem (i,i+1) and store quantities in I
!
      DO 20 I = 0, CELLS
!
         XPOS = (REAL(I) - 0.5)*DX
!ciao
         DL = D(1, I)
         UL = U(1, I)
         PL = P(1, I)
         CL = C(1, I)
!
         DR = D(1, I + 1)
         UR = U(1, I + 1)
         PR = P(1, I + 1)
         CR = C(1, I + 1)
!
!        Calculate estimates for wave speeds using adaptive
!        approximate-state Riemann solvers
!
         CALL ESTIME(SL, SM, SR)
!
         IF(SL.GE.0.0)THEN
!
!           Right-going supersonic flow
!
            DO K = 1, 3
               FI(1, K, I) = FD(K, I)
            ENDDO
!
         ENDIF
!
         IF(SL.LE.0.0.AND.SR.GE.0.0)THEN
!
!           Subsonic flow
!
            DO K = 1, 3
               HLLFLUX = SR*FD(K, I) - SL*FD(K, I + 1)
               HLLFLUX = HLLFLUX + SL*SR*(CS(1, K, I + 1) - CS(1, K, I))
               FI(1, K,I) = HLLFLUX/(SR - SL)
            ENDDO
!
         ENDIF
!
         IF(SR.LE.0.0)THEN
!
!           Left-going supersonic flow
!
            DO K = 1, 3
               FI(1, K, I) = FD(K, I + 1)
            ENDDO
!
         ENDIF
 
 20   CONTINUE
!
      END
!----------------------------------------------------------------------*
 
!======================================================================*
      SUBROUTINE ESTIME(SL, SM, SR)
!======================================================================*
!
!     Purpose: to compute wave speed estimates for the HLLC Riemann
!              solver using and adaptive approximate-state Riemann
!              solver including the PVRS, TRRS and TSRS solvers.
!              See Chap. 9, Ref. 1
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      REAL    DL, UL, PL, CL, DR, UR, PR, CR,                           &
           &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8 ,             &
           &        CUP, GEL, GER, PM, PMAX, PMIN, PPV, PQ,             &
           &        PTL, PTR, QMAX, QUSER, SL, SM, SR, UM  
!
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR
!
      QUSER = 2.0
!
!     Compute guess pressure from PVRS Riemann solver
!
      CUP  = 0.25*(DL + DR)*(CL + CR)
      PPV  = 0.5*(PL + PR) + 0.5*(UL - UR)*CUP
      PPV  = MAX(0.0, PPV)
      PMIN = MIN(PL,  PR)
      PMAX = MAX(PL,  PR)
      QMAX = PMAX/PMIN
!
      IF(QMAX.LE.QUSER.AND.(PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN
!
!        Select PRVS Riemann solver
!
         PM = PPV
         UM = 0.5*(UL + UR) + 0.5*(PL - PR)/CUP
!
      ELSE
         IF(PPV.LT.PMIN)THEN
!
!           Select Two-Rarefaction Riemann solver
!
            PQ  = (PL/PR)**G1
            UM  = (PQ*UL/CL + UR/CR + G4*(PQ - 1.0))/(PQ/CL + 1.0/CR)
            PTL = 1.0 + G7*(UL - UM)/CL
            PTR = 1.0 + G7*(UM - UR)/CR
            PM  = 0.5*(PL*PTL**G3 + PR*PTR**G3)
!
         ELSE
!
!           Use Two-Shock Riemann solver with PVRS as estimate
!
            GEL = SQRT((G5/DL)/(G6*PL + PPV))
            GER = SQRT((G5/DR)/(G6*PR + PPV))
            PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER)
            UM  = 0.5*(UL + UR) + 0.5*(GER*(PM - PR) - GEL*(PM - PL))
         ENDIF
      ENDIF
!
!     Find speeds
!
      IF(PM.LE.PL)THEN
         SL = UL - CL
      ELSE
         SL = UL - CL*SQRT(1.0 + G2*(PM/PL - 1.0))
      ENDIF
!
      SM = UM
!
      IF(PM.LE.PR)THEN
         SR = UR + CR
      ELSE
         SR = UR + CR*SQRT(1.0 + G2*(PM/PR - 1.0))
      ENDIF
!
      END
!----------------------------------------------------------------------*
 
!======================================================================*
      SUBROUTINE RPGODU(CELLS, INTFLX)
!======================================================================*
!
!     Purpose: to compute an intercell Godunov flux using
!              the exact Riemann solver
!
!     Theory is found in Ref. 1, Chaps. 4 and 9
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      INTEGER  I, CELLS, IDIM, INTFLX
!
      REAL     C, CL, CR, D, DL, DR, DSAM, ENERGS, FI, P, PL,           &
           &         PM, PR, PSAM, U, UL, UM, UR, USAM, XOVERT,         &
           &         GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
!
      PARAMETER (IDIM = 3000)
!
      DIMENSION D(2, -1:IDIM+2), U(2, -1:IDIM+2), P(2, -1:IDIM+2),      &
           &          C(2, -1:IDIM+2),                                  &
           &          FI(2, 3,-1:IDIM+2)
!
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR
      COMMON /PRIMIT/ D, U, P
      COMMON /SOUNDS/ C
      COMMON /FLUXES/ FI
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
!
!     Solve Riemann problem (i,i+1) and store quantities in i
!
!     Set value x/t = 0 (along t-axis)
!
      XOVERT = 0.0
!
      DO I = 0, CELLS
!
         DL = D(2, I)
         UL = U(2, I)
         PL = P(2, I)
         CL = C(2, I)
!
         DR = D(2, I + 1)
         UR = U(2, I + 1)
         PR = P(2, I + 1)
         CR = C(2, I + 1)
!
!        Solve Riemann problem exactly for star values
!        of pressure PM and velocity UM
!
         CALL EXACT(PM, UM)
!
!        Sample solution of Riemann problem at x/t = 0
!        to find Godunov state (DSAM, PSAM, USAM)
!
         CALL SAMPLE(PM, UM, XOVERT, DSAM, PSAM, USAM)
!
!        Compute intercell flux at Godunov state
!
         FI(2, 1,I) = DSAM*USAM
         FI(2, 2,I) = DSAM*USAM*USAM + PSAM
         ENERGS  = 0.5*USAM*USAM*DSAM + PSAM/G8
         FI(2, 3,I) = USAM*(ENERGS + PSAM)
!
      ENDDO
!
      END
!----------------------------------------------------------------------*

!======================================================================*
      SUBROUTINE EXACT(P, U)
!======================================================================*
!
      IMPLICIT NONE
!
!     Purpose: to compute PM and UM  in the Star Region using
!              the EXACT Riemann solver
!
!              We use exact relations for density and EXACT
!              solution for sonic flow in sampling routine
!
!     Declaration of variables
!
      INTEGER I, NRITER
      REAL    DL, UL, PL, CL, DR, UR, PR, CR,                           &
           &        CHANGE, FL, FLD, FR, FRD, P, POLD, PSTART, TOLPRE,  &
           &        U, UDIFF
!
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR
      DATA TOLPRE, NRITER/1.0E-05, 20/
!
!     Guessed value PSTART is computed
!
      CALL GUESSP(PSTART)
!
      POLD  = PSTART
      UDIFF = UR - UL
!
      DO 10 I = 1, NRITER
!
         CALL PREFUN(FL, FLD, POLD, DL, PL, CL)
         CALL PREFUN(FR, FRD, POLD, DR, PR, CR)
         P      = POLD - (FL + FR + UDIFF)/(FLD + FRD)
         CHANGE = 2.0*ABS((P - POLD)/(P + POLD))
         IF(CHANGE.LE.TOLPRE)GOTO 20
         IF(P.LT.0.0)P = TOLPRE
         POLD  = P
!
 10   CONTINUE
!
      WRITE(6,*)'Divergence in Newton-Raphson iteration'
!
 20   CONTINUE
!
!     Compute velocity in Star region
!
      U = 0.5*(UL + UR + FR - FL)
!
      END
!----------------------------------------------------------------------*
      
!======================================================================*
      SUBROUTINE SAMPLE(PM, UM, S, D, P, U)
!======================================================================*
!
!     Purpose: to sample the solution throughout the wave pattern
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      REAL    DL, UL, PL, CL, DR, UR, PR, CR,                           &
           &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8,              &
           &        C, CML, CMR, D, P, PM, PML, PMR,  S,                &
           &        SHL, SHR, SL, SR, STL, STR, U, UM
!
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR
!
      IF(S.LE.UM)THEN
!
!        Sampling point lies to the left of the contact discontinuity
!
         IF(PM.LE.PL)THEN
!
!           Left rarefaction
!
            SHL = UL - CL
!
            IF(S.LE.SHL)THEN
!
!              Sampled point is left data state
!
               D = DL
               U = UL
               P = PL
            ELSE
               CML = CL*(PM/PL)**G1
               STL = UM - CML
!
               IF(S.GT.STL)THEN
!
!                 Sampled point is Star Left state
!
                  D = DL*(PM/PL)**(1.0/GAMMA)
                  U = UM
                  P = PM
               ELSE
!
!                 Sampled point is inside left fan
!
                  U = G5*(CL + G7*UL + S)
                  C = G5*(CL + G7*(UL - S))
                  D = DL*(C/CL)**G4
                  P = PL*(C/CL)**G3
               ENDIF
            ENDIF
         ELSE
!
!           Left shock
!
            PML = PM/PL
            SL  = UL - CL*SQRT(G2*PML + G1)
!
            IF(S.LE.SL)THEN
!
!              Sampled point is left data state
!
               D = DL
               U = UL
               P = PL
 
            ELSE
!
!              Sampled point is Star Left state
!
               D = DL*(PML + G6)/(PML*G6 + 1.0)
               U = UM
               P = PM
            ENDIF
         ENDIF
      ELSE
!
!        Sampling point lies to the right of the contact discontinuity
!
         IF(PM.GT.PR)THEN
!
!           Right shock
!
            PMR = PM/PR
            SR  = UR + CR*SQRT(G2*PMR + G1)
!
            IF(S.GE.SR)THEN
!
!              Sampled point is right data state
!
               D = DR
               U = UR
               P = PR
            ELSE
!
!              Sampled point is Star Right state
!
               D = DR*(PMR + G6)/(PMR*G6 + 1.0)
               U = UM
               P = PM
            ENDIF
         ELSE
!
!           Right rarefaction
!
            SHR = UR + CR
!
            IF(S.GE.SHR)THEN
!
!              Sampled point is right data state
!
               D = DR
               U = UR
               P = PR
            ELSE
               CMR = CR*(PM/PR)**G1
               STR = UM + CMR
!
               IF(S.LE.STR)THEN
!
!                 Sampled point is Star Right state
!
                  D = DR*(PM/PR)**(1.0/GAMMA)
                  U = UM
                  P = PM
               ELSE
!
!                 Sampled point is inside left fan
!
                  U = G5*(-CR + G7*UR + S)
                  C = G5*(CR - G7*(UR - S))
                  D = DR*(C/CR)**G4
                  P = PR*(C/CR)**G3
               ENDIF
            ENDIF
         ENDIF
      ENDIF
!
      END
!----------------------------------------------------------------------*

!======================================================================*
      SUBROUTINE GUESSP(PM)
!======================================================================*
!
!     Purpose: to provide a guessed value for pressure
!              in the Star Region
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      REAL    DL, UL, PL, CL, DR, UR, PR, CR,                           &
           &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8,              &   
           &        CUP, GEL, GER, PM, PMAX, PMIN, PPV, PQ,             &
           &        PTL, PTR, QMAX, QUSER, UM
!
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR
!
      QUSER = 2.0
!
!     Compute guess pressure from PVRS Riemann solver
!
      CUP  = 0.25*(DL + DR)*(CL + CR)
      PPV  = 0.5*(PL + PR) + 0.5*(UL - UR)*CUP
      PPV  = MAX(0.0, PPV)
      PMIN = MIN(PL,  PR)
      PMAX = MAX(PL,  PR)
      QMAX = PMAX/PMIN
!
      IF(QMAX.LE.QUSER.AND.(PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN
!
!        Select PVRS Riemann solver
!
         PM = PPV
      ELSE
         IF(PPV.LT.PMIN)THEN
!
!           Select Two-Rarefaction Riemann solver
!
            PQ  = (PL/PR)**G1
            UM  = (PQ*UL/CL + UR/CR + G4*(PQ - 1.0))/(PQ/CL + 1.0/CR)
            PTL = 1.0 + G7*(UL - UM)/CL
            PTR = 1.0 + G7*(UM - UR)/CR
            PM  = 0.5*(PL*PTL**G3 + PR*PTR**G3)
         ELSE
!
!           Select Two-Shock Riemann solver with PVRS as estimate
!
            GEL = SQRT((G5/DL)/(G6*PL + PPV))
            GER = SQRT((G5/DR)/(G6*PR + PPV))
            PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER)
         ENDIF
      ENDIF
!
      END
!----------------------------------------------------------------------*

!======================================================================*
      SUBROUTINE PREFUN(F,FD,P,DK,PK,CK)
!======================================================================*
!
!     Purpose: to evaluate the pressure functions FL and FR
!              in exact Riemann solver and their derivatives
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      REAL    AK, BK, CK, DK, F, FD, P, PK, PRATIO, QRT,                &
           &        GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      !
      COMMON /GAMMAS/ GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
!
      IF(P.LE.PK)THEN
!
!        Rarefaction wave
!
         PRATIO = P/PK
         F      = G4*CK*(PRATIO**G1 - 1.0)
         FD     = (1.0/(DK*CK))*PRATIO**(-G2)
      ELSE
!
!        Shock wave
!
         AK  = G5/DK
         BK  = G6*PK
         QRT = SQRT(AK/(BK + P))
         F   = (P - PK)*QRT
         FD  = (1.0 - 0.5*(P - PK)/(BK + P))*QRT
      ENDIF
!
      END
!
!----------------------------------------------------------------------*
