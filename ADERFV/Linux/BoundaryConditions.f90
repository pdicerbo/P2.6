SUBROUTINE BoundaryConditions 
  USE typesDef
  ! Local variables 
  REAL :: j,k,iFace
  REAL :: Qbc(nVar),Fbc(nVar,d),Vbc(nVar)
  !
  ! Fix boundary data  
  !Vbc = (/ 1., 0., 0., 0., 1. /)    ! primitive variables     
  !CALL PDEPrim2Cons(qBC,Vbc)        ! convert into conservative variables    
  !
  DO iFace = 1, nFace
     ! Here, we need to take care of the boundary conditions 
     ! For the moment, we use either simple extrapolation (copy from inside the domain) 
     ! or impose a constant value 
     IF(Face(iFace)%Left.EQ.0) THEN
        Face(iFace)%qL(:) = Face(iFace)%qR(:)
        Face(iFace)%FL(:) = Face(iFace)%FR(:)
     ENDIF
     IF(Face(iFace)%Right.EQ.0) THEN 
        Face(iFace)%qR(:) = Face(iFace)%qL(:)
        Face(iFace)%FR(:) = Face(iFace)%FL(:)
     ENDIF
  ENDDO
  ! 
END SUBROUTINE BoundaryConditions
    
    
