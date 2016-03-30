SUBROUTINE WriteDataGnuplot
  USE typesDef   
  IMPLICIT NONE 
  INTEGER            :: i,iDim,iErr 
  REAL               :: QN(nVar),VN(nVar)
  REAL               :: x0(d), xGP(d) 
  !
  OPEN(31,file='rho.dat',status='unknown',position='append')
  OPEN(32,file='v.dat',status='unknown',position='append')
  OPEN(33,file='p.dat',status='unknown',position='append')
  WRITE(31,227) time
  WRITE(32,227) time
  WRITE(33,227) time  
  !
  DO i = 1, nElem      
     x0 = x(:,tri(1,i)) 
     QN(:) = uh(:,i)
        !
     xGP = x0 + 0.5*dx(:) 
     CALL PDECons2Prim(VN,QN,iErr)
        !
     WRITE(31,321) xGP(1), VN(1)
     WRITE(32,321) xGP(1), VN(2)
     WRITE(33,321) xGP(1), VN(5)
  ENDDO
  WRITE(31,110);WRITE(32,110);WRITE(33,110)
  CLOSE(31);CLOSE(32);CLOSE(33)  

110 FORMAT(' ')  
227 format('"Time = ',E13.6)  
321 FORMAT(1x,2(E21.12,1x))

END SUBROUTINE WriteDataGnuplot
