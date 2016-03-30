!*---------------------------------------------------------------*!
! Here are all the different time updates: ftcd, upwind, lax,... *!
!*---------------------------------------------------------------*!

!*------\
!* FTCS  >
!*------/
subroutine ftcs_update(xnpo,xn,yn,Dt,Dx,JG,bc_flag,eqtype_flag,vv,nu)

  Use nrtype
  Implicit none

  integer(I4B)                    :: JG,bc_flag,eqtype_flag,j

  real(dp), dimension(JG)         :: xnpo,xn,yn
  real(dp)                        :: Dt,Dx,vv,nu,Dt_o_Dx

  Dt_o_Dx = vv*Dt/(2.0d0*Dx)

  !* select the equation type: advection or wave *!
  if (eqtype_flag.eq.1) then
     !* advection equation *! 
     do j = 2, JG - 1
        xnpo(j) = xn(j) - Dt_o_Dx*( xn(j+1) - xn(j-1) )
     enddo
  else  if (eqtype_flag.eq.2) then
     !* wave      equation *! 
     do j = 2, JG-1
        xnpo(j) = xn(j) + Dt_o_Dx*( yn(j+1) - yn(j-1) )
     enddo
  else 
     !* sanity check *!
     write(*,*)'problems with the eqtype_flag',bc_flag
     stop
  end if


  !*--------------------\
  !* boundary conditions >
  !*--------------------/
  if (bc_flag.eq.1) then
     !* sommerfeld *!
     xnpo(1)  = xn(2)    - (Dx - Dt)/(Dx + Dt)*(xnpo(2) - xn(1)     )
     xnpo(JG) = xn(JG-1) + (Dx - Dt)/(Dx + Dt)*(xn(JG)  - xnpo(JG-1))
  else if (bc_flag.eq.2) then
     !* reflect *!
     xnpo(1)  = 0.0d0
     xnpo(JG) = 0.0d0
  else if (bc_flag.eq.3) then
     !* periodic: xn(0) = xn(JG); xn(JG+1) = xn(1) *!
     xnpo(1)  = xn(1)  - Dt_o_Dx*( xn(2) - xn(JG)   )
     xnpo(JG) = xn(JG) - Dt_o_Dx*( xn(1) - xn(JG-1) )
  else 
     !* sanity check *!
     write(*,*)'problems with the bc_flag',bc_flag
     stop
  end if

  return

end subroutine ftcs_update

!*------\
!* UPWD  >
!*------/
subroutine upwd_update(xnpo,xn,yn,Dt,Dx,JG,bc_flag,eqtype_flag,vv,nu)

  Use nrtype
  Implicit none

  integer(I4B)                    :: JG,bc_flag,eqtype_flag,j

  real(dp), dimension(JG)         :: xnpo,xn,yn
  real(dp)                        :: Dt,Dx,vv,nu,Dt_o_Dx,Dt_o_Dx2

  Dt_o_Dx  = vv*Dt/Dx
  Dt_o_Dx2 = vv*Dt/Dx/Dx

  !* select the equation type: advection or wave *!
  if (eqtype_flag.eq.1) then
     !* advection equation *! 
     !* need to distinguish according to the sign of vv *!
     if (vv.gt.0.0d0) then 
        do j = 2, JG - 1
           xnpo(j) = xn(j) - Dt_o_Dx*( xn(j)   - xn(j-1) )
        enddo
     else 
        do j = 2, JG - 1
           xnpo(j) = xn(j) - Dt_o_Dx*( xn(j+1) - xn(j)  )
        enddo
     end if
  else  if (eqtype_flag.eq.2) then
     !* wave      equation *! 
     write(*,*)'Upwind is an ill-posed method for solving a wave eqn'
     stop
  else if (eqtype_flag.eq.3) then
     !* burgers   equation *! 
     !* need to distinguish according to the sign of vv *!
     
!!$     ! Conservative scheme
!!$     do j = 2, JG - 1
!!$        xnpo(j) = xn(j) - 0.5d0*Dt/Dx*( xn(j)**2 - xn(j-1)**2 )
!!$     enddo

     ! Non-conservative
     do j = 2, JG - 1
        xnpo(j) = xn(j) - Dt/Dx*xn(j)*( xn(j) - xn(j-1) )
     enddo
     
!!$     if (abs(vv).gt.0.0d0) then 
!!$        do j = 2, JG - 1
!!$           xnpo(j) = xn(j) - 0.5d0*Dt_o_Dx*( xn(j)**2 - xn(j-1)**2 ) +   &
!!$                &        nu*Dt_o_Dx2*( xn(j+1) - 2.0d0*xn(j) + xn(j-1) )
!!$        enddo
!!$     else 
!!$        do j = 2, JG - 1
!!$           xnpo(j) = xn(j) - 0.5d0*Dt_o_Dx*( xn(j+1)**2 - xn(j)**2 ) +   &
!!$                &        nu*Dt_o_Dx2*( xn(j+1) - 2.0d0*xn(j) + xn(j-1) )
!!$        enddo
!!$     end if
  else 
     !* sanity check *!
     write(*,*)'problems with the eqtype_flag',bc_flag
     stop
  end if

  !*--------------------\
  !* boundary conditions >
  !*--------------------/
  if (bc_flag.eq.1) then
     !* sommerfeld *!
     xnpo(1)  = xn(2)    - (Dx - Dt)/(Dx + Dt)*(xnpo(2) - xn(1)     )
     xnpo(JG) = xn(JG-1) + (Dx - Dt)/(Dx + Dt)*(xn(JG)  - xnpo(JG-1))
  else if (bc_flag.eq.2) then
     !* reflect *!
     xnpo(1)  = 0.0d0
     xnpo(JG) = 0.0d0
  else if (bc_flag.eq.3) then
     if (eqtype_flag.eq.1) then
        !* periodic: xn(0) = xn(JG); xn(JG+1) = xn(1) *!
        if (abs(vv).gt.0.0d0) then 
           xnpo(1)  = xn(1)  - Dt_o_Dx*( xn(1)  - xn(JG)   )
           xnpo(JG) = xn(JG) - Dt_o_Dx*( xn(JG) - xn(JG-1) )
        else 
           xnpo(1)  = xn(1)  + Dt_o_Dx*( xn(2)  - xn(1) )
           xnpo(JG) = xn(JG) + Dt_o_Dx*( xn(JG) - xn(1) )
        end if
     else if (eqtype_flag.eq.3) then
        !* burgers   equation *! 
        !* need to distinguish according to the sign of vv *!
        if (abs(vv).gt.0.0d0) then 
           xnpo(1) = xn(1)   - 0.5d0*Dt_o_Dx*( xn(1)**2 - xn(JG)**2 ) +  &
                &        nu*Dt_o_Dx2*( xn(2) - 2.0d0*xn(1) + xn(JG) )
           xnpo(JG) = xn(JG) - 0.5d0*Dt_o_Dx*( xn(JG)**2 - xn(JG-1)**2 )+&
                &        nu*Dt_o_Dx2*( xn(1) - 2.0d0*xn(JG) + xn(JG-1) )
        else if (abs(vv).lt.0.0d0) then 
           write(*,*)'not implemented yet';stop
        end if
     end if
  else 
     !* sanity check *!
     write(*,*)'problems with the bc_flag',bc_flag
     stop
  end if

  return

end subroutine upwd_update

!*------\
!* LAXF  >
!*------/
subroutine laxf_update(xnpo,xn,yn,Dt,Dx,JG,bc_flag,eqtype_flag,vv,nu)

  Use nrtype
  Implicit none

  integer(I4B)                    :: JG,bc_flag,eqtype_flag,j

  real(dp), dimension(JG)         :: xnpo,xn,yn
  real(dp)                        :: Dt,Dx,vv,nu,Dt_o_Dx,Dt_o_Dx2

  Dt_o_Dx  = vv*Dt/(2.0d0*Dx)
  Dt_o_Dx2 = vv*Dt/(Dx**2)

  !* select the equation type: advection or wave *!
  if (eqtype_flag.eq.1) then
     !* advection equation *! 
     do j = 2, JG - 1
        xnpo(j) = 0.5d0*(xn(j+1) + xn(j-1)) - Dt_o_Dx*( xn(j+1) - xn(j-1) )
     enddo
  else  if (eqtype_flag.eq.2) then
     !* wave      equation *! 
     do j = 2, JG - 1
        xnpo(j) = 0.5d0*(xn(j+1) + xn(j-1)) + Dt_o_Dx*( yn(j+1) - yn(j-1) )
     enddo
  else  if (eqtype_flag.eq.3) then
     !* burgers   equation *! 
     do j = 2, JG - 1
        xnpo(j) = 0.5d0*(xn(j+1) + xn(j-1)) -                            &
             &    0.5d0*Dt_o_Dx*( xn(j+1)**2 - xn(j-1)**2 ) +            &
             &       nu*Dt_o_Dx2*( xn(j+1) - 2.0d0*xn(j) + xn(j-1) )

        ! Conservative scheme for the Burger equation
!!$
!!$        xnpo ( j) = xn ( j) - vv*Dt/(2.0d0*Dx) * ( xn ( j)**2 - xn ( j-1)**2)
     enddo
  else 
     !* sanity check *!
     write(*,*)'problems with the eqtype_flag',bc_flag
     stop
  end if

  !*--------------------\
  !* boundary conditions >
  !*--------------------/
  if (bc_flag.eq.1) then
     !* sommerfeld *!
     xnpo(1)  = xn(2)    - (Dx - Dt)/(Dx + Dt)*(xnpo(2) - xn(1)     )
     xnpo(JG) = xn(JG-1) + (Dx - Dt)/(Dx + Dt)*(xn(JG)  - xnpo(JG-1))
  else if (bc_flag.eq.2) then
     !* reflect *!
     xnpo(1)  = 0.0d0
     xnpo(JG) = 0.0d0
  else if (bc_flag.eq.3) then
     !* periodic: xn(0) = xn(JG); xn(JG+1) = xn(1) *!
     if (eqtype_flag.eq.1) then
        !* advection equation *! 
        xnpo(1)  = 0.5d0*(xn(2) + xn(JG)  ) - Dt_o_Dx*( xn(2) - xn(JG)   )
        xnpo(JG) = 0.5d0*(xn(1) + xn(JG-1)) - Dt_o_Dx*( xn(1) - xn(JG-1) )
     else if (eqtype_flag.eq.2) then
        !* wave      equation *! 
        xnpo(1)  = 0.5d0*(xn(2) + xn(JG)  ) + Dt_o_Dx*( yn(2) - yn(JG)   )
        xnpo(JG) = 0.5d0*(xn(1) + xn(JG-1)) + Dt_o_Dx*( yn(1) - yn(JG-1) )
     else if (eqtype_flag.eq.3) then
        !* burgers   equation *! 
        xnpo(1)  = 0.5d0*(xn(2) + xn(JG)  ) -                            &
             &     0.5d0*Dt_o_Dx*( xn(2)**2 - xn(JG)**2   )
        xnpo(JG) = 0.5d0*(xn(1) + xn(JG-1)) -                            &
             &     0.5d0*Dt_o_Dx*( xn(1)**2 - xn(JG-1)**2 )
     end if
  else 
     !* sanity check *!
     write(*,*)'problems with the bc_flag',bc_flag
     stop
  end if

  return

end subroutine laxf_update

!*------\
!* LAXW  >
!*------/
subroutine laxw_update(xnpo,xn,yn,Dt,Dx,JG,bc_flag,eqtype_flag,vv,nu)

  Use nrtype
  Implicit none

  integer(I4B)                    :: JG,bc_flag,eqtype_flag,j

  real(dp), dimension(JG)         :: xnpo,xn,yn
  real(dp)                        :: Dt,Dx,vv,nu,Dt_o_Dx

  Dt_o_Dx = vv*Dt/Dx

  !* select the equation type: advection or wave *!
  if (eqtype_flag.eq.1) then
     !* advection equation *! 
     do j = 2, JG-1
        xnpo(j) = xn(j) - Dt_o_Dx*0.5d0*( (xn(j+1) - xn(j-1))            &
             &          - Dt_o_Dx*(xn(j+1) - 2.0d0*xn(j) + xn(j-1)) )
     enddo
  else  if (eqtype_flag.eq.2) then
     !* wave      equation *! 
     do j = 2, JG-1
        xnpo(j) = xn(j) + Dt_o_Dx*0.5d0*( (yn(j+1) - yn(j-1))            &
             &          + Dt_o_Dx*(xn(j+1) - 2.0d0*xn(j) + xn(j-1)) )
     enddo
  else  if (eqtype_flag.eq.3) then
     !* burgers   equation *! 
  else 
     !* sanity check *!
     write(*,*)'problems with the eqtype_flag',bc_flag
     stop
  end if


  !*--------------------\
  !* boundary conditions >
  !*--------------------/
  if (bc_flag.eq.1) then
     !* sommerfeld *!
     xnpo(1)  = xn(2)    - (Dx - Dt)/(Dx + Dt)*(xnpo(2) - xn(1)     )
     xnpo(JG) = xn(JG-1) + (Dx - Dt)/(Dx + Dt)*(xn(JG)  - xnpo(JG-1))
  else if (bc_flag.eq.2) then
     !* reflect *!
     xnpo(1)  = 0.0d0
     xnpo(JG) = 0.0d0
  else if (bc_flag.eq.3) then
     !* periodic: xn(0) = xn(JG); xn(JG+1) = xn(1) *!
     if (eqtype_flag.eq.1) then
        xnpo(1) = xn(1)   - Dt_o_Dx*0.5d0*( (xn(2) - xn(JG))             &
             &            - Dt_o_Dx*(xn(2) - 2.0d0*xn(1)  + xn(JG)) )
        xnpo(JG) = xn(JG) - Dt_o_Dx*0.5d0*( (xn(1) - xn(JG-1))           &
             &            - Dt_o_Dx*(xn(1) - 2.0d0*xn(JG) + xn(JG-1)) )
     else if (eqtype_flag.eq.2) then
        xnpo(1)  = xn(1)  + Dt_o_Dx*0.5d0*( (yn(2) - yn(JG))             &
             &            + Dt_o_Dx*(xn(2) - 2.0d0*xn(1)  + xn(JG)) )
        xnpo(JG) = xn(JG) + Dt_o_Dx*0.5d0*( (yn(1) - yn(JG-1))           &
             &            + Dt_o_Dx*(xn(1) - 2.0d0*xn(JG) + xn(JG-1)) )
     end if
  else 
     !* sanity check *!
     write(*,*)'problems with the bc_flag',bc_flag
     stop
  end if

  return

end subroutine laxw_update

!*------\
!* LEAP  >
!*------/
subroutine leap_update(xnpo,xn,xnmo,Dt,Dx,JG,bc_flag,eqtype_flag,vv,nu)

  Use nrtype
  Implicit none

  integer(I4B)                    :: JG,bc_flag,eqtype_flag,j

  real(dp), dimension(JG)         :: xnpo,xn,xnmo
  real(dp)                        :: Dt,Dx,vv,nu,Dt_o_Dx

  Dt_o_Dx = vv*Dt/Dx

  !* select the equation type: advection or wave *!
  if (eqtype_flag.eq.1) then
     !* advection equation *! 
     do j = 2, JG-1
        xnpo(j) = xnmo(j) - Dt_o_Dx*( xn(j+1) - xn(j-1) ) 
     enddo
  else  if (eqtype_flag.eq.2) then
     !* wave      equation *! 
     do j = 2, JG-1
        xnpo(j) = Dt_o_Dx**2*xn(j+1) + 2.0d0*xn(j)*(1.0d0 - Dt_o_Dx**2) + &
             &    Dt_o_Dx**2*xn(j-1) - xnmo(j)
     enddo
  else  if (eqtype_flag.eq.3) then
     !* burgers   equation *! 
     do j = 2, JG-1
        xnpo(j) = xnmo(j) - 0.5d0*Dt_o_Dx*( xn(j+1)**2 - xn(j-1)**2 )
     enddo
  else 
     !* sanity check *!
     write(*,*)'problems with the eqtype_flag',bc_flag
     stop
  end if

  !*--------------------\
  !* boundary conditions >
  !*--------------------/
  if (bc_flag.eq.1) then
     !* sommerfeld *!
     xnpo(1)  = xn(2)    - (Dx - Dt)/(Dx + Dt)*(xnpo(2) - xn(1)     )
     xnpo(JG) = xn(JG-1) + (Dx - Dt)/(Dx + Dt)*(xn(JG)  - xnpo(JG-1))
  else if (bc_flag.eq.2) then
     !* reflect *!
     xnpo(1)  = 0.0d0
     xnpo(JG) = 0.0d0
  else if (bc_flag.eq.3) then
     if (eqtype_flag.eq.1) then
        !* periodic: xn(0) = xn(JG); xn(JG+1) = xn(1) *!
        xnpo(1)  = xnmo(1)  - Dt_o_Dx*( xn(2) - xn(JG)   )
        xnpo(JG) = xnmo(JG) - Dt_o_Dx*( xn(1) - xn(JG-1) )
     else if (eqtype_flag.eq.2) then
        xnpo(1)  = Dt_o_Dx**2*xn(2) + 2.0d0*xn(1)* (1.0d0 - Dt_o_Dx**2) + &
             &     Dt_o_Dx**2*xn(JG)   - xnmo(1)
        xnpo(JG) = Dt_o_Dx**2*xn(1) + 2.0d0*xn(JG)*(1.0d0 - Dt_o_Dx**2) + &
             &     Dt_o_Dx**2*xn(JG-1) - xnmo(JG)
     else if (eqtype_flag.eq.3) then
        xnpo(1)  = xnmo(1)  - 0.5d0*Dt_o_Dx*( xn(2)**2 - xn(JG)**2   )
        xnpo(JG) = xnmo(JG) - 0.5d0*Dt_o_Dx*( xn(1)**2 - xn(JG-1)**2 )
     end if
  else 
     !* sanity check *!
     write(*,*)'problems with the bc_flag',bc_flag
     stop
  end if

  return

end subroutine leap_update

!*------\
!* ITCN  >
!*------/
subroutine itcn_update(xnpo,xn,yn,Dt,Dx,JG,bc_flag,eqtype_flag,vv,nu)

  Use nrtype
  Implicit none

  integer(I4B)                    :: JG,bc_flag,eqtype_flag,j

  real(dp), dimension(JG)         :: xnpo,xn,yn
  real(dp)                        :: Dt,Dx,vv,nu,Dt_o_Dx

  Dt_o_Dx = vv*Dt/2.0d0/Dx

  !* select the equation type: advection or wave *!
  if (eqtype_flag.eq.1) then
     !* advection equation *! 
     do j = 2, JG-1
        xnpo(j) = xn(j) - Dt_o_Dx*( yn(j+1) - yn(j-1) )
     enddo
  else  if (eqtype_flag.eq.2) then
     !* wave      equation *! 
     do j = 2, JG-1
        xnpo(j) = xn(j) + Dt_o_Dx*( yn(j+1) - yn(j-1) )
     enddo
  else 
     !* sanity check *!
     write(*,*)'problems with the eqtype_flag',bc_flag
     stop
  end if

  !*--------------------\
  !* boundary conditions >
  !*--------------------/
  if (bc_flag.eq.1) then
     !* sommerfeld *!
     xnpo(1)  = xn(2)    - (Dx - Dt)/(Dx + Dt)*(xnpo(2) - xn(1)     )
     xnpo(JG) = xn(JG-1) + (Dx - Dt)/(Dx + Dt)*(xn(JG)  - xnpo(JG-1))
  else if (bc_flag.eq.2) then
     !* reflect *!
     xnpo(1)  = 0.0d0
     xnpo(JG) = 0.0d0
  else if (bc_flag.eq.3) then
     if (eqtype_flag.eq.1) then
        !* periodic: xn(0) = xn(JG); xn(JG+1) = xn(1) *!
        xnpo(1)  = xn(1)  - Dt_o_Dx*( yn(2) - yn(JG)   )
        xnpo(JG) = xn(JG) - Dt_o_Dx*( yn(1) - yn(JG-1) )
     else if (eqtype_flag.eq.2) then
        xnpo(1)  = xn(1)  + Dt_o_Dx*( yn(2) - yn(JG)   )
        xnpo(JG) = xn(JG) + Dt_o_Dx*( yn(1) - yn(JG-1) )
     end if
  else 
     !* sanity check *!
     write(*,*)'problems with the bc_flag',bc_flag
     stop
  end if

  return

end subroutine itcn_update



