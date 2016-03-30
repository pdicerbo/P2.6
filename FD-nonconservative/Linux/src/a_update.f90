!* routine for the update of the analytic solution *!
subroutine analytic_update(xg,Npts,xmax,id_flag,eqtype_flag,bc_flag,         &
     &                 no_gaussians,vv,nu,tt,xpeak,xpeak0,sigma,         &
     &                 xstep01,xstep02,un_anal,wave_no)

  Use nrtype

  Implicit none

  integer(I4B)                    :: Npts,id_flag,eqtype_flag,bc_flag
  integer(I4B)                    :: no_gaussians
  integer(I4B)                    :: nx,i
  real(dp), dimension(1:Npts)     :: un_anal,xg
  real(dp)                        :: xstep1,xstep2,xstep01,xstep02
  real(dp)                        :: tt,vv,nu,wave_no,xmax,tmp
  real(dp)                        :: xpeak,xpeak0,sigma

  Intent(IN)  :: xg,Npts,xmax,id_flag,eqtype_flag,tt,                    &
       &         xpeak0,sigma,wave_no
  Intent(OUT) :: un_anal

  if (bc_flag.eq.1) then !* sommerfeld bcs *!
     no_gaussians = 0
  end if

  !* set the shock speed to 0.5 in the case of a burgers *!
  !* equation with square wave                           *!
  if (id_flag.eq.2.and.eqtype_flag.eq.3) then
     vv = 0.5d0
  end if

  !* update the analytic solution *!
  if (eqtype_flag.eq.1.or.eqtype_flag.eq.3) then
     !*--------------------\ *!
     !* advection equation  >*!
     !*--------------------/ *!
     if (id_flag.eq.1) then 
        !---------------!
        !* gaussian    *!
        !---------------!
        do nx = 1, Npts  
           un_anal(nx) = 0.0d0
           do i = -no_gaussians, no_gaussians
              !* set the peak of the different gaussians *!
              xpeak = xpeak0 + xmax*i
              !* the analytic is just a train of no_gaussians gaussians *!
              un_anal(nx) = un_anal(nx) +                                &
                   &        exp( -((xg(nx) - vv*tt - xpeak)/sigma)**2 )  
           end do
        end do
     else if (id_flag.eq.2) then 
        !---------------!
        !* square wave *!
        !---------------!
        !* advect the edges *!
        xstep1 = xstep01 + vv*tt
        xstep2 = xstep02 + vv*tt

        !* rescale them if out of the box *!
100     if(xstep1.gt.1.0d0.and.bc_flag.ne.1) then
           xstep1 = xstep1 - xmax
           go to 100
        end if
110     if(xstep2.gt.1.0d0.and.bc_flag.ne.1) then
           xstep2 = xstep2 - xmax
           go to 110 
        end if

        do nx = 1, Npts  
           if (xstep1.lt.xstep2) then
              if   (xg(nx).lt.xstep1.or.xg(nx).gt.xstep2)  then
                 !* set to 1 if between xstep1 and xstep 2 *!
                 un_anal(nx) = 0.0d0
              else 
                 un_anal(nx) = 1.0d0
              end if
           else 
              if   (xg(nx).lt.xstep1.and.xg(nx).gt.xstep2)  then
                 !* set to 1 if between xstep1 and xstep 2 *!
                 un_anal(nx) = 0.0d0
              else 
                 un_anal(nx) = 1.0d0
              end if
           end if
        end do

     else if (id_flag.eq.3) then 
        !---------------!
        !* wave packet *!
        !---------------!
        do nx = 1, Npts  

           un_anal(nx) = 0.0d0
           do i = -no_gaussians, no_gaussians
              !* set the peak of the different gaussians *!
              xpeak = xpeak0 + xmax*i
              !* the analytic is just a train of no_gaussians gaussians *!
              un_anal(nx) = un_anal(nx) +                               &
                   &        exp( -((xg(nx) - vv*tt - xpeak)/sigma)**2 )  
           end do
           un_anal(nx) = un_anal(nx)*sin( wave_no*(xg(nx) - vv*tt) )
        end do
     else if (id_flag.eq.4) then 
        !---------------!
        !* sine wave   *!
        !---------------!
        do nx = 1, Npts  
           tmp = wave_no*(xg(nx) - vv*tt)
           if (xg(nx).gt.tmp.and.bc_flag.eq.1) then
              un_anal(nx) = 0.0d0
           else 
              un_anal(nx) = sin( wave_no*(xg(nx) - vv*tt) )
           end if
        end do
     else 
        !----------------!
        !* sanity check *!
        !----------------!
        write(*,*)'problems with the id_flag',id_flag
        stop
     end if

  else if (eqtype_flag.eq.2) then
     !*--------------------\ *!
     !* wave      equation  >*!
     !*--------------------/ *!
     if (id_flag.eq.1) then 
        !---------------!
        !* gaussian    *!
        !---------------!
        do nx = 1, Npts  
           un_anal(nx) = 0.0d0
           do i = -no_gaussians, no_gaussians
              !* set the peak of the different gaussians *!
              xpeak = xpeak0 + xmax*i
              !* the analytic is just a train of no_gaussians gaussians *!
              un_anal(nx) = un_anal(nx) +                                &
                   &   0.5d0*exp( -((xg(nx) - vv*tt - xpeak)/sigma)**2 ) &  
                   &         +                                           & 
                   &   0.5d0*exp( -((xg(nx) + vv*tt - xpeak)/sigma)**2 )  
           end do
        end do
     else if (id_flag.eq.2) then 
        !---------------!
        !* square wave *!
        !---------------!
        ! to be implemented... !
     else if (id_flag.eq.3) then 
        !---------------!
        !* wave packet *!
        !---------------!
        ! to be implemented... !
     else if (id_flag.eq.4) then 
        !---------------!
        !* sine wave   *!
        !---------------!
        do nx = 1, Npts  
           un_anal(nx) = 0.5d0*sin( wave_no*(xg(nx) - vv*tt) ) +         &
                &        0.5d0*sin( wave_no*(xg(nx) + vv*tt) )
        end do
     end if
  end if

end subroutine analytic_update
