program Driver

  !-----------------------------------------------------------------------!
  !
  ! This routine solves the equation u_t = -+ u_x using a number of 
  ! different algorithm. The initial data consists of a Gaussian, a 
  ! square wave or a wave packet 
  !
  !-----------------------------------------------------------------------!

  Use nrtype
  Use InterfaceModule
  Use a_update
  Use i_data

  Implicit none

  real(dp), pointer, dimension(:) :: xg 
  real(dp), pointer, dimension(:) :: un,unpo,unmo,un_anal
  real(dp), pointer, dimension(:) :: rn,rnpo,rnmo
  real(dp), pointer, dimension(:) :: sn,snpo,snmo
  real(dp), pointer, dimension(:) :: unrhs,rnrhs,snrhs
  real(dp), dimension(20)         :: params
  real(dp)                        :: xstep1,xstep2
  real(dp)                        :: T_max,xmax
  real(dp)                        :: tt,t0,x0,norm
  real(dp)                        :: Dt,Dx
  real(dp)                        :: cfl,vv,nu,wave_no
  real(dp)                        :: xpeak,xpeak0,no_periods
  real(dp), parameter             :: sigma      = 0.1d0

  integer(I4B)                    :: Npts,nt_stps,nx,nt,i,j,k
  integer(I4B)                    :: dump_frqncy
  integer(I4B)                    :: id_fl,update_fl,bc_fl,eqtype_fl
  integer(I4B)                    :: no_gaussians

  character(len=18)               :: fname_u
  character(len=21)               :: fname_norm
  character(len=19)               :: fname_ua


  call readin(params)

  ! setting the various parameters
  Npts        = params(1)
  no_periods  = params(2)
  cfl         = params(3)
  dump_frqncy = params(4)
  vv          = params(5)
  nu          = params(6)
  wave_no     = params(7)*2.0d0*PI
  xstep1      = params(8)
  xstep2      = params(9)
  id_fl       = params(10)
  update_fl   = params(11)
  bc_fl       = params(12)
  eqtype_fl   = params(13)

  allocate( xg(     1:Npts) )
  allocate( unpo(   1:Npts) )
  allocate( un(     1:Npts) )
  allocate( unmo(   1:Npts) )
  allocate( un_anal(1:Npts) )
  allocate( rn(     1:Npts) )
  allocate( rnmo(   1:Npts) )
  allocate( rnpo(   1:Npts) )
  allocate( sn(     1:Npts) )
  allocate( snmo(   1:Npts) )
  allocate( snpo(   1:Npts) )
  allocate( unrhs(  1:Npts) )
  allocate( rnrhs(  1:Npts) )
  allocate( snrhs(  1:Npts) )

  xmax         = 1.0d0  
  x0           = 0.0d0
  Dx           = xmax/(Npts - 1)
  Dt           = cfl*Dx
  T_max        = no_periods*xmax*abs(vv)
  no_gaussians = int(T_max) + 1
  nt_stps      = T_max/Dt       !* note that this is an integer *!
  Dt           = T_max/nt_stps  !* modify Dt to stop code at T_max *!
  xpeak        = xmax/2.0d0
  xpeak0       = xpeak
  nt           = 0
  tt           = 0.0d0
  dump_frqncy  = min(dump_frqncy,nt_stps/20)
  if (dump_frqncy.eq.0) dump_frqncy = 1

  write(*,*)'' 
  write(*,*)'T_max=',real(T_max),'Nt_steps=',real(nt_stps),              &
       &     'DT=',real(Dt)
  write(*,*)'Effective CFL=',Dt/Dx
  write(*,*)'' 


  !* Initial data *!
  norm = 0.0

  call initialdata(xg,x0,xmax,Dx,Dt,Npts,xpeak0,xpeak,sigma,wave_no,     &
       &           no_gaussians,vv,nu,un,unpo,unmo,un_anal,rn,sn,        &
       &           xstep1,xstep2,id_fl,update_fl,eqtype_fl,bc_fl,        &
       &           fname_u,fname_norm,fname_ua)

  open(unit=33,file=fname_u   ,status='unknown')
  open(unit=34,file=fname_norm,status='unknown')
  open(unit=35,file=fname_ua  ,status='unknown')

  !* dump the initial data *!
  call writeout(xg,Npts,tt,nt,unpo,un_anal)

  !*----------------------*!
  !* start the time steps *!
  !*----------------------*!
  do nt = 1, nt_stps  

     tt = t0 + Dt*nt  !* increment time *!

     if (update_fl.eq.0) then 
        !*------\   *! 
        !* FTCS  >  *!
        !*------/   *! 
        if (eqtype_fl.eq.1) then 
           !* advection equation *!
           call ftcs_update(unpo,un,unmo,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
        else if (eqtype_fl.eq.2) then 
           !* wave      equation *!
           call ftcs_update(rnpo,rn,sn,  Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
           call ftcs_update(snpo,sn,rn,  Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
           unpo(:) = un(:) + Dt*sn(:)
        end if
     else if (update_fl.eq.1) then 
        !*------\   *! 
        !* UPWD  >  *!
        !*------/   *! 
        if (eqtype_fl.eq.1) then 
           !* advection equation *!
           call upwd_update(unpo,un,unmo,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
        else if (eqtype_fl.eq.2) then 
           write(*,*)'Upwind is an ill-posed method for solving a wave eqn'
           stop
        else if (eqtype_fl.eq.3) then 
           !* burgers   equation *!
           call upwd_update(unpo,un,unmo,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
        end if
     else if (update_fl.eq.2) then 
        !*------\   *! 
        !* LAXF  >  *!
        !*------/   *! 
        if (eqtype_fl.eq.1) then 
           call laxf_update(unpo,un,unmo,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
        else if (eqtype_fl.eq.2) then 
           call laxf_update(rnpo,rn,sn,  Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
           call laxf_update(snpo,sn,rn,  Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
           unpo(:) = un(:) + 0.5d0*Dt*( snpo(:) + sn(:) )
        else if (eqtype_fl.eq.3) then 
           call laxf_update(unpo,un,unmo,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
        end if
     else if (update_fl.eq.3) then 
        !*------\   *! 
        !* LAXW  >  *!
        !*------/   *! 
        if (eqtype_fl.eq.1) then 
           call laxw_update(unpo,un,unmo,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
        else if (eqtype_fl.eq.2) then 
           call laxw_update(rnpo,rn,sn,  Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
           call laxw_update(snpo,sn,rn,  Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
           unpo(:) = un(:) + 0.5d0*Dt*( snpo(:) + sn(:) )
        end if
     else if (update_fl.eq.4) then 
        !*------\   *! 
        !* LEAP  >  *!
        !*------/   *! 
        !* Note that there is no difference here between      *!
        !* the advection, the wave and the burgers equations  *!
        call leap_update(unpo,un,unmo,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu) 
     else if (update_fl.eq.5) then 
        !*------\   *! 
        !* ITCN  >  *!
        !*------/   *! 
        if (eqtype_fl.eq.1) then 
           !======================================== 1st step =======!
           ! update u to get (1)^{\tilde u}^{n+1}                    !
           unrhs(:) = un(:)                                          !
           call itcn_update(unpo,un,unrhs,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu)
           ! average now to get (1)^{\bar u}^{n+1/2}                 !
           unrhs(:) = 0.5d0*( unpo(:) + un(:) )                      !
           !======================================== 1st step =======!

           !======================================== 1st iteration ==!
           ! update u to get (2)^{\tilde u}^{n+1}                    !
           call itcn_update(unpo,un,unrhs,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu)
           ! average now to get (2)^{\bar u}^{n+1/2}                 !
           unrhs(:) = 0.5d0*( unpo(:) + un(:) )                      !
           !======================================== 1st iteration ==!

           !======================================== 2nd iteration ==!
           ! update u to get u^{n+1} !
           call itcn_update(unpo,un,unrhs,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu)
           !======================================== 2nd iteration ==!
        else if (eqtype_fl.eq.2) then 

            snrhs(:) = sn(:)
            rnrhs(:) = rn(:)
            call itcn_update(rnpo,rn,snrhs,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu)
            call itcn_update(snpo,sn,rnrhs,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu)

            snrhs(:) = 0.5d0*( snpo(:) + sn(:) )
            rnrhs(:) = 0.5d0*( rnpo(:) + rn(:) )
            call itcn_update(rnpo,rn,snrhs,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu)
            call itcn_update(snpo,sn,rnrhs,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu)

            snrhs(:) = 0.5d0*( snpo(:) + sn(:) )
            rnrhs(:) = 0.5d0*( rnpo(:) + rn(:) )
            call itcn_update(rnpo,rn,snrhs,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu)
            call itcn_update(snpo,sn,rnrhs,Dt,Dx,Npts,bc_fl,eqtype_fl,vv,nu)

            unpo(:) = un(:) + 0.5d0*Dt*( snpo(:) + sn(:) )
        end if
        else 
           !*------\   *! 
           !* check >  *!
           !*------/   *! 
           write(*,*)'problems with the update_fl',update_fl
           stop
        end if

        !* dump data at dump_frqncy *!
        if ( mod(nt,dump_frqncy) == 0 ) then 

           !* update the analytic *!
           call analytic_update(xg,Npts,xmax,id_fl,eqtype_fl,bc_fl,           &
                &           no_gaussians,vv,nu,tt,xpeak,xpeak0,sigma,     &
                &           xstep1,xstep2,un_anal,wave_no)

           call writeout(xg,Npts,tt,nt,unpo,un_anal)
        end if

        !* update the solution *!
        unmo (:) = un   (:)
        un   (:) = unpo (:)
        rn   (:) = rnpo (:)
        sn   (:) = snpo (:)     

     end do

     nt = nt - 1
     !* A last dump of the data *!
     call analytic_update(xg,Npts,xmax,id_fl,eqtype_fl,bc_fl,                 &
          &           no_gaussians,vv,nu,tt,xpeak,xpeak0,sigma,           &
          &           xstep1,xstep2,un_anal,wave_no)
     call writeout(xg,Npts,tt,nt,unpo,un_anal)

     deallocate( xg     )
     deallocate( unpo   )
     deallocate( un     )
     deallocate( unmo   )
     deallocate( un_anal)
     deallocate( rn     )
     deallocate( rnmo   )
     deallocate( rnpo   )
     deallocate( sn     )
     deallocate( snmo   )
     deallocate( snpo   )
     deallocate( unrhs  )
     deallocate( rnrhs  )
     deallocate( snrhs  )

     write(*,*)''
     write(*,*)' So long! '

     close(33)
     close(34)
     close(35)

   end program Driver


