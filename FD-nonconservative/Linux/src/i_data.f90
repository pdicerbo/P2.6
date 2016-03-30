subroutine initialdata(xg,x0,xmax,Dx,Dt,Npts,xpeak0,xpeak,sigma,wave_no,     &
       &               no_gaussians,vv,nu,un,unpo,unmo,un_anal,rn,sn,        &
       &               xstep1,xstep2,id_flag,update_flag,eqtype_flag,bc_flag,&
       &               fname_u,fname_norm,fname_ua)

  Use nrtype
  Use a_update
  Implicit none

  integer(I4B)                    :: Npts,nx
  integer(I4B)                    :: id_flag,update_flag,eqtype_flag,       &
       &                             bc_flag
  integer(I4B)                    :: no_gaussians

  real(dp), dimension(Npts)       :: xg 
  real(dp), dimension(Npts)       :: un,unpo,unmo,un_anal,rn,sn
  real(dp)                        :: xstep1,xstep2
  real(dp)                        :: Dx,Dt,wave_no,x0,xmax,vv,nu
  real(dp)                        :: xpeak0,xpeak,sigma,tt

  character(len=18)               :: fname_u
  character(len=21)               :: fname_norm
  character(len=19)               :: fname_ua

  Intent(IN)  :: x0,xmax,Dx,Dt,Npts,xpeak0,xpeak,sigma,wave_no,        &
       &         no_gaussians,vv,nu,xstep1,xstep2,id_flag,update_flag, &
       &         eqtype_flag

  Intent(OUT) :: xg,un,unpo,unmo,un_anal,rn,sn,fname_u,fname_norm,fname_ua

  if (update_flag.eq.0) then
     if      (eqtype_flag.eq.1) then
        write(fname_u,   300)'ftcs_adv'
        write(fname_norm,310)'ftcs_adv'
     else if (eqtype_flag.eq.2) then
        write(fname_u,   300)'ftcs_wav'
        write(fname_norm,310)'ftcs_wav'
     else if (eqtype_flag.eq.3) then
        write(fname_u,   300)'ftcs_brg'
        write(fname_norm,310)'ftcs_brg'
     end if
  else if (update_flag.eq.1) then
     if      (eqtype_flag.eq.1) then
        write(fname_u,   300)'upwd_adv'
        write(fname_norm,310)'upwd_adv'
     else if (eqtype_flag.eq.2) then
        write(fname_u,   300)'upwd_wav'
        write(fname_norm,310)'upwd_wav'
     else if (eqtype_flag.eq.3) then
        write(fname_u,   300)'upwd_brg'
        write(fname_norm,310)'upwd_brg'
     end if
  else if (update_flag.eq.2) then
     if      (eqtype_flag.eq.1) then
        write(fname_u,   300)'laxf_adv'
        write(fname_norm,310)'laxf_adv'
     else if (eqtype_flag.eq.2) then
        write(fname_u,   300)'laxf_wav'
        write(fname_norm,310)'laxf_wav'
     else if (eqtype_flag.eq.3) then
        write(fname_u,   300)'laxf_brg'
        write(fname_norm,310)'laxf_brg'
     end if
  else if (update_flag.eq.3) then
     if      (eqtype_flag.eq.1) then
        write(fname_u,   300)'laxw_adv'
        write(fname_norm,310)'laxw_adv'
     else if (eqtype_flag.eq.2) then
        write(fname_u,   300)'laxw_wav'
        write(fname_norm,310)'laxw_wav'
     else if (eqtype_flag.eq.3) then
        write(fname_u,   300)'laxw_brg'
        write(fname_norm,310)'laxw_brg'
     end if
  else if (update_flag.eq.4) then
     if      (eqtype_flag.eq.1) then
        write(fname_u,   300)'leap_adv'
        write(fname_norm,310)'leap_adv'
     else if (eqtype_flag.eq.2) then
        write(fname_u,   300)'leap_wav'
        write(fname_norm,310)'leap_wav'
     else if (eqtype_flag.eq.3) then
        write(fname_u,   300)'leap_brg'
        write(fname_norm,310)'leap_brg'
     end if
  else if (update_flag.eq.5) then
     if      (eqtype_flag.eq.1) then
        write(fname_u,   300)'itcn_adv'
        write(fname_norm,310)'itcn_adv'
     else if (eqtype_flag.eq.2) then
        write(fname_u,   300)'itcn_wav'
        write(fname_norm,310)'itcn_wav'
     else if (eqtype_flag.eq.3) then
        write(fname_u,   300)'itcn_brg'
        write(fname_norm,310)'itcn_brg'
     end if
  else 
     !* sanity check *!
     write(*,*)'problems with the update_flag',update_flag
     stop
  end if

  if      (eqtype_flag.eq.1) then
     write(fname_ua,  320)'adv'
  else if (eqtype_flag.eq.2) then
     write(fname_ua,  320)'wav'
  else if (eqtype_flag.eq.3) then
     write(fname_ua,  320)'brg'
  end if

  !* set the initial data *!

  do nx = 1, Npts  
     xg(nx)  = x0 + Dx*(nx-1)

     if (id_flag.eq.1) then 
        !* gaussian    *!
        un(nx) = exp( -((xg(nx) - xpeak)/sigma)**2 )  
        rn(nx) = -2.0d0*(xg(nx) - xpeak)*                                &
             &   exp(-( (xg(nx) - xpeak)/sigma)**2)/sigma**2
        sn(nx) = 0.0d0

     else if (id_flag.eq.2) then 
        !* square wave *!
        !* these are the extrema of the steps; the first component *!
        !* is the initial one.
        if (xg(nx).lt.xstep1.or.xg(nx).gt.xstep2) then
           un(nx) = 0.0d0
        else 
           un(nx) = 1.0d0
        end if
        rn(nx) = 0.0d0
        sn(nx) = 0.0d0

     else if (id_flag.eq.3) then 
        !* wave packet *!
        un(nx) = exp( -((xg(nx) - xpeak)/sigma)**2 )*                    &
             &   sin(wave_no*xg(nx))  
        rn(nx) = exp( -((xg(nx) - xpeak)/sigma)**2 )*(                   &
             &  wave_no*cos(wave_no*xg(nx)) -                            &
             &    2.0d0*sin(wave_no*xg(nx))*( (xg(nx) - xpeak)/sigma )**2)
        sn(nx) = 0.0d0
     else if (id_flag.eq.4) then 
        !* sine wave *!
        un(nx) = sin(wave_no*xg(nx))  
        rn(nx) = cos(wave_no*xg(nx))*wave_no  
        sn(nx) = 0.0d0
     else 
        !* sanity check *!
        write(*,*)'problems with the id_flag',id_flag
        stop
     end if
  end do

  if (  (update_flag.eq.4.and.eqtype_flag.eq.1).or.                        &
       & update_flag.eq.4.and.eqtype_flag.eq.3      ) then 
     !* Leapfrog needs special care for the -1 time-level *!
     !* in the case of the advection equation; in this    *!
     !* it is convenient to go back in time with the      *!
     !* analytic  solution                                *!
     tt = - Dt  
     call analytic_update(xg,Npts,xmax,id_flag,eqtype_flag,bc_flag,           &
          &           no_gaussians,vv,nu,tt,xpeak,xpeak0,sigma,           &
          &           xstep1,xstep2,unmo,wave_no)
     write(*,*)'porc!'
  else 
     unmo(:)    = un(:)
  end if

  un_anal(:) = un(:)
  unpo(:)    = un(:)


300 format( 'out/u_'   ,g8.8,'.dat' )
310 format( 'out/norm_',g8.8,'.dat' )
320 format( 'out/ua_'  ,g3.3,'.dat' )

end subroutine initialdata
