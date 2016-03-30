Module a_update

  Interface
     subroutine analytic_update(xg,Npts,xmax,id_fl,eqtype_fl,bc_fl,           &
          &                 no_gaussians,vv,nu,tt,xpeak,xpeak0,sigma,     &
          &                 xstep1,xstep2,un_anal,wave_no)

       Use nrtype

       Implicit none

       integer(I4B)                    :: Npts,id_fl,eqtype_fl,bc_fl
       integer(I4B)                    :: no_gaussians
       integer(I4B)                    :: nx,i
       real(dp), dimension(1:Npts)     :: un_anal,xg
       real(dp)                        :: xstep1,xstep2
       real(dp)                        :: tt,vv,nu,wave_no,xmax,tmp
       real(dp)                        :: xpeak,xpeak0,sigma

       Intent(IN)  :: xg,Npts,xmax,id_fl,eqtype_fl,nu,tt,                &
            &         xpeak0,sigma,wave_no
       Intent(OUT) :: un_anal

     end subroutine analytic_update
  end Interface

end Module a_update



Module InterfaceModule 


  Interface
     subroutine writeout(xg,nx_pnts,tt,nt,unpo,un_anal)

       Use nrtype

       Implicit none

       integer(I4B)                    :: nx_pnts,nt
       integer(I4B)                    :: nx
       real(dp), dimension(1:nx_pnts)  :: unpo,un_anal,xg
       real(dp)                        :: tt,norm
     end subroutine writeout
  end Interface

  Interface 
     subroutine readin(params)

       Use nrtype
       Implicit none

       integer(I4B)               :: i,no_params
       real(dp), dimension(20)    :: params

       character(len=25)          :: title
       character(len=25)          :: descr_param(20)

       INTENT(OUT)                :: params


     end subroutine readin
  end Interface

  Interface 
     subroutine ftcs_update(xnpo,xn,xnmo,dt,h,JG,bc_fl,eqtype_fl,vv,nu)

       Use nrtype
       Implicit none

       integer(I4B)                    :: JG,bc_fl,eqtype_fl,j

       real(dp), dimension(JG)         :: xnpo,xn,xnmo
       real(dp)                        :: dt,h,vv,nu,dt_o_dx

     end subroutine ftcs_update
  end Interface

  Interface 
     subroutine upwd_update(xnpo,xn,xnmo,dt,h,JG,bc_fl,eqtype_fl,vv,nu)

       Use nrtype
       Implicit none

       integer(I4B)                    :: JG,bc_fl,eqtype_fl,j

       real(dp), dimension(JG)         :: xnpo,xn,xnmo
       real(dp)                        :: dt,h,vv,nu,dt_o_dx

     end subroutine upwd_update
  end Interface

  Interface 
     subroutine laxf_update(xnpo,xn,xnmo,dt,h,JG,bc_fl,eqtype_fl,vv,nu)

       Use nrtype
       Implicit none

       integer(I4B)                    :: JG,bc_fl,eqtype_fl,j

       real(dp), dimension(JG)         :: xnpo,xn,xnmo
       real(dp)                        :: dt,h,vv,nu,dt_o_dx

     end subroutine laxf_update
  end Interface

  Interface 
     subroutine laxw_update(xnpo,xn,xnmo,dt,h,JG,bc_fl,eqtype_fl,vv,nu)

       Use nrtype
       Implicit none

       integer(I4B)                    :: JG,bc_fl,eqtype_fl,j

       real(dp), dimension(JG)         :: xnpo,xn,xnmo
       real(dp)                        :: dt,h,vv,nu,dt_o_dx

     end subroutine laxw_update
  end Interface

  Interface
     subroutine leap_update(xnpo,xn,xnmo,dt,h,JG,bc_fl,eqtype_fl,vv,nu)

       Use nrtype
       Implicit none

       integer(I4B)                    :: JG,bc_fl,eqtype_fl,j

       real(dp), dimension(JG)         :: xnpo,xn,xnmo
       real(dp)                        :: dt,h,vv,nu,dt_o_dx

     end subroutine leap_update
  end Interface

  Interface
     subroutine itcn_update(xnpo,xn,xnmo,dt,h,JG,bc_fl,eqtype_fl,vv,nu)

       Use nrtype
       Implicit none

       integer(I4B)                    :: JG,bc_fl,eqtype_fl,j

       real(dp), dimension(JG)         :: xnpo,xn,xnmo
       real(dp)                        :: dt,h,vv,nu,dt_o_dx
     end subroutine itcn_update
  end Interface

end Module InterfaceModule

Module i_data

  Interface
     subroutine initialdata(xg,x0,xmax,Dx,Dt,Npts,xpeak0,xpeak,sigma,wave_no, &
          &               no_gaussians,vv,nu,un,unpo,unmo,un_anal,rn,sn,      &
          &               xstep1,xstep2,                                      &
          &               id_flag,update_flag,eqtype_flag,bc_flag,            &
          &               fname_u,fname_norm,fname_ua)

       Use nrtype
       Use a_update
       Implicit none

       integer(I4B)                    :: Npts,nx
       integer(I4B)                    :: id_flag,update_flag,eqtype_flag,   &
            &                             bc_flag
       integer(I4B)                    :: no_gaussians

       real(dp), dimension(Npts)       :: xg 
       real(dp), dimension(Npts)       :: un,unpo,unmo,un_anal,rn,sn
       real(dp)                        :: xstep1,xstep2
       real(dp)                        :: Dx,Dt,wave_no,x0,xmax,vv,nu
       real(dp)                        :: xpeak0,xpeak,sigma,tt

       character(len=14)               :: fname_u
       character(len=17)               :: fname_norm
       character(len=15)               :: fname_ua

       Intent(IN)  :: x0,xmax,Dx,Dt,Npts,xpeak0,xpeak,sigma,wave_no,        &
            &         no_gaussians,vv,nu,xstep1,xstep2,id_flag,update_flag, &
            &         eqtype_flag

       Intent(OUT) :: xg,un,unpo,unmo,un_anal,rn,sn,fname_u,fname_norm,fname_ua

     end subroutine initialdata
  end Interface
  
end Module i_data
