subroutine writeout(xg,Npts,tt,nt,unpo,un_anal)

  Use nrtype

  Implicit none

  integer(I4B)                    :: Npts,nt
  integer(I4B)                    :: nx
  real(dp), dimension(1:Npts)     :: unpo,un_anal,xg
  real(dp)                        :: tt,norm

  write(33,120)tt
  write(35,120)tt
  do nx = 1, Npts
     write(33,100)xg(nx),unpo(nx)
     write(35,100)xg(nx),un_anal(nx)
  end do
  write(33,110) !* separation line *!
  write(35,110) !* separation line *!

  !* calculate the norm *!
  norm = 0.0
  do nx = 1, Npts
     norm = norm + (unpo(nx))**2
  end do
  write(34,*)tt,log(norm),norm

  write(*,*)' Taking step:',nt,', time',tt

100 format(1x,2(f13.6,1x))
110 format('')
120 format('"Time =',e13.6)

end subroutine writeout
