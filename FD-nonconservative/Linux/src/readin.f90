subroutine readin(params)

  Use nrtype
  Implicit none

  integer(I4B)               :: i,no_params
  real(dp), dimension(20)    :: params

  character(len=25)          :: title
  character(len=25)          :: descr_param(20)
  
  INTENT(OUT)                :: params


  !* Read the input file *!
  open (unit=50, file='inputfile.par',status='old')
  !* read how many params *!
  read (50,240) title, no_params
  write(*,*) title, no_params

  do i = 1, no_params
     !* read the actual params *!
     read (50,280) descr_param(i),params(i)
     write(*,*)descr_param(i),params(i),i
  enddo

240 format(a25,i2.2)
280 format(a25,f10.3)

  return
end subroutine readin
