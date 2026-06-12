program dkjf
  implicit none

  integer :: i,j, tmp,t
  character(len=100) :: cha
  integer, parameter ::nx=362,ny=292,nt=2
  real :: dum
  real, dimension(nx,ny) :: qsrmix, palbi, alb, ziceld,picefr,zqsr,qsr_ice,qsr_tot
  open(11,file='./iceflx.txt',form='formatted')
!  open(12,file='./qsrmix.txt',form='formatted')
  open(51,file='./output3.gdat',form='unformatted')
  open(52,file='./output3.ctl',form='formatted')

  do t = 1, nt
    do j = 1, ny
      do i = 1, nx
        read(11,*) cha,cha,tmp,tmp,qsr_tot(i,j),qsr_ice(i,j),&
                   palbi(i,j),alb(i,j),ziceld(i,j),picefr(i,j)
!      read(12,*) qsrmix(i,j) 

!      zqsr(i,j)=qsrmix(i,j)*(1.-palbi(i,j))/(1.-(alb(i,j)*ziceld(i,j)+palbi(i,j)*picefr(i,j)))

      enddo
    enddo
    write(51) qsr_ice(:,:) 
    write(51) qsr_tot(:,:) 
  enddo

  write(52,*) "dset ^output3.gdat"
  write(52,*) "undef -9.99e08"
  write(52,*) "title dfj"
  write(52,*) "options sequential"
  write(52,*) "xdef 362 linear 0 0.1"
  write(52,*) "ydef 292 linear 0 0.1"
  write(52,*) "zdef 1 linear 0 1"
  write(52,*) "tdef 2 linear jan2010 1mo"
  write(52,*) "vars 2"
  write(52,*) "qsr_ice 0 99 z"
  write(52,*) "qsr_tot 0 99 z"
  write(52,*) "endvars"
end
