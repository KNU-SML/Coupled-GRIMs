#!/bin/csh -f

cat test_2020_yng/fcst.o* | grep "soojin data" > z.txt
cat > fortr.f90 << EOF
program dkjf
  implicit none

  integer :: i,j, tmp
  character(len=100) :: cha
  integer, parameter ::nx=362,ny=292
  real, dimension(nx,ny) :: qnsoce, qnsice, qsrmix, rain, snow, ievp, tevp, qsr, emp

  open(11,file='./z.txt',form='formatted')
  open(51,file='./output.gdat',form='unformatted')
  open(52,file='./output.ctl',form='formatted')

  do j = 1, ny
    do i = 1, nx
      read(11,*) cha,cha,tmp,tmp,qnsoce(i,j), qnsice(i,j), qsrmix(i,j), &
                 rain(i,j), snow(i,j), ievp(i,j), tevp(i,j), emp(i,j), qsr(i,j)
    enddo
  enddo

  write(51) qnsoce
  write(51) qnsice
  write(51) qsrmix
  write(51) rain
  write(51) snow  
  write(51) ievp  
  write(51) tevp  
  write(51) emp   
  write(51) qsr   

  write(52,*) "dset ^output.gdat"
  write(52,*) "undef -9.99e08"
  write(52,*) "title dfj"
  write(52,*) "options sequential"
  write(52,*) "xdef 362 linear 0 0.1"
  write(52,*) "ydef 292 linear 0 0.1"
  write(52,*) "zdef 1 linear 0 1"
  write(52,*) "tdef 1 linear jan2010 1mo"
  write(52,*) "vars 9"
  write(52,*) "qnsoce 0 99 z"
  write(52,*) "qnsice 0 99 z"
  write(52,*) "qsrmix 0 99 z"
  write(52,*) "rain   0 99 z"
  write(52,*) "snow   0 99 z"
  write(52,*) "ievp   0 99 z"
  write(52,*) "tevp   0 99 z"
  write(52,*) "emp    0 99 z"
  write(52,*) "qsr    0 99 z"
  write(52,*) "endvars"
end
EOF

ifort fortr.f90
./a.out

