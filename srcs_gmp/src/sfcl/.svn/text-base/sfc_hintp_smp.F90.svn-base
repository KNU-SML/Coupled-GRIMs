!
   subroutine sfc_hintp_smp(im,jm,ylin,grid,xlong,ylat,gval)
!-------------------------------------------------------------------------------
   real                 ::  grid(im,jm)
   real                 ::  gval
   real                 ::  ylin(jm)
   real                 ::  dum
!-------------------------------------------------------------------------------
   dxin = 360./float(im)
!
   alamd = xlong
   i1=alamd/dxin+1.001
   iindx1 = i1
   i2=i1+1
   if(i2.gt.im) i2=1
   iindx2 = i2
   dum =(i2-1.001)*dxin
!   print*,'i,i1,i2',i,i1,i2,alamd,dum
   ddx = (alamd-float(i1-1)*dxin)/dxin
!
   j2=1
   aphi = ylat
   do 50 jj=1,jm
     if (aphi.lt.ylin(jj)) go to 50
     j2=jj
     !
     go to 42
     !
50 continue
42 continue
   if(j2.gt.2) go to 43
   j1=1
!
   j2=2
   !
   go to 44
   !
43 continue
   if(j2.le.jm) go to 45
   j1=jm-1
   j2=jm
   !
   go to 44
   !
45 continue
   j1=j2-1
44 continue
!  print *,'j,j1,j2',j,j1,j2
   jindx1 = j1
   jindx2 = j2
   ddy = (aphi-ylin(j1))/(ylin(j2)-ylin(j1))
!
   y = ddy
   j1 = jindx1
   j2 = jindx2
   x = ddx
   i1 = iindx1
   i2 = iindx2
   gval = (1.-x)*(1.-y)*grid(i1,j1)+(1.-y)*x*grid(i2,j1)+                      &
              (1.-x)*y*grid(i1,j2)+x*y*grid(i2,j2)
!
!    print*,'hintp'
!    write(*, 9) i1,j1,i2,j2,grid(i1,j1),grid(i2,j1),grid(i1,j2),grid(i2,j2)
!
9  format(4i5,4e13.5)
!
   return
   end subroutine sfc_hintp_smp
!-------------------------------------------------------------------------------
