!
   subroutine chgr_comp_maxmin(f,imax,jmax,kmax,title)
!-------------------------------------------------------------------------------
   real               ::  f(imax,jmax,kmax)
   character(len=8)   ::  title
!-------------------------------------------------------------------------------
   print 99, title
   99 format(2x,'title=',a8)
!
   do k = 1,kmax
!
     fmax=f(1,1,k)
     iimax=1
     jjmax=1
     fmin=f(1,1,k)
     iimin=1
     jjmin=1
!
     do j = 1,jmax
       do i = 1,imax
         if(fmax.lt.f(i,j,k)) then
           fmax=f(i,j,k)
           iimax=i
           jjmax=j
         endif
         if(fmin.gt.f(i,j,k)) then
           fmin=f(i,j,k)
           iimin=i
           jjmin=j
         endif
       enddo
     enddo
!
     print 100, k,fmax,iimax,jjmax,fmin,iimin,jjmin
100  format(2x,'level=',i3,' max=',e12.4,' at i=',i5,' j=',i5,                 &
             ' min=',e12.4,' at i=',i5,' j=',i5)
!
   enddo
!
   return
   end subroutine chgr_comp_maxmin
!-------------------------------------------------------------------------------
