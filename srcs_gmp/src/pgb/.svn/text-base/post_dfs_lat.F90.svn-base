!
   subroutine post_dfs_lat(k,a,w)
!-------------------------------------------------------------------------------
   use constant, only : pi_
!-------------------------------------------------------------------------------
   save
   real a(k),w(k)
!-------------------------------------------------------------------------------
   kh=k/2
   dlt=pi_/k
!
   do j = 1,kh
     a(j)=cos((j-0.5)*dlt)
     a(k+1-j)=-a(j)
   enddo
!
   sindlt=2.*sin(dlt*0.5)
#ifdef ORO
   w(1)=1.-cos(dlt*0.5)
   w(k)=w(1)
!
   do j = 2,kh
#else
!
   do j = 1,kh
#endif
     w(j)=sin((j-0.5)*dlt)*sindlt
     w(k+1-j)=w(j)
   enddo
!
   if(k.ne.kh*2) then
     a(kh+1)=0.
     w(kh+1)=sindlt*0.5
   endif
!
   return
   end 
