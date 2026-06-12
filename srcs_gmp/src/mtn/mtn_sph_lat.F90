#include <define.h>
   subroutine mtn_sph_lat(k,a,w)                                                
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:  mtn_sph_lat   compute equally-spaced latitudes 
!                                                                               
! abstract: computes cosines of colatitude and gaussian weights                 
!   on equally-spaced latitudes.                                                
!                                                                               
! usage:    call mtn_sph_lat(k,a,w)                         
!                                                                               
!   input argument list:                                                        
!     k        - number of latitudes.                                           
!                                                                               
!   output argument list:                                                       
!     a        - real (k) cosines of colatitude.                                
!     w        - real (k) gaussian weights.                                     
!                                                                               
!-------------------------------------------------------------------------------
   save
   real  ::  a(k),w(k)
!
   pi=4.*atan(1.)
   kh=k/2
   dlt=pi/(k-1)
!
   do j = 1,kh
     a(j)=cos((j-1)*dlt)
     a(k+1-j)=-a(j)
   enddo
!
   w(1)=1.-cos(dlt*0.5)
   w(k)=w(1)
   sindlt=2.*sin(dlt*0.5)
#ifdef CRAY_THREAD
!dir$ ivdep
#endif
   do j = 2,kh
     w(j)=sin((j-1)*dlt)*sindlt
     w(k+1-j)=w(j)
   enddo
!
   if(k.ne.kh*2) then
     a(kh+1)=0.
     w(kh+1)=sindlt*0.5
   endif
!
   return
   end subroutine mtn_sph_lat
!-------------------------------------------------------------------------------
