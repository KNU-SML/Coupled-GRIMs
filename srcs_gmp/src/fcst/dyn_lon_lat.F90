#include "define.h"
   subroutine dyn_lon_lat(xlon,xlat,colrad,lon,lat)                             
!-------------------------------------------------------------------------------
   use constant, only   :   pi_
!-------------------------------------------------------------------------------
!
! abstract : get normal lon,lat in radians   
!
!-------------------------------------------------------------------------------
   integer, intent(in)  ::  lon,lat
   real, intent(out)    ::  xlon(lon,lat),xlat(lon,lat)
   real, intent(in)     ::  colrad(lat/2)                    
   real,parameter       ::  tpi=2.e0*pi_,hpi=0.5e0*pi_
!-------------------------------------------------------------------------------
   bphi = tpi / lon                                                          
   do j = 1,lat
     do i = 1,lon                                                              
       xlon(i,j) = (i-1) * bphi                                                
     enddo
   enddo
   do j = 1,lat/2
     do i = 1,lon
       xlat(i,j) = hpi - colrad(j)                                             
     enddo
#ifdef DBG
     write(6,'(A,I4,F10.5)')'hoon:j,xlat=',j,xlat(1,j)*180/pi_
#endif
   enddo
   do j=lat/2+1,lat
     jj=lat+1-j
     do i = 1,lon
       xlat(i,j) = - xlat(i,jj)
     enddo
   enddo
   return                                                                    
   end subroutine dyn_lon_lat 
