#include <define.h>
   subroutine post_equall_lat(jh,slat,clat,wlat)
!-------------------------------------------------------------------------------
!
! subprogram documentation block                                            
!                                                                               
! subprogram:    post_equall_lat   compute equally-spaced latitude functions    
!                                                                               
! abstract: computes sines and cosines and gaussian weights                     
!           of equally-spaced latitudes from pole to equator.                   
!           the weights are computed based on ellsaesser (jam,1966).            
!                                                                               
! program history log:                                                          
!   91-10-31  mark iredell                                                      
!   93-12-28  iredell  modified weights based on ellsaesser                     
!                                                                               
! usage:    call post_equall_lat(jh,slat,clat,wlat)             
!                                                                               
!   input argument list:                                                        
!     jh       - integer number of latitudes in a hemisphere                    
!                                                                               
!   output argument list:                                                       
!     slat     - real (jh) sines of latitude                                    
!     clat     - real (jh) cosines of latitude                                  
!     wlat     - real (jh) gaussian weights                                     
!                                                                               
! subprograms called:                                                           
!   minv         solves full matrix problem                                     
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only : jo_, jjh=>johf_
!-------------------------------------------------------------------------------
   real, parameter      ::  pi=3.14159265358979
   real                 ::  slat(jh),clat(jh),wlat(jh)
   real                 ::  awork(jjh,jjh+1)
#define DEFAULT
#ifdef MINV
#undef DEFAULT
   real                 ::  bwork(jjh*2)
#endif
#ifdef DEFAULT
   integer              ::  iwork(jjh*2)
#endif
!
   dlat=0.5*pi/(jh-1)                                                        
   slat(1)=1.                                                                
   clat(1)=0.                                                                
!
   do j = 2,jh-1                                                               
     slat(j)=cos((j-1)*dlat)                                                 
     clat(j)=sin((j-1)*dlat)                                                 
   enddo                                                                     
!
   slat(jh)=0.                                                               
   clat(jh)=1.                                                               
!
   do js = 1,jh                                                                
     do j = 1,jh                                                               
       awork(js,j)=cos(2*(js-1)*(j-1)*dlat)                                  
     enddo                                                                   
   enddo                                                                     
!
#define DEFAULT
#ifdef MINV
#undef DEFAULT
   do js = 1,jh                                                                
      awork(js,jh+1)=-1./(4*(js-1)**2-1)                                      
   enddo                                                                     
!
   call minv (awork,jh,jh,bwork,da,1.e-12,1,0)                               
!
   do j = 1,jh                                                                 
      wlat(j)=awork(j,jh+1)                                                   
   enddo                                                                     
#endif
#ifdef DEFAULT
   call sph_semi_inverse_matrix(awork,jh,det,iwork(1),iwork(jh+1))                             
!
   do j = 1,jh                                                                 
     wlat(j)=0.                                                              
   enddo                                                                     
!
   do j = 1,jh                                                                 
     do jj = 1,jh                                                              
       wlat(j)=wlat(j)+awork(jj,j)*(-1./(4.*(jj-1)**2-1))                    
     enddo                                                                   
   enddo                                                                     
#endif
!
   return                                                                    
   end subroutine post_equall_lat
!-------------------------------------------------------------------------------
