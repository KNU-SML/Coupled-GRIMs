#include <define.h>
   subroutine post_sph_lat
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ... 
!
!    [post_sph_lat] * ----- [post_equall_colat] *
!                       |-- [post_equall_lat] *
!
!-------------------------------------------------------------------------------
   end subroutine post_sph_lat
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_equall_colat(k,a,w)
!-------------------------------------------------------------------------------
!
! subprogram:  post_equall_colat   compute equally-spaced latitudes
!
! abstract: computes cosines of colatitude and gaussian weights
!   on equally-spaced latitudes.
!
! program history log:
!    1992-10-31  iredell                development
!    2000-03-09  songyou hong           cvs verion setup
!    2009-10-01  jung-eun kim           f90 format with standard physics modules
!    2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call post_equall_colat(k,a,w)
!
!   input argument list:
!     k        - number of latitudes.
!
!   output argument list:
!     a        - real (k) cosines of colatitude.
!     w        - real (k) gaussian weights.
!
!   remarks: fortran 9x extensions are used.
!
!-------------------------------------------------------------------------------
   save
   integer              ::  k
   real                 ::  a(k),w(k)
!-------------------------------------------------------------------------------
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
!
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
   end subroutine post_equall_colat
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
#ifdef DFS
   subroutine post_equall_lat(jmax,slat,clat,ioff)
#else
   subroutine post_equall_lat(jh,slat,clat,wlat)
#endif
!-------------------------------------------------------------------------------
   use constant, only : pi=>pi_
   use paramodel, only : jo_, jjh=>johf_
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:    post_equall_lat  compute equally-spaced latitude functions 
!                                                                               
! abstract: computes sines and cosines and gaussian weights                     
!           of equally-spaced latitudes from pole to equator.                   
!           the weights are computed based on ellsaesser (jam,1966).            
!                                                                               
! program history log:                                                          
!   1991-10-31  mark iredell                                                      
!   1993-12-28  iredell  modified weights based on ellsaesser                     
!   1995-01-01  kanamitsu modified to run on workstations                         
!   2000-03-09  songyou hong           cvs verion setup, prognostic clouds
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
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
#ifdef DFS
   real                ::  slat(jmax),clat(jmax,3)
   integer             ::  ioff
#else
   real                ::  slat(jh),clat(jh),wlat(jh)
   real                ::  awork(jjh,jjh+1)
#define DEFAULT
#ifdef MINV
#undef DEFAULT
   real                ::  bwork(jjh*2)
#endif
#ifdef DEFAULT
   integer             ::  iwork(jjh*2)
#endif
#endif          /* not DFS */
!-------------------------------------------------------------------------------
!
#ifdef DFS
!
! no pole and equator
!  
   if (ioff.eq.1) then
     dlat=pi/jmax
     do j = 1,jmax
       rlat=0.5*pi-(j-0.5)*dlat
       slat(j)=sin(rlat)
       clat(j,1)=cos(rlat)
       clat(j,2)=1./clat(j,1)
       clat(j,3)=clat(j,2)*clat(j,2)
     enddo
   else
     dlat=pi/(jmax-1)                                                        
     do j = 1,jmax
       rlat=0.5*pi-(j-1)*dlat
       slat(j)=sin(rlat)                                                 
       clat(j,1)=cos(rlat)
     enddo                                                                     
   endif
#else
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
   call iminv(awork,jh,det,iwork(1),iwork(jh+1))                             
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
#endif          /* not DFS */
!
   return                                                                    
   end subroutine post_equall_lat
