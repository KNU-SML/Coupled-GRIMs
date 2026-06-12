#include <define.h>
   subroutine chgr_divor2wind(di,ze,uu,vv,e)             
!-------------------------------------------------------------------------------
!
! subprogram:    chgr_divor2wind      computes wind from vorticity, divergence. 
!                                                                               
! abstract: calculates   spectral representation of cosine-weighted             
!   wind components from spectral representation of vorticity and               
!   divergence.                                                                 
!                                                                               
! program history log:                                                          
!   1988-05-13  joseph sela                             
!   2000-01-01  song-you hong          usgs topo, kagwd options
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call chgr_divor2wind (di, ze, uu, vv, e)   
!   input argument list:                                                        
!     di       - triangular spectral divergence at one level.                   
!     ze       - triangular spectral vorticity  at one level.                   
!     e        - array of constants computed in sph_poly_epsilon1.  
!                used to calculate cosine-weighted latitudinal derivatives   
!                of stream function or velocity potential.                      
!                                                                               
!   output argument list:                                                       
!     uu       - triangular spectral zonal pseudo-wind at one level.            
!     vv       - triangular spectral meridional pseudo-wind                     
!                at one level.                                                  
!
!-------------------------------------------------------------------------------
   use paramodel, only : jcap_,jcap1_,jcap2_,lnt2_,lnut2_,rerth_
   use comchgr, only   : d, z, u, v
   save
!-------------------------------------------------------------------------------
   real  ::  e(jcap1_,jcap2_)                                                
   real  ::  di(lnt2_),ze(lnt2_),uu(lnut2_),vv(lnut2_)                       
!-------------------------------------------------------------------------------
!       complex cil                                                             
!
   rjcap1=jcap1_                                                             
!                                                                               
   call chgr_sph_transp(di,d,-1)
   call chgr_sph_transp(ze,z,-1)
!  l=0                                                                     
   xn = 0.0e0                                                               
!  u(x,1,1) =  e(2,1)*z(y,2,1)                                             
   u(1,1,1) =  e(1,2)*z(1,1,2)                                             
   u(2,1,1) =  e(1,2)*z(2,1,2)                                             
!  v(x,1,1) = -e(2,1)*d(y,2,1)                                             
   v(1,1,1) = -e(1,2)*d(1,1,2)                                             
   v(2,1,1) = -e(1,2)*d(2,1,2)                                             
   do i = 2,jcap_                                                            
     xn = xn + 1.0e0                                                          
     r = 1.0e0/xn                                                             
     r1 = 1.0e0/(xn+1.0e0)                                                   
!    u(x,i,1)=-e(i,1)*z(y,i-1,1)*r+e(i+1,1)*z(y,i+1,1)*r1                      
     u(1,1,i)=-e(1,i)*z(1,1,i-1)*r+e(1,i+1)*z(1,1,i+1)*r1                      
     u(2,1,i)=-e(1,i)*z(2,1,i-1)*r+e(1,i+1)*z(2,1,i+1)*r1                      
!    v(x,i,1)= e(i,1)*d(y,i-1,1)*r-e(i+1,1)*d(y,i+1,1)*r1                      
     v(1,1,i)= e(1,i)*d(1,1,i-1)*r-e(1,i+1)*d(1,1,i+1)*r1                      
     v(2,1,i)= e(1,i)*d(2,1,i-1)*r-e(1,i+1)*d(2,1,i+1)*r1                      
   enddo
   xn = xn + 1.0e0                                                          
   r = 1.0e0/xn                                                             
!  u(x,jcap1_,1) = -e(jcap1_,1)*z(y,jcap_,1)*r                             
   u(1,1,jcap1_) = -e(1,jcap1_)*z(1,1,jcap_)*r                             
   u(2,1,jcap1_) = -e(1,jcap1_)*z(2,1,jcap_)*r                             
!  v(x,jcap1_,1) =  e(jcap1_,1)*d(y,jcap_,1)*r                             
   v(1,1,jcap1_) =  e(1,jcap1_)*d(1,1,jcap_)*r                             
   v(2,1,jcap1_) =  e(1,jcap1_)*d(2,1,jcap_)*r                             
   xn = xn + 1.0e0                                                          
   r = 1.0e0/xn                                                             
!  u(x,jcap2_,1) = -e(jcap2_,1)*z(y,jcap1_)*r                              
   u(1,1,jcap2_) = -e(1,jcap2_)*z(1,1,jcap1_)*r                            
   u(2,1,jcap2_) = -e(1,jcap2_)*z(2,1,jcap1_)*r                            
!  v(x,jcap2_,1) =  e(jcap2_,1)*d(y,jcap1_)*r                              
   v(1,1,jcap2_) =  e(1,jcap2_)*d(1,1,jcap1_)*r                            
   v(2,1,jcap2_) =  e(1,jcap2_)*d(2,1,jcap1_)*r                            
!  fin l=0                                                               
!
   do 1000 l = 2,jcap1_                                                        
     xl = l                                                                    
     xn = xl-1.0e0                                                            
     xll = xn                                                                  
     r = 1.0e0/xn                                                             
     r1 = 1.0e0/xl                                                            
!    cil=cmplx(0.,1.)*(l-1)                                                  
!    u(x,1,l)=-cil*d(y,1,l)*r*r1 + e(2,l)*z(y,2,l)*r1                        
     u(1,l,1)= d(2,l,1)*r*r1*xll+e(l,2)*z(1,l,2)*r1                          
     u(2,l,1)=-d(1,l,1)*r*r1*xll+e(l,2)*z(2,l,2)*r1                          
!    v(x,1,l)=-cil*z(y,1,l)*r*r1 - e(2,l)*d(y,2,l)*r1                        
     v(1,l,1)= z(2,l,1)*r*r1*xll-e(l,2)*d(1,l,2)*r1                          
     v(2,l,1)=-z(1,l,1)*r*r1*xll-e(l,2)*d(2,l,2)*r1                          
!
     do 2 i=2,jcap_                                                            
       xn = xn+1.0e0                                                            
!
       if(xn.gt.rjcap1) go to 2                        
!
       r = 1.0e0/xn                                                             
       r1 = 1.0e0/(xn+1.0e0)                                                   
       r2 = r*r1                                      
       rl = r2*xll                                   
!      u(x,i,l)=-e(i,l)*z(y,i-1,l)*r-cil*d(y,i,l)*r2                           
!                   + e(i+1,l)*z(y,i+1,l)*r1                                 
       u(1,l,i)=-e(l,i)*z(1,l,i-1)*r+d(2,l,i)*rl+e(l,i+1)*z(1,l,i+1)*r1
       u(2,l,i)=-e(l,i)*z(2,l,i-1)*r-d(1,l,i)*rl+e(l,i+1)*z(2,l,i+1)*r1
!      v(x,i,l)= e(i,l)*d(y,i-1,l)*r-cil*z(y,i,l)*r2                           
!                   - e(i+1,l)*d(y,i+1,l)*r1                                 
       v(1,l,i)= e(l,i)*d(1,l,i-1)*r+z(2,l,i)*rl-e(l,i+1)*d(1,l,i+1)*r1
       v(2,l,i)= e(l,i)*d(2,l,i-1)*r-z(1,l,i)*rl-e(l,i+1)*d(2,l,i+1)*r1
     2 continue                                                     
!
     xn = xn +1.0e0                                                 
!
     if(xn.gt.rjcap1) go to 1000                                   
!
     r = 1.0e0/xn                                                             
     r1 = r/(xn+1.0e0)                                                        
!    u(x,jcap1_,l)=-e(jcap1_,l)*z(y,jcap_,l)*r-cil*d(y,jcap1_,l)*r1          
     u(1,l,jcap1_)=-e(l,jcap1_)*z(1,l,jcap_)*r+d(2,l,jcap1_)*r1*xll            
     u(2,l,jcap1_)=-e(l,jcap1_)*z(2,l,jcap_)*r-d(1,l,jcap1_)*r1*xll            
!    v(x,jcap1_,l)= e(jcap1_,l)*d(y,jcap_,l)*r-cil*z(y,jcap1_,l)*r1          
     v(1,l,jcap1_)= e(l,jcap1_)*d(1,l,jcap_)*r+z(2,l,jcap1_)*r1*xll            
     v(2,l,jcap1_)= e(l,jcap1_)*d(2,l,jcap_)*r-z(1,l,jcap1_)*r1*xll            
     xn = xn + 1.0e0                                                          
     r = 1.0e0/xn                                                             
!    u(x,jcap2_,l) = -e(jcap2_,l)*z(y,jcap1_,l)*r                            
     u(1,l,jcap2_) = -e(l,jcap2_)*z(1,l,jcap1_)*r                            
     u(2,l,jcap2_) = -e(l,jcap2_)*z(2,l,jcap1_)*r                            
!    v(x,jcap2_,l) =  e(jcap2_,l)*d(y,jcap1_,l)*r                            
     v(1,l,jcap2_) =  e(l,jcap2_)*d(1,l,jcap1_)*r                            
     v(2,l,jcap2_) =  e(l,jcap2_)*d(2,l,jcap1_)*r                            
   1000 continue                                                  
!....                                                                           
   call chgr_sph_transp_wind(uu,u,1)
   call chgr_sph_transp_wind(vv,v,1)
   do j = 1,lnut2_                                                           
     uu(j)=uu(j)*rerth_                                                        
     vv(j)=vv(j)*rerth_                                                        
   enddo
!
   return                                                                    
   end                                                                       
!-------------------------------------------------------------------------------
