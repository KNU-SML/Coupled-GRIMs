#include <define.h>
   subroutine chgr_sph_funct (qlnt, qlnv, colrad, lat)
!-------------------------------------------------------------------------------
!
! subprogram:    chgr_sph_funct       evaluates associated legendre functions. 
!                                                                               
! abstract: evaluates the required values of the normalized                     
!   associated legendre function at a prescribed colatitude.                    
!   a standard recursion relation is used with real arithmetic.                 
!                                                                               
! program history log:                                                          
!   1988-10-25  joseph sela    
!   2000-01-01  song-you hong          cvs
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call chgr_sph_funct (qlnt, qlnv, colrad, lat)
!   input argument list:                                                        
!     colrad   - half precision colatitudes in radians for which                
!                the associated legendre functions are to be                    
!                computed.                                                      
!     lat      - index which indicates the current latitude.                    
!                                                                               
!   output argument list:                                                       
!     qlnt     - doubled scalar triangle of                                     
!                half precision associated legendre functions.                  
!     qlnv     - doubled vector triangle of                                     
!                half precision associated legendre functions.                  
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only : lnt2_,lnut2_,latg2_,jcap1_,jcap2_,twoj1_,lnut_
   use comchgr, only   : deps, rdeps, dx, y
!-------------------------------------------------------------------------------
   save
!                                                                               
   real   ::   qlnt(lnt2_)                                             
   real   ::   qlnv(lnut2_)                                            
   real   ::   colrad(latg2_)                                            
!                                                                               
   real   ::   x(jcap1_)                                            
   real   ::   dpln(lnut2_)                                            
   real   ::   colr                                                    
   real   ::   sinlat                                                    
   real   ::   cos2                                                    
   real   ::   prod                                                    
!                                                                              
   data ifir /0/                                                             
!-------------------------------------------------------------------------------
!
   if  (ifir .eq. 1)  go to 500                                              
!
     ifir = 1                                                             
     do ll = 1,jcap1_                                                        
       rdeps(ll) = 0.0                                                    
     enddo
     lplus = jcap1_                                                            
     len   = jcap1_                                                            
     do inde=2,jcap2_                                                      
       do ll = 1,len                                                         
         l = ll - 1                                                         
         n = l + inde - 1                                                   
         rdeps(ll+lplus) = (n*n - l*l) / (4.0 * n*n - 1.0)                  
       enddo
       lplus = lplus + len 
       len = len - 1
     enddo
     do i = jcap2_,lnut_                                                     
       rdeps(i) = sqrt(rdeps(i))                                          
     enddo
     do i = 1,lnut_                                                          
       deps(2*i-1) = rdeps(i)                                             
       deps(2*i  ) = rdeps(i)                                             
     enddo
     ibegin = twoj1_ + 1                                                       
     do i = ibegin,lnut2_                                                    
       rdeps(i) = 1.0/deps(i)                                             
     enddo
     do ll = 1,jcap1_                                                        
       x(ll) = ll*2+1                                                     
     enddo
     do ll = 1,jcap1_                                                        
       y(ll) = x(ll)/(x(ll)-1.)                                           
     enddo
     do ll = 1,jcap1_                                                        
       x(ll) = sqrt(x(ll))                                                
     enddo
     do ll = 1,jcap1_                                                        
       dx(2*ll-1) = x(ll)                                                 
       dx(2*ll  ) = x(ll)                                                 
     enddo
!
   500 continue
!                                                                              
   colr   = colrad(lat)
   sinlat = cos(colr)                                                        
   cos2   = 1.0 - sinlat * sinlat                                            
   prod   = 1.0                                                              
!
   do ll = 1,jcap1_                                                        
     x(ll) = 0.5*prod                                                   
!    if (prod .lt. 1.0e-75)  prod=0.0                                   
     prod = prod*cos2*y(ll)                                             
   enddo
!
   do ll = 1,jcap1_                                                        
     x(ll) = sqrt(x(ll))                                                
   enddo
!
   do ll = 1,jcap1_                                                        
     dpln(2*ll-1) = x(ll)                                               
     dpln(2*ll  ) = x(ll)                                               
   enddo
!
   lplus = twoj1_                                                            
   do ll = 1,twoj1_                                                        
     dpln(ll+lplus) = dx(ll) * sinlat * dpln(ll)                        
   enddo
!
   lp2 = 0                                                                   
   lp1 =     twoj1_                                                          
   lp0 = 2 * twoj1_                                                          
   len =     twoj1_ - 2                                                      
!
   do n = 3,jcap2_                                                        
!dir$ ivdep                                                                     
     do ll = 1,len                                                           
       dpln(ll+lp0) = (sinlat * dpln(ll+lp1)                                   &
                      - deps(ll+lp1) * dpln(ll+lp2)) * rdeps(ll+lp0)
     enddo
     lp2 = lp1                                                                 
     lp1 = lp0                                                                 
     lp0 = lp0 + len                                                           
     len = len - 2                                                             
   enddo
!                                                                              
   do i = 1,lnut2_                                                         
     qlnv(i) = dpln(i)                                                  
   enddo
!                                                                              
   len = 2 * twoj1_ - 2                                                      
   do ll = 1,len                                                           
     qlnt(ll) = qlnv(ll)                                                
   enddo
!
   lpv = 2 * twoj1_                                                          
   lpt = lpv - 2                                                             
   len = twoj1_ - 4                                                          
   do n=3,jcap1_                                                         
     do ll = 1,len                                                           
       qlnt(ll+lpt) = qlnv(ll+lpv)                                        
     enddo
     lpv = lpv + len + 2                                                       
     lpt = lpt + len                                                           
     len = len - 2                                                             
   enddo
!                                                                              
   return                                                                    
   end                                                                       
!-------------------------------------------------------------------------------

