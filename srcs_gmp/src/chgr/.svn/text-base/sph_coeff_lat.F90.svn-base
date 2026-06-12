#include <define.h>
   subroutine sph_coeff_lat(fp,fm,fln,qln,n)                                    
!-------------------------------------------------------------------------------
!
! subprogram: sph_coeff_lat   accumulates spectral from fourier input.          
!                                                                               
! abstract: accumulates spectral coefficients from fourier input.               
!   call to subroutine accumulates the contribution at the given                
!   latitude only.                                                              
!                                                                               
! program history log:                                                          
!   1988-04-04  joseph sela                                                    
!   2000-01-01  song-you hong          usgs topo, kagwd options
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call sph_coeff_lat(fp, fm, fln, qln, n)                           
!   input argument list:                                                        
!     fp       - array of fourier  coefficients for  symmetric part.            
!     fm       - array of fourier  coefficients for asymmetric part.            
!     fln      - array of spectral coefficients to be summed into               
!                at the current latitude.                                       
!     qln      - array of associated legendre functions.                        
!     n        - number of fields in fp, fm, fln arrays.                        
!                                                                               
!   output argument list:                                                       
!     fln      - array of spectral coefficients with contribution               
!                from the current latitude.                                     
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_,lonf2_,lnt2_,jcap1_,twoj1_
!-------------------------------------------------------------------------------
   save
!-------------------------------------------------------------------------------
   real  ::  fp(lonf2_,n), fm(lonf_,n),                                        &
             qln(lnt2_), fln(lnt2_,n)                                                 
   real  ::  s(lnt2_)                                                  
!-------------------------------------------------------------------------------
!
   npair = (jcap1_-3)/2                                                      
!
   do k = 1,n                                                                
     do i = 1,twoj1_                                                        
       s(i) = fp(i,k) * qln(i)                                          
     enddo
     len = twoj1_ - 2                                                         
     do i = 1,len                                                           
       s(i+twoj1_) = fm(i,k) * qln(i+twoj1_)                            
     enddo
     iplus = twoj1_*2 - 2                                                     
     len   = twoj1_ - 4                                                      
     do j = 1,npair                                                         
       do i = 1,len                                                       
         s(i+iplus) = fp(i,k) * qln(i+iplus)                          
       enddo
       iplus = iplus + len                                          
       len = len - 2                                               
       do i = 1,len                                                 
         s(i+iplus) = fm(i,k) * qln(i+iplus)                   
       enddo
       iplus = iplus + len                                     
       len = len - 2                                          
     enddo
     do i = 1,len                                              
       s(i+iplus) = fp(i,k) * qln(i+iplus)                 
     enddo
     do i = 1,lnt2_                                         
       fln(i,k) = fln(i,k) + s(i)                       
     enddo
   enddo
!
   return                                                                    
   end                                                                       
!-------------------------------------------------------------------------------
