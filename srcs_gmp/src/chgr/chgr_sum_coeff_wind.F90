#include <define.h>
   subroutine chgr_sum_coeff_wind (vln,ap,qln,n)
!-------------------------------------------------------------------------------
!
! subprogram: chgr_sum_coeff_wind    produces vector fourier coefficients. 
!                                                                               
! abstract: uses spherical harmonics of a pseudo vector to produce              
!   the fourier coefficients of the fields on the northern                      
!   and southern hemisphere latitude circles for which                          
!   the input legendre functions apply.                                         
!                                                                               
! program history log:                                                          
!   1988-04-04  joseph sela 
!   2000-01-01  song-you hong          usgs topo, kagwd options
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    chgr_sum_coeff_wind (vln,ap,qln,n)  
!   input argument list:                                                        
!     vln      - array of spherical harmonics.                                  
!     qln      - array of associated legendre functions.                        
!     n        - number of fields in vln and ap arrays.                         
!                                                                               
!   output argument list:                                                       
!     ap       - array of fourier representations of fields on                  
!                northern and southern hemisphere latitude circles.             
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_,lonf2_,lnut2_,jcap2_,twoj1_
   save
!-------------------------------------------------------------------------------
   real                 ::  ap(lonf2_,n),                                      &
                            qln(lnut2_), vln(lnut2_,n) 
   real                 ::  s(lnut2_), sev(twoj1_), sod(twoj1_) 
!-------------------------------------------------------------------------------
!                                                                              
   npair = (jcap2_-2)/2 
!
   do k = 1,n         
     do i = 1,lnut2_
       s(i) = qln(i) * vln(i,k) 
     enddo
     do i = 1,twoj1_            
       sev(i) = s(i)          
       sod(i) = s(i+twoj1_)  
     enddo
     iplus = twoj1_*2        
     len   = twoj1_-2       
!                                                                              
     do j=1,npair          
       do i = 1,len      
         sev(i) = sev(i) + s(i+iplus) 
       enddo
       iplus = iplus + len            
       len = len - 2                 
!                                                                              
       do i = 1,len                 
         sod(i) = sod(i) + s(i+iplus) 
       enddo
       iplus = iplus + len            
       len = len - 2                 
     enddo
!                                                                              
     do i = 1,twoj1_                                                         
       ap(i      ,k) = sev(i) + sod(i)
     enddo
     do i = 1,twoj1_                                                         
       ap(i+lonf_,k) = sev(i) - sod(i)
     enddo
   enddo
!
   return                                                                    
   end                                                                       
!-------------------------------------------------------------------------------
