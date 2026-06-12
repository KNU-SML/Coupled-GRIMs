#include <define.h>
   subroutine chgr_sum_coeff (fln,ap,qln,n)                                   
!-------------------------------------------------------------------------------
!
! subprogram: chgr_sum_coeff   produces scalar fourier coefficients.              
!                                                                               
! abstract: uses spherical harmonics of a scalar to produce the                 
!   fourier coefficients of the scalar on the northern                          
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
! usage:    call chgr_sum_coeff (fln,ap,qln,n)      
!   input argument list:                                                       
!     fln      - array of spherical harmonics.                                  
!     qln      - array of associated legendre functions.                        
!     n        - number of fields in fln and ap arrays.                         
!                                                                               
!   output argument list:                                                       
!     ap       - array of fourier representations of fields on                  
!                northern and southern hemisphere latitude circles.             
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only : jcap1_,twoj1_,lnt2_,lonf_, lonf2_
!-------------------------------------------------------------------------------
   save
!-------------------------------------------------------------------------------
   real                 ::  ap(lonf2_,n),                                      &
                            qln(lnt2_), fln(lnt2_,n) 
   real                 ::  s(lnt2_), sev(twoj1_), sod(twoj1_)
!-------------------------------------------------------------------------------
!                                                                              
   npair = (jcap1_-3)/2                                                      
!
   do k = 1,n                                                                
     do i = 1,lnt2_                                                          
       s(i) = qln(i) * fln(i,k)                            
     enddo
     do i = 1,twoj1_                                                         
       sev(i) = s(i)   
     enddo
     len = twoj1_ - 2                                                          
     do i = 1,len                                                            
       sod(i) = s(i+twoj1_)  
     enddo
     sod(twoj1_-1) = 0.0e0                                                    
     sod(twoj1_  ) = 0.0e0                                                    
     iplus = twoj1_*2 - 2                                                      
     len   = twoj1_ - 4                                                        
!                                                                              
     do j = 1,npair                                                          
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
     do i = 1,len   
       sev(i) = sev(i) + s(i+iplus) 
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
