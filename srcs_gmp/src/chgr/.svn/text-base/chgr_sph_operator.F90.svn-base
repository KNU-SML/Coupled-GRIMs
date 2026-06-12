#include <define.h>
   subroutine chgr_sph_operator(am,ap,bm,bp,fln,qlnwcs,qder,n)     
!-------------------------------------------------------------------------------
!
! subprogram: chgr_sph_operator                             
!                                                                               
! abstract: computes curl operator for                                          
!   the terms aln and bln from eq 28 (sela 1980, 1982).                         
!   input arrays am, ap, bm, bp consist of combinations                         
!   of northern and southern hemisphere fourier coefficients                    
!   (at a given latitude) which have been                                       
!   previously combined in subroutine symns.                                    
!   subroutine psu22 is analogous to subroutine msu22.                          
!                                                                               
! program history log:                                                          
!   1988-04-11  joseph sela 
!   2000-01-01  song-you hong          cvs
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call chgr_sph_operator(am,ap,bm,bp,fln,qlnwcs,qder,n) 
!   input argument list:                                                        
!     am       - array of multipliers of antisymmetric pln.                     
!                difference output from symns for a.                            
!     ap       - array of multipliers of     symmetric pln.                     
!                sum        output from symns for a.                            
!     bm       - array of multipliers of antisymmetric pln.                     
!                difference output from symns for b.                            
!     bp       - array of multipliers of     symmetric pln.                     
!                sum        output from symns for b.                            
!     fln      - array of coefficients which already have                       
!                contributions from other latitudes.                            
!                each call to psu22 adds a contribution to fln                  
!                for the current latitude.                                      
!     qlnwcs   - array of legendres times gaussian weights.                     
!     qder     - array of derivatives of legendres.                             
!     n        - number of fields in arrays am, ap, bm, bp, fln.                
!                                                                               
!   output argument list:                                                       
!     fln      - array of coefficients modified by the contribution             
!                at the current latitude.                                       
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_,lonf2_,twoj1_,lnt2_,jcap1_
!
   save
   real  ::   am(lonf_,n), ap(lonf2_,n),                                       &
              bm(lonf_,n), bp(lonf2_,n),                                       &
              qlnwcs(lnt2_), qder(lnt2_), fln(lnt2_,n) 
   real  ::   s(lnt2_), sev(twoj1_), sod(twoj1_)                        
!-------------------------------------------------------------------------------
   npair = (jcap1_-3)/2                                                      
!
   do k = 1,n                                                                
!
     do i = 1,twoj1_                                                         
       s(i) = am(i,k) * qder(i) 
     enddo
     len = twoj1_ - 2                                                          
     do i = 1,len                                                            
       s(i+twoj1_) = ap(i,k) * qder(i+twoj1_) 
     enddo
     iplus = twoj1_*2 - 2  
     len   = twoj1_ - 4                                                        
!                                                                              
     do j = 1,npair  
       do i = 1,len                                                            
         s(i+iplus) = am(i,k) * qder(i+iplus) 
       enddo
       iplus = iplus + len
       len = len - 2 
!                                                                              
       do i = 1,len                                                            
         s(i+iplus) = ap(i,k) * qder(i+iplus)
       enddo
       iplus = iplus + len
       len = len - 2 
     enddo
!                                                                              
     do i = 1,len                                                            
       s(i+iplus) = am(i,k) * qder(i+iplus)
     enddo
!                                                                              
     do i = 1,lnt2_                                                          
       fln(i,k) = fln(i,k) + s(i)
     enddo
!                                                                              
     do l=1,twoj1_,2                                                       
       sev(l  ) = -bp(l+1,k)
       sev(l+1) =  bp(l  ,k)
     enddo
     len = twoj1_ - 2                                                          
     do l=1,len,2                                                          
       sod(l  ) = -bm(l+1,k) 
       sod(l+1) =  bm(l  ,k)
     enddo
!                                                                              
     do i = 1,twoj1_                                                         
       s(i) = sev(i) * qlnwcs(i)
     enddo
     len = twoj1_ - 2                                                          
     do i = 1,len                                                            
       s(i+twoj1_) = sod(i) * qlnwcs(i+twoj1_)
     enddo
     iplus = twoj1_*2 - 2                                                      
     len   = twoj1_ - 4                                                        
!                                                                              
     do j = 1,npair                                                          
       do i = 1,len                                                            
         s(i+iplus) = sev(i) * qlnwcs(i+iplus)
       enddo
       iplus = iplus + len
       len = len - 2
!                                                                              
       do i = 1,len                                                            
         s(i+iplus) = sod(i) * qlnwcs(i+iplus)
       enddo
       iplus = iplus + len
       len = len - 2
     enddo
!                                                                              
     do i = 1,len                                                            
       s(i+iplus) = sev(i) * qlnwcs(i+iplus) 
     enddo
!                                                                              
     do i = 1,lnt2_
       fln(i,k) = fln(i,k) + s(i) 
     enddo
!
   enddo
!
   return                                                                    
   end                                                                       
!-------------------------------------------------------------------------------
