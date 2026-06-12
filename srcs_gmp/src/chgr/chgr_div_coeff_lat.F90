#include <define.h>
   subroutine chgr_div_coeff_lat(am,ap,bm,bp,fln,qlnwcs,qder,n)     
!-------------------------------------------------------------------------------
!
! subprogram: chgr_div_coeff_lat          
!                                                                               
! abstract: computes negative of divergence operator for                        
!   the terms a and b from eq 27 (sela 1980, 1982).                             
!   input arrays am, ap, bm, bp consist of combinations                         
!   of northern and southern hemisphere fourier coefficients                    
!   previously computed in symns.                                               
!                                                                               
! usage:    call chgr_div_coeff_lat(am,ap,bm,bp,fln,qlnwcs,qder,n)   
!   input argument list:                                                        
!     am       - array of multipliers of antisymmetric pln.                     
!     ap       - array of multipliers of     symmetric pln.                     
!     bm       - array of multipliers of antisymmetric pln.                     
!     bp       - array of multipliers of     symmetric pln.                     
!     fln      - array of coefficients which already have                       
!                contributions from other latitudes.                            
!     qlnwcs   - array of legendres times weights.                              
!     qder     - array of derivatives of legendres.                             
!     n        - number of fields in arrays am, ap, bm, bp, fln.                
!                                                                               
!   output argument list:                                                       
!     fln      - array of coefficients modified by the contribution             
!                at the current latitude.                                       
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only : jcap_,jcap1_,lnt2_,twoj1_,lonf_,lonf2_
!
   save
   real                 ::  am(lonf_,n), ap(lonf2_,n),                         &
                            bm(lonf_,n), bp(lonf2_,n),                         &
                            qlnwcs(lnt2_), qder(lnt2_), fln(lnt2_,n)     
   real                 ::  s(lnt2_), sev(twoj1_), sod(twoj1_)            
!-------------------------------------------------------------------------------
!                                                                              
   npair = (jcap1_-3)/2                                                      
!
   do k = 1,n                                                                
     do i = 1,twoj1_                                                         
       s(i) = bm(i,k) * qder(i)             
     enddo
      len = twoj1_ - 2                                                          
     do i = 1,len                                                            
       s(i+twoj1_) = bp(i,k) * qder(i+twoj1_) 
     enddo
     iplus = twoj1_*2 - 2                                                      
     len   = twoj1_ - 4                                                        
!                                                                              
     do j = 1,npair                                                          
       do i = 1,len                                                            
         s(i+iplus) = bm(i,k) * qder(i+iplus) 
       enddo
       iplus = iplus + len                   
       len = len - 2                        
!                                                                              
       do i = 1,len                                                            
         s(i+iplus) = bp(i,k) * qder(i+iplus) 
       enddo
       iplus = iplus + len                   
       len = len - 2                        
     enddo
!                                                                              
     do i = 1,len                                                            
       s(i+iplus) = bm(i,k) * qder(i+iplus)   
     enddo
!                                                                              
     do i = 1,lnt2_                                                          
       fln(i,k) = fln(i,k) + s(i)            
     enddo
!                                                                              
     do l=1,twoj1_,2                        
       sev(l  ) =  ap(l+1,k)               
       sev(l+1) = -ap(l  ,k)  
     enddo
     len = twoj1_ - 2                                                          
     do l=1,len,2                                                          
       sod(l  ) =  am(l+1,k) 
       sod(l+1) = -am(l  ,k)
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
   enddo
!
   return                                                                    
   end                                                                       
!-------------------------------------------------------------------------------
