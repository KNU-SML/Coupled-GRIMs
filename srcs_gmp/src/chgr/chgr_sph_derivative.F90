#include <define.h>
   subroutine chgr_sph_derivative (qlnt,qlnv,qdert,eps,lat,qlnwct,rcs2,wgtl) 
!-------------------------------------------------------------------------------
!
! subprogram:    chgr_sph_derivative      computes derivatives of legendres.  
!                                                                               
! abstract: computes derivatives of associated legendre functions               
!   and for convenience, other required products of                             
!   legendres and factors involving wave number and latitude.                   
!   the resulting arrays are required for the application                       
!   of divergence and curl operators in msu22 and psu22.                        
!                                                                               
! program history log:                                                          
!   1988-11-02  joseph sela
!   2000-01-01  song-you hong          cvs 
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call chgr_sph_derivative (qlnt,qlnv,qdert,eps,lat,qlnwct,rcs2,wgtl)
!   input argument list:                                                        
!     qlnt     - doubled scalar triangular                                      
!                array of associated legendre functions at                      
!                a given latitude.                                              
!                on input, values of qlnt are a subset of qlnv.                 
!     qlnv     - doubled vector triangular                                      
!                array of associated legendre functions at                      
!                a given latitude.                                              
!     eps      - array of function of wave number computed in sph_poly_epsilon1
!                eps is used only during first call to chgr_sph_derivative.     
!     lat      - latitude index.                                                
!     rcs2     - array of constants computed in chgr_sph_gaussian 
!                (1/sin(lat)**2).          
!     wgtl     - weight at gaussian latitude.                                   
!                multiplier of output arrays qlnt, qdert, qlnwct.               
!                                                                               
!   output argument list:                                                       
!     qlnt     - doubled scalar triangular                                      
!                array of qlnt*n*(n+1)*1/a**2   times wgtl.                     
!     qdert    - doubled scalar triangular                                      
!                array of legendre derivatives  times wgtl.                     
!     qlnwct   - doubled scalar triangular                                      
!                array of qlnt*l*rcs2(lat)/a    times wgtl.                     
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only : lnut2_,lnt2_,jcap_,jcap1_,jcap2_,twoj1_,latg2_
   use constant, only  : rerth_
   use comchgr, only   : dxa, dxb, dxc, dxd
!-------------------------------------------------------------------------------
   save
!-------------------------------------------------------------------------------
   real   ::   qlnt(lnt2_)                                            
   real   ::   qlnv(lnut2_)                                           
   real   ::   qdert(lnt2_)                                            
   real   ::   eps(jcap1_,jcap2_)                                    
   real   ::   qlnwct(lnt2_)                                            
   real   ::   rcs2(latg2_)                                           
!                                                                              
   real   ::   dxint(lnt2_)                                            
   real   ::   dx(twoj1_,jcap2_)                                    
   real   ::   deps(twoj1_,jcap2_)                                    
!                                                                              
   data ifir /0/                                                             
!-------------------------------------------------------------------------------
!
   if (ifir .eq. 1)  go to 500                                              
!
     ifir = 1                                                             
     do ll = 1,twoj1_                                                        
       dxint(2*ll-1) = ll                                                 
       dxint(2*ll  ) = ll                                                 
     enddo
     lp = 0                                                                    
     do i = 1,jcap2_                                                        
       do ll = 1,twoj1_
         dx(ll,i) = dxint(ll+lp)                                            
       enddo
       lp = lp + 2
     enddo
     do i = 1,jcap2_                                                        
       do ll = 1,jcap1_
         deps(2*ll-1,i) = eps(ll,i)                                         
         deps(2*ll  ,i) = eps(ll,i)                                         
       enddo
      enddo
      lp1 = twoj1_                                                              
      len = twoj1_ - 2                                                          
      do i = 1,jcap_                                                         
        do ll = 1,len                                                           
          dxa(ll+lp1) =  dx(ll,i+1) * deps(ll,i+1)                           
          dxb(ll+lp1) = -dx(ll,i  ) * deps(ll,i+2)                           
        enddo
        lp1 = lp1 + len
        len = len - 2
     enddo
     do i = 1,jcap2_                                                        
       do ll = 1,twoj1_                                                        
         dx(ll,i) = dx(ll,i) - 1.e0                                        
       enddo
     enddo
     do ll = 1,twoj1_                                                        
       dxb(ll) = -dx(ll,1) * deps(ll,2)                                   
     enddo
     lp  = 0                                                                   
     len = twoj1_                                                              
     do i = 1,jcap1_                                                        
       do ll = 1,len                                                           
         dxc(ll+lp) = dx(ll,1)                                              
         dxd(ll+lp) = dx(ll,i) * dx(ll,i+1)                                 
       enddo
       lp  = lp  + len 
       len = len - 2
      enddo
   500 continue
!                                                                              
!    compute pln derivatives                                                   
!
   wcsa=rcs2(lat)*wgtl/rerth_                                                
   raa=wgtl/( rerth_ * rerth_ )                                              
!                                                                              
   lp0 = 0                                                                   
   lp1 = twoj1_                                                              
   len = twoj1_ - 2                                                          
!
   do i = 1,jcap_                                                         
     do ll = 1,len                                                           
       qdert(ll+lp1) = qlnt(ll+lp0) * dxa(ll+lp1)                         
     enddo
     lp0 = lp1                                                                 
     lp1 = lp1 + len                                                           
     len = len - 2                                                             
   enddo
!
   lbegin = twoj1_ + 1                                             
!
   do ll = lbegin,lnt2_                                                    
     qdert(ll) = qdert(ll) + qlnv(ll+twoj1_) * dxb(ll)                  
   enddo
!
   do ll = 1,twoj1_                                                        
     qdert(ll) = qlnv(ll+twoj1_) * dxb(ll)                              
   enddo
!
   do ll = 1,lnt2_                                                         
     qdert(ll) = qdert(ll) * wcsa                                       
   enddo
!                                                                              
   do ll = 1,lnt2_                                                         
     qlnwct(ll) = qlnt(ll) * dxc(ll) * wcsa                             
   enddo
!                                                                              
   do ll = 1,lnt2_                                                         
     qlnt(ll) = qlnt(ll) * dxd(ll) * raa                                
   enddo
!                                                                              
   return                                                                    
   end subroutine chgr_sph_derivative
!-------------------------------------------------------------------------------
