#include <define.h>
   subroutine post_sph_fft_coeff(m,im,ix,nc,nctop,km,wgt,clat,pln,plntop,mp,   &
                                 f,spc,spctop)                    
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:    post_sph_fft_coeff      analyze spectral from fourier   
!                                                                               
! abstract: analyzes spectral coefficients from fourier coefficients            
!           for a latitude pair (northern and southern hemispheres).            
!           vector components are multiplied by cosine of latitude.             
!                                                                               
! program history log:                                                          
!    1992-10-31  iredell                development
!    2000-03-09  songyou hong           cvs verion setup, prognostic clouds
!    2009-10-01  jung-eun kim           f90 format with standard physics modules
!    2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call post_sph_fft_coeff(m,im,ix,nc,nctop,km,wgt,clat,pln,plntop,mp,
!    &                  f,spc,spctop)                                           
!                                                                               
!   input argument list:                                                        
!     m        - integer spectral truncation                                    
!     im       - integer even number of fourier coefficients                    
!     ix       - integer dimension of fourier coefficients (ix>=im+2)           
!     nc       - integer dimension of spectral coefficients                     
!                (nc>=(m+1)*(m+2))                                              
!     nctop    - integer dimension of spectral coefficients over top            
!                (nctop>=2*(m+1))                                               
!     km       - integer number of fields                                       
!     wgt      - real gaussian weight                                           
!     clat     - real cosine of latitude                                        
!     pln      - real ((m+1)*(m+2)/2) legendre polynomials                      
!     plntop   - real (m+1) legendre polynomial over top                        
!     mp       - integer (km) identifiers (0 for scalar, 1 for vector)          
!     f        - real (ix,2,km) fourier coefficients combined                   
!     spc      - real (nc,km) spectral coefficients                             
!     spctop   - real (nctop,km) spectral coefficients over top                 
!                                                                               
!   output argument list:                                                       
!     spc      - real (nc,km) spectral coefficients                             
!     spctop   - real (nctop,km) spectral coefficients over top                 
!                                                                               
!-------------------------------------------------------------------------------
   integer  ::  mp(km)                                                            
   real     ::  pln((m+1)*(m+2)/2),plntop(m+1)                                       
   real     ::  f(ix,2,km)                                                           
   real     ::  spc(nc,km),spctop(nctop,km)                                          
   real     ::  fw(2,2,km)                                                           
!-------------------------------------------------------------------------------
!
!  for each zonal wavenumber, analyze terms over total wavenumber.              
!  analyze even and odd polynomials separately.                                 
!  commented code replaced by library calls.                                    
!
   lx=min(m,im/2)                                                            
   do l = 0,lx                                                                 
     nt=mod(m+1-l,2)+1                                                       
!
     do k = 1,km                                                               
       if(mp(k).eq.0) then                                                   
         fw(1,1,k)=wgt*(f(2*l+1,1,k)+f(2*l+1,2,k))                           
         fw(2,1,k)=wgt*(f(2*l+2,1,k)+f(2*l+2,2,k))                           
         fw(1,2,k)=wgt*(f(2*l+1,1,k)-f(2*l+1,2,k))                           
         fw(2,2,k)=wgt*(f(2*l+2,1,k)-f(2*l+2,2,k))                           
       else                                                                  
         fw(1,1,k)=wgt*clat*(f(2*l+1,1,k)+f(2*l+1,2,k))                      
         fw(2,1,k)=wgt*clat*(f(2*l+2,1,k)+f(2*l+2,2,k))                      
         fw(1,2,k)=wgt*clat*(f(2*l+1,1,k)-f(2*l+1,2,k))                      
         fw(2,2,k)=wgt*clat*(f(2*l+2,1,k)-f(2*l+2,2,k))                      
         spctop(2*l+1,k)=spctop(2*l+1,k)+plntop(l+1)*fw(1,nt,k)              
         spctop(2*l+2,k)=spctop(2*l+2,k)+plntop(l+1)*fw(2,nt,k)              
       endif                                                                 
     enddo                                                                   
     is=l*(2*m+1-l)                                                          
     ip=is/2+1                                                               
!
     do n = l,m,2                                                              
       do k = 1,km                                                             
         spc(is+2*n+1,k)=spc(is+2*n+1,k)+pln(ip+n)*fw(1,1,k)                 
         spc(is+2*n+2,k)=spc(is+2*n+2,k)+pln(ip+n)*fw(2,1,k)                 
       enddo                                                                 
     enddo                                                                   
!
     do n = l+1,m,2                                                            
       do k = 1,km                                                             
         spc(is+2*n+1,k)=spc(is+2*n+1,k)+pln(ip+n)*fw(1,2,k)                 
         spc(is+2*n+2,k)=spc(is+2*n+2,k)+pln(ip+n)*fw(2,2,k)                 
       enddo                                                                 
     enddo                                                                   
!
   enddo                                                                     
!
   return                                                                    
   end subroutine post_sph_fft_coeff                                                                       
