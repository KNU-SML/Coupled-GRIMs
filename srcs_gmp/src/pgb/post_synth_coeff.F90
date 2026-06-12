#include <define.h>
   subroutine post_synth_coeff
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [post_synth_coeff]
!           |
!           |-- [post_synth_coeff1] *
!           |-- [post_synth_coeff2] *
!
!-------------------------------------------------------------------------------
   end subroutine post_synth_coeff
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_synth_coeff1(m,im,nc,nctop,km,pln,plntop,spc,spctop,f) 
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:    post_synth_coeff1      synthesize fourier from spectral 
!                                                                               
! abstract: synthesizes fourier coefficients from spectral coefficients         
!           for a latitude pair (northern and southern hemispheres).            
!                                                                               
! program history log:                                                          
!    1992-10-31  iredell                development
!    2000-03-09  songyou hong           cvs verion setup
!    2009-10-01  jung-eun kim           f90 format with standard physics modules
!    2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call post_synth_coeff1(m,im,nc,nctop,km,pln,plntop,spc,spctop,f)   
!                                                                               
!   input argument list:                                                        
!     m        - integer spectral truncation                                    
!     im       - integer dimension of fourier coefficients (im>=2*(m+1))        
!     nc       - integer dimension of spectral coefficients                     
!                (nc>=(m+1)*(m+2))                                              
!     nctop    - integer dimension of spectral coefficients over top            
!                (nctop>=2*(m+1))                                               
!     km       - integer number of fields                                       
!     pln      - real ((m+1)*(m+2)/2) legendre polynomial                       
!     plntop   - real (m+1) legendre polynomial over top                        
!     spc      - real (nc,km) spectral coefficients                             
!     spctop   - real (nctop,km) spectral coefficients over top                 
!                                                                               
!   output argument list:                                                       
!     f        - real (im,2,km) fourier coefficients for latitude pair          
!                                                                               
! subprograms called:                                                           
!   sgemvx1      cray library matrix times vector                               
!                                                                               
!-------------------------------------------------------------------------------
   real  ::  pln((m+1)*(m+2)/2),plntop(m+1)                                       
   real  ::  spc(nc,km),spctop(nctop,km)                                          
   real  ::  f(im,2,km)                                                           
!-------------------------------------------------------------------------------
   mmax=min(m,(im-2)/2)
!
!  initialize fourier coefficients with terms over top of the spectrum.         
!  initialize even and odd polynomials separately.                              
!
   ltope=mod(mmax+1,2)
   ltopo=1-ltope                                                             
   do k = 1,km                                                                 
     do l = ltope,mmax,2
       f(2*l+1,1,k)=plntop(l+1)*spctop(2*l+1,k)                              
       f(2*l+2,1,k)=plntop(l+1)*spctop(2*l+2,k)                              
       f(2*l+1,2,k)=0.                                                       
       f(2*l+2,2,k)=0.                                                       
     enddo                                                                   
     do l = ltopo,mmax,2
       f(2*l+1,1,k)=0.                                                       
       f(2*l+2,1,k)=0.                                                       
       f(2*l+1,2,k)=plntop(l+1)*spctop(2*l+1,k)                              
       f(2*l+2,2,k)=plntop(l+1)*spctop(2*l+2,k)                              
     enddo                                                                   
   enddo                                                                     
!
!  for each zonal wavenumber, synthesize terms over total wavenumber.           
!  synthesize even and odd polynomials separately.                              
!  commented code replaced by library calls.                                    
!
   do l = 0,mmax                                                               
     is=l*(2*m+1-l)
     ip=is/2+1                                                               
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
     call sgemvx1(km,(m+2-l)/2,1.,spc(is+2*l+1,1),nc,4,pln(ip+l),2,            &
                  1.,f(2*l+1,1,1),im*2)                                      
     call sgemvx1(km,(m+2-l)/2,1.,spc(is+2*l+2,1),nc,4,pln(ip+l),2,            &
                  1.,f(2*l+2,1,1),im*2)                                      
     call sgemvx1(km,(m+1-l)/2,1.,spc(is+2*l+3,1),nc,4,pln(ip+l+1),2,          &
                  1.,f(2*l+1,2,1),im*2)                                      
     call sgemvx1(km,(m+1-l)/2,1.,spc(is+2*l+4,1),nc,4,pln(ip+l+1),2,          &
                  1.,f(2*l+2,2,1),im*2)                                      
#endif
#ifdef DEFAULT
     do n = l,m,2
       do k = 1,km                                                             
         f(2*l+1,1,k)=f(2*l+1,1,k)+pln(ip+n)*spc(is+2*n+1,k)                 
         f(2*l+2,1,k)=f(2*l+2,1,k)+pln(ip+n)*spc(is+2*n+2,k)                 
       enddo                                                                 
     enddo                                                                   
     do n = l+1,m,2
       do k = 1,km                                                             
         f(2*l+1,2,k)=f(2*l+1,2,k)+pln(ip+n)*spc(is+2*n+1,k)                 
         f(2*l+2,2,k)=f(2*l+2,2,k)+pln(ip+n)*spc(is+2*n+2,k)                 
       enddo                                                                 
     enddo                                                                   
#endif
   enddo                                                                     
!
!  separate fourier coefficients from each hemisphere.                          
!  odd polynomials contribute negatively to the southern hemisphere.            
!
   do k = 1,km   
     do l = 0,mmax
       f1r=f(2*l+1,1,k)                                                      
       f1i=f(2*l+2,1,k)                                                      
       f(2*l+1,1,k)=f1r+f(2*l+1,2,k)                                         
       f(2*l+2,1,k)=f1i+f(2*l+2,2,k)                                         
       f(2*l+1,2,k)=f1r-f(2*l+1,2,k)                                         
       f(2*l+2,2,k)=f1i-f(2*l+2,2,k)                                         
     enddo                                                                   
   enddo                                                                     
! 
!  zero out fourier waves outside of spectrum                                   
!
   if(2*m+3.gt.im) return
   do l2 = 2*m+3,im                                                            
     do k = 1,km                                                               
       f(l2,1,k)=0.                                                          
       f(l2,2,k)=0.                                                          
     enddo                                                                   
   enddo                                                                     
!
   return                                                                    
   end subroutine post_synth_coeff1                                                                      
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_synth_coeff2(m,im,ix,nc,nctop,km,clat,pln,plntop,mp,        &
                                spc,spctop,f)                          
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:    post_synth_coeff2     synthesize fourier from spectral    
!                                                                               
! abstract: synthesizes fourier coefficients from spectral coefficients         
!           for a latitude pair (northern and southern hemispheres).            
!           vector components are divided by cosine of latitude.                
!                                                                               
! program history log:                                                          
!   91-10-31  mark iredell                                                      
!                                                                               
! usage:    call post_synth_coeff2(m,im,ix,nc,nctop,km,clat,pln,plntop,mp,
!    &                  spc,spctop,f)                                           
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
!     clat     - real cosine of latitude                                        
!     pln      - real ((m+1)*(m+2)/2) legendre polynomial                       
!     plntop   - real (m+1) legendre polynomial over top                        
!     spc      - real (nc,km) spectral coefficients                             
!     spctop   - real (nctop,km) spectral coefficients over top                 
!     mp       - integer (km) identifiers (0 for scalar, 1 for vector)          
!                                                                               
!   output argument list:                                                       
!     f        - real (ix,2,km) fourier coefficients for latitude pair          
!                                                                               
! subprograms called:                                                           
!   sgemvx1      cray library matrix times vector                               
!                                                                               
!-------------------------------------------------------------------------------
   integer  ::  mp(km)                                                            
   real     ::  pln((m+1)*(m+2)/2),plntop(m+1)                                       
   real     ::  spc(nc,km),spctop(nctop,km)                                          
   real     ::  f(ix,2,km)                                                           
!-------------------------------------------------------------------------------
!
!  synthesis over pole.                                                         
!  zero out fourier waves.                                                      
!
   if(clat.eq.0) then                                                        
     do k = 1,km                                                               
       do l = 0,im/2                                                           
         f(2*l+1,1,k)=0.                                                     
         f(2*l+2,1,k)=0.                                                     
         f(2*l+1,2,k)=0.                                                     
         f(2*l+2,2,k)=0.                                                     
       enddo                                                                 
     enddo                                                                   
!
!  initialize fourier coefficients with terms over top of the spectrum.         
!  initialize even and odd polynomials separately.                              
!
     ltope=mod(m+1,2)                                                        
     do k = 1,km                                                               
       l=mp(k)                                                               
       if(l.eq.ltope) then                                                   
         f(2*l+1,1,k)=plntop(l+1)*spctop(2*l+1,k)                            
         f(2*l+2,1,k)=plntop(l+1)*spctop(2*l+2,k)                            
         f(2*l+1,2,k)=0.                                                     
         f(2*l+2,2,k)=0.                                                     
       else                                                                  
         f(2*l+1,1,k)=0.                                                     
         f(2*l+2,1,k)=0.                                                     
         f(2*l+1,2,k)=plntop(l+1)*spctop(2*l+1,k)                            
         f(2*l+2,2,k)=plntop(l+1)*spctop(2*l+2,k)                            
       endif                                                                 
     enddo                                                                   
!
!  for each zonal wavenumber, synthesize terms over total wavenumber.           
!  synthesize even and odd polynomials separately.                              
!
     do k = 1,km                                                               
       l=mp(k)                                                               
       is=l*(2*m+1-l)                                                        
       ip=is/2+1                                                             
       do n = l,m,2                                                            
         f(2*l+1,1,k)=f(2*l+1,1,k)+pln(ip+n)*spc(is+2*n+1,k)                 
         f(2*l+2,1,k)=f(2*l+2,1,k)+pln(ip+n)*spc(is+2*n+2,k)                 
       enddo                                                                 
       do n = l+1,m,2                                                          
         f(2*l+1,2,k)=f(2*l+1,2,k)+pln(ip+n)*spc(is+2*n+1,k)                 
         f(2*l+2,2,k)=f(2*l+2,2,k)+pln(ip+n)*spc(is+2*n+2,k)                 
       enddo                                                                 
     enddo                                                                   
! 
!  separate fourier coefficients from each hemisphere.                          
!  odd polynomials contribute negatively to the southern hemisphere.            
!
#ifdef CRAY_THREAD
!dir$ ivdep                                                                     
#endif
     do k = 1,km                                                               
       l=mp(k)                                                               
       f1r=f(2*l+1,1,k)                                                      
       f1i=f(2*l+2,1,k)                                                      
       f(2*l+1,1,k)=f1r+f(2*l+1,2,k)                                         
       f(2*l+2,1,k)=f1i+f(2*l+2,2,k)                                         
       f(2*l+1,2,k)=f1r-f(2*l+1,2,k)                                         
       f(2*l+2,2,k)=f1i-f(2*l+2,2,k)                                         
     enddo                                                                   
!                                                                               
!  synthesis over finite latitude.                                              
!  initialize fourier coefficients with terms over top of the spectrum.         
!  initialize even and odd polynomials separately.                              
!                                                                               
   else                                                                      
     lx=min(m,im/2)                                                          
     ltope=mod(m+1,2)                                                        
     ltopo=1-ltope                                                           
     do k = 1,km                                                               
       if(mp(k).eq.0) then                                                   
         do l = ltope,lx,2                                                     
           f(2*l+1,1,k)=0.                                                   
           f(2*l+2,1,k)=0.                                                   
           f(2*l+1,2,k)=0.                                                   
           f(2*l+2,2,k)=0.                                                   
         enddo                                                               
         do l = ltopo,lx,2                                                     
           f(2*l+1,1,k)=0.                                                   
           f(2*l+2,1,k)=0.                                                   
           f(2*l+1,2,k)=0.                                                   
           f(2*l+2,2,k)=0.                                                   
         enddo                                                               
       else                                                                  
         do l = ltope,lx,2                                                     
           f(2*l+1,1,k)=plntop(l+1)*spctop(2*l+1,k)                          
           f(2*l+2,1,k)=plntop(l+1)*spctop(2*l+2,k)                          
           f(2*l+1,2,k)=0.                                                   
           f(2*l+2,2,k)=0.                                                   
         enddo                                                               
         do l = ltopo,lx,2                                                     
           f(2*l+1,1,k)=0.                                                   
           f(2*l+2,1,k)=0.                                                   
           f(2*l+1,2,k)=plntop(l+1)*spctop(2*l+1,k)                          
           f(2*l+2,2,k)=plntop(l+1)*spctop(2*l+2,k)                          
         enddo                                                               
       endif                                                                 
     enddo                                                                   
! 
!  for each zonal wavenumber, synthesize terms over total wavenumber.           
!  synthesize even and odd polynomials separately.                              
!  commented code replaced by library calls.                                    
!
     do l = 0,lx                                                               
       is=l*(2*m+1-l)                                                        
       ip=is/2+1                                                             
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
       call sgemvx1(km,(m+2-l)/2,1.,spc(is+2*l+1,1),nc,4,                      &
                    pln(ip+l),2,1.,f(2*l+1,1,1),ix*2)                        
       call sgemvx1(km,(m+2-l)/2,1.,spc(is+2*l+2,1),nc,4,                      &
                    pln(ip+l),2,1.,f(2*l+2,1,1),ix*2)                        
       call sgemvx1(km,(m+1-l)/2,1.,spc(is+2*l+3,1),nc,4,                      &
                    pln(ip+l+1),2,1.,f(2*l+1,2,1),ix*2)                      
       call sgemvx1(km,(m+1-l)/2,1.,spc(is+2*l+4,1),nc,4,                      &
                    pln(ip+l+1),2,1.,f(2*l+2,2,1),ix*2)                      
#endif
#ifdef DEFAULT
       do n = l,m,2                                                            
         do k = 1,km                                                           
           f(2*l+1,1,k)=f(2*l+1,1,k)+pln(ip+n)*spc(is+2*n+1,k)               
           f(2*l+2,1,k)=f(2*l+2,1,k)+pln(ip+n)*spc(is+2*n+2,k)               
         enddo                                                               
       enddo                                                                 
       do n = l+1,m,2                                                          
         do k = 1,km                                                           
           f(2*l+1,2,k)=f(2*l+1,2,k)+pln(ip+n)*spc(is+2*n+1,k)               
           f(2*l+2,2,k)=f(2*l+2,2,k)+pln(ip+n)*spc(is+2*n+2,k)               
         enddo                                                               
       enddo                                                                 
#endif
     enddo                                                                   
!
!  separate fourier coefficients from each hemisphere.                          
!  odd polynomials contribute negatively to the southern hemisphere.            
!  divide vector components by cosine latitude.                                 
!
     do k = 1,km                                                               
       do l = 0,lx                                                             
         f1r=f(2*l+1,1,k)                                                    
         f1i=f(2*l+2,1,k)                                                    
         f(2*l+1,1,k)=f1r+f(2*l+1,2,k)                                       
         f(2*l+2,1,k)=f1i+f(2*l+2,2,k)                                       
         f(2*l+1,2,k)=f1r-f(2*l+1,2,k)                                       
         f(2*l+2,2,k)=f1i-f(2*l+2,2,k)                                       
       enddo                                                                 
       if(mp(k).eq.1) then                                                   
         do l = 0,lx                                                           
           f(2*l+1,1,k)=f(2*l+1,1,k)/clat                                    
           f(2*l+2,1,k)=f(2*l+2,1,k)/clat                                    
           f(2*l+1,2,k)=f(2*l+1,2,k)/clat                                    
           f(2*l+2,2,k)=f(2*l+2,2,k)/clat                                    
         enddo                                                               
       endif                                                                 
     enddo                                                                   
!
!  zero out fourier waves outside of spectrum                                   
!
     do l = lx+1,im/2                                                          
       do k = 1,km                                                             
         f(2*l+1,1,k)=0.                                                     
         f(2*l+2,1,k)=0.                                                     
         f(2*l+1,2,k)=0.                                                     
         f(2*l+2,2,k)=0.                                                     
       enddo                                                                 
     enddo                                                                   
   endif                                                                     
! 
   return                                                                    
   end subroutine post_synth_coeff2                                                                       
