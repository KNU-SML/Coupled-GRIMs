#include <define.h>
   subroutine post_legendre_poly
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [post_legendre_poly]
!           |
!           |-- [post_sph_poly] *
!           |-- [post_legendre] *
!
!-------------------------------------------------------------------------------
   end subroutine post_legendre_poly
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_sph_poly(m,slat,clat,eps,epstop,pln,plntop)    
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:    post_sph_poly        compute legendre polynomials 
!                                                                               
! abstract: evaluates the orthonormal associated legendre polynomials           
!           in the spectral triangle at a given latitude.                       
!           subprogram post_spectral_field should be called already.    
!           if l is the zonal wavenumber, n is the total wavenumber,            
!           and eps(l,n)=sqrt((n**2-l**2)/(4*n**2-1)) then                      
!           the following bootstrapping formulas are used:                      
!           pln(0,0)=sqrt(0.5)                                                  
!           pln(l,l)=pln(l-1,l-1)*clat*sqrt(float(2*l+1)/float(2*l))            
!           pln(l,n)=(slat*pln(l,n-1)-eps(l,n-1)*pln(l,n-2))/eps(l,n)           
!           synthesis at the pole needs only two zonal wavenumbers.             
!           scalar fields are synthesized with zonal wavenumber 0 while         
!           vector fields are synthesized with zonal wavenumber 1.              
!           (thus polar vector fields are implicitly divided by clat.)          
!           the following bootstrapping formulas are used at the pole:          
!           pln(0,0)=sqrt(0.5)                                                  
!           pln(1,1)=sqrt(0.75)                                                 
!           pln(l,n)=(pln(l,n-1)-eps(l,n-1)*pln(l,n-2))/eps(l,n)                
!                                                                               
! usage:    call post_sph_poly(m,slat,clat,eps,epstop,pln,plntop)   
!                                                                               
!   input argument list:                                                        
!     m        - integer spectral truncation                                    
!     slat     - real sine of latitude                                          
!     clat     - real cosine of latitude                                        
!     eps      - real ((m+1)*(m+2)/2) sqrt((n**2-l**2)/(4*n**2-1))              
!     epstop   - real (m+1) sqrt((n**2-l**2)/(4*n**2-1)) over top               
!                                                                               
!   output argument list:                                                       
!     pln      - real ((m+1)*(m+2)/2) legendre polynomial                       
!     plntop   - real (m+1) legendre polynomial over top                        
!                                                                               
!-------------------------------------------------------------------------------
   real eps((m+1)*(m+2)/2),epstop(m+1)                                       
   real pln((m+1)*(m+2)/2),plntop(m+1)                                       
!-------------------------------------------------------------------------------
!
!  iteratively compute pln within spectral triangle at pole                     
!
   if(clat.eq.0.) then                                                       
     pln(1)=sqrt(0.5)                                                        
     pln(m+2)=sqrt(0.75)                                                     
     pln(2)=pln(1)/eps(2)                                                    
     pln(m+3)=pln(m+2)/eps(m+3)                                              
     pln(3)=(pln(2)-eps(2)*pln(1))/eps(3)                                    
     do n = 3,m                                                                
       i=n+1                                                                 
       pln(i)=(pln(i-1)-eps(i-1)*pln(i-2))/eps(i)                            
       i=n+m+1                                                               
       pln(i)=(pln(i-1)-eps(i-1)*pln(i-2))/eps(i)                            
     enddo                                                                   
     do i = 2*m+2,(m+1)*(m+2)/2                                                
       pln(i)=0.                                                             
     enddo                                                                   
! 
!  compute polynomials over top of spectral triangle                            
!
     i=m+2                                                                   
     plntop(1)=(pln(i-1)-eps(i-1)*pln(i-2))/epstop(1)                        
     i=2*m+2                                                                 
     plntop(2)=(pln(i-1)-eps(i-1)*pln(i-2))/epstop(2)                        
     do l = 2,m                                                                
       plntop(l+1)=0.                                                        
     enddo                                                                   
!
!  iteratively compute pln(l,l) (bottom hypotenuse of triangle)                 
!
   else                                                                      
     nml=0                                                                   
     i=1                                                                     
     pln(i)=sqrt(0.5)                                                        
     do l = 1,m-nml                                                            
       plni=pln(i)                                                           
       i=l*(2*m+3-l)/2+(nml+1)                                               
       pln(i)=plni*clat*sqrt(float(2*l+1)/float(2*l))                        
     enddo                                                                   
!
!  compute pln(l,l+1) (diagonal next to bottom hypotenuse of triangle)          
!
     nml=1                                                                   
#ifdef CRAY_THREAD
!dir$ ivdep                                                                     
#endif
     do l = 0,m-nml                                                            
       i=l*(2*m+3-l)/2+(nml+1)                                               
       pln(i)=slat*pln(i-1)/eps(i)                                           
     enddo                                                                   
!
!  compute remaining pln in spectral triangle                                   
!
     do nml = 2,m                                                              
#ifdef CRAY_THREAD
!dir$ ivdep                                                                     
#endif
       do l = 0,m-nml                                                          
         i=l*(2*m+3-l)/2+(nml+1)                                             
         pln(i)=(slat*pln(i-1)-eps(i-1)*pln(i-2))/eps(i)                     
       enddo                                                                 
     enddo                                                                   
!
!  compute polynomials over top of spectral triangle                            
!
     do l = 0,m                                                                
       nml=m+1-l                                                             
       i=l*(2*m+3-l)/2+(nml+1)                                               
       plntop(l+1)=(slat*pln(i-1)-eps(i-1)*pln(i-2))/epstop(l+1)             
     enddo                                                                   
   endif                                                                     
!
   return                                                                    
   end subroutine post_sph_poly                                                                      
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_legendre(ipd,cosclt,maxwv,iromb,ex,px,pnm)
!-------------------------------------------------------------------------------
!
! abstract: 
!   - compute legendre polynomials
!   - evaluates the normalized associated legendre polynomials
!   at a given cosine of colatitude.
!
! usage:    call post_legendre(ipd,cosclt,maxwv,iromb,ex,px,pnm)
!
!   input argument list:
!     ipd      - ipd = 0 to evaluate polynomial itself,
!                ipd = 1 to evaluate polynomial derivative.
!     cosclt   - cosine of colatitude at which to evalutate.
!     maxwv    - spectral truncation
!     iromb    - iromb = 0 for triangular truncation
!                iromb = 1 for rhomboidal truncation
!
!   work argument list:
!     ex       - real (maxwv+2) work area
!     px       - real (maxwv+3) work area
!
!   output argument list:
!     pnm      - real (kmax) polynomial values,
!                where kmax=(maxwv+1)*(iromb+1)*maxwv+2)/2.
!
!-------------------------------------------------------------------------------
   save
   integer              ::  ipd,iromb,maxwv
   real                 ::  pnm((maxwv+1)*((iromb+1)*maxwv+2)/2)
   real                 ::  ex(0:maxwv+1),px(-1:maxwv+1)
!-------------------------------------------------------------------------------
   sinclt=sqrt(1.-cosclt**2)
   k=0
   ex(0)=0.
   px(-1)=0.
   px(0)=sqrt(0.5)
   do m = 0,maxwv
     if(m.gt.0) px(0)=px(0)*sinclt/(ex(1)*sqrt(float(2*m)))
     do n = m+1,maxwv+iromb*m+1
       ex(n-m)=sqrt(float(n**2-m**2)/float(4*n**2-1))
     enddo
     do n = m+1,maxwv+iromb*m+1
       px(n-m)=(cosclt*px(n-m-1)-ex(n-m-1)*px(n-m-2))/ex(n-m)
     enddo
     if(ipd.eq.0) then
       do n = m,maxwv+iromb*m
         k=k+1
         pnm(k)=px(n-m)
       enddo
     else
       do n = m,maxwv+iromb*m
         k=k+1
         pnm(k)=n*ex(n-m+1)*px(n-m+1)-(n+1)*ex(n-m)*px(n-m-1)
       enddo
     endif
   enddo
!
   return
   end subroutine post_legendre
