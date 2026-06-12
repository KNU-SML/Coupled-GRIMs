#include <define.h>
   subroutine post_wind_transform
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [post_wind_transform]
!           |
!           |-- [post_wind2divor] *
!           |-- [post_divor2wind] *
!
!-------------------------------------------------------------------------------
   end subroutine post_wind_transform
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_wind2divor(m,enn1,elonn1,eon,eontop,u,v,utop,vtop,d,z)   
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:  post_wind2divor   compute divergence and vorticity from winds 
!                                                                               
! abstract: computes the divergence and vorticity from wind components          
!           in spectral space. 
!           subprogram post_spectral_field should be called already.        
!           if l is the zonal wavenumber, n is the total wavenumber,            
!           eps(l,n)=sqrt((n**2-l**2)/(4*n**2-1)) and a is earth radius,        
!           then the divergence d is computed as                                
!             d(l,n)=i*l*a*u(l,n)                                               
!                    +eps(l,n+1)*n*a*v(l,n+1)-eps(l,n)*(n+1)*a*v(l,n-1)         
!           and the vorticity z is computed as                                  
!             z(l,n)=i*l*a*v(l,n)                                               
!                    -eps(l,n+1)*n*a*u(l,n+1)+eps(l,n)*(n+1)*a*u(l,n-1)         
!           where u is the zonal wind and v is the meridional wind.             
!           u and v are weighted by the secant of latitude.                     
!           extra terms are used over top of the spectral triangle.             
!           advantage is taken of the fact that eps(l,l)=0                      
!           in order to vectorize over the entire spectral triangle.            
!                                                                               
! program history log:                                                          
!    1992-10-31  iredell                development
!    2000-03-09  songyou hong           cvs verion setup
!    2009-10-01  jung-eun kim           f90 format with standard physics modules
!    2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call post_wind2divor(m,enn1,elonn1,eon,eontop,u,v,utop,vtop,d,z)   
!                                                                               
!   input argument list:                                                        
!     m        - integer spectral truncation                                    
!     enn1     - real ((m+1)*(m+2)/2) n*(n+1)/a**2                              
!     elonn1   - real ((m+1)*(m+2)/2) l/(n*(n+1))*a                             
!     eon      - real ((m+1)*(m+2)/2) epsilon/n*a                               
!     eontop   - real (m+1) epsilon/n*a over top                                
!     u        - real ((m+1)*(m+2)) zonal wind (over coslat)                    
!     v        - real ((m+1)*(m+2)) merid wind (over coslat)                    
!     utop     - real (2*(m+1)) zonal wind (over coslat) over top               
!     vtop     - real (2*(m+1)) merid wind (over coslat) over top               
!                                                                               
!   output argument list:                                                       
!     d        - real ((m+1)*(m+2)) divergence                                  
!     z        - real ((m+1)*(m+2)) vorticity                                   
!                                                                               
!-------------------------------------------------------------------------------
   real  ::  enn1((m+1)*(m+2)/2),elonn1((m+1)*(m+2)/2)                            
   real  ::  eon((m+1)*(m+2)/2),eontop(m+1)                                       
   real  ::  u((m+1)*(m+2)),v((m+1)*(m+2)),utop(2*(m+1)),vtop(2*(m+1))            
   real  ::  d((m+1)*(m+2)),z((m+1)*(m+2))                                        
!-------------------------------------------------------------------------------
!
!  compute terms from the spectral triangle                                     
!
   i=1                                                                       
   d(2*i-1)=0.                                                               
   d(2*i)=0.                                                                 
   z(2*i-1)=0.                                                               
   z(2*i)=0.                                                                 
!
   do i = 2,(m+1)*(m+2)/2-1                                                    
     d(2*i-1)=-elonn1(i)*u(2*i)+eon(i+1)*v(2*i+1)-eon(i)*v(2*i-3)            
     d(2*i)=elonn1(i)*u(2*i-1)+eon(i+1)*v(2*i+2)-eon(i)*v(2*i-2)             
     z(2*i-1)=-elonn1(i)*v(2*i)-eon(i+1)*u(2*i+1)+eon(i)*u(2*i-3)            
     z(2*i)=elonn1(i)*v(2*i-1)-eon(i+1)*u(2*i+2)+eon(i)*u(2*i-2)             
   enddo                                                                     
!
   i=(m+1)*(m+2)/2                                                           
   d(2*i-1)=-elonn1(i)*u(2*i)-eon(i)*v(2*i-3)                                
   d(2*i)=elonn1(i)*u(2*i-1)-eon(i)*v(2*i-2)                                 
   z(2*i-1)=-elonn1(i)*v(2*i)+eon(i)*u(2*i-3)                                
   z(2*i)=elonn1(i)*v(2*i-1)+eon(i)*u(2*i-2)                                 
! 
!  compute terms from over top of the spectral triangle                         
!
   do l = 0,m                                                                  
     i=l*(2*m+1-l)/2+m+1                                                     
     d(2*i-1)=d(2*i-1)+eontop(l+1)*vtop(2*l+1)                               
     d(2*i)=d(2*i)+eontop(l+1)*vtop(2*l+2)                                   
     z(2*i-1)=z(2*i-1)-eontop(l+1)*utop(2*l+1)                               
     z(2*i)=z(2*i)-eontop(l+1)*utop(2*l+2)                                   
   enddo                                                                     
!
!  multiply by laplacian term                                                   
!
   do i = 2,(m+1)*(m+2)/2                                                      
     d(2*i-1)=d(2*i-1)*enn1(i)                                               
     d(2*i)=d(2*i)*enn1(i)                                                   
     z(2*i-1)=z(2*i-1)*enn1(i)                                               
     z(2*i)=z(2*i)*enn1(i)                                                   
   enddo                                                                     
!
   return                                                                    
   end subroutine post_wind2divor                                                                       
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_divor2wind(m,enn1,elonn1,eon,eontop,d,z,u,v,utop,vtop)
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:    post_divor2wind  compute winds from divergence and vorticity        
!                                                                               
! abstract: computes the wind components from divergence and vorticity          
!           in spectral space. 
!           subprogram post_spectral_field should be called already.        
!           if l is the zonal wavenumber, n is the total wavenumber,            
!           eps(l,n)=sqrt((n**2-l**2)/(4*n**2-1)) and a is earth radius,        
!           then the zonal wind component u is computed as                      
!             u(l,n)=-i*l/(n*(n+1))*a*d(l,n)                                    
!                    +eps(l,n+1)/(n+1)*a*z(l,n+1)-eps(l,n)/n*a*z(l,n-1)         
!           and the meridional wind component v is computed as                  
!             v(l,n)=-i*l/(n*(n+1))*a*z(l,n)                                    
!                    -eps(l,n+1)/(n+1)*a*d(l,n+1)+eps(l,n)/n*a*d(l,n-1)         
!           where d is divergence and z is vorticity.                           
!           extra terms are computed over top of the spectral triangle.         
!           advantage is taken of the fact that eps(l,l)=0                      
!           in order to vectorize over the entire spectral triangle.            
!                                                                               
! program history log:                                                          
!   1991-10-31  mark iredell                                                      
!   2000-01-01  hann-ming henry juang  mpi
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!                                                                               
! usage:    call post_divor2wind(m,enn1,elonn1,eon,eontop,d,z,         
!    &                 u,v,utop,vtop)                                           
!                                                                               
!   input argument list:                                                        
!     m        - integer spectral truncation                                    
!     enn1     - real ((m+1)*(m+2)/2) n*(n+1)/a**2                              
!     elonn1   - real ((m+1)*(m+2)/2) l/(n*(n+1))*a                             
!     eon      - real ((m+1)*(m+2)/2) epsilon/n*a                               
!     eontop   - real (m+1) epsilon/n*a over top                                
!     d        - real ((m+1)*(m+2)) divergence                                  
!     z        - real ((m+1)*(m+2)) vorticity                                   
!                                                                               
!   output argument list:                                                       
!     u        - real ((m+1)*(m+2)) zonal wind (times coslat)                   
!     v        - real ((m+1)*(m+2)) merid wind (times coslat)                   
!     utop     - real (2*(m+1)) zonal wind (times coslat) over top              
!     vtop     - real (2*(m+1)) merid wind (times coslat) over top              
!                                                                               
!-------------------------------------------------------------------------------
   real  ::  enn1((m+1)*(m+2)/2),elonn1((m+1)*(m+2)/2)                            
   real  ::  eon((m+1)*(m+2)/2),eontop(m+1)                                       
   real  ::  d((m+1)*(m+2)),z((m+1)*(m+2))                                        
   real  ::  u((m+1)*(m+2)),v((m+1)*(m+2)),utop(2*(m+1)),vtop(2*(m+1))            
!-------------------------------------------------------------------------------
!
!  compute winds in the spectral triangle                                       
!
   i=1                                                                       
   u(2*i-1)=eon(i+1)*z(2*i+1)                                                
   u(2*i)=eon(i+1)*z(2*i+2)                                                  
   v(2*i-1)=-eon(i+1)*d(2*i+1)                                               
   v(2*i)=-eon(i+1)*d(2*i+2)                                                 
!
   do i = 2,(m+1)*(m+2)/2-1                                                    
     u(2*i-1)=elonn1(i)*d(2*i)+eon(i+1)*z(2*i+1)-eon(i)*z(2*i-3)             
     u(2*i)=-elonn1(i)*d(2*i-1)+eon(i+1)*z(2*i+2)-eon(i)*z(2*i-2)            
     v(2*i-1)=elonn1(i)*z(2*i)-eon(i+1)*d(2*i+1)+eon(i)*d(2*i-3)             
     v(2*i)=-elonn1(i)*z(2*i-1)-eon(i+1)*d(2*i+2)+eon(i)*d(2*i-2)            
   enddo                                                                     
!
   i=(m+1)*(m+2)/2                                                           
   u(2*i-1)=elonn1(i)*d(2*i)-eon(i)*z(2*i-3)                                 
   u(2*i)=-elonn1(i)*d(2*i-1)-eon(i)*z(2*i-2)                                
   v(2*i-1)=elonn1(i)*z(2*i)+eon(i)*d(2*i-3)                                 
   v(2*i)=-elonn1(i)*z(2*i-1)+eon(i)*d(2*i-2)                                
!
!  compute winds over top of the spectral triangle                              
!
   do l = 0,m                                                                  
     i=l*(2*m+1-l)/2+m+1                                                     
     utop(2*l+1)=-eontop(l+1)*z(2*i-1)                                       
     utop(2*l+2)=-eontop(l+1)*z(2*i)                                         
     vtop(2*l+1)=eontop(l+1)*d(2*i-1)                                        
     vtop(2*l+2)=eontop(l+1)*d(2*i)                                          
   enddo                                                                     
!
   return                                                                    
   end subroutine post_divor2wind  
!-------------------------------------------------------------------------------                                                                    
