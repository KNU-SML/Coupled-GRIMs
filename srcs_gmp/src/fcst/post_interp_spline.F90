#include <define.h>
   subroutine post_interp_spline
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ... 
!
!   [post_max_wind] or [rpost_max_wind]
!         |
!    (post_interp_spline)
!         |
!         |-- [post_spline_derivative] * 
!         |     : compute 2nd derivatives for cubic spline
!         |-- [post_spline_max] * 
!               : determine maximum value of cubic spline
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   end subroutine post_interp_spline
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_spline_derivative(l,n,x,f,s)                                              
!-------------------------------------------------------------------------------
!
! subprogram: post_spline_derivative  compute 2nd derivatives for cubic splines 
!                                                                               
! abstract: compute the second derivatives of cubic spline profiles             
!   in preparation for cubic spline interpolations.                             
!   cubic splines are piecewise cubic polynomials fitting the data              
!   with continuous first and second derivatives at interior points             
!   and second derivatives set to zero at and beyond the end points.            
!   the computations are done by marching up then down the profiles.            
!   note the inner dimension of the data is the number of profiles.             
!                                                                               
! program history log:                                                          
!   1992-10-31  iredell                development                                                   
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call post_spline_derivative(l,n,x,f,s)                                              
!                                                                               
!   input argument list:                                                        
!     l        - integer number of profiles                                     
!     n        - integer number of points in each profile                       
!     x        - real (n) monotonically increasing abscissa values              
!     f        - real (l,n) data values                                         
!                                                                               
!   output argument list:                                                       
!     s        - real (l,n) 2nd derivative of f with respect to x               
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only : levs_
!-------------------------------------------------------------------------------
   real  ::  x(n),f(l,n),s(l,n)                                              
   real  ::  rho(levs_-1)
!-------------------------------------------------------------------------------
!
!  initialize end points                                                        
! 
   rho(1)=0.                                                                 
   do i = 1,l                                                                  
     s(i,1)=0.                                                               
     s(i,n)=0.                                                               
   enddo                                                                     
! 
!  march up the profiles                                                        
!
   do k = 2,n-1                                                                
     hm1=x(k)-x(k-1)                                                         
     rh=1./(x(k+1)-x(k))                                                     
     rho(k)=-1./(hm1*(rho(k-1)+2.)*rh+2.)                                    
     do i = 1,l                                                                
       d=6.*((f(i,k+1)-f(i,k))*rh-(f(i,k)-f(i,k-1))/hm1)*rh                  
       s(i,k)=(hm1*s(i,k-1)*rh-d)*rho(k)                                     
     enddo                                                                   
   enddo                                                                     
! 
!  march down the profiles                                                      
!
   do k = n-1,2,-1                                                             
     do i = 1,l                                                                
       s(i,k)=rho(k)*s(i,k+1)+s(i,k)                                         
     enddo                                                                   
   enddo                                                                     
!
   return                                                                    
   end subroutine post_spline_derivative 
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine post_spline_max(l,n,x,f,s,p,xp,fp)                                      
!-------------------------------------------------------------------------------
!
! subprogram:    post_spline_max   find maximum value using cubic splines    
!                                                                               
! abstract: compute the maximum data value of cubic spline profiles.            
!   cubic splines are piecewise cubic polynomials fitting the data              
!   with continuous first and second derivatives at interior points             
!   and second derivatives set to zero at and beyond the end points.            
!   subprogram post_spline_derivative must be already called to 
!   compute 2nd derivatives.        
!   note the inner dimension of the data is the number of profiles.             
!                                                                               
! program history log:                                                          
!   1992-10-31  iredell                                                           
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!                                                                               
! usage:    call post_spline_max(l,n,x,f,s,p,xp,fp)            
!                                                                               
!   input argument list:                                                        
!     l        - integer number of profiles                                     
!     n        - integer number of points in each profile                       
!     x        - real (n) monotonically increasing abscissa values              
!     f        - real (l,n) data values                                         
!     s        - real (l,n) 2nd derivative of f 
!                (from subprogram post_spline_derivative)        
!                                                                               
!   output argument list:                                                       
!     p        - real (l) point number                                          
!     xp       - real (l) abscissa values of maximum value                      
!     fp       - real (l) maximum data values                                   
!                                                                               
!-------------------------------------------------------------------------------
   real  ::  x(n),f(l,n),s(l,n),p(l),xp(l),fp(l)                             
!-------------------------------------------------------------------------------
!
!  find maximum given value                                                     
!
   do i = 1,l                                                                  
     p(i)=1                                                                  
     fp(i)=f(i,1)                                                            
   enddo                                                                     
   do k = 2,n                                                                  
     do i = 1,l                                                                
       if(f(i,k).gt.fp(i)) then                                              
         p(i)=k                                                              
         fp(i)=f(i,k)                                                        
       endif                                                                 
     enddo                                                                   
   enddo                                                                     
!
!  determine maximum value of cubic spline                                      
!
   do i = 1,l                                                                  
     k1=nint(p(i))                                                           
     kt=k1+sign(1,n+1-2*k1)                                                  
     dx=x(k1)-x(kt)                                                          
     df=f(i,k1)-f(i,kt)                                                      
     s1=s(i,k1)                                                              
     st=s(i,kt)                                                              
     dp=df/dx+dx*(2*s1+st)/6                                                 
     k2=k1+sign(1.,dp)                                                       
     if(k2.ge.1.and.k2.le.n) then                                            
       x1=x(k1)                                                              
       x2=x(k2)                                                              
       xm=(x2+x1)/2                                                          
       dx=x2-x1                                                              
       f1=f(i,k1)                                                            
       f2=f(i,k2)                                                            
       df=f2-f1                                                              
       s1=s(i,k1)                                                            
       s2=s(i,k2)                                                            
       sm=(s2+s1)/2                                                          
       ds=s2-s1                                                              
       if(ds.ne.0.) then                                                     
         xpa=xm-sm*dx/ds                                                     
         xpb=sqrt((dx**2*(4*sm**2-s1*s2)/(3*ds)-2*df)/ds)                    
         xp(i)=xpa+xpb                                                       
         sp=s1+ds*(xp(i)-x(k1))/dx                                           
         if(sp.gt.0.) xp(i)=xpa-xpb                                          
       elseif(s1.lt.0.) then                                                 
         xp(i)=xm-df/(dx*s1)                                                 
       else                                                                  
         xp(i)=x1                                                            
       endif                                                                 
       dxp=xp(i)-x1                                                          
       p(i)=k1+dxp/dx                                                        
       fp(i)=f1+dxp/dx*(df-(dx-dxp)*(dxp*ds+dx*(2*s1+s2))/6)                 
     else                                                                    
       p(i)=k1                                                               
       xp(i)=x(k1)                                                           
       fp(i)=f(i,k1)                                                         
     endif                                                                   
   enddo                                                                     
!
   return
   end subroutine post_spline_max
