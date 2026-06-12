#include <define.h>
   subroutine post_gauss2latlon(gauin,imxin,jmxin,xlonw,xlatn,                 &
                                dxlon,dxlat,regout,imxout,                     &
                                jmxout,undef,luptr)    
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_,latg_
!-------------------------------------------------------------------------------
!                                                                               
!  abstract: interpolation from lat/lon grid to other lat/lon grid                        
!
!  ::: structure ::: This file contains ... 
!
!     [post_guass2latlon] * ----- [chgr_gaussian_lat] *
!                             |-- [dyn_bessel_function] *
!
!-------------------------------------------------------------------------------
   real                      ::  gauin (imxin,jmxin)                                             
   real                      ::  regout(imxout,jmxout)                                           
   real                      ::  gaul(latg_),regl(latg_)                                             
   real,allocatable,save     ::  ddx(:)
   real,allocatable,save     ::  ddy(:)
   integer,allocatable,save  ::  iindx1(:)
   integer,allocatable,save  ::  iindx2(:)
   integer,allocatable,save  ::  jindx1(:)
   integer,allocatable,save  ::  jindx2(:)
!-------------------------------------------------------------------------------
   data  ifp/0/                                                               
   save ifp
!-------------------------------------------------------------------------------
   allocate(iindx1(lonf_))
   allocate(iindx2(lonf_))
   allocate(jindx1(latg_))
   allocate(jindx2(latg_))
   allocate(ddx(lonf_))
   allocate(ddy(latg_))
!                                                                               
   if(ifp.ne.0) go to 111                                                    
   ifp=1                                                                     
!                                                                               
   write(luptr,*) 'imxin=',imxin,' jmxin=',jmxin                             
   write(luptr,*) 'xlatn=',xlatn,' xlonw=',xlonw,                              &
                  ' dxlat=',dxlat,' dxlon=',dxlon                            
   write(luptr,*) 'imxout=',imxout,' jmxout=',jmxout                         
!                                                                               
   call chgr_gaussian_lat(gaul,jmxin)   
!                                                                               
   do 20 j = 1,jmxout                                                          
     regl(j)=xlatn-float(j-1)*dxlat                                            
   20 continue                                                                  
!                                                                               
   dxin =360./float(imxin)                                                   
!                                                                               
   do i = 1,imxout                                                          
     alamd=xlonw+float(i-1)*dxlon                                              
     if(alamd.lt.0.) alamd=360.+alamd                                          
     i1=alamd/dxin+1.001                                                       
     if(i1.gt.imxin) i1=1                                                      
     iindx1(i)=i1                                                              
     i2=i1+1                                                                   
     if(i2.gt.imxin) i2=1                                                      
     iindx2(i)=i2                                                              
     ddx(i)=(alamd-float(i1-1)*dxin)/dxin                                      
   enddo 
!                                                                               
   do j = 1,jmxout                                                          
     aphi=regl(j)                                                              
     do 50 jj = 1,jmxin                                                          
       if(aphi.lt.gaul(jj)) go to 50                                             
       j2=jj                                                                     
       !
       go to 42                                                                  
       !
     50 continue                                                                  
     j2=jmxin                                                                  
     42 continue                                                                  
     if(j2.gt.2) go to 43                                                      
     j1=1                                                                      
     j2=2                                                                      
     !
     go to 44                                                                  
     !
     43 continue                                                                  
     if(j2.le.jmxin) go to 45                                                  
     j1=jmxin-1                                                                
     j2=jmxin                                                                  
     !
     go to 44                                                                  
     !
     45 continue
     j1=j2-1                                                                   
     44 continue
     jindx1(j)=j1                                                              
     jindx2(j)=j2                                                              
     ddy(j)=(aphi-gaul(j1))/(gaul(j2)-gaul(j1))                                
   enddo 
!                                                                               
 111 continue                                                                  
!                                                                               
!     write(luptr,*) 'iindx1'                                                   
!     write(luptr,*) (iindx1(n),n=1,imxout)                                     
!     write(luptr,*) 'iindx2'                                                   
!     write(luptr,*) (iindx2(n),n=1,imxout)                                     
!     write(luptr,*) 'jindx1'                                                   
!     write(luptr,*) (jindx1(n),n=1,jmxout)                                     
!     write(luptr,*) 'jindx2'                                                   
!     write(luptr,*) (jindx2(n),n=1,jmxout)                                     
!     write(luptr,*) 'ddy'                                                      
!     write(luptr,*) (ddy(n),n=1,jmxout)                                        
!     write(luptr,*) 'ddx'                                                      
!     write(luptr,*) (ddx(n),n=1,jmxout)                                        
!                                                                               
   do j = 1,jmxout                                                          
     y=ddy(j)                                                                  
     j1=jindx1(j)                                                              
     j2=jindx2(j)                                                              
     do i = 1,imxout                                                          
       x=ddx(i)                                                                  
       i1=iindx1(i)                                                              
       i2=iindx2(i)                                                              
       if(gauin(i1,j1).eq.undef.or.gauin(i2,j1).eq.undef.or.                   &
         gauin(i1,j2).eq.undef.or.gauin(i2,j2).eq.undef) then               
         regout(i,j)=undef                                                       
       else                                                                      
         regout(i,j)=(1.-x)*(1.-y)*gauin(i1,j1)+(1.-y)*x*gauin(i2,j1)+         &
                     (1.-x)*y*gauin(i1,j2)+x*y*gauin(i2,j2)                         
       endif                                                                     
     enddo    
   enddo     
!                                                                               
   sum1=0.                                                                   
   sum2=0.                                                                   
   do i = 1,imxin                                                           
     sum1=sum1+gauin(i,1)                                                      
     sum2=sum2+gauin(i,jmxin)                                                  
   enddo 
   sum1=sum1/float(imxin)                                                    
   sum2=sum2/float(imxin)                                                    
   do i = 1,imxin                                                              
     if(gauin(i,1).eq.undef) sum1=undef                                      
     if(gauin(i,jmxin).eq.undef) sum2=undef                                  
   enddo                                                                     
!                                                                               
   do i = 1,imxout                                                          
     if(abs(regl(1)).eq.90.) then                                              
       regout(i,     1)=sum1                                                   
     endif                                                                     
     if(abs(regl(jmxout)).eq.90.) then                                         
       regout(i,jmxout)=sum2                                                   
     endif                                                                     
   enddo 
!                                                                               
   return                                                                    
   end subroutine post_gauss2latlon                                                                       
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_gaussian_lat(gaul,k)       
!-------------------------------------------------------------------------------
   real, allocatable :: a(:)
   real              :: gaul(1)                                                         
!-------------------------------------------------------------------------------
!
   allocate(a(k))
   esp=1.e-6                                                                 
   c=(1.e0-(2.e0/3.14159265358979e0)**2)*0.25e0                              
   fk=k                                                                      
   kk=k/2                                                                    
!
   call dyn_bessel_function(a,kk)
!
   do is = 1,kk                                                             
     xz=cos(a(is)/sqrt((fk+0.5e0)**2+c))                                       
     iter=0                                                                    
     10 pkm2=1.e0                                                                 
     pkm1=xz                                                                   
     iter=iter+1                                                               
     if(iter.gt.10) go to 70                                                   
     do 20 n = 2,k                                                               
       fn=n                                                                      
       pk=((2.e0*fn-1.e0)*xz*pkm1-(fn-1.e0)*pkm2)/fn                             
       pkm2=pkm1                                                                 
     20 pkm1=pk                                                                   
     pkm1=pkm2                                                                 
     pkmrk=(fk*(pkm1-xz*pk))/(1.e0-xz**2)                                      
     sp=pk/pkmrk                                                               
     xz=xz-sp                                                                  
     avsp=abs(sp)                                                              
     if(avsp.gt.esp) go to 10                                                  
     a(is)=xz                                                                  
   enddo
!
   if(k.eq.kk*2) go to 50                                                    
!
   a(kk+1)=0.e0                                                              
   pk=2.e0/fk**2                                                             
   do n = 2,k,2                                                             
     fn=n                                                                      
     pk=pk*fn**2/(fn-1.e0)**2                                                  
   enddo
!
   50 continue                                                                  
!
   do n = 1,kk                                                              
     l=k+1-n                                                                   
     a(l)=-a(n)                                                                
   enddo 
!                                                                               
   radi=180./(4.*atan(1.))                                                   
   do n = 1,k                                                              
     gaul(n)=90.-acos(a(n))*radi                                               
   enddo 
!                                                                               
   print *,'gaussian lat (deg) for jmax=',k                                  
   print *,(gaul(n),n=1,k)                                                   
!                                                                               
   deallocate(a)
   return                                                                    
   70 write(6,6000)                                                             
   6000 format(5x,14herror in gauaw)                                          
   stop                                                                      
!
   end subroutine chgr_gaussian_lat
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine dyn_bessel_function(bes,n)
!-------------------------------------------------------------------------------
   real  ::  bes(n)                                                          
   real  ::  bz(50)                                                          
!                                                                               
   data pi /3.14159265358979e0/                                               
   data bz / 2.4048255577e0, 5.5200781103e0,                                   & 
        8.6537279129e0,11.7915344391e0,14.9309177086e0,18.0710639679e0,        & 
        21.2116366299e0,24.3524715308e0,27.4934791320e0,30.6346064684e0,       &  
        33.7758202136e0,36.9170983537e0,40.0584257646e0,43.1997917132e0,       &  
        46.3411883717e0,49.4826098974e0,52.6240518411e0,55.7655107550e0,       &  
        58.9069839261e0,62.0484691902e0,65.1899648002e0,68.3314693299e0,       &  
        71.4729816036e0,74.6145006437e0,77.7560256304e0,80.8975558711e0,       &  
        84.0390907769e0,87.1806298436e0,90.3221726372e0,93.4637187819e0,       &  
        96.6052679510e0,99.7468198587e0,102.888374254e0,106.029930916e0,       &  
        109.171489649e0,112.313050280e0,115.454612653e0,118.596176630e0,       &  
        121.737742088e0,124.879308913e0,128.020877005e0,131.162446275e0,       &  
        134.304016638e0,137.445588020e0,140.587160352e0,143.728733573e0,       &  
        146.870307625e0,150.011882457e0,153.153458019e0,156.295034268e0/         
!-------------------------------------------------------------------------------
!
   nn=n                                                                      
   if(n.le.50) go to 12
!
   bes(50)=bz(50)
   do j = 51,n
     bes(j)=bes(j-1)+pi
   enddo
   nn=49
!
   12 do 15 j=1,nn
   15 bes(j)=bz(j)
!
   return           
   end subroutine dyn_bessel_function 
!-------------------------------------------------------------------------------
