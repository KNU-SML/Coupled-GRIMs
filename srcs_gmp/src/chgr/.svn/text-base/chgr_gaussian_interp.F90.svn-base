#include <define.h>
   subroutine chgr_gaussian_interp
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [chgr_gaussian_interp]
!      |
!      |--- [chgr_interp_gaussian_sph] *
!      |--- [chgr_interp_gaussian_dfs] *
!      |--- [chgr_gaussian_lat] *
!      |--- [dyn_bessel_function] *
!
!-------------------------------------------------------------------------------
   end subroutine chgr_gaussian_interp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_interp_gaussian_sph(gauin,imxin,jmxin,gauout,               &
                                       imxout,jmxout)  
!-------------------------------------------------------------------------------
!
! abstract: interpolation from gaussian grid to other gaussian grid   
!
!-------------------------------------------------------------------------------
   integer  ::  imxin,jmxin,imxout,jmxout
!
   real     ::  gauin (imxin ,jmxin )                                           
   real     ::  gauout(imxout,jmxout)                                           
!
   real,    allocatable, dimension(:), save  :: gaulin,gaulou,ddx,ddy
   integer, allocatable, dimension(:), save  :: iindx1,iindx2,jindx1,jindx2
!                                                                               
   data ifp/0/                                                               
!-------------------------------------------------------------------------------                                                         
   if(imxin.eq.imxout.and.jmxin.eq.jmxout) then                              
     do j = 1,jmxin                                                          
       do i = 1,imxin                                                          
            gauout(i,j)=gauin(i,j)                                             
       enddo
     enddo
     return                                                                  
   endif                                                                     
!                                                                             
   if(ifp.eq.0) then
     allocate(ddx(imxout),ddy(jmxout))
     allocate(iindx1(imxout),iindx2(imxout))
     allocate(jindx1(jmxout),jindx2(jmxout))
     allocate(gaulin(jmxin),gaulou(jmxout))
!
     ifp=1                                                                   
!                                                                           
     call chgr_gaussian_lat(gaulin,jmxin)  
     do j = 1,jmxin                                                        
       gaulin(j)=90.-gaulin(j)                                            
     enddo
!                                                                               
     call chgr_gaussian_lat(gaulou,jmxout) 
     do j = 1,jmxout                      
       gaulou(j)=90.-gaulou(j)           
     enddo
!                                       
     dxin =360./float(imxin )          
     dxout=360./float(imxout)         
!                                                                               
     do i = 1,imxout                 
       alamd=float(i-1)*dxout       
       i1=alamd/dxin+1.001         
       iindx1(i)=i1               
       i2=i1+1                   
       if(i2.gt.imxin) i2=1     
       iindx2(i)=i2            
!      print *,'i,i1,i2',i,i1,i2 
       ddx(i)=(alamd-float(i1-1)*dxin)/dxin  
     enddo
!                                                                               
     j2=1                                   
     do j = 1,jmxout                       
       aphi=gaulou(j)                     
       do 50 jj=1,jmxin                  
!
         if(aphi.lt.gaulin(jj)) go to 50
!
         j2=jj                         
!
         go to 42                     
!
       50 continue                   
       42 continue                  
!
       if(j2.gt.2) go to 43        
!
       j1=1                       
       j2=2                      
!
       go to 44                 
!
       43 continue             
!
       if(j2.le.jmxin) go to 45 
!
       j1=jmxin-1              
       j2=jmxin               
!
       go to 44              
!
       45 continue          
       j1=j2-1             
       44 continue        
!      print *,'j,j1,j2',j,j1,j2     
       jindx1(j)=j1      
       jindx2(j)=j2     
       ddy(j)=(aphi-gaulin(j1))/(gaulin(j2)-gaulin(j1))  
     enddo
!                                                                               
!    print *,'iindx1'                                                          
!    print *,(iindx1(n),n=1,imxout)                                            
!    print *,'iindx2'                                                          
!    print *,(iindx2(n),n=1,imxout)                                            
!    print *,'jindx1'                                                          
!    print *,(jindx1(n),n=1,jmxout)                                            
!    print *,'jindx2'                                                          
!    print *,(jindx2(n),n=1,jmxout)                                            
!    print *,'ddy'                                                             
!    print *,(ddy(n),n=1,jmxout)                                               
!    print *,'ddx'                                                             
!    print *,(ddx(n),n=1,jmxout)                                               
!                                                                               
   endif         ! ifp = 0
!                                                                               
   if(imxin.gt.imxout.or.jmxin.gt.jmxout) then
     do j = 1,jmxout
       j1=jindx1(j)
       if (ddy(j).gt.0.5) j1=jindx2(j)
       do i = 1,imxout
         i1=iindx1(i)
         if (ddx(i).gt.0.5) i1=iindx2(i)
         gauout(i,j)=gauin(i1,j1)
       enddo
     enddo
   else
     do j = 1,jmxout                                    
       y=ddy(j)                                        
       j1=jindx1(j)                                   
       j2=jindx2(j)                                  
       do i = 1,imxout                             
         x=ddx(i)                                   
         i1=iindx1(i)                            
         i2=iindx2(i)                             
         gauout(i,j)=(1.-x)*(1.-y)*gauin(i1,j1)+(1.-y)*x*gauin(i2,j1)+         &
                     (1.-x)*y*gauin(i1,j2)+x*y*gauin(i2,j2)  
       enddo
     enddo
   endif
!
   return                                                                    
   end subroutine chgr_interp_gaussian_sph                  
!
!-------------------------------------------------------------------------------
   subroutine chgr_interp_gaussian_dfs(gauin,imxin,jmxin,gauout,imxout,jmxout)
!-------------------------------------------------------------------------------
!
! abstract: interpolation from gaussian grid to other gaussian grid  
!
!-------------------------------------------------------------------------------
   real, allocatable, dimension(:), save    :: ddx,ddy
   integer, allocatable, dimension(:), save :: iindx1,iindx2,jindx1,jindx2
!
   real     ::  gauin (imxin ,jmxin )                                           
   real     ::  gauout(imxout,jmxout)
   real     ::  gaulin(jmxin),gaulou(jmxout)                                   
   integer  ::  imxin ,jmxin, imxout,jmxout
!                                                                               
   data ifp/0/                                                               
!-------------------------------------------------------------------------------     
   if (ifp.eq.0) then
     ifp=1                                                                     
     allocate(ddx(imxout),ddy(jmxout))
     allocate(iindx1(imxout),iindx2(imxout))
     allocate(jindx1(jmxout),jindx2(jmxout))
!                                                                               
     call chgr_gaussian_lat(gaulin,jmxin)                                     
     do j = 1,jmxin                                                           
       gaulin(j)=90.-gaulin(j)      ! real latitude
     enddo
!
!    grid start from N.H
!
     dlat=180./jmxout
     do j = 1,jmxout                                                          
       gaulou(j)=90.-(j-0.5)*dlat
     enddo
!                                                                               
     dxin =360./float(imxin )                                                  
     dxout=360./float(imxout)                                                  
!                                                                               
     do i = 1,imxout                                                         
       alamd=float(i-1)*dxout                                               
       i1=alamd/dxin+1.001                                                 
       iindx1(i)=i1                                                       
       i2=i1+1                                                           
       if(i2.gt.imxin) i2=1                                             
       iindx2(i)=i2                                                    
       ddx(i)=(alamd-float(i1-1)*dxin)/dxin                           
     enddo
!                                                                               
     j2=1
     do j = 1,jmxout                                                 
       aphi=gaulou(j)                                               
       do 50 jj=1,jmxin                                            
!
         if(aphi.lt.gaulin(jj)) go to 50                          
!
         j2=jj                                                   
         exit
50     continue                                                 
       if(j2.gt.2) go to 43                                    
       j1=1                                                   
       j2=2                                                  
!
       go to 44                                             
!
43     continue                                            
       if(j2.le.jmxin) go to 45                           
       j1=jmxin-1                                        
       j2=jmxin                                         
!
       go to 44                                        
!
45     continue                                       
       j1=j2-1                                       
44     continue                                     
!      print *,'j,j1,j2',j,j1,j2                   
       jindx1(j)=j1                               
       jindx2(j)=j2                              
       ddy(j)=(aphi-gaulin(j1))/(gaulin(j2)-gaulin(j1))
     enddo
!                                                                               
!    print *,'iindx1'                                                          
!    print *,(iindx1(n),n=1,imxout)                                            
!    print *,'iindx2'                                                          
!    print *,(iindx2(n),n=1,imxout)                                            
!    print *,'jindx1'                                                          
!    print *,(jindx1(n),n=1,jmxout)                                            
!    print *,'jindx2'                                                          
!    print *,(jindx2(n),n=1,jmxout)                                            
!    print *,'ddy'                                                             
!    print *,(ddy(n),n=1,jmxout)                                               
!    print *,'ddx'                                                             
!    print *,(ddx(n),n=1,jmxout)                                               
!
   endif         ! first call
!
   do j = 1,jmxout                                                          
     y=ddy(j)                                                                  
     j1=jindx1(j)                                                              
     j2=jindx2(j)                                                              
     do i = 1,imxout                                                          
       x=ddx(i)                                 
       i1=iindx1(i)                            
       i2=iindx2(i)                           
       gauout(i,j)=(1.-x)*(1.-y)*gauin(i1,j1)+(1.-y)*x*gauin(i2,j1)+           &
                     (1.-x)*y*gauin(i1,j2)+x*y*gauin(i2,j2)
     enddo
   enddo
!
   return                                                                    
   end subroutine chgr_interp_gaussian_dfs
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_gaussian_lat(gaul,k)                    
!-------------------------------------------------------------------------------
   save                                                                      
!                                                                               
   real, allocatable  ::  a(:)
   real               ::  gaul(k)            
!-------------------------------------------------------------------------------
   allocate(a(k))
!
   esp=1.e-14                                                                
   c=(1.e0-(2.e0/3.14159265358979e0)**2)*0.25e0                              
   fk=k                                                                      
   kk=k/2                                                                    
   call dyn_bessel_function(a,kk)   
   do is = 1,kk                                                             
     xz=cos(a(is)/sqrt((fk+0.5e0)**2+c))                                       
     iter=0                                                                    
     10 pkm2=1.e0                           
     pkm1=xz                                                                   
     iter=iter+1                                                               
!
     if(iter.gt.10) go to 70                                                   
!
     do n = 2,k                                                               
       fn=n                                
       pk=((2.e0*fn-1.e0)*xz*pkm1-(fn-1.e0)*pkm2)/fn 
       pkm2=pkm1                                    
       pkm1=pk                                     
     enddo
!
     pkm1=pkm2                                                                 
     pkmrk=(fk*(pkm1-xz*pk))/(1.e0-xz**2)                                      
     sp=pk/pkmrk                                                               
     xz=xz-sp                                                                  
     avsp=abs(sp)                                                              
!
     if(avsp.gt.esp) go to 10                                                  
!
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
   50 continue                                                                  
!
   do n = 1,kk                                                              
     l=k+1-n                                                                   
     a(l)=-a(n)                                                                
   enddo
!                                                                               
   radi=180.e0/(4.e0*atan(1.e0))                                             
   do n = 1,k                                                              
     gaul(n)=acos(a(n))*radi                                                   
   enddo
   print *,'gaussian lat (deg) for jmax=',k                                  
   print *,(gaul(n),n=1,k)                                                   
!                                                                               
   deallocate(a)
   return                                                                    
   70 write(6,6000)                                                             
   6000 format(5x,14herror in gauaw)                                          
!
   stop
   end subroutine chgr_gaussian_lat               
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine dyn_bessel_function(bes,n)  
!-------------------------------------------------------------------------------
   save                                                                      
!
   real  ::  bes(n)                                                          
   real  ::  bz(50)                                                          
!                                                                               
   data pi/3.14159265358979e0/                                               
   data bz         / 2.4048255577e0, 5.5200781103e0,                           &
       8.6537279129e0,11.7915344391e0,14.9309177086e0,18.0710639679e0,         &
      21.2116366299e0,24.3524715308e0,27.4934791320e0,30.6346064684e0,         &
      33.7758202136e0,36.9170983537e0,40.0584257646e0,43.1997917132e0,         &
      46.3411883717e0,49.4826098974e0,52.6240518411e0,55.7655107550e0,         &
      58.9069839261e0,62.0484691902e0,65.1899648002e0,68.3314693299e0,         &
      71.4729816036e0,74.6145006437e0,77.7560256304e0,80.8975558711e0,         &
      84.0390907769e0,87.1806298436e0,90.3221726372e0,93.4637187819e0,         &
      96.6052679510e0,99.7468198587e0,102.888374254e0,106.029930916e0,         & 
      109.171489649e0,112.313050280e0,115.454612653e0,118.596176630e0,         &
      121.737742088e0,124.879308913e0,128.020877005e0,131.162446275e0,         &
      134.304016638e0,137.445588020e0,140.587160352e0,143.728733573e0,         &
      146.870307625e0,150.011882457e0,153.153458019e0,156.295034268e0/         
!-------------------------------------------------------------------------------
   nn=n                                                                      
   if(n.le.50) go to 12                                                      
   bes(50)=bz(50)                                                            
   do j = 51,n                                                               
     bes(j)=bes(j-1)+pi                                                        
   enddo
   nn=49                                                                     
   12 do j = 1,nn                                                              
     bes(j)=bz(j)                                                              
   enddo
!
   return                                                                    
   end subroutine dyn_bessel_function
!-------------------------------------------------------------------------------
