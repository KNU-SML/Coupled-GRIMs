#include <define.h>
   subroutine post_pgb_2d
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ... 
!
!    [post_pgb_2d] *  1D (vs. in fcst : 2D) 
!         |
!         |-- [post_lifting_index] * 
!         |-- [post_max_wind] *
!                  |
!          (post_interp_spline) 
!                  |-- [post_spline_derivative]
!                  |     : compute 2nd derivatives for cubic spline
!                  |-- [post_spline_max]
!                        : determine maximum value of cubic spline
!
!-------------------------------------------------------------------------------
   end subroutine post_pgb_2d
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_lifting_index(io,im,ix,km,sl,ps,t,q,tm,tli)    
!-------------------------------------------------------------------------------
   use paramodel, only : levs_
!-------------------------------------------------------------------------------
!
!  subprogram:    post_lifting_index      compute best lifted index from sigma              
!                                                                               
! abstract: computes the best lifted index from profiles in sigma.              
!   the best lifted index is here computed by finding the parcel                
!   below sigma 0.80 with the warmest equivalent potential temperature,         
!   then raising it to 500 mb and subtracting its parcel temperature            
!   from the environment temperature.                                           
!                                                                               
! program history log:                                                          
!    1992-10-31  iredell                development
!    2000-03-09  songyou hong           cvs verion setup, prognostic clouds
!    2009-10-01  jung-eun kim           f90 format with standard physics modules
!    2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call post_lifting_index(im,ix,km,sl,ps,t,q,tm,tli)  
!                                                                               
!   input argument list:                                                        
!     im       - integer number of points                                       
!     ix       - integer first dimension of upper air data                      
!     km       - integer number of sigma levels                                 
!     sl       - real (km) sigma values                                         
!     ps       - real (im) surface pressure in kpa                              
!     t        - real (ix,km) temperature in k                                  
!     q        - real (ix,km) specific humidity in kg/kg                        
!     tm       - real (im) 500 mb temperature in k                              
!                                                                               
!   output argument list:                                                       
!     tli      - real (ix,km) best lifted index in k                            
!                                                                               
! subprograms called:                                                           
!   (fpkap)   - function to compute pressure to the kappa                       
!   (ftdp)    - function to compute dewpoint temperature                        
!   (ftlcl)   - function to compute lifting condensation level                  
!   (fthe)    - function to compute equivalent potential temperature            
!   (ftma)    - function to compute moist adiabat temperature                   
!                                                                               
!-------------------------------------------------------------------------------
   integer          ::  io
   real             ::  sl(km),ps(im),t(ix,km),q(ix,km),tm(im),tli(im)                  
   real, parameter  ::  cp= 1.0046e+3 ,rd= 2.8705e+2 ,rv= 4.6150e+2
   real, parameter  ::  rk=rd/cp,eps=rd/rv,epsm1=rd/rv-1.
   real, parameter  ::  slift=0.80,plift=50.
   real             ::  slk(levs_),psk(2*io),slkma(2*io),thema(2*io)                 
!-------------------------------------------------------------------------------
! 
!  initialize variables                                                         
!
   pliftk=(plift/100.)**rk                                                   
!
   do k = 1,km                                                                 
     slk(k)=sl(k)**rk                                                        
   enddo                                                                     
!
   do i = 1,im                                                                 
     psk(i)=fpkap(ps(i))                                                     
     slkma(i)=0.                                                             
     thema(i)=0.                                                             
   enddo                                                                     
!
!  select the warmest equivalent potential temperature                          
!  between the surface and sigma slift1                                         
!
   k=1                                                                       
   do while(sl(k).gt.slift)                                                   
     do i = 1,im                                                               
       p=sl(k)*ps(i)                                                         
       pv=p*q(i,k)/(eps-epsm1*q(i,k))                                        
       tdpd=max(t(i,k)-ftdp(pv),0.)                                          
       tlcl=ftlcl(t(i,k),tdpd)                                               
       slklcl=slk(k)*tlcl/t(i,k)                                             
       thelcl=fthe(tlcl,slklcl*psk(i))                                       
       if(thelcl.gt.thema(i)) then                                           
         slkma(i)=slklcl                                                     
         thema(i)=thelcl                                                     
       endif                                                                 
     enddo                                                                   
   k=k+1                                                                   
   enddo                                                                     
!
!  lift the parcel to 500 mb along a dry adiabat below the lcl                  
!  or along a moist adiabat above the lcl.                                      
!  the lifted index is the environment minus parcel temperature.                
!
   do i = 1,im                                                                 
     if(ps(i).gt.plift.and.thema(i).gt.0.) then                              
       slkp=pliftk/psk(i)                                                    
       slkc=min(slkp,slkma(i))                                               
       tlift=slkp/slkc*ftma(thema(i),slkc*psk(i),qma)                        
       tli(i)=tm(i)-tlift                                                    
     else                                                                    
       tli(i)=0.                                                             
     endif                                                                   
   enddo                                                                     
!
   return                                                                    
   end subroutine post_lifting_index                                                                      
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_max_wind(io,im,ix,km,sl,                                    &
                     ps,u,v,                                                   &
                     spdmw,umw,vmw,pmw,smw)                                  
!-------------------------------------------------------------------------------
   use paramodel, only : levs_
!-------------------------------------------------------------------------------
!                                                                               
!  subprogram:    post_max_wind      sigma to maxwind interpolation            
!                                                                               
! abstract: locates the maximum wind speed level (maxwind level) and            
!   returns the wind speed, components and pressure at that level.              
!   the maxwind level is restricted to be between 50kpa and 7kpa.               
!   the maxwind level is identified by cubic spline interpolation               
!   of the wind speeds in log pressure.                                         
!                                                                               
! program history log:                                                          
!   92-10-31  iredell                                                           
!                                                                               
! usage:    call post_max_wind(im,ix,km,sl,                    
!    &                  ps,u,v,                                                 
!    &                  spdmw,umw,vmw,pmw,smw)                                  
!                                                                               
!   input argument list:                                                        
!     im       - integer number of points                                       
!     ix       - integer first dimension of upper air data                      
!     km       - integer number of sigma levels                                 
!     sl       - real (km) sigma values                                         
!     ps       - real (im) surface pressure in kpa                              
!     u        - real (ix,km) zonal wind in m/s                                 
!     v        - real (ix,km) merid wind in m/s                                 
!                                                                               
!   output argument list:                                                       
!     spdmw    - real (im) maxwind wind speed in m/s                            
!     umw      - real (im) maxwind zonal wind in m/s                            
!     vmw      - real (im) maxwind merid wind in m/s                            
!     pmw      - real (im) maxwind pressure in kpa                              
!     smw      - real (im) maxwind sigma layer number                           
!                                                                               
! subprograms called:                                                           
!   post_spline_derivative    compute second derivatives for cubic spline    
!   post_spline_max           determine maximum value of cubic spline    
!                                                                               
!-------------------------------------------------------------------------------
   real             ::  sl(km),ps(im)                                                   
   real             ::  u(ix,km),v(ix,km)                                               
   real             ::  spdmw(im),umw(im),vmw(im),pmw(im),smw(im)                       
   real             ::  s(levs_),spd(2*io,levs_),d2spd(2*io,levs_)                    
   real, parameter  ::  pmwbot=500.e-1,pmwtop=70.e-1                                   
!-------------------------------------------------------------------------------
!
!  fix vertical coordinate proportional to log pressure                         
!  and calculate wind speeds between pmwbot and pmwtop                          
!
   do k = 1,km                                                                 
     s(k)=-log(sl(k))                                                        
   enddo                                                                     
!
   do k = 1,km                                                                 
     do i = 1,im                                                               
       p=sl(k)*ps(i)                                                         
       if(p.le.pmwbot.and.p.ge.pmwtop) then                                  
         spd(i,k)=sqrt(u(i,k)**2+v(i,k)**2)                                  
       else                                                                  
         spd(i,k)=0.                                                         
       endif                                                                 
     enddo                                                                   
   enddo                                                                     
!
!  use spline routines to determine maxwind level and wind speed                
!
   call post_spline_derivative(im,km,s,spd,d2spd)                   
   call post_spline_max(im,km,s,spd,d2spd,smw,pmw,spdmw)        
!
!  compute maxwind pressure and wind components                                 
!
   do i = 1,im                                                                 
     pmw(i)=exp(-pmw(i))*ps(i)                                               
     k=int(smw(i))                                                           
     if(float(k).eq.smw(i)) then                                             
       ub=u(i,k)                                                             
       vb=v(i,k)                                                             
     else                                                                    
       ub=(k+1-smw(i))*u(i,k)+(smw(i)-k)*u(i,k+1)                            
       vb=(k+1-smw(i))*v(i,k)+(smw(i)-k)*v(i,k+1)                            
     endif                                                                   
     spdb=sqrt(ub**2+vb**2)                                                  
     umw(i)=ub*spdmw(i)/spdb                                                 
     vmw(i)=vb*spdmw(i)/spdb                                                 
   enddo                                                                     
!
   return                                                                    
   end subroutine post_max_wind                                                                      
!-------------------------------------------------------------------------------
