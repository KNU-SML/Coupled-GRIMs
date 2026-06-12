#include <define.h>
   subroutine post_interp_sigma
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [post_interp_sigma]
!           |
!           |-- [post_sigma2pressure] *
!           |-- [post_sigma2tropopause] *
!
!-------------------------------------------------------------------------------
   end subroutine post_interp_sigma
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_sigma2pressure(io,jo,im,ix,km,si,sl,                        &
                    ps,us,vs,os,zs,zi,ts,rs,qs,qcs,qrs,qis,qss,qgs,            &
                    nccns,ncs,nrs,                                             &
                    o3s,tkes,ko,po,                                            &
                    lpqc,lpqr,lpqi,lpqs,lpqg,lpo3,lptke,lpnccn,lpnc,lpnr,      &
#ifdef DCMIP
                    q1s,q2s,q3s,q4s,                                           &
                    lpq1,lpq2,lpq3,lpq4,                                       &
                    q1p,q2p,q3p,q4p,                                           &
#endif
#ifdef DFS
                    vos,vop,                                                   &
#endif
                    up,vp,op,zp,tp,rp,qp,qcp,qrp,qip,qsp,qgp,                  &
                    nccnp,ncp,nrp,o3p,tkep)
!-------------------------------------------------------------------------------
   use paramodel, only : ko_,levs_
!-------------------------------------------------------------------------------
!
!  subprogram:    post_sigma2pressure       sigma to pressure interpolation
!
! abstract: interpolates winds, omega, height, temperature and humidity
!   from the sigma coordinate system to the mandatory pressure levels.
!   assumes that relative humidity, temperature, geopotential heights,
!   wind components and vertical velocity vary linearly in the vertical
!   with the log of pressure.  underground heights are obtained using
!   the shuell method and underground temperatures are obtained using
!   a constant moist adiabatic lapse rate.  heights above the top sigma
!   level are integrated hydrostatically.  otherwise fields are held
!   constant outside the sigma structure and no extrapolation is done.
!
! program history log:
!   1992-10-31  sela,newell,gerrity,ballish,deaven,iredell
!   2000-03-09  songyou hong           cvs verion setup
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call post_sigma2pressure(im,ix,km,si,sl,
!    &                 ps,us,vs,os,zs,zi,ts,rs,qs,
!    &                 ko,po,up,vp,op,zp,tp,rp,sp)
!
!   input argument list:
!     im       - integer number of points
!     ix       - integer first dimension of upper air data
!     km       - integer number of sigma levels
!     si       - real (km+1) sigma interface values
!     sl       - real (km) sigma values
!     ps       - real (im) surface pressure in kpa
!     us       - real (ix,km) zonal wind in m/s
!     vs       - real (ix,km) merid wind in m/s
!     os       - real (ix,km) vertical velocity in pa/s
!     zs       - real (ix,km) heights on the full levels in m
!     zi       - real (ix,km) heights on the interfaces in m
!     ts       - real (ix,km) temperature in k
!     rs       - real (ix,km) relative humidity in percent
!     qs       - real (ix,km) specific humidity in kg/kg
!     ko       - integer number of pressure levels
!     po       - real (ko) mandatory pressures in kpa
!
!   output argument list:
!     up       - real (ix,ko) zonal wind in m/s
!     vp       - real (ix,ko) merid wind in m/s
!     op       - real (ix,ko) vertical velocity in pa/s
!     zp       - real (ix,ko) heights in m
!     tp       - real (ix,ko) temperature in k
!     rp       - real (ix,ko) relative humidity in percent
!
! subprograms called:
!   isrchfltx - find first value in an array less than target value
!
!-------------------------------------------------------------------------------
!
   real             ::  si(ix,km+1),sl(ix,km),ps(im)
   real             ::  us(ix,km),vs(ix,km),os(ix,km)
   real             ::  zs(ix,km),zi(ix,km),ts(ix,km),rs(ix,km),qs(ix,km)
   real             ::  o3s(ix,km),qcs(ix,km),qrs(ix,km),qis(ix,km),           &
                        qss(ix,km),qgs(ix,km),tkes(ix,km),                     &
                        nccns(ix,km),ncs(ix,km),nrs(ix,km)
   real             ::  po(ko)
   logical lpo3, lptke, lpqc, lpqr, lpqi, lpqs, lpqg
#ifdef DCMIP
   real             ::  q1s(ix,km),q2s(ix,km),q3s(ix,km),q4s(ix,km)
   logical lpq1,lpq2,lpq3,lpq4
   real             ::  q1p(ix,ko),q2p(ix,ko),q3p(ix,ko),q4p(ix,ko)
#endif
   logical lpnccn, lpnc, lpnr
#ifdef DFS
   real             ::  vos(ix,km),vop(ix,ko)
#endif
   real             ::  up(ix,ko),vp(ix,ko),op(ix,ko)
   real             ::  zp(ix,ko),tp(ix,ko),rp(ix,ko)
   real             ::  qp(ix,ko),o3p(ix,ko),tkep(ix,ko)
   real             ::  qcp(ix,ko),qrp(ix,ko),qip(ix,ko),qsp(ix,ko),qgp(ix,ko)
   real             ::  nccnp(ix,ko),ncp(ix,ko),nrp(ix,ko)
   real             ::  sp(2*io+6,ko_)
   real             ::  asi(levs_),asl(levs_),apo(ko_),aps(2*io)
   real, parameter  ::  g= 9.8000e+0 ,rd= 2.8705e+2 ,rv= 4.6150e+2
   real, parameter  ::  rog=rd/g,fvirt=rv/rd-1.
   real, parameter  ::  gammam=-6.5e-3,zshul=75.,tvshul=290.66
!-------------------------------------------------------------------------------
!
!  compute log pressures for interpolation
!
   do i = 1,im
     do k = 1,km
       asi(k)=log(si(i,k))
       asl(k)=log(sl(i,k))
     enddo
     do k = 1,ko
       apo(k)=log(po(k))
     enddo
     aps(i)=log(ps(i))
!
!  determine sigma layers bracketing pressure layer
!  and interpolate to obtain real sigma layer number
!
     kd=1
     do k = 1,ko
       ask=apo(k)-aps(i)
       kd=kd+isrchfltx(km-kd-1,asl(kd+1),1,ask)-1
       sp(i,k)=kd+(asl(kd)-ask)/(asl(kd)-asl(kd+1))
     enddo
!
!  interpolate sigma to pressure
!
     do k = 1,ko
       ask=apo(k)-aps(i)
!
!  below ground use shuell method to obtain height, constant lapse rate
!  to obtain temperature, and hold other fields constant
!
       if(ask.gt.0.) then
#ifdef DFS
         vop(i,k)=vos(i,1)
#endif
         up(i,k)=us(i,1)
         vp(i,k)=vs(i,1)
         op(i,k)=os(i,1)
         tvsf=ts(i,1)*(1.+fvirt*qs(i,1))-gammam*(zs(i,1)-zi(i,1))
         if(zi(i,1).gt.zshul) then
           tvsl=tvsf-gammam*zi(i,1)
           if(tvsl.gt.tvshul) then
             if(tvsf.gt.tvshul) then
               tvsl=tvshul-5.e-3*(tvsf-tvshul)**2
             else
               tvsl=tvshul
             endif
           endif
           gammas=(tvsf-tvsl)/zi(i,1)
         else
           gammas=0.
         endif
         part=rog*ask
         zp(i,k)=zi(i,1)-tvsf*part/(1.+0.5*gammas*part)
         tp(i,k)=ts(i,1)+gammam*(zp(i,k)-zs(i,1))
         rp(i,k)=rs(i,1)
         qp(i,k)=qs(i,1)
         if(lpqc) qcp(i,k)=qcs(i,1)
         if(lpqr) qrp(i,k)=qrs(i,1)
         if(lpqi) qip(i,k)=qis(i,1)
         if(lpqs) qsp(i,k)=qss(i,1)
         if(lpqg) qgp(i,k)=qgs(i,1)
         if(lpnccn) nccnp(i,k)=nccns(i,1)
         if(lpnc) ncp(i,k)=ncs(i,1)
         if(lpnr) nrp(i,k)=nrs(i,1)
         if(lpo3) o3p(i,k)=o3s(i,1)
         if(lptke) tkep(i,k)=tkes(i,1)
#ifdef DCMIP
         if(lpq1) q1p(i,k)=q1s(i,1)
         if(lpq2) q2p(i,k)=q2s(i,1)
         if(lpq3) q3p(i,k)=q3s(i,1)
         if(lpq4) q4p(i,k)=q4s(i,1)
#endif
!
!  above top sigma ground integrate height hydrostatically
!  and hold other fields constant
!
       elseif(sp(i,k).ge.km) then ! ask
#ifdef DFS
         vop(i,k)=vos(i,km)
#endif
         up(i,k)=us(i,km)
         vp(i,k)=vs(i,km)
         op(i,k)=os(i,km)
         tvkm=ts(i,km)*(1.+fvirt*qs(i,km))
         zp(i,k)=zs(i,km)+rog*tvkm*(asl(km)-ask)
         tp(i,k)=ts(i,km)
         rp(i,k)=rs(i,km)
         qp(i,k)=qs(i,km)
         if(lpqc) qcp(i,k)=qcs(i,km)
         if(lpqr) qrp(i,k)=qrs(i,km)
         if(lpqi) qip(i,k)=qis(i,km)
         if(lpqs) qsp(i,k)=qss(i,km)
         if(lpqg) qgp(i,k)=qgs(i,km)
         if(lpnccn) nccnp(i,k)=nccns(i,km)
         if(lpnc) ncp(i,k)=ncs(i,km)
         if(lpnr) nrp(i,k)=nrs(i,km)
         if(lpo3) o3p(i,k)=o3s(i,km)
         if(lptke) tkep(i,k)=tkes(i,km)
#ifdef DCMIP
         if(lpq1) q1p(i,k)=q1s(i,km)
         if(lpq2) q2p(i,k)=q2s(i,km)
         if(lpq3) q3p(i,k)=q3s(i,km)
         if(lpq4) q4p(i,k)=q4s(i,km)
#endif
!
!  within sigma structure, interpolate fields linearly in log pressure
!  between bracketing full sigma layers except heights are interpolated
!  between the nearest full sigma layer and the nearest sigma interface
!
       else !ask
         kd=max(int(sp(i,k)),1)
         ku=kd+1
         wu=sp(i,k)-kd
         wd=1.-wu
#ifdef DFS
         vop(i,k)=wu*vos(i,ku)+wd*vos(i,kd)
#endif
         up(i,k)=wu*us(i,ku)+wd*us(i,kd)
         vp(i,k)=wu*vs(i,ku)+wd*vs(i,kd)
         op(i,k)=wu*os(i,ku)+wd*os(i,kd)
         ki=int(sp(i,k))+1
         di=asi(ki)-ask
         kl=nint(ki-0.5+sign(0.5,di))
         wl=di/(asi(ki)-asl(kl))
         wi=1.-wl
         zp(i,k)=wi*zi(i,ki)+wl*zs(i,kl)
         tp(i,k)=wu*ts(i,ku)+wd*ts(i,kd)
         rp(i,k)=wu*rs(i,ku)+wd*rs(i,kd)
         qp(i,k)=wu*qs(i,ku)+wd*qs(i,kd)
         if(lpqc) qcp(i,k)=wu*qcs(i,ku)+wd*qcs(i,kd)
         if(lpqr) qrp(i,k)=wu*qrs(i,ku)+wd*qrs(i,kd)
         if(lpqi) qip(i,k)=wu*qis(i,ku)+wd*qis(i,kd)
         if(lpqs) qsp(i,k)=wu*qss(i,ku)+wd*qss(i,kd)
         if(lpqg) qgp(i,k)=wu*qgs(i,ku)+wd*qgs(i,kd)
         if(lpnccn) nccnp(i,k)=wu*nccns(i,ku)+wd*nccns(i,kd)
         if(lpnc) ncp(i,k)=wu*ncs(i,ku)+wd*ncs(i,kd)
         if(lpnr) nrp(i,k)=wu*nrs(i,ku)+wd*nrs(i,kd)
         if(lpo3) o3p(i,k)=wu*o3s(i,ku)+wd*o3s(i,kd)
         if(lptke) tkep(i,k)=wu*tkes(i,ku)+wd*tkes(i,kd)
#ifdef DCMIP
         if(lpq1) q1p(i,k)=wu*q1s(i,ku)+wd*q1s(i,kd)
         if(lpq2) q2p(i,k)=wu*q2s(i,ku)+wd*q2s(i,kd)
         if(lpq3) q3p(i,k)=wu*q3s(i,ku)+wd*q3s(i,kd)
         if(lpq4) q4p(i,k)=wu*q4s(i,ku)+wd*q4s(i,kd)
#endif
       endif
     enddo ! k
   enddo ! i
! 
   return
   end subroutine post_sigma2pressure
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_sigma2tropopause(im,ix,km,sl,ps,u,v,t,                      &
                                    ptp,utp,vtp,ttp,shtp,stp)   
!-------------------------------------------------------------------------------
   use paramodel, only : levs_
!-------------------------------------------------------------------------------
!
!  subprogram:    post_sigma2tropopause      sigma to tropopause interpolation                 
!                                                                               
! abstract: locates the tropopause pressure level and interpolates              
!   the winds and temperature and wind shear to the tropopause.                 
!   the tropopause is identified by the lowest level above 450 mb               
!   where the temperature lapse rate -dt/dz becomes less than 2 k/km.           
!   the tropopause is not allowed higher than 85 mb.                            
!   interpolations are done linearly in log of pressure.                        
!                                                                               
! program history log:                                                          
!   1992-10-31  mccalla,iredell                                                   
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!                                                                               
! usage:    call post_sigma2tropopause(im,ix,km,sl,   
!    &                  ps,u,v,t,                                               
!    &                  ptp,utp,vtp,ttp,shtp,stp)                               
!                                                                               
!   input argument list:                                                        
!     im       - integer number of points                                       
!     ix       - integer first dimension of upper air data                      
!     km       - integer number of sigma levels                                 
!     sl       - real (km) sigma values                                         
!     ps       - real (im) surface pressure in kpa                              
!     u        - real (ix,km) zonal wind in m/s                                 
!     v        - real (ix,km) merid wind in m/s                                 
!     t        - real (ix,km) temperature in k                                  
!                                                                               
!   output argument list:                                                       
!     ptp      - real (im) tropopause pressure in kpa                           
!     utp      - real (im) tropopause zonal wind in m/s                         
!     vtp      - real (im) tropopause merid wind in m/s                         
!     ttp      - real (im) tropopause temperature in k                          
!     shtp     - real (im) tropopause wind speed shear in (m/s)/m               
!     stp      - real (im) tropopause sigma layer number                        
!                                                                               
!-------------------------------------------------------------------------------
   real             ::  sl(km),ps(im)                                                   
   real             ::  u(ix,km),v(ix,km),t(ix,km)                                      
   real             ::  ptp(im),utp(im),vtp(im),ttp(im),shtp(im),stp(im)                
   real             ::  asl(levs_)                                                      
   real, parameter  ::  g= 9.8000e+0 ,rd= 2.8705e+2                                    
   real, parameter  ::  rog=rd/g                                                       
   real, parameter  ::  ptbot=450.e-1,pttop=85.e-1,gamt=2.e-3                          
!-------------------------------------------------------------------------------
!
   fgamma(k)=(t(i,k-1)-t(i,k+1))/(rog*t(i,k)*(asl(k-1)-asl(k+1)))            
!
!  identify tropopause as first layer above ptbot but below pttop               
!  where the temperature lapse rate drops below gamt                            
!    stp is real interpolated sigma layer number of tropopause                  
!
   do k = 1,km                                                                 
     asl(k)=log(sl(k))                                                         
   enddo                                                                     
!
   do i = 1,im                                                                 
     k=3                                                                     
     pu=ps(i)*sl(k)                                                          
     do while(k.lt.km-1.and.pu.gt.ptbot)                                      
       k=k+1                                                                 
       pu=ps(i)*sl(k)                                                        
     enddo                                                                   
     gamd=fgamma(k-1)                                                        
     gamd=max(gamd,gamt)                                                     
     gamu=fgamma(k)                                                          
     do while(k.lt.km-1.and.pu.gt.pttop.and.gamu.gt.gamt)                     
       k=k+1                                                                 
       pu=ps(i)*sl(k)                                                        
       gamd=gamu                                                             
       gamu=fgamma(k)                                                        
     enddo                                                                   
     gamu=min(gamu,gamt)                                                     
     stp(i)=k-(gamt-gamu)/(gamd-gamu)                                        
   enddo                                                                     
! 
!  interpolate tropopause pressure, temperature, winds and wind shear           
!  tropopause pressure is constrained to be between ptbot and pttop             
!
   do i = 1,im                                                                 
     kd=stp(i)                                                               
     ku=kd+1                                                                 
     wu=stp(i)-kd                                                            
     dlp=asl(ku)-asl(kd)                                                     
     ptp(i)=ps(i)*sl(kd)*exp(wu*dlp)                                         
     if(ptp(i).gt.ptbot) then                                                
       wu=wu+log(ptbot/ptp(i))/dlp                                           
       ptp(i)=ptbot                                                          
     elseif(ptp(i).lt.pttop) then                                            
       wu=wu+log(pttop/ptp(i))/dlp                                           
       ptp(i)=pttop                                                          
     endif                                                                   
     ttp(i)=t(i,kd)+wu*(t(i,ku)-t(i,kd))                                     
     utp(i)=u(i,kd)+wu*(u(i,ku)-u(i,kd))                                     
     vtp(i)=v(i,kd)+wu*(v(i,ku)-v(i,kd))                                     
     spdd=sqrt(u(i,kd)**2+v(i,kd)**2)                                        
     spdu=sqrt(u(i,ku)**2+v(i,ku)**2)                                        
     shtp(i)=(spdu-spdd)/(rog*0.5*(t(i,ku)+t(i,kd))*dlp)                     
   enddo                                                                     
! 
   return                                                                    
   end subroutine post_sigma2tropopause
