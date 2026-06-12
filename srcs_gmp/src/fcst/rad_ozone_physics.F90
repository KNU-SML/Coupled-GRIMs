#include <define.h>
   subroutine rad_ozone_physics (im,ix,kx,deltim,ozi,ozo,prsi,prsl,xlat,       &
#ifndef O3CHEM
                               pis,lat)
#else
                               t1,q1,czmn,lat,inistp)
#endif
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
   use constant, only  : rd_,rv_,pi_,fv_
#ifndef O3CHEM
   use comfcst, only   : jo3, ko3, prdin, disin
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! variables : 
!     im       - number of profiles to compute    
!     t1       - (LONS2,levs_) temperature in k  
!     q1       - (LONS2,levh_) specific humidity in kg/kg
!     pstar    - (LONS2) surface pressure in kpa        
!     sl       - (levs_) p/psfc at middle of layer
!     deltim   - time step in secs
!     czmn     - mean cosine solar zenith angle for all gaussian lats
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-10-07  yeon-joo lee           add o3chem
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
#undef DEBUG
   integer              ::  i,k,n,l,im,ix,kx,j1,j2
#ifdef O3CHEM 
   integer              :: m, inistp
#endif 
   integer              ::  lat,kmin,kmax
   real                 ::  dt,deltim,dp,dphii,elat,xdeg,tem,temp,ozib
   real                 ::  pmin,pmax,psmin,psmax
#ifdef DEBUG
   real                 ::  lond,latd
#endif
#ifndef O3CHEM
   real,parameter       ::  blat=-85.0,dphi=10.0
#else
   real,parameter       ::  xavog=6.022e23 ! Avogadro's number 6.022e23/mole
   real,parameter       ::  weight_air=28.964
   real,parameter       ::  percnt_o2=0.23133
   integer,parameter    ::  weight_o3=48,weight_o2=32 
#endif
   real                 ::  prsi(ix,kx+1),prsl(ix,kx)
#ifndef O3CHEM
   real                 ::  ozi(ix,kx),ozo(ix,kx),xlat(ix),pis(ix)
#else
   real                 ::  ozi(ix,kx),ozo(ix,kx),xlat(ix)
   real                 ::  t1(ix,kx),q1(ix,kx),czmn(ix)
#endif
#ifndef O3CHEM
   real                 ::  po3(ko3), slog(ILOTS,levs_),                       &
                            prdout(ILOTS,ko3),                                 &
                            disout(ILOTS,ko3)
   real                 ::  wk1(ILOTS), wk2(ILOTS),                            &
                            wk3(ILOTS),                                        &
                            wk4(ILOTS),                                        &
                            wk5(ILOTS), wkp(ILOTS,2),                          &
                            wkd(ILOTS,2),                                      &
                            prod(ILOTS),dist(ILOTS),                           &
                            ddy(ILOTS)
  integer               ::  jindx1(ILOTS),jindx2(ILOTS)
#else
!
!---- chemistry pack : Aeronomy of the Middle Atm.,guy brasseur, 1986
!
   real                 ::  coeffj1a(9),coeffj1b(9),coeffj1c(9)
   real                 ::  coeffj2a(9),coeffj2b(9),coeffj2c(9)
   real                 ::  coeffj1t(9),coeffj2t(9)
   real                 ::  pressa(9)
!
   real                 ::  press(levs_),vtj(ILOTS,levs_),rho(ILOTS,levs_)
   real                 ::  prodrate(ILOTS,levs_),dissrate(ILOTS,levs_)
   real                 ::  coeffj1(levs_),coeffj2(levs_)
   real                 ::  o1,o2,o3,air,xk2,xk3,o3_out,o_out,dpress,dcoeffj
!
   data coeffj1a/5.07e-18, 8.21e-15, 5.25e-13, 7.42e-12, 4.93e-11,             &
                 7.48e-10, 9.02e-10, 9.13e-10, 9.14e-10/
!
   data coeffj1b/0.0, 0.0, 0.0, 0.0, 0.0,                                      &
             0.0, 0.0, 3.67e-20, 7.40e-8/
!
   data coeffj1c/1.78e-18, 5.60e-15, 4.14e-13, 5.04e-12, 2.47e-12,             &
             1.85e-10, 6.48e-10, 6.99e-9, 4.34e-8/
!
   data coeffj2a/1.44e-5, 1.22e-5, 1.96e-5, 4.26e-5, 1.26e-4,                  &
              4.44e-3, 7.52e-3, 7.66e-3, 7.67e-3/
!
   data coeffj2b/5.58e-5, 5.96e-5, 6.42e-5, 7.16e-5, 7.83e-5,                  &
              7.47e-5, 7.82e-5, 1.13e-4, 1.02e-4/
!
   data coeffj2c/3.04e-4, 3.40e-4, 3.40e-4, 3.40e-4, 3.40e-4,                  &
              3.40e-4, 3.40e-4, 3.40e-4, 3.40e-4/
!
!-------reference : ch.4 table 4.3
!
      data pressa/ 269.0, 122.0,  55.0,  25.0, 11.5,                           &
                           1.40, 0.20, 0.0089, 0.00064 /
#endif
!-------------------------------------------------------------------------------
!
   dt=deltim * 2
!
#ifdef O3CHEM
!
! chemistry pack ; j1,j2 iteration
!
   do i = 1,im
     do k = 1,levs_
       press(k) = prsl(i,k)*10.0   ! mb
     enddo
     do m = 1,9
       coeffj1t(m) = coeffj1a(m)+coeffj1b(m)+coeffj1c(m)
       coeffj2t(m) = coeffj2a(m)*0.3+coeffj2b(m)*0.7
     enddo
     do k = 1,levs_
       do m = 1,8
       if(press(k).gt.pressa(1)) then
         coeffj1(k) = coeffj1t(1)*0.01
         coeffj2(k) = coeffj2t(1)*0.01
       endif
       if(press(k).le.pressa(m).and.press(k).ge.pressa(m+1)) then
          dpress    = pressa(m)-pressa(m+1)
         dcoeffj   = coeffj1t(m+1)-coeffj1t(m)
         coeffj1(k)= dcoeffj/dpress*(pressa(m)-press(k))+coeffj1t(m)
         dcoeffj   = coeffj2t(m+1)-coeffj2t(m)
         coeffj2(k)= dcoeffj/dpress*(pressa(m)-press(k))+coeffj2t(m)
#ifdef DEBUG
         lond = 0
         latd = 42
         if(lat.eq.latd.and.j.eq.lond ) then
           print *,'press(k)=',press(k),k
           print *,'pressa(m)=',pressa(m),m
           print *,'pressa(m-m+1)=',pressa(m)-pressa(m+1)
           print *,'coeffj1,coeffj2=',coeffj1(k),coeffj2(k)
         endif
#endif
       endif
       if(press(k).lt.pressa(9)) then
         coeffj1(k) = coeffj1t(9)
         coeffj2(k) = coeffj2t(9)
       endif 
     enddo
   enddo
!
!--- prod, diss rate calculate ----------
!
   do k = 1,levs_   ! modified loop (ozonepack1) ; density
     vtj(i,k)  = t1(i,k) * (1.+fv_*q1(i,k))
     rho(i,k)   = (prsl(i,k)/rd_) / vtj(i,k) ! density tons/m**3=g/cm**3
     air = rho(i,k)/weight_air*xavog !  concentration of [M] cm-3
     o2  = rho(i,k)/weight_o2*xavog*percnt_o2  ! concentration /cm**3
     xk2 = 6.0e-34 * (300.0/t1(i,k))**(2.3) 
     xk3 = 8.0e-12 * exp(-2060.0/t1(i,k)) 
!
! ref.: http://en.wikipedia.org/wiki/Atmospheric#Mass_percentage_calculations
! for the first time step, equilibrium state o3 is used to prod,diss
! and from the next time step, former ozone value is used
!   
     if(inistp.eq.1 ) then
       o3 = o2*sqrt((xk2*air*coeffj1(k))/(coeffj2(k)*xk3))
       ozi(i,k) = o3/rho(i,k)*weight_o3/xavog
#ifdef DEBUG
       lond = 0
       latd = 42
       if ( lat.eq.latd.and.i.eq.lond ) then
         print *,'o3 .steady. calculation -------'
         print '("k2,k3,air,j2,j3= ",5e11.3)',xk2,xk3,air,                     &
         coeffj1(k),coeffj2(k)
         print '("(k2*air*coeffj1(k))/(coeffj2(k)*k3) =",e11.3)',              &
            ((xk2*air*coeffj1(k))/(coeffj2(k)*xk3))
         print '("sqrt()=",e11.3)',                                            &
           sqrt((xk2*air*coeffj1(k))/(coeffj2(k)*xk3))
         print '("o3 = ",e11.3)',                                              &
           o2*sqrt((xk2*air*coeffj1(k))/(coeffj2(k)*xk3))
         print *,'--------------------------------'
       endif
#endif
     else
       o3 = ozi(i,k)*rho(i,k)/weight_o3*xavog
#ifdef DEBUG
       if ( lat.eq.latd.and.j.eq.lond ) then
         print *,'o3 calculation ------------------'
         print '("ozi,ro = ",2e11.3)',ozi(i,k),rho(i,k)
         print '("weight_o3*xavog = ",e11.3)',weight_o3*xavog
         print '("ozi(i,k)*rho(i,k)= ",e11.3)'                                 &
                ,ozi(i,k)*rho(i,k)
         print '("o3 = ",e11.3)',
                ozi(i,k)*rho(i,k)/weight_o3*xavog
         print *,'--------------------------------'
       endif
#endif
     endif
!
!  for night time ------------------
!  include extinguish daytime and nighttime - using solar zenith angle   
!
     if(czmn(i).le.0.0 ) then
       coeffj1(k) = 0.0
       coeffj2(k) = 0.0
     else
       coeffj1(k) = coeffj1(k)*czmn(i)
       coeffj2(k) = coeffj2(k)*czmn(i)
     endif
       o1 = (coeffj2(k)*o3/(xk2*o2*air))
       prodrate(i,k) = xk2*o1*o2*air
       dissrate(i,k) = coeffj2(k) + xk3*o1
#ifdef DEBUG
       lond = 0
       latd = 42
       if ( lat.eq.latd.and.i.eq.lond ) then
         print '("j2,j3,cosz= ",3e11.2)',                                      &
            coeffj1(k),coeffj2(k),czmn(i)
         print '("p,air,o2,o3,o=",f6.1,4e11.2)',press(k),                      &
            air,o2,o3,o
       endif
#endif
!
! change unit from mixing ratio to concentration
!
       ozib = ozi(i,k)*rho(i,k)/weight_o3*xavog
!
! ozone prediction
!
       o3 = (ozib+prodrate(i,k)*dt) / (1.0+dissrate(i,k)*dt)
       ozo(i,k) = o3*weight_o3 / (xavog*rho(i,k))
#ifdef DEBUG
       lond = 0
       latd = 42
       if (lat.eq.latd.and.i.eq.lond ) then
         print '("ozo,ozi,ozo-ozi=",3e13.3)',                                  &
            ozo(i,k),ozi(i,k),ozo(i,k)-ozi(i,k)
         print*,' prodrate(i,k)*dt ',prodrate(i,k)*dt,                         &
            dissrate(i,k)*dt,dt
         print*,' ozib o3 ',ozib,o3
         print *,'1. ozi[kg/kg], 2. ozi[/cm3], 3. ozi[kg/kg]=',                &
            ozi(i,k),ozi(i,k)*rho(i,k)/weight_o3*xavog,                        &
            ozib*weight_o3/(xavog*rho(i,k))
         print '("prod,diss= ",2e11.2)',prodrate(i,k),dissrate(i,k)
       endif 
#endif
     enddo
   enddo
#else 
   dp = 0.2844
   do n = 1,ko3
     po3(n) = 101.30 * exp(-dp*(float(n)-0.5))
     po3(n) = log(po3(n))
   enddo
   do k = 1,levs_
     do i = 1,im
       slog(i,k) = log(prsl(i,k)/prsi(i,1))
     enddo
   enddo
   dphii = 1.0 / dphi
   elat  = blat + (jo3-1)*dphi
!
   do i = 1,im
     xdeg = xlat(i)*180./pi_
     ddy(i) = (xdeg - blat) * dphii + 1.0
     jindx1(i)  = ddy(i)
     jindx2(i)  = jindx1(i) + 1
     ddy(i)     = ddy(i) - jindx1(i)
     if(xdeg .le. blat) then
       jindx1(i) = 1
       jindx2(i) = 1
     endif
     if (xdeg .ge. elat) then
       jindx1(i) = jo3
       jindx2(i) = jo3
     endif
     jindx1(i) = min(max (jindx1(i),1),jo3)
     jindx2(i) = min(max (jindx2(i),1),jo3)
!         print*,' rad_ozone_physics i xdeg ddy j1 j2 ',i,xdeg,ddy(i),
!    1             jindx1(i),jindx2(i)
   enddo
! 
   do l = 1,ko3
     do i = 1,im
       j1 = jindx1(i)
       j2 = jindx2(i)
       tem= 1.0 - ddy(i)
       prdout(i,l) = (tem * prdin(j1,l) + ddy(i) * prdin(j2,l))                &
                   * 1.655
       disout(i,l) = tem * disin(j1,l) + ddy(i) * disin(j2,l)
     enddo
   enddo
!
   do i = 1,im
     wkp(i,1)     = prdout(i,1)
     wkp(i,2)     = prdout(i,ko3)
     wkd(i,1)     = disout(i,1)
     wkd(i,2)     = disout(i,ko3)
   enddo
!
   psmin =  1.0e10
   psmax = -1.0e10
   do i = 1,im
     psmin = min(psmin,pis(i))
     psmax = max(psmax,pis(i))
   enddo
!
   do l = 1,levs_
     do i = 1,im
       pmin = psmin + slog(i,l)
       pmax = psmax + slog(i,l)
       kmax = 1
       kmin = 1
       do k = 1,ko3-1
         if (pmin .lt. po3(k)) kmax = k
         if (pmax .lt. po3(k)) kmin = k
       enddo
       kmax = min(kmax+1,ko3-1)
!        print *,' pmin=',pmin,' pmax=',pmax,' kmin=',kmin
!    *, ' kmax=',kmax,' po3=',po3(kmin),po3(kmax),' l=',l
     enddo
     do i = 1,im
       wk1(i) = pis(i) + slog(i,l)
       prod(i) = 0.0
       dist(i) = 0.0
     enddo
     tem = 10.0 * exp(wk1(1))
     do k = kmin,kmax
       do i = 1,im
         wk2(i)     = prdout(i,k)
         wk3(i)     = prdout(i,k+1)
         wk4(i)     = disout(i,k)
         wk5(i)     = disout(i,k+1)
       enddo
       temp = 1.0 / (po3(k) - po3(k+1))
       do i = 1,im
         if (wk1(i) .lt. po3(k) .and. wk1(i) .ge. po3(k+1)) then
           tem      = (wk1(i) - po3(k+1)) * temp
           prod(i)  = tem * wk2(i) + (1.0-tem) * wk3(i)
           dist(i)  = tem * wk4(i) + (1.0-tem) * wk5(i)
         endif
       enddo
     enddo
!
     do i = 1,im
       if (wk1(i) .lt. po3(ko3)) then
         prod(i) = wkp(i,2)
         dist(i) = wkd(i,2)
       endif
       if (wk1(i) .ge. po3(1)) then
         prod(i) = wkp(i,1)
         dist(i) = wkd(i,1)
       endif
       ozib = ozi(i,l)  ! no filling
       ozo(i,l) = (ozib + prod(i)*dt) / (1.0 + dist(i)*dt)
#ifdef DEBUG
       lond = 0
       latd = 42
       if ( lat.eq.latd.and.j.eq.lond ) then
         print '("ozo,ozi,ozo-ozi=",3e13.3)',                                  &
               ozo(i,l),ozi(i,l),ozo(i,l)-ozi(i,l)
         print '("prod,diss= ",2e11.2)',prod(i),dist(i)
       endif
#endif      
     enddo
   enddo
#endif
!
   return
   end subroutine rad_ozone_physics
