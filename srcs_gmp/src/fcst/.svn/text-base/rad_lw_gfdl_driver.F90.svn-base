#include <define.h>
   subroutine rad_lw_gfdl_driver(ipts,heatra,grnflx,topflx,                    &
#ifdef CLR
                    grnfx0,topfx0,                                             &
#endif
#ifdef NIM
                    prssi,                                                     &
#endif
                    press,temp,rh2o,qo3,cldfac,                                &
                    camt,nclds,ktop,kbtm)
!-------------------------------------------------------------------------------
!
! subroutine: rad_lw_gfdl_driver
!
! abstract: subroutine rad_lw_gfdl_driver computes temperature-corrected co2 
!           transmission functions and also computes the pressure grid and 
!           layer optical paths.
!
!          inputs:                (common blocks)
!      cldfac                          cldcom
!      press,temp,rh2o,qo3             radisw
!      camt,nclds,ktop,kbtm            radisw
!      co251,co258,cdt51,cdt58         co2bd3
!      c2d51,c2d58,co2m51,co2m58       co2bd3
!      cdtm51,cdtm58,c2dm51,c2dm58     co2bd3
!      stemp,gtemp                     co2bd3
!      co231,co238,cdt31,cdt38         co2bd2
!      c2d31,c2d38                     co2bd2
!      co271,co278,cdt71,cdt78         co2bd4
!      c2d71,c2d78                     co2bd4
!      betinw                          bdwide
!          outputs:
!      heatra,grnflx,topflx            lwout
!          called by:
!      radmn or input routine of model
!          calls:
!      rad_lw
!
! program history log:
!   1988-05-06  kenneth campana        development
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!  ::: structure :::
!
!    [rad_lw_gfdl_driver] - [ rad_lw ]
!
!-------------------------------------------------------------------------------
   use paramodel
   use hcon
   use rdparm
   use co2dta
   use rnddta
!
#ifdef NIM
   real                 ::  prssi(imbx,lp1)
#endif
   real                 ::  press(imbx,lp1),temp(imbx,lp1),rh2o(imbx,l)
   real                 ::  qo3(imbx,l)
   real                 ::  cldfac(imbx,lp1,lp1),camt(imbx,lp1)
   integer              ::  nclds(imax),ktop(imbx,lp1),kbtm(imbx,lp1)
   real                 ::  heatra(imbx,l),grnflx(imax),topflx(imax)
#ifdef CLR
   real                 ::  grnfx0(imax),topfx0(imax)
#endif
   real                 ::  delp2(imbx,l)
!
   real                 ::  qh2o(imbx,l),t(imbx,lp1)
   real                 ::  p(imbx,lp1),delp(imbx,l)
   real                 ::  co21(imbx,lp1,lp1),co2nbl(imbx,l)
   real                 ::  co2sp1(imbx,lp1),co2sp2(imbx,lp1)
   real                 ::  var1(imbx,l),var2(imbx,l)
   real                 ::  var3(imbx,l),var4(imbx,l)
   real                 ::  cntval(imbx,lp1)
   real                 ::  toto3(imbx,lp1),tphio3(imbx,lp1),totphi(imbx,lp1)
   real                 ::  totvo2(imbx,lp1),emx1(imax),emx2(imax),            &
                            empl(imbx,llp1)
!
   real                 ::  co2r(imbx,lp1),dift(imbx,lp1)
   real                 ::  co2r1(imbx,lp1),dco2d1(imbx,lp1)
   real                 ::  d2cd21(imbx,lp1),d2cd22(imbx,lp1)
   real                 ::  co2r2(imbx,lp1),dco2d2(imbx,lp1)
   real                 ::  co2mr(imbx,l),co2md(imbx,l),co2m2d(imbx,l)
   real                 ::  tdav(imbx,lp1),tstdav(imbx,lp1)
   real, target         ::  vv(imbx,l)
   real                 ::  vsum1(imax),vsum2(imax)
   real, pointer        ::  vsum3(:,:)
   real                 ::  a1(imax),a2(imax)
   real                 ::  dco2dt(imbx,lp1),d2cdt2(imbx,lp1)
!
   real, target         ::  texpsl(imbx,lp1)
   real, pointer        ::  tlsqu(:,:)
   real, pointer        ::  vsum4(:,:)
!
   vsum3=>texpsl ; tlsqu=>texpsl
   vsum4=>vv
!
!****compute flux pressures (p) and differences (delp2,delp)
!****compute flux level temperatures (t) and continuum temperature
!    corrections (texpsl)
!
   do k = 2,l
     do i = 1,ipts
#ifdef NIM
       p(i,k)=prssi(i,k)*1.0e3
#else
       p(i,k)=haf*(press(i,k-1)+press(i,k))
#endif
       t(i,k)=haf*(temp(i,k-1)+temp(i,k))
     enddo
   enddo
   do i = 1,ipts
#ifdef NIM
     p(i,1)=prssi(i,1)*1.0e3
#else
     p(i,1)=zero
#endif
     p(i,lp1)=press(i,lp1)
     t(i,1)=temp(i,1)
     t(i,lp1)=temp(i,lp1)
   enddo
   do k = 1,l
     do i = 1,ipts
       delp2(i,k)=p(i,k+1)-p(i,k)
       delp(i,k)=one/delp2(i,k)
     enddo
   enddo
!
!****compute argument for cont.temp.coeff.
!    (this is 1800.(1./temp-1./296.))
!mk   if( temp(ipts,lp1).lt.1. ) close(90)
!
   do k = 1,lp1
     do i = 1,ipts
       texpsl(i,k)=h18e3/temp(i,k)-h6p08108
!
!...then take exponential
!
       texpsl(i,k)=exp(texpsl(i,k))
     enddo
   enddo
!
!***compute optical paths for h2o and o3, using the diffusivity
!   approximation for the angular integration (1.66). obtain the
!   unweighted values(var1,var3) and the weighted values(var2,var4).
!   the quantities h3m4(.0003) and h3m3(.003) appearing in the var2 and
!   var4 expressions are the approximate voigt corrections for h2o and
!   o3,respectively.
!
   do k = 1,l
     do i = 1,ipts
       qh2o(i,k)=rh2o(i,k)*diffctr
!
!---vv is the layer-mean pressure (in atm),which is not the same as
!   the level pressure (press)
!
       vv(i,k)=haf*(p(i,k+1)+p(i,k))*p0inv
       var1(i,k)=delp2(i,k)*qh2o(i,k)*ginv
       var3(i,k)=delp2(i,k)*qo3(i,k)*diffctr*ginv
       var2(i,k)=var1(i,k)*(vv(i,k)+h3m4)
       var4(i,k)=var3(i,k)*(vv(i,k)+h3m3)
!
!  compute optical path for the h2o continuum, using roberts coeffs.
!  (betinw),and temp. correction (texpsl). the diffusivity factor
!  (which cancels out in this expression) is assumed to be 1.66. the
!  use of the diffusivity factor has been shown to be a significant
!  source of error in the continuum calcs.,but the time penalty of
!  an angular integration is severe.
!
       cntval(i,k)=texpsl(i,k)*rh2o(i,k)*var2(i,k)*betinw/                   &
                (rh2o(i,k)+rath2omw)
     enddo
   enddo
!
!   compute summed optical paths for h2o,o3 and continuum
!
   do i = 1,ipts
     totphi(i,1)=zero
     toto3(i,1)=zero
     tphio3(i,1)=zero
     totvo2(i,1)=zero
   enddo
   do k = 2,lp1
     do i = 1,ipts
       totphi(i,k)=totphi(i,k-1)+var2(i,k-1)
       toto3(i,k)=toto3(i,k-1)+var3(i,k-1)
       tphio3(i,k)=tphio3(i,k-1)+var4(i,k-1)
       totvo2(i,k)=totvo2(i,k-1)+cntval(i,k-1)
     enddo
   enddo
!
!---emx1 is the additional pressure-scaled mass from press(l) to
!   p(l). it is used in nearby layer and emiss calculations.
!---emx2 is the additional pressure-scaled mass from press(l) to
!   p(lp1). it is used in calculations between flux levels l and lp1.
!
   do i = 1,ipts
     emx1(i)=qh2o(i,l)*press(i,l)*(press(i,l)-p(i,l))*gp0inv
     emx2(i)=qh2o(i,l)*press(i,l)*(p(i,lp1)-press(i,l))*gp0inv
   enddo
!
!---empl is the pressure scaled mass from p(k) to press(k) (index 2-lp1)
!   or to press(k+1) (index lp2-ll)
!
   do k = 1,l
     do i = 1,ipts
       empl(i,k+1)=qh2o(i,k)*p(i,k+1)*(p(i,k+1)-press(i,k))*gp0inv
     enddo
   enddo
   do k = 1,lm1
     do i = 1,ipts
       empl(i,k+lp1)=qh2o(i,k+1)*p(i,k+1)*(press(i,k+1)-p(i,k+1))*gp0inv
     enddo
   enddo
   do i = 1,ipts
     empl(i,1)=var2(i,l)
     empl(i,llp1)=empl(i,ll)
   enddo
!
!***compute weighted temperature (tdav) and pressure (tstdav) integrals
!   for use in obtaining temp. difference bet. sounding and std.
!   temp. sounding (dift)
!
   do i = 1,ipts
     tstdav(i,1)=zero
     tdav(i,1)=zero
   enddo
   do k = 1,lp1
     do i = 1,ipts
       vsum3(i,k)=temp(i,k)-stemp(k)
     enddo
   enddo
   do k = 1,l
     do i = 1,ipts
       vsum2(i)=gtemp(k)*delp2(i,k)
       vsum1(i)=vsum2(i)*vsum3(i,k)
       tstdav(i,k+1)=tstdav(i,k)+vsum2(i)
       tdav(i,k+1)=tdav(i,k)+vsum1(i)
     enddo
   enddo
!
!
!****evaluate coefficients for co2 pressure interpolation (a1,a2)
   do i = 1,ipts
     a1(i)=(press(i,lp1)-p0xzp8)/p0xzp2
     a2(i)=(p0-press(i,lp1))/p0xzp2
   enddo
!
!***perform co2 pressure interpolation on all inputted transmission
!   functions and temp. derivatives
!---successively computing co2r,dco2dt and d2cdt2 is done to save
!   storage (at a slight loss in computation time)
!
   do k = 1,lp1
     do i = 1,ipts
       co2r1(i,k)=a1(i)*co231(k)+a2(i)*co238(k)
       d2cd21(i,k)=h1m3*(a1(i)*c2d31(k)+a2(i)*c2d38(k))
       dco2d1(i,k)=h1m2*(a1(i)*cdt31(k)+a2(i)*cdt38(k))
       co2r2(i,k)=a1(i)*co271(k)+a2(i)*co278(k)
       d2cd22(i,k)=h1m3*(a1(i)*c2d71(k)+a2(i)*c2d78(k))
       dco2d2(i,k)=h1m2*(a1(i)*cdt71(k)+a2(i)*cdt78(k))
     enddo
   enddo
   do k = 1,l
     do i = 1,ipts
       co2mr(i,k)=a1(i)*co2m51(k)+a2(i)*co2m58(k)
       co2md(i,k)=h1m2*(a1(i)*cdtm51(k)+a2(i)*cdtm58(k))
       co2m2d(i,k)=h1m3*(a1(i)*c2dm51(k)+a2(i)*c2dm58(k))
     enddo
   enddo
!
!***compute co2 temperature interpolations for all bands,using dift
!
!   the case where k=1 is handled first. we are now replacing
!   3-dimensional arrays by 2-d arrays, to save space. thus this
!   calculation is for (i,kp,1)
!
   do kp = 2,lp1
     do i = 1,ipts
       dift(i,kp)=tdav(i,kp)/tstdav(i,kp)
     enddo
   enddo
   do i = 1,ipts
     co21(i,1,1)=1.0
     co2sp1(i,1)=1.0
     co2sp2(i,1)=1.0
   enddo
   do kp = 2,lp1
     do i = 1,ipts
!
!---calculations for kp>1 for k=1
!
       co2r(i,kp)=a1(i)*co251(kp,1)+a2(i)*co258(kp,1)
       dco2dt(i,kp)=h1m2*(a1(i)*cdt51(kp,1)+a2(i)*cdt58(kp,1))
       d2cdt2(i,kp)=h1m3*(a1(i)*c2d51(kp,1)+a2(i)*c2d58(kp,1))
       co21(i,kp,1)=co2r(i,kp)+dift(i,kp)*(dco2dt(i,kp)+                       &
                haf*dift(i,kp)*d2cdt2(i,kp))
!
!---calculations for (effectively) kp=1,k>kp. these use the
!   same value of dift due to symmetry
!
       co2r(i,kp)=a1(i)*co251(1,kp)+a2(i)*co258(1,kp)
       dco2dt(i,kp)=h1m2*(a1(i)*cdt51(1,kp)+a2(i)*cdt58(1,kp))
       d2cdt2(i,kp)=h1m3*(a1(i)*c2d51(1,kp)+a2(i)*c2d58(1,kp))
       co21(i,1,kp)=co2r(i,kp)+dift(i,kp)*(dco2dt(i,kp)+                       &
                haf*dift(i,kp)*d2cdt2(i,kp))
     enddo
   enddo
!
!   the transmission functions used in rad_lw_spa88 may be computed now.
!---(in the 250 loop,dift really should be (i,1,k), but dift is
!    invariant with respect to k,kp,and so (i,1,k)=(i,k,1))
!
   do k = 2,lp1
     do i = 1,ipts
       co2sp1(i,k)=co2r1(i,k)+dift(i,k)*(dco2d1(i,k)+haf*dift(i,k)*           &
                  d2cd21(i,k))
       co2sp2(i,k)=co2r2(i,k)+dift(i,k)*(dco2d2(i,k)+haf*dift(i,k)*           &
                  d2cd22(i,k))
     enddo
   enddo
!
!   next the case when k=2...l
!
   do k = 2,l
!DIR$ unroll 50
     do kp= k+1,lp1
       do i = 1,ipts
         dift(i,kp)=(tdav(i,kp)-tdav(i,k))/                                    &
                 (tstdav(i,kp)-tstdav(i,k))
         co2r(i,kp)=a1(i)*co251(kp,k)+a2(i)*co258(kp,k)
         dco2dt(i,kp)=h1m2*(a1(i)*cdt51(kp,k)+a2(i)*cdt58(kp,k))
         d2cdt2(i,kp)=h1m3*(a1(i)*c2d51(kp,k)+a2(i)*c2d58(kp,k))
         co21(i,kp,k)=co2r(i,kp)+dift(i,kp)*(dco2dt(i,kp)+                     &
                haf*dift(i,kp)*d2cdt2(i,kp))
         co2r(i,kp)=a1(i)*co251(k,kp)+a2(i)*co258(k,kp)
         dco2dt(i,kp)=h1m2*(a1(i)*cdt51(k,kp)+a2(i)*cdt58(k,kp))
         d2cdt2(i,kp)=h1m3*(a1(i)*c2d51(k,kp)+a2(i)*c2d58(k,kp))
         co21(i,k,kp)=co2r(i,kp)+dift(i,kp)*(dco2dt(i,kp)+                     &
                haf*dift(i,kp)*d2cdt2(i,kp))
       enddo
     enddo
   enddo
!
!   finally the case when k=kp,k=2..lp1
!
   do k = 2,lp1
     do i = 1,ipts
       dift(i,k)=haf*(vsum3(i,k)+vsum3(i,k-1))
       co2r(i,k)=a1(i)*co251(k,k)+a2(i)*co258(k,k)
       dco2dt(i,k)=h1m2*(a1(i)*cdt51(k,k)+a2(i)*cdt58(k,k))
       d2cdt2(i,k)=h1m3*(a1(i)*c2d51(k,k)+a2(i)*c2d58(k,k))
       co21(i,k,k)=co2r(i,k)+dift(i,k)*(dco2dt(i,k)+                           &
                haf*dift(i,k)*d2cdt2(i,k))
     enddo
   enddo
!
!--- we arent doing nbl tfs on the 100 cm-1 bands .
!
   do k = 1,l
     do i = 1,ipts
       co2nbl(i,k)=co2mr(i,k)+vsum3(i,k)*(co2md(i,k)+haf*                      &
               vsum3(i,k)*co2m2d(i,k))
     enddo
   enddo
!
!***compute temp. coefficient based on t(k) (see ref.2)
!
   do k = 1,lp1
     do i = 1,ipts
       if (t(i,k).le.h25e2) then
         tlsqu(i,k)=b0+(t(i,k)-h25e2)*                                         &
                  (b1+(t(i,k)-h25e2)*                                          &
                  (b2+b3*(t(i,k)-h25e2)))
       else
         tlsqu(i,k)=b0
       endif
     enddo
   enddo
!
!***apply to all co2 tfs
!
   do k = 1,lp1
     do kp= 1,lp1
       do i = 1,ipts
         co21(i,kp,k)=co21(i,kp,k)*(one-tlsqu(i,kp))+tlsqu(i,kp)
       enddo
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       co2sp1(i,k)=co2sp1(i,k)*(one-tlsqu(i,1))+tlsqu(i,1)
       co2sp2(i,k)=co2sp2(i,k)*(one-tlsqu(i,1))+tlsqu(i,1)
     enddo
   enddo
   do k = 1,l
     do i = 1,ipts
       co2nbl(i,k)=co2nbl(i,k)*(one-tlsqu(i,k))+tlsqu(i,k)
     enddo
   enddo
   call rad_lw(ipts,heatra,grnflx,topflx,                                      &
#ifdef CLR
              grnfx0,topfx0,                                                   &
#endif
              qh2o,press,p,delp,delp2,temp,t,                                  &
              cldfac,nclds,ktop,kbtm,camt,                                     &
              co21,co2nbl,co2sp1,co2sp2,                                       &
              var1,var2,var3,var4,cntval,                                      &
              toto3,tphio3,totphi,totvo2,                                      &
              emx1,emx2,empl)
!
   return
   end subroutine rad_lw_gfdl_driver
!-------------------------------------------------------------------------------
#include <define.h>
!-------------------------------------------------------------------------------
   subroutine rad_lw(ipts,heatra,grnflx,topflx,                                &
#ifdef CLR
                    grnfx0,topfx0,                                             &
#endif
                    qh2o,press,p,delp,delp2,temp,t,                            &
                    cldfac,nclds,ktop,kbtm,camt,                               &
                    co21,co2nbl,co2sp1,co2sp2,                                 &
                    var1,var2,var3,var4,cntval,                                &
                    toto3,tphio3,totphi,totvo2,                                &
                    emx1,emx2,empl)
!-------------------------------------------------------------------------------
   use paramodel
   use hcon
   use rdparm
   use rnddta
   use funct_rad_lw, only : rad_lw_ele290, rad_lw_e290, rad_lw_e2spec,         &
                            rad_lw_e3v88,  rad_lw_spa88
!-------------------------------------------------------------------------------
#include <tabcom.h>
!
   real                 ::  qh2o(imbx,l),press(imbx,lp1)
   real                 ::  p(imbx,lp1),delp(imbx,l),delp2(imbx,l)
   real                 ::  temp(imbx,lp1)
   real                 ::  t(imbx,lp1),cldfac(imbx,lp1,lp1)
   real                 ::  camt(imbx,lp1)
   integer              ::  nclds(imax),ktop(imbx,lp1),kbtm(imbx,lp1)
   real                 ::  co21(imbx,lp1,lp1),co2nbl(imbx,l)
   real                 ::  co2sp1(imbx,lp1),co2sp2(imbx,lp1)
   real                 ::  var1(imbx,l),var2(imbx,l),var3(imbx,l)
   real                 ::  var4(imbx,l)
   real                 ::  cntval(imbx,lp1)
   real                 ::  heatra(imbx,l),grnflx(imax),topflx(imax)
!
#ifdef CLR
   real                 ::  heatr0(imbx,l),flxnt0(imbx,lp1)
   real                 ::  grnfx0(imax),topfx0(imax),gxcts0(imax)
   real                 ::  flx1e0(imax)
   real                 ::  excts0(imbx,l),ctso30(imbx,l),cts0(imbx,l)
   real                 ::  flx0(imbx,lp1)
#endif
!
   real                 ::  gxcts(imax),flx1e1(imax)
   real                 ::  avephi(imbx,lp1),emiss(imbx,lp1),emissb(imbx,lp1)
!
   real                 ::  toto3(imbx,lp1),tphio3(imbx,lp1),totphi(imbx,lp1)
   real                 ::  totvo2(imbx,lp1),emx1(imax),emx2(imax)
   real                 ::  empl(imbx,llp1)
!
   real                 ::  excts(imbx,l),ctso3(imbx,l),cts(imbx,l)
   real                 ::  e1flx(imbx,lp1)
   real                 ::  co2sp(imbx,lp1),to3spc(imbx,l),to3sp(imbx,lp1)
   real                 ::  oss(imbx,lp1),css(imbx,lp1),ss1(imbx,lp1)
   real                 ::  ss2(imbx,lp1)
   real                 ::  tc(imbx,lp1),dtc(imbx,lp1)
   real                 ::  sorc(imbx,lp1,nbly),csour(imbx,lp1)
!
!cc
!
   real, target         ::  avvo2(imbx,lp1)
   real, target         ::  avmo3(imbx,lp1)
   real, target         ::  avpho3(imbx,lp1)
   real, pointer        ::  heatem(:,:)
   real, pointer        ::  over1d(:,:)
   real, pointer        ::  c(:,:)
   real, pointer        ::  c2(:,:)
   integer              ::  itop(imax),ibot(imax),indtc(imax)
   real                 ::  to31d(imbx,lp1),cont1d(imbx,lp1)
   real                 ::  delptc(imax),ptop(imax),pbot(imax),ftop(imax)
   real                 ::  fbot(imax) ,emspec(imbx,2)
!
!---dimension of variables equivalenced to those in vtemp---
!
   real, target         ::  dsorc(imbx,lp1)
   real, target         ::  alp(imbx,llp1)
   real, target         ::  delpr1(imbx,lp1)
   real                 ::  vtmp3(imbx,lp1)
   real, pointer        ::  csub(:,:)
   real, target         ::  csub2(imbx,llp1)
   real, pointer        ::  fac1(:,:)
   real, pointer        ::  delpr2(:,:)
   real, pointer        ::  emisdg(:,:)
   real, pointer        ::  contdg(:,:)
   real, pointer        ::  to3dg(:,:)
   real, pointer        ::  flxnet(:,:)
   integer              ::  ixo(imbx,lp1)
   real, pointer        ::  vsum1(:,:)
   real, pointer        ::  flxthk(:,:)
   real, pointer        ::  z1(:,:)
!
!---dimension of variables passed to other subroutines---
!   (and not found in common blocks)
!
   real, target         ::  emd(imbx,llp1)
   real, target         ::  tpl(imbx,llp1)
   real, pointer        ::  e1cts1(:,:),e1cts2(:,:)
   real, pointer        ::  e1ctw1(:,:),e1ctw2(:,:)
!
!   it is possible to equivalence emd,tpl to the above variables
!   as they get called at different times
!
   real                 ::  fxo(imbx,lp1),dt(imbx,lp1)
   real                 ::  fxoe2(imbx,lp1),dte2(imbx,lp1)
   real                 ::  fxosp(imbx,2),dtsp(imbx,2)
!
!     dimension of local variables
!
   real                 ::  rlog(imbx,l),flx(imbx,lp1)
   real                 ::  totevv(imbx,lp1),cnttau(imbx,lp1)
!
   c=>alp   ;  csub=>alp
   c2=>csub2
   fac1=>dsorc ; over1d=>dsorc ; delpr2=>dsorc ; flxnet=>dsorc
   heatem=>delpr1
   flxthk=>avvo2 ; to3dg=>avvo2
   z1=>avmo3 ; contdg=>avmo3
   emisdg=>avpho3 ; vsum1=>avpho3
   e1cts1=>emd(:,1:lp1)
   e1cts2=>emd(:,lp2:llp1)
   e1ctw1=>tpl(:,1:lp1)
   e1ctw2=>tpl(:,lp2:llp1)
!
!          first section is table lookup for source function and
!     derivative (b and db/dt).also,the nlte co2 source function
!     is obtained
!
!---in calcs. below, decrementing the index by 9
!   accounts for the tables beginning at t=100k.
!   at t=100k.
!
   do k = 1,lp1
     do i = 1,ipts
!
!---temp. indices for e1,source
!
       vtmp3(i,k)=aint(temp(i,k)*hp1)
       fxo(i,k)=vtmp3(i,k)-9.
       dt(i,k)=temp(i,k)-ten*vtmp3(i,k)
!
!---integer index for source (used immediately)
!
       ixo(i,k)=fxo(i,k)
     enddo
   enddo

   do k = 1,l
     do i = 1,ipts
!
!---temp. indices for e2 (kp=1 layer not used in flux calculations)
!
       vtmp3(i,k)=aint(t(i,k+1)*hp1)
       fxoe2(i,k)=vtmp3(i,k)-9.
       dte2(i,k)=t(i,k+1)-ten*vtmp3(i,k)
     enddo
   enddo
!
!---special case to handle kp=lp1 layer and special e2 calcs.
!
   do i = 1,ipts
     fxoe2(i,lp1)=fxo(i,l)
     dte2(i,lp1)=dt(i,l)
     fxosp(i,1)=fxoe2(i,lm1)
     fxosp(i,2)=fxo(i,lm1)
     dtsp(i,1)=dte2(i,lm1)
     dtsp(i,2)=dt(i,lm1)
   enddo
!
!---source function for combined band 1
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),1)
       dsorc(i,k)=dsrce(ixo(i,k),1)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,1)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for combined band 2
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),2)
       dsorc(i,k)=dsrce(ixo(i,k),2)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,2)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for combined band 3
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),3)
       dsorc(i,k)=dsrce(ixo(i,k),3)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,3)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for combined band 4
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),4)
       dsorc(i,k)=dsrce(ixo(i,k),4)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,4)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for combined band 5
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),5)
       dsorc(i,k)=dsrce(ixo(i,k),5)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,5)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for combined band 6
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),6)
       dsorc(i,k)=dsrce(ixo(i,k),6)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,6)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for combined band 7
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),7)
       dsorc(i,k)=dsrce(ixo(i,k),7)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,7)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for combined band 8
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),8)
       dsorc(i,k)=dsrce(ixo(i,k),8)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,8)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for band 9 (560-670 cm-1)
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),9)
       dsorc(i,k)=dsrce(ixo(i,k),9)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,9)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for band 10 (670-800 cm-1)
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),10)
       dsorc(i,k)=dsrce(ixo(i,k),10)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,10)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for band 11 (800-900 cm-1)
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),11)
       dsorc(i,k)=dsrce(ixo(i,k),11)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,11)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for band 12 (900-990 cm-1)
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),12)
       dsorc(i,k)=dsrce(ixo(i,k),12)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,12)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for band 13 (990-1070 cm-1)
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),13)
       dsorc(i,k)=dsrce(ixo(i,k),13)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,13)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!---source function for band 14 (1070-1200 cm-1)
!
   do i = 1,ipts
     do k = 1,lp1
       vtmp3(i,k)=source(ixo(i,k),14)
       dsorc(i,k)=dsrce(ixo(i,k),14)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       sorc(i,k,14)=vtmp3(i,k)+dt(i,k)*dsorc(i,k)
     enddo
   enddo
!
!        the following subroutine obtains nlte source function for co2
!
!
!     call nlte
!
!
!---obtain special source functions for the 15 um band (csour)
!   and the window region (ss1)
!
   do k = 1,lp1
     do i = 1,ipts
       ss1(i,k)=sorc(i,k,11)+sorc(i,k,12)+sorc(i,k,14)
     enddo
   enddo
   do k = 1,lp1
     do i = 1,ipts
       csour(i,k)=sorc(i,k,9)+sorc(i,k,10)
     enddo
   enddo
!
!---compute temp**4 (tc) and vertical temperature differences
!   (oss,css,ss2,dtc). all these will be used later in flux computa-
!   tions.
!
   do k = 1,lp1
     do i = 1,ipts
       tc(i,k)=(temp(i,k)*temp(i,k))**2
     enddo
   enddo
   do k = 1,l
     do i = 1,ipts
       oss(i,k+1)=sorc(i,k+1,13)-sorc(i,k,13)
       css(i,k+1)=csour(i,k+1)-csour(i,k)
       dtc(i,k+1)=tc(i,k+1)-tc(i,k)
       ss2(i,k+1)=ss1(i,k+1)-ss1(i,k)
     enddo
   enddo
!
!
!---the followimg is a drastic rewrite of the radiation code to
!    (largely) eliminate three-dimensional arrays. the code works
!    on the following principles:
!
!          let k = fixed flux level, kp = varying flux level
!          then flux(k)=sum over kp : (deltab(kp)*tau(kp,k))
!               over all kps, from 1 to lp1.
!
!          we can break down the calculations for all ks as follows:
!
!          for all ks k=1 to lp1:
!              flux(k)=sum over kp : (deltab(kp)*tau(kp,k))  (1)
!                      over all kps, from k+1 to lp1
!          and
!              for kp from k+1 to lp1:
!                 flux(kp) = deltab(k)*tau(k,kp)              (2)
!
!          now if tau(k,kp)=tau(kp,k) (symmetrical arrays)
!          we can compute a 1-dimensional array tau1d(kp) from
!          k+1 to lp1, each time k is incremented.
!          equations (1) and (2) then become:
!
!             tau1d(kp) = (values for tau(kp,k) at the particular k)
!             flux(k) = sum over kp : (deltab(kp)*tau1d(kp))   (3)
!             flux(kp) = deltab(k)*tau1d(kp)                   (4)
!
!         the terms for tau (k,k) and other special terms (for
!         nearby layers) must, of course, be handled separately, and
!         with care.
!
!      compute "upper triangle" transmission functions for
!      the 9.6 um band (to3sp) and the 15 um band (over1d). also,
!      the
!      stage 1...compute o3 ,over transmission fctns and avephi
!---do k = 1 calculation (from flux layer kk to the top) separately
!   as vectorization is improved,and ozone cts transmissivity
!   may be extracted here.
!
   do k = 1,l
     do i = 1,ipts
       avephi(i,k)=totphi(i,k+1)
     enddo
   enddo
!
!---in order to properly evaluate emiss integrated over the (lp1)
!   layer, a special evaluation of emiss is done. this requires
!   a special computation of avephi, and it is stored in the
!   (otherwise vacant) lp1th position
!
   do i = 1,ipts
     avephi(i,lp1)=avephi(i,lm1)+emx1(i)
   enddo
!
!   compute fluxes for k=1
!
   call rad_lw_ele290(ipts,e1cts1,e1cts2,e1flx,e1ctw1,e1ctw2,emiss,            &
                      fxo,dt,fxoe2,dte2,avephi,temp,t)
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=bo3rnd(2)*tphio3(i,k+1)/toto3(i,k+1)
       to3spc(i,k)=haf*(fac1(i,k)*                                             &
            (sqrt(one+(four*ao3rnd(2)*toto3(i,k+1))/fac1(i,k))-one))
!
!   for k=1, to3sp is used instead of to31d (they are equal in this
!   case); to3sp is passed to spa90, while to31d is a work-array.
!
       to3sp(i,k)=exp(hm1ez*(to3spc(i,k)+sko3r*totvo2(i,k+1)))
       over1d(i,k)=exp(hm1ez*(sqrt(ab15wd*totphi(i,k+1))+                      &
               skc1r*totvo2(i,k+1)))
!
!---because all continuum transmissivities are obtained from the
!  2-d quantity cnttau (and its reciprocal totevv) we store both
!  of these here. for k=1, cont1d equals cnttau
!
       cnttau(i,k)=exp(hm1ez*totvo2(i,k+1))
       totevv(i,k)=1./cnttau(i,k)
     enddo
   enddo
   do k = 1,l
     do i = 1,ipts
       co2sp(i,k+1)=over1d(i,k)*co21(i,1,k+1)
     enddo
   enddo
   do k = 1,l
     do i = 1,ipts
       co21(i,k+1,1)=co21(i,k+1,1)*over1d(i,k)
     enddo
   enddo
!
!---rlog is the nbl amount for the 15 um band calculation
!
   do i = 1,ipts
     rlog(i,1)=over1d(i,1)*co2nbl(i,1)
   enddo
!
!---the terms when kp=1 for all k are the photon exchange with
!   the top of the atmosphere, and are obtained differently than
!   the other calculations
!
   do k = 2,lp1
     do i = 1,ipts
       flx(i,k)= (tc(i,1)*e1flx(i,k)                                           &
             +ss1(i,1)*cnttau(i,k-1)                                           &
             +sorc(i,1,13)*to3sp(i,k-1)                                        &
             +csour(i,1)*co2sp(i,k))                                           &
             *cldfac(i,1,k)
     enddo
   enddo

   do i = 1,ipts
     flx(i,1)= tc(i,1)*e1flx(i,1)+ss1(i,1)+sorc(i,1,13)                        &
             +csour(i,1)
   enddo
!
!---the kp terms for k=1...
!
   do kp = 2,lp1
     do i = 1,ipts
       flx(i,1)=flx(i,1)+(oss(i,kp)*to3sp(i,kp-1)                              &
                     +ss2(i,kp)*cnttau(i,kp-1)                                 &
                     +css(i,kp)*co21(i,kp,1)                                   &
                     +dtc(i,kp)*emiss(i,kp-1))*cldfac(i,kp,1)
     enddo
   enddo
!
!...    ditto for clear sky....
!
#ifdef CLR
   do k = 2,lp1
     do i = 1,ipts
       flx0(i,k)=  tc(i,1)*e1flx(i,k)                                          &
             +ss1(i,1)*cnttau(i,k-1)                                           &
             +sorc(i,1,13)*to3sp(i,k-1)                                        &
             +csour(i,1)*co2sp(i,k)
     enddo
   enddo
   !
   do i = 1,ipts
     flx0(i,1)=tc(i,1)*e1flx(i,1)+ss1(i,1)+sorc(i,1,13)                        &
             +csour(i,1)
   enddo
#endif
!
!---the kp terms for k=1...
!
#ifdef CLR
   do kp = 2,lp1
     do i = 1,ipts
       flx0(i,1)=flx0(i,1)+ oss(i,kp)*to3sp(i,kp-1)                            &
                     +ss2(i,kp)*cnttau(i,kp-1)                                 &
                     +css(i,kp)*co21(i,kp,1)                                   &
                     +dtc(i,kp)*emiss(i,kp-1)
     enddo
   enddo
#endif
!
!          subroutine rad_lw_spa88 is called to obtain exact cts for water
!     co2 and o3, and approximate cts co2 and o3 calculations.
!
   call rad_lw_spa88(ipts,excts,ctso3,gxcts,sorc,csour,                        &
#ifdef CLR
              excts0,ctso30,gxcts0,                                            &
#endif
              cldfac,temp,press,var1,var2,                                     &
              p,delp,delp2,totvo2,to3sp,to3spc,                                &
              co2sp1,co2sp2,co2sp)
!
!    this section computes the emissivity cts heating rates for 2
!    emissivity bands: the 0-160,1200-2200 cm-1 band and the 800-
!    990,1070-1200 cm-1 band. the remaining cts comtributions are
!    contained in ctso3, computed in rad_lw_spa88.
!
   do i = 1,ipts
     vtmp3(i,1)=1.
   enddo
   do k = 1,l
     do i = 1,ipts
       vtmp3(i,k+1)=cnttau(i,k)*cldfac(i,k+1,1)
     enddo
   enddo
   do k = 1,l
     do i = 1,ipts
       cts(i,k)=radcon*delp(i,k)*(tc(i,k)*                                     &
          (e1ctw2(i,k)*cldfac(i,k+1,1)-e1ctw1(i,k)*cldfac(i,k,1)) +            &
             ss1(i,k)*(vtmp3(i,k+1)-vtmp3(i,k)))
     enddo
   enddo
!
   do k = 1,l
     do i = 1,ipts
       vtmp3(i,k)=tc(i,k)*(cldfac(i,k,1)*(e1cts1(i,k)-e1ctw1(i,k)) -           &
                       cldfac(i,k+1,1)*(e1cts2(i,k)-e1ctw2(i,k)))
     enddo
   enddo
   do i = 1,ipts
     flx1e1(i)=tc(i,lp1)*cldfac(i,lp1,1)*(e1cts1(i,lp1)-e1ctw1(i,lp1))
   enddo
   do k = 1,l
     do i = 1,ipts
       flx1e1(i)=flx1e1(i)+vtmp3(i,k)
     enddo
   enddo
!
!  ... ditto for clear sky  ...
!
#ifdef CLR
   do i = 1,ipts
     vtmp3(i,1)=1.
   enddo
!
   do k = 1,l
     do i = 1,ipts
       vtmp3(i,k+1)=cnttau(i,k)
     enddo
   enddo
!
   do k = 1,l
     do i = 1,ipts
       cts0(i,k)=radcon*delp(i,k)*(tc(i,k)*                                    &
                  (e1ctw2(i,k)-e1ctw1(i,k)) +                                  &
                  ss1(i,k)*(vtmp3(i,k+1)-vtmp3(i,k)))
     enddo
   enddo
!
   do k = 1,l
     do i = 1,ipts
       vtmp3(i,k)=tc(i,k)*(e1cts1(i,k)-e1ctw1(i,k)  -                          &
                       (e1cts2(i,k)-e1ctw2(i,k)))
     enddo
   enddo
!
   do i = 1,ipts
     flx1e0(i)=tc(i,lp1)*(e1cts1(i,lp1)-e1ctw1(i,lp1))
   enddo
!
   do k = 1,l
     do i = 1,ipts
       flx1e0(i)=flx1e0(i)+vtmp3(i,k)
     enddo
   enddo
#endif
!
!---now repeat flux calculations for the k=2..lm1  cases.
!   calculations for flux level l and lp1 are done separately, as all
!   emissivity and co2 calculations are special cases or nearby layers.
!
   do k = 2,lm1
     klen=k
!
     do kk = 1,lp1-k
       do i = 1,ipts
         avephi(i,kk+k-1)=totphi(i,kk+k)-totphi(i,k)
       enddo
     enddo
     do i = 1,ipts
       avephi(i,lp1)=avephi(i,lm1)+emx1(i)
     enddo
!
!---compute emissivity fluxes (e2) for this case. note that
!   we have omitted the nearby later case (emiss(i,k,k)) as well
!   as all cases with k=l or lp1. but these cases have always
!   been handled as special cases, so we may as well compute
!    their fluxes separastely.
!
     call rad_lw_e290(ipts,emissb,emiss,avephi,klen,fxoe2,dte2)
     do kk = 1,lp1-k
       do i = 1,ipts
         avmo3(i,kk+k-1)=toto3(i,kk+k)-toto3(i,k)
         avpho3(i,kk+k-1)=tphio3(i,kk+k)-tphio3(i,k)
         avvo2(i,kk+k-1)=totvo2(i,kk+k)-totvo2(i,k)
         cont1d(i,kk+k-1)=cnttau(i,kk+k-1)*totevv(i,k-1)
       enddo
     enddo
!
     do kk = 1,lp1-k
       do i = 1,ipts
         fac1(i,kk+k-1)=bo3rnd(2)*avpho3(i,kk+k-1)/avmo3(i,kk+k-1)
         vtmp3(i,kk+k-1)=haf*(fac1(i,kk+k-1)*                                  &
                           (sqrt(one+(four*ao3rnd(2)*avmo3(i,kk+k-1))/         &
                           fac1(i,kk+k-1))-one))
         to31d(i,kk+k-1)=exp(hm1ez*(vtmp3(i,kk+k-1)+sko3r*avvo2(i,kk+k-1)))
         over1d(i,kk+k-1)=exp(hm1ez*(sqrt(ab15wd*avephi(i,kk+k-1))+            &
                              skc1r*avvo2(i,kk+k-1)))
         co21(i,kk+k,k)=over1d(i,kk+k-1)*co21(i,kk+k,k)
       enddo
     enddo
     do kp = k+1,lp1
       do i = 1,ipts
            co21(i,k,kp)=over1d(i,kp-1)*co21(i,k,kp)
       enddo
     enddo
!
!---rlog is the nbl amount for the 15 um band calculation
!
     do i = 1,ipts
       rlog(i,k)=over1d(i,k)*co2nbl(i,k)
     enddo
!
!---the kp terms for arbirrary k..
!
     do kp = k+1,lp1
        do i = 1,ipts
          flx(i,k)=flx(i,k)+(oss(i,kp)*to31d(i,kp-1)                           &
                   +ss2(i,kp)*cont1d(i,kp-1)                                   &
                   +css(i,kp)*co21(i,kp,k)                                     &
                   +dtc(i,kp)*emiss(i,kp-1))*cldfac(i,kp,k)
        enddo
     enddo
     do kp = k+1,lp1
       do i = 1,ipts
         flx(i,kp)=flx(i,kp)+(oss(i,k)*to31d(i,kp-1)                           &
                   +ss2(i,k)*cont1d(i,kp-1)                                    &
                   +css(i,k)*co21(i,k,kp)                                      &
                   +dtc(i,k)*emissb(i,kp-1))*cldfac(i,k,kp)
       enddo
     enddo
!
! ....   ditto for clear sky .. cldfac=1.
!
#ifdef CLR
     do kp = k+1,lp1
       do i = 1,ipts
         flx0(i,k)=flx0(i,k)+ oss(i,kp)*to31d(i,kp-1)                          &
                   +ss2(i,kp)*cont1d(i,kp-1)                                   &
                   +css(i,kp)*co21(i,kp,k)                                     &
                   +dtc(i,kp)*emiss(i,kp-1)
       enddo
     enddo
!
     do kp = k+1,lp1
       do i = 1,ipts
         flx0(i,kp)=flx0(i,kp)+ oss(i,k)*to31d(i,kp-1)                         &
                    +ss2(i,k)*cont1d(i,kp-1)                                   &
                    +css(i,k)*co21(i,k,kp)                                     &
                    +dtc(i,k)*emissb(i,kp-1)
       enddo
     enddo
#endif
   enddo
!
!   now do k=l case. since the kp loop is length 1, many simplifi-
!   cations occur. also, the co2 quantities (as well as the emiss
!  quantities) are computed in the nbl sedction; therefore, we want
!  only over,to3 and cont1d (over(i,l),to31d(i,l) and cont1d(i,l)
!  according to the notation. thus no call is made to the rad_lw_e290
!  subroutine.
!         the third section calculates boundary layer and nearby layer
!     corrections to the transmission functions obtained above. methods
!     are given in ref. (4).
!          the following ratios are used in various nbl calculations:
!
!   the remaining calculations are for :
!                        1) the (k,k) terms, k=2,lm1;
!                        2) the (l,l) term
!                        3) the (l,lp1) term
!                        4) the (lp1,l) term
!                        5) the (lp1,lp1) term.
!     each is uniquely handled; different flux terms are computed
!     differently
!
!
!          fourth section obtains water transmission functions
!     used in q(approx) calculations and also makes nbl corrections:
!     1) emiss (i,j) is the transmission function matrix obtained
!     by calling subroutine e1e288;
!     2) "nearby layer" corrections (emiss(i,i)) are obtained
!     using subroutine rad_lw_e3v88;
!     3) special values at the surface (emiss(l,lp1),emiss(lp1,l),
!     emiss(lp1,lp1)) are calculated.
!
!
!      obtain arguments for e1e288 and rad_lw_e3v88:
!
   do i = 1,ipts
     tpl(i,1)=temp(i,l)
     tpl(i,lp1)=haf*(t(i,lp1)+temp(i,l))
     tpl(i,llp1)=haf*(t(i,l)+temp(i,l))
   enddo
   do k = 2,l
     do i = 1,ipts
       tpl(i,k)=t(i,k)
       tpl(i,k+l)=t(i,k)
     enddo
   enddo
!
!---e2 functions are required in the nbl calculations for 2 cases,
!   denoted (in old code) as (l,lp1) and (lp1,lp1)
!
   do i = 1,ipts
     avephi(i,1)=var2(i,l)
     avephi(i,2)=var2(i,l)+empl(i,l)
   enddo
   !
   call rad_lw_e2spec(ipts,emiss,avephi,fxosp,dtsp)
!
!  call rad_lw_e3v88 for nbl h2o transmissivities
   call rad_lw_e3v88(ipts,emd,tpl,empl)
!
!   compute nearby layer and special-case transmissivities for emiss
!    using methods for h2o given in ref. (4)
!
   do k = 2,l
     do i = 1,ipts
       emisdg(i,k)=emd(i,k+l)+emd(i,k)
     enddo
   enddo
!
!   note that emx1/2 (pressure scaled paths) are now computed in
!   rad_lw_gfdl_driver
!
   do i = 1,ipts
     emspec(i,1)=(emd(i,1)*empl(i,1)-emd(i,lp1)*empl(i,lp1))/                  &
                  emx1(i) + quartr*(emiss(i,1)+emiss(i,2))
     emisdg(i,lp1)=two*emd(i,lp1)
     emspec(i,2)=two*(emd(i,1)*empl(i,1)-emd(i,llp1)*empl(i,llp1))/            &
                  emx2(i)
   enddo
!
   do i = 1,ipts
     fac1(i,l)=bo3rnd(2)*var4(i,l)/var3(i,l)
     vtmp3(i,l)=haf*(fac1(i,l)*                                                &
                  (sqrt(one+(four*ao3rnd(2)*var3(i,l))/fac1(i,l))-one))
     to31d(i,l)=exp(hm1ez*(vtmp3(i,l)+sko3r*cntval(i,l)))
     over1d(i,l)=exp(hm1ez*(sqrt(ab15wd*var2(i,l))+                            &
                  skc1r*cntval(i,l)))
     cont1d(i,l)=cnttau(i,l)*totevv(i,lm1)
     rlog(i,l)=over1d(i,l)*co2nbl(i,l)
   enddo
!
   do k = 1,l
     do i = 1,ipts
       rlog(i,k)=log(rlog(i,k))
     enddo
   enddo
!
   do k = 1,lm1
     do i = 1,ipts
       delpr1(i,k+1)=delp(i,k+1)*(press(i,k+1)-p(i,k+1))
       alp(i,k+l)=-sqrt(delpr1(i,k+1))*rlog(i,k+1)
     enddo
   enddo
!
   do k = 1,l
     do i = 1,ipts
       delpr2(i,k+1)=delp(i,k)*(p(i,k+1)-press(i,k))
       alp(i,k)=-sqrt(delpr2(i,k+1))*rlog(i,k)
     enddo
   enddo
!
   do i = 1,ipts
     alp(i,ll)=-rlog(i,l)
     alp(i,llp1)=-rlog(i,l)*sqrt(delp(i,l)*(p(i,lp1)-press(i,lm1)))
   enddo
!
!        the first computation is for the 15 um band,with the
!     for the combined h2o and co2 transmission function.
!
!       perform nbl computations for the 15 um band
!***the statement function sf in prev. versions is now explicitly
!   evaluated.
!
   do k = 1,llp1
     do i = 1,ipts
       c(i,k)=alp(i,k)*(hmp66667+alp(i,k)*(quartr+alp(i,k)*hm6666m2))
     enddo
   enddo
!
   do i = 1,ipts
     co21(i,lp1,lp1)=one+c(i,l)
     co21(i,lp1,l)=one+(delp2(i,l)*c(i,ll)-(press(i,l)-p(i,l))*                &
                 c(i,llm1))/(p(i,lp1)-press(i,l))
     co21(i,l,lp1)=one+((p(i,lp1)-press(i,lm1))*c(i,llp1)-                     &
                 (p(i,lp1)-press(i,l))*c(i,l))/(press(i,l)-press(i,lm1))
   enddo
!
   do k = 2,l
     do i = 1,ipts
       co21(i,k,k)=one+haf*(c(i,lm1+k)+c(i,k-1))
     enddo
   enddo
!
!    compute nearby-layer transmissivities for the o3 band and for the
!    one-band continuum band (to3 and emiss2). the sf2 function is
!    used. the method is the same as described for co2 in ref (4).
!
   do k = 1,lm1
     do i = 1,ipts
       csub(i,k+1)=cntval(i,k+1)*delpr1(i,k+1)
       csub(i,k+l)=cntval(i,k)*delpr2(i,k+1)
     enddo
   enddo
!
!---the sf2 function in prev. versions is now explicitly evaluated
!
   do k = 1,llm2
     do i = 1,ipts
       csub2(i,k+1)=sko3r*csub(i,k+1)
       c(i,k+1)=csub(i,k+1)*(hmp5+csub(i,k+1)*                                 &
              (hp166666-csub(i,k+1)*h41666m2))
       c2(i,k+1)=csub2(i,k+1)*(hmp5+csub2(i,k+1)*                              &
               (hp166666-csub2(i,k+1)*h41666m2))
     enddo
   enddo
!
   do i = 1,ipts
     contdg(i,lp1)=1.+c(i,llm1)
     to3dg(i,lp1)=1.+c2(i,llm1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       contdg(i,k)=one+haf*(c(i,k)+c(i,lm1+k))
       to3dg(i,k)=one+haf*(c2(i,k)+c2(i,lm1+k))
     enddo
   enddo
!
!---now obtain fluxes
!
!    for the diagonal terms...
!
   do k = 2,lp1
     do i = 1,ipts
       flx(i,k)=flx(i,k)+(dtc(i,k)*emisdg(i,k)                                 &
                    +ss2(i,k)*contdg(i,k)                                      &
                    +oss(i,k)*to3dg(i,k)                                       &
                    +css(i,k)*co21(i,k,k))*cldfac(i,k,k)
     enddo
   enddo
!
!     for the two off-diagonal terms...
!
   do i = 1,ipts
     flx(i,l)=flx(i,l)+(css(i,lp1)*co21(i,lp1,l)                               &
                     +dtc(i,lp1)*emspec(i,2)                                   &
                     +oss(i,lp1)*to31d(i,l)                                    &
                     +ss2(i,lp1)*cont1d(i,l))*cldfac(i,lp1,l)
     flx(i,lp1)=flx(i,lp1)+(css(i,l)*co21(i,l,lp1)                             &
                         +oss(i,l)*to31d(i,l)                                  &
                         +ss2(i,l)*cont1d(i,l)                                 &
                         +dtc(i,l)*emspec(i,1))*cldfac(i,l,lp1)
   enddo
!
! ...   ditto for clear sky ... cldfac =1.
!
#ifdef CLR
   do k = 2,lp1
     do i = 1,ipts
       flx0(i,k)=flx0(i,k)+ dtc(i,k)*emisdg(i,k)                               &
                    +ss2(i,k)*contdg(i,k)                                      &
                    +oss(i,k)*to3dg(i,k)                                       &
                    +css(i,k)*co21(i,k,k)
     enddo
   enddo
#endif
!
!     for the two off-diagonal terms...
!
#ifdef CLR
   do i = 1,ipts
     flx0(i,l)=flx0(i,l)+ css(i,lp1)*co21(i,lp1,l)                             &
                     +dtc(i,lp1)*emspec(i,2)                                   &
                     +oss(i,lp1)*to31d(i,l)                                    &
                     +ss2(i,lp1)*cont1d(i,l)
     flx0(i,lp1)=flx0(i,lp1)+ css(i,l)*co21(i,l,lp1)                           &
                         +oss(i,l)*to31d(i,l)                                  &
                         +ss2(i,l)*cont1d(i,l)                                 &
                         +dtc(i,l)*emspec(i,1)
   enddo
#endif
!
!     final section obtains emissivity heating rates,
!     total heating rates and the flux at the ground
!
!     .....calculate the emissivity heating rates
!
   do k = 1,l
     do i = 1,ipts
       heatem(i,k)=radcon*(flx(i,k+1)-flx(i,k))*delp(i,k)
     enddo
   enddo
!
!     .....calculate the total heating rates
!
   do k = 1,l
     do i = 1,ipts
       heatra(i,k)=heatem(i,k)-cts(i,k)-ctso3(i,k)+excts(i,k)
     enddo
   enddo
!
!     .....calculate the flux at each flux level using the flux at the
!    top (flx1e1+gxcts) and the integral of the heating rates (vsum1)
!
   do k = 1,l
     do i = 1,ipts
       vsum1(i,k)=heatra(i,k)*delp2(i,k)*radcon1
     enddo
   enddo
!
   do i = 1,ipts
     topflx(i)=flx1e1(i)+gxcts(i)
     flxnet(i,1)=topflx(i)
   enddo
!
!---only the surface value of flux (grnflx) is needed unless
!    the thick cloud section is invoked.
!
   do k = 2,lp1
     do i = 1,ipts
       flxnet(i,k)=flxnet(i,k-1)+vsum1(i,k-1)
     enddo
   enddo
!
   do i = 1,ipts
     grnflx(i)=flxnet(i,lp1)
   enddo
!
! ...   ditto for clear sky .. cldfac=1.
!
#ifdef CLR
   do k = 1,l
     do i = 1,ipts
       heatem(i,k)=radcon*(flx0(i,k+1)-flx0(i,k))*delp(i,k)
     enddo
   enddo
#endif
!
!     .....calculate the total heating rates
!
#ifdef CLR
   do k = 1,l
     do i = 1,ipts
       heatr0(i,k)=heatem(i,k)-cts0(i,k)-ctso30(i,k)+excts0(i,k)
     enddo
   enddo
#endif
!
!     .....calculate the flux at each flux level using the flux at the
!    top (flx1e1+gxcts) and the integral of the heating rates (vsum1)
!
#ifdef CLR
   do k = 1,l
     do i = 1,ipts
       vsum1(i,k)=heatr0(i,k)*delp2(i,k)*radcon1
     enddo
   enddo
!
   do i = 1,ipts
     topfx0(i)=flx1e0(i)+gxcts0(i)
     flxnt0(i,1)=topfx0(i)
   enddo
!
#endif
!
!---only the surface value of flux (grnflx) is needed unless
!    the thick cloud section is invoked.
!
#ifdef CLR
   do k = 2,lp1
     do i = 1,ipts
       flxnt0(i,k)=flxnt0(i,k-1)+vsum1(i,k-1)
     enddo
   enddo
!
   do i = 1,ipts
     grnfx0(i)=flxnt0(i,lp1)
   enddo
#endif
!
!     this is the thick cloud section.optionally,if thick cloud
!     fluxes are to be "convectively adjusted",ie,df/dp is constant,
!     for cloudy part of grid point, the following code is executed.
!***first,count the number of clouds along the lat. row. skip the
!   entire thick cloud computation of there are no clouds.
!     icnt=0
!     do 1301 i=1,ipts
!     icnt=icnt+nclds(i)
!301  continue
!     if (icnt.eq.0) go to 6999
!---find the maximum number of clouds in the latitude row
!     kclds=nclds(1)
!     do 2106 i=2,ipts
!     kclds=max(nclds(i),kclds)
!106  continue
!
!
!***obtain the pressures and fluxes of the top and bottom of
!   the ncth cloud (it is assumed that all ktop and kbtms have
!   been defined!).
!     do 1361 kk=1,kclds
!     kmin=lp1
!     kmax=0
!     do 1362 i=1,ipts
!       j1=ktop(i,kk+1)
!       if (j1.eq.1) go to 1362
!       j3=kbtm(i,kk+1)
!       if (j3.gt.j1) then
!         ptop(i)=p(i,j1)
!         pbot(i)=p(i,j3+1)
!         ftop(i)=flxnet(i,j1)
!         fbot(i)=flxnet(i,j3+1)
!***obtain the "flux derivative" df/dp (delptc)
!         delptc(i)=(ftop(i)-fbot(i))/(ptop(i)-pbot(i))
!         kmin=min(kmin,j1)
!         kmax=max(kmax,j3)
!       endif
!362  continue
!     kmin=kmin+1
!***calculate the tot. flux chg. from the top of the cloud, for
!   all levels.
!     do 1365 k=kmin,kmax
!     do 1363 i=1,ipts
!       if (ktop(i,kk+1).eq.1) go to 1363
!       if(ktop(i,kk+1).lt.k .and. k.le.kbtm(i,kk+1)) then
!         z1(i,k)=(p(i,k)-ptop(i))*delptc(i)+ftop(i)
!original flxnet(i,k)=flxnet(i,k)*(one-camt(i,kk+1)) +
!original1            z1(i,k)*camt(i,kk+1)
!         flxnet(i,k)=z1(i,k)
!       endif
!363  continue
!365  continue
!361  continue
!***using this flux chg. in the cloudy part of the grid box, obtain
!   the new fluxes, weighting the clear and cloudy fluxes:again, only
!    the fluxes in thick-cloud levels will eventually be used.
!     do 6051 k=1,lp1
!     do 6051 i=1,ipts
!     flxnet(i,k)=flxnet(i,k)*(one-camt(i,nc)) +
!    1            z1(i,k)*camt(i,nc)
!051  continue
!***merge flxthk into flxnet for appropriate levels.
!     do 1401 k=1,lp1
!     do 1401 i=1,ipts
!     if (k.gt.itop(i) .and. k.le.ibot(i)
!    1  .and.  (nc-1).le.nclds(i))  then
!          flxnet(i,k)=flxthk(i,k)
!     endif
!401  continue
!
!******end of cloud loop****!
!6001  continue
!6999  continue
!***the final step is to recompute the heating rates based on the
!   revised fluxes:
!      do 6101 k=1,l
!      do 6101 i=1,ipts
!      heatra(i,k)=radcon*(flxnet(i,k+1)-flxnet(i,k))*delp(i,k)
!6101  continue
!     the thick cloud section ends here.
!
!
   return
   end subroutine rad_lw
