#include <define.h>
   subroutine phys_cps_ccm3(imx2,imx22,kmx,                                    &
                 jcap,delt,nstep,prsi,prsl,                                    &
                 q1,t1,cldwrk,rn,kbot,ktop,icps,spd,lat,slimsk,dot,            &
                 hpbl,gamt,gamq,delx)
#ifdef CCMCNV
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
   use comadj
   use phys_ccm_module, only : ccm3_initialize,       ccm3_reset,              &
                               ccm3_virtual_temp,     ccm3_height8level,       &
                               ccm3_height8interface, ccm3_cloud_main
!-------------------------------------------------------------------------------
!
! abstract :  ncar ccm3 deep convection
!
!  ::: structure :::
!
!    [phys_cps_ccm3]
!      |
!      |---[phys_ccm_module]
!               |
!               |--- [ccm3_initailize] ----- [ccm3_constant_setup]  
!               |                       |--- [esinti]  
!               |                               |- [gestbl]  
!               |--- [ccm3_reset]               |- [gffgch]  
!               |--- [ccm3_virtual_temp]  
!               |--- [ccm3_height8level]  
!               |--- [ccm3_height8interface]  
!               |--- [ccm3_cloud_main] ----- [ccm3_cloud_buoyancy]  
!                                       |--- [ccm3_cloud_trigger]  
!                                       |--- [ccm3_cloud_property]  
!                                       |--- [ccm3_cloud_closure]  
!                                       |--- [ccm3_cloud_trans]  
!                                       |--- [ccm3_cloud_q1q2]  
!
! program history log:
!   1994-01-01  g. zhang               initial development
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
!
! input arguments
!
   integer  ::  nstep      ! current time step for ccm2_saturation_driver diag.
   integer  ::  lat                       ! latitude index (s->n)
   real     ::  tv(ILOTS,levs_)           ! virtual temperature array
   real     ::  zi(ILOTS,levs_+1)         ! height above sfc in interface
   real     ::  zm(ILOTS,levs_)           ! height above sfc in mid-layer
   real     ::  dlf(ILOTS,levs_)          ! detraining cld h20 from convection
   real     ::  pflx(ILOTS,levs_)         ! conv rain flux thru out btm of lev
   real     ::  phis(ILOTS)               ! surface geopotential
   real     ::  pmid(ILOTS,levs_)         ! pressures at model levels
   real     ::  pint(ILOTS,levs_+1)       ! pressure at model interfaces
   real     ::  pdel(ILOTS,levs_)         ! delta-p (pressure depth)
   real     ::  rpdel(ILOTS,levs_)        ! 1./pdel
   real     ::  pmln(ILOTS,levs_)         ! log(pmid)
   real     ::  piln(ILOTS,levs_+1)       ! log(pint)
   real     ::  prsi(ILOTS,levs_+1)       
   real     ::  prsl(ILOTS,levs_)          
!
! input/output arguments
!
   real     ::  t(ILOTS,levs_)            ! temperature (k)
   real     ::  q(ILOTS,levs_,1)    ! constituent mixing ratio field
   real     ::  ts(ILOTS)                 ! surface temperature
!
! output arguments
!
   real     ::  precc(ILOTS)              ! convective-scale preciptn rate
   real     ::  cnt(ILOTS)                ! top level of convective activity
   real     ::  cnb(ILOTS)                ! bottom level of convective activity
   real     ::  cmfdt(ILOTS,levs_)        ! dt/dt due to moist convection
   real     ::  cmfdq(ILOTS,levs_)        ! dq/dt due to moist convection
   real     ::  zmdt(ILOTS,levs_)         ! zhang convective temperature tendency
   real     ::  zmdq(ILOTS,levs_)         ! zhang convective moisture tendency
   real     ::  cmfdqr(ILOTS,levs_)       ! dq/dt due to moist convective rainout 
   real     ::  cmfmc(ILOTS,levs_)        ! moist convection cloud mass flux
   real     ::  cmfsl(ILOTS,levs_)        ! moist convection lw stat energy flux
!
   real     ::  q1(imx22,kmx),t1(imx22,kmx),rn(imx22)
   real     ::  spd(imx22,kmx),slimsk(imx22),dot(imx22,kmx),cldwrk(imx22)
   real     ::  hpbl(imx22),gamt(imx22),gamq(imx22)
#ifdef EXPLICIT_CLOUDINESS
   real     ::  qci(imx22,kmx),qrs(imx22,kmx)
#endif
   integer  ::  kbot(imx22),ktop(imx22),icps(imx22)
!
   real     ::  pmin
   data pmin /1.e-20/   ! minimum precip rate (ms-1)
!
!---------------------------local workspace-----------------------------
!
! most variables with appended 2 are a second copy of similar quantities
! declared above to accommodate calls to two convection schemes
!
   integer  ::  i,k                ! lon, lev, constituent indices
!-------------------------------------------------------------------------------
!
! define thermodynamic constansts
!
   call ccm3_initialize
!
! zero out precip and convective fields before accumulating terms
!
   call ccm3_reset(precc,ILOTS,0.)
   call ccm3_reset(cmfdt ,ILOTS*levs_,0.)
   call ccm3_reset(cmfdq ,ILOTS*levs_,0.)
   call ccm3_reset(zmdt  ,ILOTS*levs_,0.)
   call ccm3_reset(zmdq  ,ILOTS*levs_,0.)
   call ccm3_reset(cmfdqr,ILOTS*levs_,0.)
   call ccm3_reset(cmfmc ,ILOTS*levs_,0.)
   call ccm3_reset(cmfsl ,ILOTS*levs_,0.)
!
!!!! note the ccmcnv is top down, whereas ncep is bottom up 
!
! convert surface pressure to pascal from cb
!
   do i = 1,ILOTS
     pint(i,levs_+1) = 1000. * prsi(i,1)
   enddo
!
   do k = 1,levs_
     kk =  levs_ - k + 2
     do i = 1, ILOTS
       pint(i,k) = 1000. * prsi(i,kk)
     enddo
   enddo
!
   pmid = 0.  
   do k = 1,levs_
     do i = 1, ILOTS
       pmid(i,k) = 0.5*(pint(i,k)+pint(i,k+1))
     enddo
   enddo
!
   do k = 1,levs_
     kk = levs_ - k + 1
     do i = 1, ILOTS
       t(i,k) = t1(i,kk)
       q(i,k,1) = q1(i,kk)
     enddo
   enddo
!
! set layer thicknesses
!
   do k = 1,levs_
     do i = 1,ILOTS
       pdel(i,k) = pint(i,k+1) - pint(i,k)
     end do
   end do
!
   do k = 1,levs_
     do i = 1,ILOTS
       pmln(i,k) = log(pmid(i,k))
       rpdel(i,k) = 1./pdel(i,k)
       piln(i,k) = log(pint(i,k))
     end do
   end do
!
   do i = 1,ILOTS
     piln(i,levs_+1) = log(pint(i,levs_+1))
   end do
!
   do i = 1,ILOTS
     phis(i) = 0.  ! assume zero terrain
     ts(i) = 400.  ! no check for snow or rain
   end do
!
! calculate geopotential height 
!
   call ccm3_virtual_temp(t       ,q       ,zvir    ,tv      )
   call ccm3_height8level(piln(1,levs_+1),pmln    ,rair    ,                   &
                          gravit   ,tv  , zm      )
   call ccm3_height8interface(piln    ,pmln    ,rair    ,                      &
                              gravit  ,tv      ,zm      ,zi      ) 
   call ccm3_cloud_main(t      ,q        ,precc   ,cnt     ,cnb    ,           &
                   hpbl   ,zm       ,phis    ,zi      ,zmdq   ,                &
                  zmdt   ,pmid     ,pint    ,pdel    ,ts     ,                 &
                  delt   ,cmfmc    ,cmfdqr  ,nstep   ,lat    ,                 &
                  gamt   ,gamq    ,dlf     ,pflx     ,cldwrk)
!
! output
!
   do i = 1,ILOTS
     kbot(i) = levs_ - nint(cnb(i)) + 1
     ktop(i) = levs_ - nint(cnt(i)) + 1
   end do
!
   do i = 1,ILOTS
     if(precc(i).gt.pmin.and.kbot(i).lt.ktop(i)) then
       rn(i) = precc(i) * delt * 2.0   ! m
       icps(i) = 1
     else
       rn(i) = 0.
       icps(i) = 0
     endif
   end do
!
   do k = 1,levs_
     kk = levs_ - k + 1
     do i = 1, ILOTS
       if(icps(i).eq.1) then
         t1(i,k) = t(i,kk)
         q1(i,k) = q(i,kk,1)
       endif
     enddo
   enddo
#ifdef DBG
   i = ILOTS/4
   write(6,*)
   write(6,*) ' cape precip/hr ', cldwrk(i),precc(i) * 3600.
   write(6,*)
   write(6,*)' kbot ktop icps ', kbot(i),ktop(i),icps(i)
   write(6,*)
   write(6,*)
   write(6,145)
 145  format('lev' ,8x,'pres',4x,' t  ',7x,'  q  ',5x,'mass flx',              &
        3x,' t  tend',2x,' q  tend')
   write(6,146)
 146  format('   ' ,8x,' pa ',4x,' k  ',7x,'kg/kg',5x,' kg/m^2s',              &
        3x,' deg/s  ',2x,' kg/m^2s')
   write(6,*)
   do k = 1,levs_
     write(6,155)k,pmid(i,k),t(i,k),q(i,k,1),cmfmc(i,k),                      &
           zmdt(i,k),zmdq(i,k)
   enddo
 155  format(i2,5x,f9.2,1x,f8.2,5(2x,g10.3),3x,a)
#endif
!
   return
#endif /* CCMCNV end */
   end subroutine phys_cps_ccm3
