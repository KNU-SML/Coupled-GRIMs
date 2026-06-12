#include <define.h>
   subroutine lsm_vic1_main(ims2,imx2,msl,lat,                                 &
                  t1,q1,snoweq,tskin,qsurf,                                    &
                  smcld,stcld,dm,sigmaf,vegtype,canopy,                        &
                  dlwflx,dswflx,snowmt,snowev,delt,z0,tg3,                     &
                  gflux,zsoil,cm,ch,                                           &
                  prsi1,prsl1,prsik1,prslk1,zl1,slimsk,                        &
                  drain,evap,hflx,ep,wind,                                     &
                  nsl,month,binf,ds,dsm,ws,cef,                                &
                  expld,kstld,dphld,bubld,qrtld,bkdld,sldld,wcrld,             &
                  wpwld,smrld,smxld,sicld,dpnld,sxnld,epnld,bbnld,             &
                  apnld,btnld,gmnld,flaild,vrtld,lstsn,                        &
                  silz,snwz,csno,rsno,tsfld,tpkld,sfwld,pkwld,gheat,           &
#ifdef SMP_RA2SFC
                  slrad, snowfl, runoff, precip)
#else
                  snowfl, runoff, precip)
#endif
!-------------------------------------------------------------------------------
!
!  notes: vic1 uses one vegetation cover. we will consider multi vegcov
!         in vic2. to simplify coupling vic1/vic2 to gmp/rmp, we use the
!         following local variables to convert the related variables to
!         the uwvic variables. ecpc/ji
!
!  the following definition is for driving vic
!
!  ::: structure :::
!
!    [lsm_driver] --- [lsm_vic1_main] * ----- [phys_lsm_vic1] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use constant, only  :  cal_,rvordm1=>rvordm1_,rd=>rd_,eps=>rdorv_,          &
                          epsm1=>rdorvm1_, cp=>cp_,g=>g_,elocp=>elocp_,        &
                          sigma=>sbc_,hvap=>hvap_,convrad_,                    &
                          qmin8=>qmin8_,cb2pa=>cb2pa_
!
   integer            ::  lat
!
! input variables
!
! model basic parameter
!
   integer  ::  ims2, imx2
   integer  ::  msl         ! number of soil layers
   integer  ::  nsl         ! number of soil nodes
!
!  msub             ::  maximum number of subgrid points
!  md               ::  md (1no/2with) precipiation distribution
!
   integer, parameter  ::  msub=2  ! for viclsm1 (=2:veg,noveg)
   integer, parameter  ::  md=1    ! for viclsm1
!
   real     ::  dtime       ! vic time step (second)
   integer  ::  month       ! current month
!
! -- atmosphere data
!
   real     ::  prc(imx2)      ! precipitation rate (m h2o/time step)
   real     ::  pgcm(imx2)     ! atm bottom level pressure (pa)
   real     ::  wind(imx2)     ! atm bottom level wind speed (m/s)
   real     ::  tgcm(imx2)     ! atm bottom level temperature (kelvin)
   real     ::  qgcm(imx2)     ! atm bottom level specific humidity (kg/kg)
   real     ::  zgcm(imx2)     ! atm bottom level height above surface (m)
   real     ::  flwds(imx2)    ! downward longwave rad onto surface (w/m2)
   real     ::  sols(imx2)     ! solar rad onto surface (w/m2)
#ifdef SMP_RA2SFC
   real     ::  slrad(imx2)     ! solar rad onto surface (radiation uint)
   real     ::  netsw(imx2)     ! solar rad onto surface (w/m2)
#endif
!
! -- soil parameters
!
   real     ::  binf(imx2)       ! variable infiltration curve parameter (n/a)
   real     ::  ds(imx2)         ! fract of dsm nonlinear baseflow begins (fract)
   real     ::  dsm(imx2)        ! maximum velocity of baseflow (mm/day)
   real     ::  ws(imx2)         ! fract maxi sm nonlinear baseflow occurs (fract)
   real     ::  cef(imx2)        ! exponent used in infiltration curve (n/a)
   real     ::  silz(imx2)       ! surface roughness of bare soil (m)
   real     ::  snwz(imx2)       ! surface roughness of snowpack (m)

   real     ::  expt(msl,imx2)   ! para the vari of ksat with sm (n/a)
   real     ::  kst(msl,imx2)    ! saturated hydrologic conductivity (mm/day)
   real     ::  dph(msl,imx2)    ! thickness of soil layer (m)
   real     ::  bub(msl,imx2)    ! bubbling pressure of soil layer (cm)
   real     ::  qrt(msl,imx2)    ! quartz content of soil layer (fraction)
   real     ::  bkd(msl,imx2)    ! bulk density of soil layer (kg/m3)
   real     ::  sld(msl,imx2)    ! soil density of soil layer (kg/m3)
   real     ::  wcr(msl,imx2)    ! sm content at the critical point (mm)
   real     ::  wpw(msl,imx2)    ! sm content wilting point (mm)
   real     ::  smr(msl,imx2)    ! soil moisture residual moisture (mm)
   real     ::  smx(msl,imx2)    ! maximum soil moisture (mm)
   real     ::  dphnd(nsl,imx2)  ! thickness of soil node (m)
   real     ::  smxnd(nsl,imx2)  ! maximum soil moisture at soil node (m3/m3)
   real     ::  expnd(nsl,imx2)  ! para the vari of ksat at soil node (n/a)
   real     ::  bubnd(nsl,imx2)  ! bubbling pressure at soil node (cm)
   real     ::  alpnd(nsl,imx2)  ! para alpha at soil node
   real     ::  betnd(nsl,imx2)  ! para beta at soil node
   real     ::  gamnd(nsl,imx2)  ! para gamma at soil node
!
! -- vegetation parameters
!
   integer  ::  nveg(imx2)      ! number of vegetation type in a grid cell
   integer  ::  nvt(msub,imx2)  ! vegetation type number
   real     ::  wt(msub,imx2)      ! fraction of grid cell covered by veg
   real     ::  rt(msl,msub,imx2)  ! fraction of root in the soil layer
   real     ::  flai(msub,imx2)    ! leaf area index
!
! modified variables 
!
   real     ::  cwt(md,msub,imx2)     ! canopy water (m h2o)
   real     ::  csn(md,msub,imx2)     ! canopy snow  (m h2o)
   real     ::  smc(msl,md,msub,imx2) ! soil moisture content (mm)
   real     ::  sic(msl,md,msub,imx2) ! soil ice content (mm)
   real     ::  tnd(nsl,md,msub,imx2) ! temperature at soil nodes (k)
   real     ::  swq(md,msub,imx2)     ! snow water equivalent (m h2o)
   real     ::  rsn(md,msub,imx2)     ! snow density (kg/m^3)
   real     ::  tsf(md,msub,imx2)     ! snow surface temperature (k)
   real     ::  tpk(md,msub,imx2)     ! snow pack temperature (k)
   real     ::  sfw(md,msub,imx2)     ! surface snow water equivalent (m h2o)
   real     ::  pkw(md,msub,imx2)     ! snow pack snow water equivalent (m h2o)
   real     ::  fmu(imx2)             ! precipitation fraction
   integer  ::  lstsn(imx2)        ! number of time step since last snow
!
! output variables
!
   real     ::  sh(imx2)       ! sensible heat flux (w/m**2) [+ to ground]
   real     ::  lh(imx2)       ! latent heat flux (w/m**2) [+ to ground]
   real     ::  gh(imx2)       ! ground heat flux (w/m**2) [+ to ground]
   real     ::  ts(imx2)       ! surface radiative temperature (kelvin)
   real     ::  albd(imx2)     ! albedo (fraction)
   real     ::  ovflw(imx2)    ! overland flow (mm/time step)
   real     ::  bsflw(imx2)    ! base flow (mm/time step)
   real     ::  snowmt(imx2)   ! snow melt (m/s)
   real     ::  snowev(imx2)   ! evaporation over snow surface (w/m2)
   real     ::  snowfl(imx2)   ! snow fall (m/time step)
   real     ::  gdh(imx2)      ! ground heat storage (W/m**2) [+ to ground]
!
! passing array
!
   integer  ::   vegtype(imx2)

   real     ::   prsi1(imx2),prsl1(imx2),prsik1(imx2),prslk1(imx2)
   real     ::   zl1(imx2),t1(imx2),q1(imx2),snoweq(imx2)
   real     ::   tskin(imx2),qsurf(imx2),dm(imx2),smcld(imx2,msl)
   real     ::   stcld(imx2,nsl),sigmaf(imx2)
   real     ::   canopy(imx2),dlwflx(imx2),dswflx(imx2)
   real     ::   z0(imx2),tg3(imx2),gflux(imx2)
   real     ::   zsoil(imx2,msl),cm(imx2),ch(imx2),slimsk(imx2)
   real     ::   drain(imx2),evap(imx2),hflx(imx2),ep(imx2)
   real     ::   expld(imx2,msl),kstld(imx2,msl),dphld(imx2,msl)
   real     ::   bubld(imx2,msl),qrtld(imx2,msl),bkdld(imx2,msl)
   real     ::   sldld(imx2,msl),wcrld(imx2,msl),wpwld(imx2,msl)
   real     ::   smrld(imx2,msl),smxld(imx2,msl),sicld(imx2,msl)
   real     ::   dpnld(imx2,nsl),sxnld(imx2,nsl),epnld(imx2,nsl)
   real     ::   bbnld(imx2,nsl),apnld(imx2,nsl),btnld(imx2,nsl)
   real     ::   gmnld(imx2,nsl),flaild(imx2),vrtld(imx2,msl)
   real     ::   csno(imx2),rsno(imx2)
   real     ::   tsfld(imx2),tpkld(imx2),sfwld(imx2),pkwld(imx2)
   real     ::   sfall(imx2),runoff(imx2),precip(imx2)
   real     ::   gheat(imx2)
!
! local array
!
   real     ::   psurf(imx2),q0(imx2),theta1(imx2)
   real     ::   tv1(imx2),rho(imx2),qs1(imx2),rch(imx2)
!
!  surface energy/water balance over land and seaice
!
!  initialization
!
   im = ims2
!
!  initialize variables. all units are supposedly m.k.s. unless specifie
!  psurf is in pascals
!  wind is wind speed, theta1 is adiabatic surface temp from level 1
!  rho is density, qs1 is sat. hum. at level1 and qss is sat. hum. at surface
!  surface roughness length is converted to m from cm
!
#ifdef DBG
   print*,' enter lsm_vic1_main '
#endif
!
   do i = 1,im
     if(slimsk(i).eq.1.) then       
       q0(i) = max(q1(i),qmin8)
       theta1(i) = t1(i) / prslk1(i) * prsik1(i)
       tv1(i) = t1(i) * (1. + rvordm1 * q0(i))
       rho(i) = (cb2pa * prsl1(i)) / (rd * tv1(i))
#ifdef ICE
       qs1(i) = cb2pa * fpvs(t1(i))
#else 
       qs1(i) = cb2pa * fpvs0(t1(i))
#endif
       qs1(i) = eps*qs1(i)/(cb2pa*prsl1(i)+epsm1*qs1(i))
       qs1(i) = max(qs1(i), qmin8)
       q0(i) = min(qs1(i),q0(i))
!
!  rcp = rho cp ch v
!
       rch(i) = rho(i) * cp * ch(i) * wind(i)
     endif
   enddo
!
!     gmp/rmp -> vic: prepare variables to run vic lsm
!
! 1. model basic parameters
!
! 1.1  msl      number of soil layers
! 1.2  nsl      number of soil nodes
! 1.3  msub     number of subgrids
! 1.4  md       number of wet and dry 
!
! 1.5  dtime    vic lsm time step (sec)
!
   dtime = delt
!
! 1.6  mm       current month
!      
   mm = month
!
! 2. atmosphere forcings
!
   do i = 1,im
     if(slimsk(i).eq.1.) then       
!
! 2.1  prc      precipitation rate (m h2o/time step)
!     
       prc(i) = precip(i)
!
! 2.2  pgcm     atm bottom level pressure (pa)
!
       pgcm(i) = cb2pa*prsl1(i)
!
! 2.3  wind     atm bottom level wind speed (m/s)
!      wind(i) = wind(i)
!
! 2.4  tgcm     atm bottom level temperature (kelvin)
!
       tgcm(i) = t1(i)
!
! 2.5  qgcm     atm bottom level specific humidity (kg/kg)
!
       qgcm(i) = q0(i)
!
! 2.6  zgcm     atm bottom level height above surface (m)
!
       zgcm(i) = -rd * tv1(i) * log(prsl1(i)/prsi1(i)) / g
!
! 2.7  flwds    downward longwave rad onto surface (w/m2)
!               (positive, not net longwave)
!
       flwds(i) = dlwflx(i)
!
! 2.8  sols     solar downward radiation (w m-2; positive, not net solar)
!
       sols(i) = dswflx(i)
#ifdef SMP_RA2SFC
!    
! 2.9  netsw    net solar radiation at the surface (w m-2; positive)
!               convert unit from radiation uint to w m-2
!
       slrad(i) = slrad(i)*convrad_*(-1.)
       netsw(i) = slrad(i)-dlwflx(i)
#endif
     end if
   end do
!
! 3. soil parameters
!
!      binf ,ds   ,dsm  ,ws   ,cef  ,expt ,kst  ,dph  ,bub  ,qrt  ,
!      bkd  ,sld  ,wcr  ,wpw  ,silz ,snwz ,smr  ,smx  ,dphnd,smxnd,
!      expnd,bubnd,alpnd,betnd,gamnd
!
   do i = 1,im
     if(slimsk(i).eq.1.) then       
       do k = 1,msl
         expt(k,i)= expld(i,k)
         kst(k,i) = kstld(i,k)
         dph(k,i) = dphld(i,k)
         bub(k,i) = bubld(i,k)
         qrt(k,i) = qrtld(i,k)
         bkd(k,i) = bkdld(i,k)
         sld(k,i) = sldld(i,k)
         wcr(k,i) = wcrld(i,k)
         wpw(k,i) = wpwld(i,k)
         smr(k,i) = smrld(i,k)
         smx(k,i) = smxld(i,k)
       end do
       do k = 1,nsl
         dphnd(k,i) = dpnld(i,k)
         smxnd(k,i) = sxnld(i,k)
         expnd(k,i) = epnld(i,k)
         bubnd(k,i) = bbnld(i,k)
         alpnd(k,i) = apnld(i,k)
         betnd(k,i) = btnld(i,k)
         gamnd(k,i) = gmnld(i,k)
       end do
     endif
   enddo
!
! 4. vegetation parameters
!
! 4.1  nveg(imx2)              number of vegetation type in a grid cell
! 4.2  nvt(msub,imx2)          vegetation type
! 4.3  wt(msub,imx2)           fraction of grid cell covered by veg
! 4.4  rt(msl,msub,imx2)       fraction root in the soil layer
! 4.5  flai(msub,imx2)         leaf area index
!
   do i = 1,im
     if(slimsk(i).eq.1.) then       
       nveg(i) = 1
!
       if(sigmaf(i).gt.1.e-5.and.sigmaf(i).lt.(1.0-1.e-5))then
         if(vegtype(i).gt.0 .and. vegtype(i).lt.12) then
           nveg(i) = 2
         end if
       end if
       if(nveg(i).eq.2) then
         nvt(1,i) = vegtype(i)
         nvt(2,i) = 12                ! vic bare soil
         wt(1,i)  = sigmaf(i)
         wt(2,i)  = 1.0 - sigmaf(i)
         do k = 1,msl
           rt(k,1,i) = vrtld(i,k)
           rt(k,2,i) = 0.
         end do
         flai(1,i) = flaild(i)
         flai(2,i) = 0.0
       else
         if(vegtype(i).lt.12) then
           nvt(1,i) = vegtype(i)
         else
           nvt(1,i) = 12
         endif
         wt(1,i)  = 1.0
         do k = 1,msl
           rt(k,1,i) = vrtld(i,k)
         end do
         flai(1,i) = flaild(i)
       end if
     end if
   end do
!
! 5. modified variables 
!
! 5.1  cwt(md,msub,imx2)      canopy water (mm h2o)
! 5.2  csn(md,msub,imx2)      canopy snow  (mm h2o)
! 5.3  smc(msl,md,msub,imx2)  soil moisture content (mm)
! 5.4  sic(msl,md,msub,imx2)  soil ice content (mm)
! 5.5  tnd(nsl,md,msub,imx2)  temperature at soil nodes (k)
! 5.6  swq(md,msub,imx2)      snow water equivalent (mm h2o)
! 5.7  rsn(md,msub,imx2)      snow density (kg/m^3)
! 5.8  tsf(md,msub,imx2)      snow surface temperature (k)
! 5.9  tpk(md,msub,imx2)      snow pack temperature (k)
! 5.a  sfw(md,msub,imx2)      surface snow water equivalent (mm h2o)
! 5.b  pkw(md,msub,imx2)      snow pack snow water equivalent (mm h2o)
! 5.c  fmu(imx2)              precipitation fraction
! 5.d  lstsn(imx2)            number of time step since last snow
!
   do i = 1,im
     if(slimsk(i).eq.1.) then 
!            do jd = 1,md
!               cwt(jd,1,i) = canopy(i)/1000.0   ! mm -> m
!               csn(jd,1,i) = csno(i)            ! m
!            end do
       do nv = 1,nveg(i)
         do jd = 1,md
           do m = 1,msl
             smc(m,jd,nv,i) = smcld(i,m)*dph(m,i)*1000.0  ! mm
             sic(m,jd,nv,i) = sicld(i,m)
           end do
           do n = 1,nsl
             tnd(n,jd,nv,i) = stcld(i,n)
           end do
           swq(jd,nv,i) = snoweq(i)/1000.0    ! mm -> m
           rsn(jd,nv,i) = rsno(i)
!
           tsf(jd,nv,i) = tsfld(i)
           tpk(jd,nv,i) = tpkld(i)
           sfw(jd,nv,i) = sfwld(i)
           pkw(jd,nv,i) = pkwld(i)
           if(nv.eq.1)then
             cwt(jd,nv,i) = canopy(i)/1000.0
             csn(jd,nv,i) = csno(i)
           else
             cwt(jd,nv,i) = 0.0
             csn(jd,nv,i) = 0.0
           endif
         end do
       end do
     end if
     fmu(i) = 1.0
   end do
!
!    call vic lsm
!
#ifdef DBG
   print*,' --- in lsm_vic1_main --- before phys_lsm_vic1'
   print*,'lat=',lat,' im =',im
#endif
   do i = 1,im
     if(slimsk(i) .eq. 1.) then
       sh(i) = 0.0
       lh(i) = 0.0
       gh(i) = 0.0
       ts(i) = 0.0
       albd(i) = 0.0
       ovflw(i) = 0.0
       bsflw(i) = 0.0
       sfall(i) = 0.0
       gdh(i) = 0.0
!
#ifdef DBG
       print*,'before phys_lsm_vic1'
       print*, i,     msl,          nsl,         msub,         md,             &
                    dtime,           mm,                                       &
                   prc(i),      pgcm(i),      wind(i),    tgcm(i),             &
                  qgcm(i),      zgcm(i),     flwds(i),    sols(i),             &
                  binf(i),        ds(i),       dsm(i),      ws(i),             &
                expt(1,i),     kst(1,i),     dph(1,i),     cef(i),             &
                 bub(1,i),     qrt(1,i),     bkd(1,i),   sld(1,i),             &
                 wcr(1,i),     wpw(1,i),      silz(i),    snwz(i),             &
                 smr(1,i),     smx(1,i),   dphnd(1,i), smxnd(1,i),             &
               expnd(1,i),   bubnd(1,i),   alpnd(1,i), betnd(1,i),             &
               gamnd(1,i),      nveg(i),     nvt(1,i),    wt(1,i),             &
                rt(1,1,i),    flai(1,i),   cwt(1,1,i), csn(1,1,i),             &
             smc(1,1,1,i), sic(1,1,1,i), tnd(1,1,1,i), swq(1,1,i),             &
               rsn(1,1,i),   tsf(1,1,i),   tpk(1,1,i), sfw(1,1,i),             &
               pkw(1,1,i),       fmu(i),     sfall(i),   lstsn(i),             &
                    sh(i),        lh(i),        gh(i),      ts(i),             &
                  albd(i),     ovflw(i),     bsflw(i),  snowmt(i),             &
                snowev(i)
#endif
!
       call phys_lsm_vic1(msl,      nsl,         msub,         md,             &
                    dtime,           mm,            i,        lat,             &
                   prc(i),      pgcm(i),      wind(i),    tgcm(i),             &
                  qgcm(i),      zgcm(i),     flwds(i),    sols(i),             &
                  binf(i),        ds(i),       dsm(i),      ws(i),             &
                expt(1,i),     kst(1,i),     dph(1,i),     cef(i),             &
                 bub(1,i),     qrt(1,i),     bkd(1,i),   sld(1,i),             &
                 wcr(1,i),     wpw(1,i),      silz(i),    snwz(i),             &
                 smr(1,i),     smx(1,i),   dphnd(1,i), smxnd(1,i),             &
               expnd(1,i),   bubnd(1,i),   alpnd(1,i), betnd(1,i),             &
               gamnd(1,i),      nveg(i),     nvt(1,i),    wt(1,i),             &
                rt(1,1,i),    flai(1,i),   cwt(1,1,i), csn(1,1,i),             &
             smc(1,1,1,i), sic(1,1,1,i), tnd(1,1,1,i), swq(1,1,i),             &
               rsn(1,1,i),   tsf(1,1,i),   tpk(1,1,i), sfw(1,1,i),             &
               pkw(1,1,i),       fmu(i),     sfall(i),   lstsn(i),             &
                    sh(i),        lh(i),        gh(i),      ts(i),             &
                  albd(i),     ovflw(i),     bsflw(i),  snowmt(i),             &
#ifdef SMP_RA2SFC
                snowev(i),       gdh(i),      z0(i),   netsw(i))
#else
                snowev(i),       gdh(i),      z0(i))
#endif
!
#ifdef DBG
       print *,'in lsm_vic1_main i=',i,' lstsn=',lstsn(i)
       print *,'   tgcm=',tgcm(i),' ts=',ts(i)
       print *,'   tsf =',tsf(1,1,i),' tnd=',tnd(1,1,1,i)
       print*,'after phys_lsm_vic1'
       print*, i,     msl,          nsl,         msub,         md,             &
                    dtime,           mm,                                       &
                   prc(i),      pgcm(i),      wind(i),    tgcm(i),             &
                  qgcm(i),      zgcm(i),     flwds(i),    sols(i),             &
                  binf(i),        ds(i),       dsm(i),      ws(i),             &
                expt(1,i),     kst(1,i),     dph(1,i),     cef(i),             &
                 bub(1,i),     qrt(1,i),     bkd(1,i),   sld(1,i),             &
                 wcr(1,i),     wpw(1,i),      silz(i),    snwz(i),             &
                 smr(1,i),     smx(1,i),   dphnd(1,i), smxnd(1,i),             &
               expnd(1,i),   bubnd(1,i),   alpnd(1,i), betnd(1,i),             &
               gamnd(1,i),      nveg(i),     nvt(1,i),    wt(1,i),             &
                rt(1,1,i),    flai(1,i),   cwt(1,1,i), csn(1,1,i),             &
             smc(1,1,1,i), sic(1,1,1,i), tnd(1,1,1,i), swq(1,1,i),             &
               rsn(1,1,i),   tsf(1,1,i),   tpk(1,1,i), sfw(1,1,i),             &
               pkw(1,1,i),       fmu(i),     sfall(i),   lstsn(i),             &
                    sh(i),        lh(i),        gh(i),      ts(i),             &
                  albd(i),     ovflw(i),     bsflw(i),  snowmt(i),             &
                snowev(i)
       print*,'- in lsm_vic1_main - after phys_lsm_vic1 ts (i) =',ts(i),sh(i)
#endif
     endif 
   enddo
!
! 6.  update the modified variables
!
   do i = 1,im
     if(slimsk(i).eq.1.) then 
       canopy(i) = 0.
       csno(i)   = 0.
       do m = 1,msl
         smcld(i,m) = 0.0
         sicld(i,m) = 0.0
       end do
       do n = 1,nsl
         stcld(i,n) = 0.0
       end do
       snoweq(i)= 0.0
       rsno(i)  = 0.0
       tsfld(i) = 0.0
       tpkld(i) = 0.0
       sfwld(i) = 0.0
       pkwld(i) = 0.0

       do jd = 1,md           ! since just one-tile
         if(jd.eq.1) then
           canopy(i) = canopy(i) + cwt(jd,1,i)*fmu(i)
           csno(i)   = csno(i) + csn(jd,1,i)*fmu(i)
         else
           canopy(i) = canopy(i) + cwt(jd,1,i)*(1-fmu(i))
           csno(i)   = csno(i) + csn(jd,1,i)*(1-fmu(i))
         endif
       enddo
!
       do nv = 1,nveg(i)
         do jd = 1,md
           if(jd.eq.1) then
             wtld = fmu(i)*wt(nv,i)
           else
             wtld = (1.0-fmu(i))*wt(nv,i)
           end if
!
!                  canopy(i) = canopy(i) + cwt(jd,nv,i)*wtld
!                  csno(i)   = csno(i) + csn(jd,nv,i)*wtld
!
           do m = 1,msl
             if(dph(m,i).gt.0) then
               smc(m,jd,nv,i)=smc(m,jd,nv,i)/dph(m,i)/1000.0
             else
               smc(m,jd,nv,i) = 0.0
             end if
             smcld(i,m) = smcld(i,m) + smc(m,jd,nv,i)*wtld
             sicld(i,m) = sicld(i,m) + sic(m,jd,nv,i)*wtld
           end do
           do n = 1,nsl
             stcld(i,n) = stcld(i,n) + tnd(n,jd,nv,i)*wtld
           end do
!
           snoweq(i)= snoweq(i) + swq(jd,nv,i)*wtld
           rsno(i)  = rsno(i) + rsn(jd,nv,i)*wtld
           tsfld(i) = tsfld(i) + tsf(jd,nv,i)*wtld
           tpkld(i) = tpkld(i) + tpk(jd,nv,i)*wtld
           sfwld(i) = sfwld(i) + sfw(jd,nv,i)*wtld
           pkwld(i) = pkwld(i) + pkw(jd,nv,i)*wtld
         end do
       end do
       snoweq(i) = snoweq(i) * 1000.0    ! m -> mm
       canopy(i) = canopy(i) * 1000.0    ! m -> mm
     end if
   end do
!
!   vic -> gmp/rmp: prepare variables for return to parent model
!
! 7. output (o):
!
!   return the following output fields to parent model
!
!  snowmt,snowev,gflux,
!  drain,evap,hflx,ep,
!  snowfl, runoff
!
   do i = 1,im
     gflux(i) = 0.0              ! initialize nonland grid cell
     snowfl(i) = 0.0
     runoff(i) = 0.0
     drain(i) = 0.0
     gheat(i) = 0.0
!
     if(slimsk(i) .eq. 1.) then
!
!   tskin  land surface temperature (k)
!
       tskin(i) = ts(i)
!
!   evap   actual latent heat flux (w m-2: positive, if upward from sfc)
!
       evap(i) = - lh(i)
!
!   hflx       sensible heat flux (w m-2: positive, if upward from sfc)
!
       hflx(i) = - sh(i)
!
!   gflux      soil heat flux (w m-2: negative if downward from surface)
!
       gflux(i) = gh(i)
!
!   gheat      ground heat storage (w m-s: negative if downward from surface)
!
       gheat(i) = gdh(i)
!
!   ep         potential evaporation (w m-2)
!
! -- temporary, use osulsm2 to compute ep ------------------------------
!
!  compute potential evaporation for land
!
       t2 = ts(i) * ts(i)
       t4 = t2 * t2
!
!  rcap = fnet - sigma t**4 + gflx - rho cp ch v (t1-theta1)
!
       rcap = (1.0-albd(i))*sols(i) - sigma * t4 + gflux(i)                  &
                -rch(i) * (ts(i) - theta1(i))
!
!  rsmall = 4 sigma t**3 / rch + 1
!
       rsmall = 4.*sigma*ts(i)*t2/rch(i) + 1.
!
!  delta = l / cp * dqs/dt
!
       delta = elocp * eps * hvap * qs1(i) / (rd * t2)
!
       ep(i) = elocp*rsmall*rch(i)*(qs1(i)-q0(i))+rcap*delta
       ep(i) = ep(i) / (rsmall + delta)
!
! -- end of computing ep -----------------------------------------------
!
!   runoff     surface runoff (mm s-1), not infiltrating the surface
!   drain    subsurface runoff (mm s-1), drainage out bottom
!
       snowfl(i) = sfall(i)                 ! m/time step
       runoff(i) = ovflw(i)/dtime
       drain(i)  = bsflw(i)/dtime
     endif
   enddo
!
!   dm       ratio of actual/potential evap (dimensionless)
!
!  compute qsurf and dm
!
   do i = 1,im
     if(slimsk(i) .eq. 1.) then
       qsurf(i) = q1(i) + evap(i) / (elocp * rch(i))
       dm(i) = 1.
     endif
   enddo
!
   return
   end subroutine lsm_vic1_main
!
!-------------------------------------------------------------------------------
#include <define.h>
   subroutine phys_lsm_vic1(  msl,     nsl,   msub,     md,                    &
                            dtime,  month,       i,    lat,                    &
                              prc,   pgcm,    wind,   tgcm,                    &
                             qgcm,   zgcm,   flwds,   sols,                    &
                             binf,     ds,     dsm,     ws,                    &
                             expt,    kst,     dph,    cef,                    &
                              bub,    qrt,     bkd,    sld,                    &
                              wcr,    wpw,    silz,   snwz,                    &
                              smr,    smx,   dphnd,  smxnd,                    &
                            expnd,  bubnd,   alpnd,  betnd,                    &
                            gamnd,   nveg,     nvt,     wt,                    &
                               rt,   flai,     cwt,    csn,                    &
                              smc,    sic,     tnd,    swq,                    &
                              rsn,    tsf,     tpk,    sfw,                    &
                              pkw,    fmu,   sfall,  lstsn,                    &
                               sh,     lh,      gh,     ts,                    &
                             albd,  ovflw,   bsflw, snowmt,                    &
#ifdef SMP_RA2SFC
                           snowev,    gdh,   ez0,  netsw)
#else
                           snowev,    gdh,   ez0)
#endif
!-------------------------------------------------------------------------------
!
! source file:       phys_lsm_vic1.f
!
! purpose:           simulate the interaction between the land surface
!                    and the atmosphere over a cell
! program history log:
!     03-11-01  ji chen
!     modified from 4.0.3 uw vic c version
!
!-------------------------------------------------------------------------------
   use vic_veglib
!-------------------------------------------------------------------------------
#include <vartyp.h>
!
! input variables 
!
! model basic parameter
!
   integer  ::  msl              ! number of soil layers
   integer  ::  nsl              ! number of soil nodes for soil temperature
   integer  ::  msub             ! maximum number of subgrid points
   integer  ::  md               ! mdist (1no/2with) precipiation distribution

   real     ::  dtime            ! vic time step (second)
   integer  ::  month            ! current month
   integer  ::  i, lat           ! long and lati indexes
!
! atmosphere data
!
   real     ::  prc              ! precipitation rate (m h2o/time step)
   real     ::  pgcm             ! atm bottom level pressure (pa)
   real     ::  wind             ! atm bottom level wind (m/s)
   real     ::  tgcm             ! atm bottom level temperature (kelvin)
   real     ::  qgcm             ! atm bottom level specific humidity (kg/kg)
   real     ::  zgcm             ! atm bottom level height above surface (m)
   real     ::  flwds            ! downward longwave radi onto surface (w/m2)
   real     ::  sols             ! solar rad onto srf (w/m2)
#ifdef SMP_RA2SFC
   real     ::  netsw            ! net solar rad onto srf (w/m2)
#endif
!
! soil parameters
!
   real     ::  binf             ! variable infiltration curve parameter (n/a)
   real     ::  ds               ! fract of dsm nonlinear baseflow begins (fract)
   real     ::  dsm              ! maximum velocity of baseflow (mm/day)
   real     ::  ws               ! fract maxi sm nonlinear baseflow occurs (fract)
   real     ::  cef              ! exponent used in infiltration curve (n/a)
   real     ::  expt(msl)        ! para the vari of ksat with sm (n/a)
   real     ::  kst(msl)         ! saturated hydrologic conductivity (mm/day)
   real     ::  dph(msl)         ! thickness of soil layer (m)
   real     ::  bub(msl)         ! bubbling pressure of soil layer (cm)
   real     ::  qrt(msl)         ! quartz content of soil layer (fraction)
   real     ::  bkd(msl)         ! bulk density of soil layer (kg/m3)
   real     ::  sld(msl)         ! soil density of soil layer (kg/m3)
   real     ::  wcr(msl)         ! fract sm content at the critical point (mm)
   real     ::  wpw(msl)         ! fractional sm content wilting point (mm)
   real     ::  silz             ! surface roughness of bare soil (m)
   real     ::  snwz             ! surface roughness of snowpack (m)
   real     ::  smr(msl)         ! soil moisture residual moisture (mm)
   real     ::  smx(msl)         ! maximum soil moisture (mm)
   real     ::  dphnd(nsl)       ! thickness of soil node (m)
   real     ::  smxnd(nsl)       ! maximum soil moisture at soil node (m3/m3)
   real     ::  expnd(nsl)       ! para the vari of ksat at soil node (n/a)
   real     ::  bubnd(nsl)       ! bubbling pressure at soil node (cm)
   real     ::  alpnd(nsl)       ! para alpha at soil node
   real     ::  betnd(nsl)       ! para beta at soil node
   real     ::  gamnd(nsl)       ! para gamma at soil node
   real     ::  vz0            ! roughness length for vegetation
!
! vegetation parameters
!
   integer  ::  nveg             ! number of vegetation type in a grid cell
   integer  ::  nvt(msub)        ! vegetation type number
   real     ::  wt(msub)         ! fraction of grid cell covered by veg
   real     ::  rt(msl,msub)     ! fraction of root in the soil layer
   real     ::  flai(msub)       ! leaf area index
!
! modified variables 
!
   real     ::  cwt(md,msub)     ! canopy water (m h2o)
   real     ::  csn(md,msub)     ! canopy snow  (m h2o)
   real     ::  smc(msl,md,msub) ! soil moisture content (mm)
   real     ::  sic(msl,md,msub) ! soil ice content (mm)
   real     ::  tnd(nsl,md,msub) ! temperature at soil nodes (k)
   real     ::  swq(md,msub)     ! snow water equivalent (m h2o)
   real     ::  rsn(md,msub)     ! snow density (kg/m3)
   real     ::  tsf(md,msub)     ! snow surface temperature (k)
   real     ::  tpk(md,msub)     ! snow pack temperature (k)
   real     ::  sfw(md,msub)     ! surface snow water equivalent (m h2o)
   real     ::  pkw(md,msub)     ! snow pack snow water equivalent (m h2o)
   real     ::  fmu              ! precipitation fraction (fraction)

   integer  ::  lstsn            ! number of time step since last snowevant
                                 ! > 0: snow covered ground, < 0: no snow
!
! output variables 
!
   real     ::  sh               ! sensible heat flux (w/m**2) [+ to ground]
   real     ::  lh               ! latent heat flux (w/m**2) [+ to ground]
   real     ::  gh               ! ground heat flux (w/m**2) [+ to ground]
   real     ::  ts               ! surface radiative temperature (kelvin)
   real     ::  albd             ! albedo (fraction)
   real     ::  ovflw            ! overland flow (mm/time step)
   real     ::  bsflw            ! base flow (mm/time step)
   real     ::  snowmt           ! snow melt (m/s)
   real     ::  snowev           ! evaporation over snow surface (w/m2)
   real     ::  sfall            ! snow fall (m h2o/time step)
   real     ::  gdh              ! ground heat storage (w/m**2) [+ to ground]
   real     ::  ez0            ! roughness length for vegetation
!
! local variables 
!
   integer  ::  vtype         ! vegetation type
   integer  ::  ndist         ! number of precipitation districts
!
   real     ::  rhoair           ! air density (kg/m3)
   real     ::  vpair            ! vapor pressure (pa)
   real     ::  vpd              ! vapor pressure deficit (pa)
   real     ::  svp              ! function name to compute vapor pressure
!
   real     ::  rfall            ! rain fall (m h2o/time step)
   real     ::  ppt              ! effective precipitation (m/time step)
!
   real     ::  u(3), ra(3)      ! unified wind and resistance for vic
!
   real     ::  fflai            ! fractional subgrid lai
   real     ::  fwt              ! fractional subgrid weight
   real     ::  fswq             ! fractional subgrid swq
   real     ::  frsn             ! fractional subgrid rsn
   real     ::  ftsf             ! fractional subgrid tsf
   real     ::  ftpk             ! fractional subgrid tpk
   real     ::  fsfw             ! fractional subgrid sfw
   real     ::  fpkw             ! fractional subgrid pkw
   real     ::  fcwt             ! fractional subgrid cwt
   real     ::  fcsn             ! fractional subgrid csn
   real     ::  fsmc(msl)        ! fractional subgrid smc
   real     ::  fsic(msl)        ! fractional subgrid sic
   real     ::  frt(msl)         ! fractional subgrid rt
   real     ::  ftnd(nsl)        ! fractional subgrid tnd
   real     ::  falbd            ! fractional albedo over subgrid
   real     ::  flh              ! fractional subgrid latent heat
   real     ::  fsh              ! fractional subgrid sensible heat
   real     ::  fgh              ! fractional subgrid ground heat
!
   real     ::  fsnowm           ! fractional subgrid snow melt (m/time step)
   real     ::  snowme           ! fractional subgrid snow melt (w/m2)
   real     ::  fsnowe           ! fractional subgrid snow evaporation
   real     ::  frtf             ! fractional subgrid rtf
   real     ::  fovflw           ! fractional subgrid overland flow
   real     ::  fbsflw           ! fractional subgrid base flow
   real     ::  fgdh             ! fractional subgrid ground heat storage
   real     ::  fz0            ! fractional subgrid roughness length
!
   real     ::  sncvfr           ! snow cover fraction (fraction)
   real     ::  tmpt             ! temporary soil temperature (k)
   real     ::  sfc_max_unfwat   ! function name
!
!  t0c           ::   ice/water mix temperature (k)
!
   real, parameter  ::  rd = 2.8705e+2, t0c=273.15

   integer  ::  m, n, nv, idist, k     ! loop indexes
!-------------------------------------------------------------------------------
!
! following is land surface scheme from vic 
!
#ifdef DBG
   print *,'--- enter phys_lsm_vic1 --- lat=',lat,' lon i=',i
#endif
   
   sh = 0.
   lh = 0.
   gh = 0.
   ts = 0.
   albd = 0.
   ovflw = 0.
   bsflw = 0.
   snowmt = 0.
   snowev = 0.
   gdh = 0.
   ez0 = 0.

!
! compute air density, vapor pressure
!

   rhoair = pgcm / (rd * tgcm)       ! air density (kg/m3)
   vpair = pgcm*qgcm/0.622
   vpd = svp(tgcm) - vpair
   if(vpd.lt.0.0) vpd = 0.0
!
! depended on the air temperature,
! the precipitation is splitted into rainfall and/or snowfall.
!

#ifdef DBG
   print *,'vic check in phys_lsm_vic1 before vic_rainorsnow'
#endif
!
   rfall = 0.0
   sfall = 0.0
!
   call vic_rainorsnow(msl,   nsl,  msub,    md,                               &
                  nveg,   prc,  tgcm,   fmu,                                   &
                 lstsn, rfall, sfall,   cwt,                                   &
                   csn,   swq,   rsn,   tsf,                                   &
                   tpk,   sfw,   pkw,   smc,                                   &
                   sic,   tnd)

!
! start land surface process simulation
!

   if(fmu.ne.1.) then
     ndist = md
   else 
     ndist = 1
   end if
!
   do nv = 1,nveg     ! nveg, number of different land covers in a cell
!
     vtype = nvt(nv)           ! vegetation type
#ifdef SMP_RA2SFC
!
! replace LAI from RA2 with LAI from vic_veglib 
!
     fflai = veg_lai(month,vtype)          ! leaf area index
#else
     fflai = flai(nv)          ! leaf area index
#endif

#ifdef DBG
     print *,'phys_lsm_vic1 vty=', vtype, fflai
#endif

     do m = 1,msl
       frt(m) = rt(m,nv)      ! root content
     end do
!
     do idist = 1,ndist
       if(idist .eq. 1) then           ! in the wet district
         fwt = fmu*wt(nv)
       else                            ! in the dry district
         fwt = (1.0 - fmu)*wt(nv)
       end if
!
       sncvfr = 0.0                    ! initialize snowcoverfraction
       fswq = swq(idist,nv)
       frsn = rsn(idist,nv)
       ftsf = tsf(idist,nv)
       ftpk = tpk(idist,nv)
       fsfw = sfw(idist,nv)
       fpkw = pkw(idist,nv)
!
       fcwt = cwt(idist,nv)
       fcsn = csn(idist,nv)
!
       do m = 1,msl
         fsmc(m) = smc(m,idist,nv)
         fsic(m) = sic(m,idist,nv)
       end do
!
       do n = 1,nsl
         ftnd(n) = tnd(n,idist,nv)
       end do

!
! compute the vic_exch resistance
!
#ifdef DBG
       print *,'vic check in phys_lsm_vic1 before vic_exch'
#endif
       call vic_exch(vtype,  month, silz, snwz, vz0,                           &
                        zgcm,   wind,    u,   ra)
! 
! update the soil ice content according to the new soil condition
! 

       do m = 1,msl
         if(m.eq.1) then
           tmpt = (ftnd(1)+ftnd(2))/2.0
         else
           tmpt = ftnd(m+1)
         endif
         if(tmpt.lt.t0c) then
           fsic(m) = fsmc(m)-                                                  &
           sfc_max_unfwat(tmpt,smx(m),bub(m),expt(m))
           if(fsic(m).lt.0.0) fsic(m) = 0.0
         else
           fsic(m) = 0.0
         endif
       enddo

! 
! solve land surface fluxes
!
#ifdef DBG
       print *,'vic check in phys_lsm_vic1 before surface_fluxes'
       print *,'vtype ',vtype,' month=',month,fsmc,fsic
#endif

       call vic_surface_flux(msl,    nsl,  dtime,  month,                      &
                     rfall,  sfall,   pgcm,   tgcm,                            &
                     flwds,   sols,  vpair,    vpd,                            &
                    rhoair,      u,     ra,   binf,                            &
                       dph,    qrt,    bkd,    sld,                            &
                       wcr,    wpw,    smr,  dphnd,                            &
                     smxnd,  expnd,  bubnd,  alpnd,                            &
                     betnd,  gamnd,  vtype,    frt,                            &
                     fflai,  lstsn,   fswq,   frsn,                            &
                      ftsf,   ftpk,   fsfw,   fpkw,                            &
                      fcwt,   fcsn,   fsmc,   fsic,                            &
                      ftnd,  falbd,    flh,    fsh,                            &
                       fgh,   frtf, fsnowm, fsnowe,                            &
#ifdef SMP_RA2SFC
                      fgdh,  netsw, sncvfr)
#else
                      fgdh, sncvfr)
#endif

#ifdef DBG
       print *,'in phys_lsm_vic1 after vic_surface_flux',fsmc,fsic
#endif

! 
! update the soil ice content after computing new soil temperature
! 

       do m = 1,msl
         if(m.eq.1) then
           tmpt = (ftnd(1)+ftnd(2))/2.0
         else
           tmpt = ftnd(m+1)
         endif
         if(tmpt.lt.t0c) then
           fsic(m) = fsmc(m)-sfc_max_unfwat(tmpt,smx(m),bub(m),expt(m))
           if(fsic(m).lt.0.0) fsic(m) = 0.0
         else
           fsic(m) = 0.0
         endif
       enddo
!
! update variables for continous run
!
       cwt(idist,nv) = fcwt
       csn(idist,nv) = fcsn
!
       swq(idist,nv) = fswq
       rsn(idist,nv) = frsn
       tsf(idist,nv) = ftsf
       tpk(idist,nv) = ftpk
       sfw(idist,nv) = fsfw
       pkw(idist,nv) = fpkw
!
       do n = 1,nsl
         tnd(n,idist,nv) = ftnd(n)
       end do
!
       if(md.ne.1) then
         if(fmu.eq.1) then
           cwt(md,nv) = fcwt
           csn(md,nv) = fcsn
!               
           swq(md,nv) = fswq
           rsn(md,nv) = frsn
           tsf(md,nv) = ftsf
           tpk(md,nv) = ftpk
           sfw(md,nv) = fsfw
           pkw(md,nv) = fpkw
!               
           do n = 1,nsl
             tnd(n,md,nv) = ftnd(n)
           end do
         end if
       end if
!
! aggragate output variables
!
       albd = albd + falbd*fwt
!
       lh = lh + flh*fwt
       sh = sh + fsh*fwt
       gh = gh + fgh*fwt
       snowmt = snowmt + fsnowm*fwt/dtime  ! m/time step-> m/s
       snowev = snowev + fsnowe*fwt
       gdh = gdh + fgdh*fwt
! roughness length
       if (vtype .eq. 12) then
         fz0=silz
       else
         fz0=vz0
       endif
       ez0 = ez0 + fwt*(sncvfr*snwz+(1-sncvfr)*fz0)
!
       ts = ts + fwt*(sncvfr*ftsf+(1.-sncvfr)*ftnd(1))
! 
! compute land surface runoff, and soil moisture movement
! 
#ifdef DBG
       print *,'vic check in phys_lsm_vic1 before runoff', fsmc,fsic
       print *,' fsnowm, frtf, sncvfr ',fsnowm,frtf,sncvfr
       print *,' others ',msl,  dtime,    ppt,  binf,                          &
                      ws,     ds,    dsm,   smx,                               &
                     smr,    kst,   expt,                                      &
                  fovflw, fbsflw
#endif
       ppt = fsnowm+(1.0-sncvfr)*frtf
!
       call vic_runoff_solver(msl,  dtime,    ppt,   binf,                     &
                                 ws,     ds,    dsm,    cef,                   &
                                smx,    smr,    kst,   expt,                   &
                                dph,    bub,                                   &
                               fsmc,   fsic, fovflw, fbsflw)

#ifdef DBG
       print *,' after runoff',fsmc,fsic,fovflw, fbsflw
#endif
!
! -- update variables for continous run
!
       do m = 1,msl
         smc(m,idist,nv) = fsmc(m)
         sic(m,idist,nv) = fsic(m)
         if (md.ne.1) then
           if (fmu.eq.1.0) then
             smc(m,md,nv) = fsmc(m)
             sic(m,md,nv) = fsic(m)
           end if
         end if
       end do
!
! -- aggragate output variables
!
       ovflw = ovflw + fovflw*fwt
       bsflw = bsflw + fbsflw*fwt
     end do
   end do
!
   return
   end subroutine phys_lsm_vic1
!
!------------------------------------------------------------------------------
