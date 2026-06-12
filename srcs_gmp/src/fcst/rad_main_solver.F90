#include <define.h>
   subroutine rad_main_solver(ix,ix2,prsi,prsl,                                &
                     lat,xlat,solc,rsin1,rcos1,rcos2,slmsk,                    &
                     qqh2o,tt,o3qo3,tsfc,coszro,cldary,cldtot,cldcnv,          &
#if defined (RRTMGSW) || defined (RRTMGLW) 
                     snoweq,vtype,tti,                                         &
#endif
#ifndef SWRMDC
                     albedo,albvb,albnb,albvd,albnd,paers,                     &
#else
                     albedo,albvb,albnb,albvd,albnd,                           &
                     kprf,idxc,cmix,denn,rhrh,                                 &
#endif
                     itimsw,itimlw,ko3,kalb,iswsrc,                            &
                     ibnd,icwp,qci,icfc,                                       &
                     htlw,slwup,slwdn,tlwup,                                   &
                     htsw,sswup,sswdn,tswup,tswdn,                             &
#ifdef NIM_DIAG
                     fnet,                                                     &
#endif
#ifdef CLR
                     tlwup0,tswup0,sswdn0,sswup0,slwdn0,                       &
#endif
#ifdef EXPLICIT_CLOUDINESS
                     qci_cps,                                                  &
                     taucld,cldwp,cldip,                                       &
#endif
                     sswfvb,sswfnb,sswfvd,sswfnd)
!-------------------------------------------------------------------------------
!
! subroutine: rad_main_solver
!
! program history log:
!   1994-09-01  yu-tai hou             chou's sw radiation
!   1998-09-01  yu-tai hou             aerosol updates
!   1999-09-01  song-you hong          add ice effects and prognostic clouds
!
! -- input preparation
!   taulw(i,l)   - cloud optical depth for each layer, k=1-->top layer
!   cld0 (i,l)   - unweighted cloud fraction, k=1-->top layer
!   pl2  (i,lp1) - interfacing pressure (mb), k=1-->model top
!
!-paramodel
!  *******************************************************************!
!  *  rewrite the gfdl code "radfs" to call m.d.chous sw radiation
!  *    and to call the gfdl lw radiation         y-t hou dec1994
!  *  updates made to opac aerosol algorithm (1998)...y.h...sep99
!
!  !
!  *  argument list:
!  *    input
!  *      ix,ix2      - first dimensions of arrays
!  *      qs          - surface pressure in cb
!  *      prsi        - model pressure interface level (k=1 at the sfc)
!  *      prsl        - model pressure layer mean value
!  *      lat,xlat    - lat index and radians
!  *      solc        - solar constant in ly/min
!  *      rsin1,rcos1,rcos2
!  *                  - sin and cos lat for ozone interpolations
!  *      slmsk       - land/sea/ice mask (0:sea.1:land,2:ice)
!  *      qqh2o       - sepcific humidity in g/g   (k=1 at the sfc)
!  *      tt          - temperature in k           (k=1 at the sfc)
!  *      o3qo3       - ozone concentration in g/g (k=1 at the sfc)
!  *      tsfc        - surface temperature in k
!  *      coszro      - cosine of zenith angle
!  *      cldary      - packed cloud array         (k=1 at the sfc)
!  *      cldtot      - stratiform cloud
!  *      cldcnv      - convective cloud
!  *      albedo      - surface albedo from climotology
!  *      albvb,albvd - vis band albedoes for beam and diff radiation
!  *      albnb,albnd - nir band albedoes for beam and diff radiation
#ifndef SWRMDC
!  *      paers       - aerosol profiles (in fractions)
#else
!  *      kprf        - tropospheric aerosol profile type index
!  *      idxc,cmix   - aerosol component index and mixing ratio
!  *      denn        - aerosol number densities of 1st and 2nd layers
!  *      rhrh        - relative humidity in fraction
#endif
!  *    control flags
!  *      itimsw,itimlw
!  *                  - sw, lw radiation calls duration in hour
!  *      ko3         - ozone data, =1 input data; =0 gfdl climotology
!  *      kalb        - sfc alb, =0 climotology, comp over oceans
!  *                    =1 input four components from calling progm
!  *      iswsrc      - flags for selection of sw absorbers
!  *                    1:aerosols, 2:o2, 3:co2, 4:h2o, 5:o3
!  *                    =0:without; =1 with
!  *    output
!  *      htlw        - lw heating rates in k/sec
!  *      slwup       - sfc upward lw flux in w/m**2
!  *      slwdn       - sfc downward lw flux in w/m**2
!  *      tlwup       - toa upward lw flux in w/m**2
!  *      htsw        - sw heating rates in k/sec
!  *      sswup       - sfc upward sw flux in w/m**2
!  *      sswdn       - sfc downward sw flux in w/m**2
!  *      tswup       - toa upward sw flux in w/m**2
!  *      tswdn       - toa downward sw in w/m**2
!  *      tlwup0      - clear sky toa upward lw flux
!         tswup0      - clear sky toa upward sw flux
!  *      sswdn0      - clear sky sfc downward sw flux
!  *      sswup0      - clear sky sfc upward sw flux
!  *      slwdn0      - clear sky sfc downward lw flux
!  *      sswfvb      - vis beam down sw flux at sfc in w/m**2
!  *      sswfvd      - vis diff down sw flux at sfc in w/m**2
!  *      sswfnb      - nir beam down sw flux at sfc in w/m**2
!  *      sswfnd      - nir diff down sw flux at sfc in w/m**2
!  *******************************************************************!
!
!-------------------------------------------------------------------------------
#if defined (RRTMGSW) || defined (RRTMGLW)
   use constant, only           : pi_,g_,cp_,rd_,qmin_
#else
   use constant, only           : pi_,g_
#endif
   use blockdata_rad_gfdl, only : albd,za,trn,dza
#ifndef SWRMDC
#ifndef ICECLOUD
   use rdparm
#else
   use rdparm8
#endif
#else
   use rdparm99
   use aerparm
#endif
#if defined(RRTMGSW) || defined(RRTMGLW) 
#ifndef DFS
   use comfver, only            : sl
#else
   use dfsvar, only             : sl
#endif
#ifdef RRTMGSW
   use module_ra_rrtmg_sw
#endif
#ifdef RRTMGLW
   use module_ra_rrtmg_lw
#endif
#ifdef RMP
   use rscomloc, only           : rcenlat
#endif
   use paramodel, only          : ngases_, ncloud_, ntotal_, kcloud_ !jhc
#endif
   use comfcst, only            : dduo3n, ddo3n2, ddo3n3, ddo3n4,              &
                                  degrad,hsigma,daysec,rco2
#ifdef DBG
#ifdef DFS
   use dfsvar, only             : iope
#else
   use comio, only              : iope
#endif
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! --- input from from calling program
!
   integer              ::  ix,ix2,lat
   real                 ::  solc,rsin1,rcos1,rcos2
   integer              ::  itimsw,itimlw,ko3,kalb
   integer              ::  ibnd,icwp,icfc
   real                 ::  prsi(ix2,lp1), prsl(ix2,l),slmsk(ix2)
#if defined (RRTMGSW) || defined (RRTMGLW)
   real                 ::  snoweq(ix2),vtype(ix2),tti(ix2,lp1)
   real                 ::  qqh2o(ix2,lh), tt(ix2,l), tsfc(ix2),   o3qo3(ix2,l)
#else
   real                 ::  qqh2o(ix2,l), tt(ix2,l), tsfc(ix2),   o3qo3(ix2,l)
#endif	
   real                 ::  cldary(ix2,l),cldtot(ix2,l),cldcnv(ix2,l)
   real                 ::  coszro(ix2),albedo(ix2),xlat(ix2)
   real                 ::  albvb(ix2),   albnb(ix2),   albvd(ix2),   albnd(ix2)
#ifndef SWRMDC
   real                 ::  paers(ix2,nae-1)
   integer              ::  iswsrc(nsrc)
   real                 ::  qci(ix2,l)
#else
   real                 ::  rhrh(ix2,l)
   real                 ::  qci(ix2,l)
   real                 ::  denn(ndn,ix2)
   real                 ::  cmix(nxc,ix2)
   integer              ::  iswsrc(nsrc)
   integer              ::  idxc(nxc,ix2)
   integer              ::  kprf(ix2)
#endif
#if defined (RRTMGSW) || defined (RRTMGLW)
   real                 ::  prsi_pa(imbx,lp1)
   real                 ::  prsl_pa(imbx,l)
   real                 ::  t8w(imbx,lp1)
   real                 ::  t3d(imbx,l)
   real                 ::  snow(imbx)
   real                 ::  vgtype(imbx)
   real                 ::  slmsk2(imbx)
#endif
!
! --- output to calling program
!
   real                 ::  htlw (ix2,l), slwup(ix2), slwdn(ix2), tlwup(ix2)
   real                 ::  htsw (ix2,l), sswup(ix2), sswdn(ix2)
   real                 ::  tswup(ix2), tswdn(ix2)
#ifdef NIM_DIAG
   real                 ::  fnet(ix2,l+1)
#endif
!
! --- four components of downward sw flux
!
   real                 ::  sswfvb(ix2), sswfnb(ix2),sswfvd(ix2),sswfnd(ix2)
#ifdef CLR
   real                 ::  tlwup0(ix2), tswup0(ix2),sswdn0(ix2),sswup0(ix2)
   real                 ::  slwdn0(ix2)
#endif
#ifdef EXPLICIT_CLOUDINESS
   real                 ::  qci_cps(ix2,l)
   real                 ::  taucld(ix2,l),cldwp(ix2,l),cldip(ix2,l),tmp2
#endif
!
! --- internal arrays
!
! --- local arrays for radiative quantities
!
   real                 ::  hlw(imbx,l), tlwuc(imax), slwnc(imax), tswdc(imax)
   real                 ::  hsw(imbx,l), tswuc(imax), sswuc(imax), sswdc(imax)
   real                 ::  tlwu0(imax), tswu0(imax), sswd0(imax), sswu0(imax)
   real                 ::  slwn0(imax)
!
! --- local arrays for clouds
!
   integer              ::  nclds(imax)
   real                 ::  cfacsw(imbx,lp1), cfaclw(imbx,lp1,lp1)
   integer              ::  ktop (imbx,lp1), kbtm  (imbx,lp1)
   real                 ::  cldlw(imbx,lp1), cldsw (imbx,l),   taucl (imbx,l)
#if defined (RRTMGSW) || defined (RRTMGLW)
   real                 ::  cldfra(imbx,l), gcldsw(imbx,l), gtaucl(imbx,l) !jhc (rrtmg)
#endif
   real                 ::  cwp  (imbx,l),   cip   (imbx,l),   fice  (imbx,l)
   real                 ::  rew  (imbx,l),   rei   (imbx,l)
   real                 ::  taulw(imbx,l),   cld0(imbx,l) , clwp(imbx,l)
!
! --- local arrays for others
!
   real                 ::  prssi(imbx,lp1), prssl(imbx,lp1), temp (imbx,lp1)
   real                 ::  sfcalb(imax)
   real                 ::  rh2o (imbx,l),   qo3  (imbx,l)  , coszen(imax)
#if defined (RRTMGSW) || defined (RRTMGLW)
   real                 ::  grh2o (imbx,l), tsk(imbx,1)
   real                 ::  gqq(imbx,lh)
   real                 ::  w1, w2
#endif	
#ifndef SWRMDC
   real                 ::  paera(imbx,nae)
#else
   real                 ::  rhrh1(imbx,l)
   real                 ::  denn1(ndn,imax)
   real                 ::  cmix1(nxc,imax)
   integer              ::  idxc1(nxc,imax)
   integer              ::  kprf1(imax)
#endif
   real                 ::  alvb (imax),     alnb (imax),     alvd (imax)
   real                 ::  alnd (imax)
   real                 ::  gdfvb(imax),     gdfnb(imax),     gdfvd(imax)
   real                 ::  gdfnd(imax)
   integer              ::  jjrow(imax)
   real                 ::  do3v (imax),     do3vp(imax),    tthan(imax)
   real                 ::  pl2  (imbx,lp1)
!
! --- surface albedo (used for the old scheme)
!
   integer              ::  nrchnk,irchnk,ibeg,iend,ipts,i,ir,j,k,iq,kk
   integer              ::  jx,ipts1
   real                 ::  ssolar,tmp1,rhodz,th2,zen,dzen
   real                 ::  alb1,alb2,albd0,alvd1
#if defined (RRTMGSW) || defined (RRTMGLW)
   real                 ::  center_lat
   logical              ::  f_qv,f_qc,f_qr,f_qi,f_qs,f_qg,f_qndrop
   logical              ::  allowed_to_read
   integer              ::  ncloud, kcloud, ktrace
   integer              ::  ids,ide, jds,jde, kds,kde,                         &
                            ims,ime, jms,jme, kms,kme,                         &
                            its,ite, jts,jte, kts,kte
#endif
#ifdef DBG
   integer              ::  latd
#endif
!-------------------------------------------------------------------------------
!
! initialize local variables
!
   hlw=0.     ; tlwuc=0.   ; slwnc=0.   ; tswdc=0.
   hsw=0.     ; tswuc=0.   ; sswuc=0.   ; sswdc=0.
   tlwu0=0.   ; tswu0=0.   ; sswd0=0.   ; sswu0=0.
   slwn0=0.   ; nclds=0
   cfacsw=0.  ; cfaclw=0.
   ktop=0     ; kbtm=0
   cldlw=0.   ; cldsw=0.   ; taucl=0.
   cwp=0.     ; cip=0.     ; fice=0.
   rew=0.     ; rei=0.
   taulw=0.   ; cld0=0.    ; clwp=0.
   prssi=0.   ; prssl=0.   ; temp=0.
   sfcalb=0.
   rh2o=0.    ; qo3=0.     ; coszen=0.
#ifndef SWRMDC
   paera=0.
#else
   rhrh1=0.   ; denn1=0.   ; cmix1=0.   ; idxc1=0  ;  kprf1=0
#endif
   alvb=0.    ; alnb=0.    ; alvd=0.    ; alnd=0.
   gdfvb=0.   ; gdfnb=0.   ; gdfvd=0.   ; gdfnd=0.
   jjrow=0    ; do3v=0.    ; do3vp=0.   ; tthan=0. ;  pl2=0.
#if defined (RRTMGSW) || defined (RRTMGLW)
!
   gcldsw=0.  ; gtaucl=0. !jhc (rrtmg)
! define the dimension constants
!  
   ids = 1
   ide = imbx
   jds = 1
   jde = 1
   kds = 1
   kde = lp1
!  
   ims = 1
   ime = imbx
   jms = 1
   jme = 1
   kms = 1
   kme = lp1
!
   its = 1
   ite = imbx
   jts = 1
   jte = 1
   kts = 1
   kte = l
!
   ncloud = ntotal_
#ifdef WSM1
   kcloud = 1    ! temporary (not calculated)
#else
!   kcloud = (ngases_+1)*kte+1
   kcloud = kcloud_ ! jhc (rrtmg)
#endif
   ktrace = kte+1
#endif
!
!===> ... begin here
!       solc,the solar constant is scaled to a more current value.
!       i.e. if solc=2.0 ly/min then ssolar=1.96 ly/min. then
!       convert unit to w/m**2
!
   ssolar = solc * 0.98e0 * 6.97667e2
!
   nrchnk=(ix-1)/imax+1
   do irchnk = 1,nrchnk
     ibeg = (irchnk - 1) * imax + 1
     iend = ibeg + imax - 1
     if (iend .gt. ix) iend = ix
     ipts = iend - ibeg + 1
!
!===> ... assign temp,pressures, cosin of zenith angle
!         note: the nmc variables are in mks, gfdl lw variables
!         are in cgs units, so prssl is in dyns/cm**2. but
!         sw prssi is in mb.
!
     do i = 1,ipts
       ir = i + ibeg - 1
       prssl(i,lp1) = 1.0e4 * prsi(ir,1)
       prssi(i,lp1) = 10.0e0 * prsi(ir,1)
       pl2  (i,lp1) = 10.0e0 * prsi(ir,1)
       temp (i,lp1) = max(100.,tsfc(ir))
#if defined (RRTMGSW) || defined (RRTMGLW)
       prsi_pa(i,lp1) = 1.0e3 * prsi(ir,lp1)
       prsl_pa(i,l) = 1.0e3 * prsl(ir,l)
       t8w(i,lp1) = tti(ir,lp1)
       t3d(i,l)   = tt(ir,l)
       snow(i)    = snoweq(ir)
       vgtype(i)  = vtype(ir)
       slmsk2(i)  = slmsk(ir)
#endif
       coszen(i)    = coszro(ir)
#ifndef SWRMDC
!
!===> ... currently the convective type of aerosol is not in use
!         to give zeroes here.
!
       paera(i,6)   = 0.0e0
#else
       kprf1(i)     = kprf(ir)
#endif
     enddo
#ifndef SWRMDC
     do k = 1,5
       do i = 1,ipts
         ir = i + ibeg - 1
         paera(i,k) = paers(ir,k)
       enddo
     enddo
#else
     do i = 1,ipts
       ir = i + ibeg - 1
       do j = 1,nxc
         idxc1(j,i) = idxc(j,ir)
         cmix1(j,i) = cmix(j,ir)
       end do
       do j = 1,ndn
         denn1(j,i) = denn(j,ir)
       end do
     enddo
#endif
!
!===> ... all radiation variables have k=1 at the top of the atmosphere.
!         while the forcast model variables have k=1 at the surface.
!
     do k = 1,l
       do i = 1,ipts
         ir = i + ibeg - 1
         temp (i,k) = max(100.,tt(ir,lp1-k))
         prssl(i,k) = 1.0e4 * prsl(ir,lp1-k)
         prssi(i,k) = 10.0e0 * prsi(ir,lp2-k)
         pl2  (i,k) = 10.0e0 * prsi(ir,lp2-k)
         rh2o (i,k) = max(3.0e-6, qqh2o(ir,lp1-k))
#if defined (RRTMGSW) || defined (RRTMGLW)
         tsk (i,1) = max(100.,tsfc(ir))
         prsi_pa(i,k) = prsi(ir,k)*1.e3
         t8w(i,k) = tti(ir,k)
         grh2o (i,k) = max(3.0e-6, qqh2o(ir,k))
#endif
#ifdef SWRMDC
         rhrh1(i,k) = rhrh(ir,lp1-k)
#endif
       enddo
     enddo
#if defined (RRTMGSW) || defined (RRTMGLW)
!    do k = 1,lp1
     do k = 1,l
       do i = 1,ipts
         ir = i + ibeg - 1
         prsl_pa(i,k) = prsl(ir,k)*1.e3
         t3d(i,k) = tt(ir,k)
       enddo
     enddo
#endif
!
     if(icwp.eq.1) then
       tmp1 = 1.0e5 / g_
       do k = 1,l
         do i = 1,ipts
           ir = i + ibeg - 1
           rhodz = (prssi(i,k+1)-prssi(i,k))*tmp1
#ifndef EXPLICIT_CLOUDINESS
           clwp (i,k) = max(qci(ir,lp1-k),0.) * rhodz          ! g/m**2
#else
#ifdef CPS_QSUM
           tmp2 = qci(ir,lp1-k) + qci_cps(ir,lp1-k)
           clwp (i,k) = max(tmp2,0.) * rhodz                   ! g/m**2
#else
           clwp (i,k) = max(qci(ir,lp1-k),0.) * rhodz          ! g/m**2
#endif
#endif
         enddo
       enddo
     endif
!***************************************************!
!         ozone section
!***************************************************!
     if (ko3 .gt. 0) then
       do k = 1,l
         do i = 1,ipts
!
! ken+moorthi: set up minimum positive value for o3
!
           qo3(i,k) = max(o3qo3(i+ibeg-1,lp1-k), 1.0e-10)
         enddo
       enddo
     else
       do i = 1,ipts
         ir = i + ibeg - 1
         th2= 0.2e0 * xlat(ir) * degrad
         jjrow(i) = 19.001e0 - th2
         tthan(i) = (19-jjrow(i)) - th2
       enddo
!
!===> ... compute climatological zonal mean ozone,
!         seasonal and spatial interpolation done below.
!
       do k = 1,l
         do i = 1,ipts
           do3v(i)  = dduo3n(jjrow(i),k) + rsin1*ddo3n2(jjrow(i),k)            &
                + rcos1*ddo3n3(jjrow(i),k)                                     &
                + rcos2*ddo3n4(jjrow(i),k)
           do3vp(i) = dduo3n(jjrow(i)+1,k) + rsin1*ddo3n2(jjrow(i)+1,k)        &
                + rcos1*ddo3n3(jjrow(i)+1,k)                                   &
                + rcos2*ddo3n4(jjrow(i)+1,k)
!
!===> ... now latitudinal interpolation, and
!         convert o3 into mass mixing ratio(original data mpy by 1.e4)
!
           qo3(i,k) = 1.0e-4 * (do3v(i)+tthan(i)*(do3vp(i)-do3v(i)))
         enddo
       enddo
     endif
!*****************************************************!
!         cloud optical properties section
!*****************************************************!
#ifdef EXPLICIT_CLOUDINESS
     do k = 1,l
       do i = 1,ipts
         ir = i + ibeg - 1
         taucld(ir,k) = 0.0
         cldwp(ir,k)  = 0.0
         cldip(ir,k)  = 0.0
       enddo
     enddo
#endif
!
#ifndef SWRMDC
#ifndef ICECLOUD
     call rad_cloud_comp(ix2,prssi,temp,cldary,ibeg,ipts,xlat                  &
               ,ktop,kbtm,nclds,cldlw,taucl,cfacsw,cldsw                       &
               ,taulw,cld0)
#else
     call rad_cld_property_ice(ix2,prsi,prsl,prssi                             &
               ,temp,cldtot,cldcnv,ibeg,ipts,lat                               &
               ,icwp,clwp,slmsk,xlat,ktop,kbtm,nclds                           &
               ,cldlw,taucl,cfacsw,cldsw                                       &
               ,cwp,cip,rew,rei,fice                                           &
               ,taulw,cld0)
#endif
#else
     call rad_cld_property_nasa(ix2,prssi,temp,cldtot,cldcnv,ibeg,ipts,lat     &
               ,icwp,clwp,slmsk,xlat,ktop,kbtm,nclds                           &
               ,cldlw,taucl,cfacsw,cldsw                                       &
               ,cwp,cip,rew,rei,fice                                           &
               ,taulw,cld0)
#endif
#ifdef EXPLICIT_CLOUDINESS
     do k = 1,l
       do i = 1,ipts
         ir = i + ibeg - 1
         kk = l - k + 1                   ! re-ordering (sfc->top)
         taucld(ir,kk) = max(taulw(i,k),0.)       ! tau_0
         if (icwp.eq.1) then
           cldwp(ir,kk) = max(cwp(i,k),0.)       ! g/m**2
           cldip(ir,kk) = max(cip(i,k),0.)       ! g/m**2
         endif
       enddo
     enddo
#endif
!*****************************************************!
!         surface albedo section
!*****************************************************!
!
     if (kalb .eq. 0) then
       iq = int(20.0e0 * 0.537e0 + 1.0e0)
       do i = 1,ipts
         ir = i + ibeg -1
         sfcalb(i) = albedo(ir)
!
!===> ... the following code gets albedo from payne,1972 tables if
!      1) open sea point (slmsk=0);   2) kalb=0
!
         if(coszen(i).gt.0.0 .and. slmsk(ir).eq.0.0) then
           zen = degrad*acos(max(coszen(i),0.0e0))
           if(zen .ge. 74.e0) jx = int(0.5e0*(90.e0-zen)+1.e0)
           if(zen .lt. 74.e0 .and. zen .ge. 50.e0)                         &
                     jx = int(0.25e0*(74.e0-zen) + 9.0e0)
           if(zen.lt.50.e0) jx = int(.1e0*(50.e0-zen)+15.0e0)
           dzen = -(zen-za(jx))/dza(jx)
           alb1 = albd(iq,jx)+dzen*(albd(iq,jx+1)-albd(iq,jx))
           alb2 = albd(iq+1,jx)+dzen*(albd(iq+1,jx+1)-albd(iq+1,jx))
           sfcalb(i) = alb1+20.e0*(alb2-alb1)*(0.537e0-trn(iq))
         endif
       enddo
       do i = 1,ipts
         alvd(i) = sfcalb(i)
         alnd(i) = sfcalb(i)
         alvb(i) = sfcalb(i)
         alnb(i) = sfcalb(i)
       enddo
!
!===> ... visible and near ir direct beam albedo,if not ocean nor snow
!         function of cosine solar zenith angle..
!
       do i = 1,ipts
         if (slmsk(i+ibeg-1).gt.0.0e0 .and. sfcalb(i).le.0.5e0) then
           albd0 = -18.0e0 * (0.5e0 - acos(coszen(i))/pi_)
           albd0 = exp (albd0)
           alvd1 = (alvd(i) - 0.054313e0) / 0.945687e0
           alvb(i) = alvd1 + (1.0e0 - alvd1) * albd0
           alnb(i) = alvb(i)
         end if
       enddo
     else
       do i = 1,ipts
         ir = i + ibeg - 1
         alvd(i) = albvd(ir)
         alnd(i) = albnd(ir)
         alvb(i) = albvb(ir)
         alnb(i) = albnb(ir)
       enddo
     end if
#if defined (RRTMGSW) || defined (RRTMGLW)
     if(ntotal_.gt.1) then
       do k = ktrace,lh
         do i = 1,ipts
           ir = i + ibeg - 1
           gqq(i,k)=max(qqh2o(ir,k),qmin_)
         enddo
       enddo
     endif
!jhc rrtmg/ : convert cldsw(top-to-bottom) to gcldsw (bottom-to-top)
     do k=1,l
       do i=1,ipts
         gcldsw(i,k)=cldsw(i,lp1-k)
         gtaucl(i,k)=taucl(i,lp1-k)
       enddo
     enddo
!/jhc rrtmg
#endif
!**********************************************!
!         final check before radiation calls
!**********************************************!
     if (ipts .lt. imax) then
       ipts1 = ipts + 1
       do i = ipts1,imax
         coszen(i) = coszen(ipts)
         alvd(i) = alvd(ipts)
         alnd(i) = alnd(ipts)
         alvb(i) = alvb(ipts)
         alnb(i) = alnb(ipts)
       enddo
       do k = 1,lp1
         do i = ipts1,imax
           prssl(i,k) = prssl(ipts,k)
           prssi(i,k) = prssi(ipts,k)
           pl2  (i,k) = pl2  (ipts,k)
           temp (i,k) = temp (ipts,k)
#if defined (RRTMGSW) || defined (RRTMGLW)
           prsi_pa(i,k) = prsi_pa(ipts,k)
           t8w(i,k) = t8w(ipts,k)
#endif
         enddo
       enddo
       do k = 1,l
         do i = ipts1,imax
#if defined (RRTMGSW) || defined (RRTMGLW)
           grh2o(i,k) = grh2o(ipts,k)
           prsl_pa(i,k) = prsl_pa(ipts,k)
           t3d(i,k) = t3d(ipts,k)
#endif
           rh2o(i,k) = rh2o(ipts,k)
           qo3 (i,k) = qo3 (ipts,k)
         enddo
       enddo
#if defined (RRTMGSW) || defined (RRTMGLW)
       do k = ktrace,lh
         do i = ipts1,imax
           gqq(i,k)=gqq(ipts,k)
         enddo
       enddo
!jhc/
       do k=1, l
         do i=ipts1,imax
           gcldsw(i,k)=gcldsw(ipts,k)
           gtaucl(i,k)=gtaucl(ipts,k)
         enddo
       enddo
!/jhc
#endif
     end if
!
     if(itimsw .eq. 0) go to 300
!************************************************!
!         calling chous sw radiation routine
!************************************************!
!
!===> ... set timer for sw rad call
!
!     t00  = second()
!     ovhd = second() - t00
!     tbef = second()
!
#if defined (RRTMGSW) || defined (RRTMGLW)
#ifdef RMP
     center_lat = rcenlat
#else
     center_lat = 0.
#endif
!       
     f_qv = .true.
#ifdef WSM1
     f_qc = .false.
     f_qr = .false.
     f_qi = .false.
     f_qs = .false.
     f_qg = .false.
     f_qndrop = .false.
#endif
#ifdef WSM3
     f_qc = .true.
     f_qr = .true.
     f_qi = .false.
     f_qs = .false.
     f_qg = .false.
     f_qndrop = .false.
#endif
#ifdef WSM5
     f_qc = .true.
     f_qr = .true.
     f_qi = .true.
     f_qs = .true.
     f_qg = .false.
     f_qndrop = .false.
#endif
#if defined (WSM6) || defined (WDM6)
     f_qc = .true.
     f_qr = .true.
     f_qi = .true.
     f_qs = .true.
     f_qg = .true.
     f_qndrop = .false.
#endif
!
#endif
!
#if defined(SWRMDC) || defined (RRTMGSW)
#ifdef SWRMDC
     call rad_sw_nasa_driver(ipts,ssolar,iswsrc,ibnd,prssi,                    &
              temp,rh2o,qo3,rco2,coszen,taucl,                                 &
              cldsw,cfacsw,icfc,icwp,cwp,cip,rew,rei,fice,                     &
              alvb,alvd,alnb,alnd,kprf1,idxc1,cmix1,denn1,rhrh1,               &
              hsw,tswuc,tswdc,sswuc,sswdc,                                     &
              tswu0,sswu0,sswd0,                                               &
              gdfvb,gdfvd,gdfnb,gdfnd)
#endif
#ifdef RRTMGSW
     call rad_sw_rrtmg_driver(                                                 &
              hsw,                                                             &
              tswuc,tswdc,sswuc,sswdc,                                         &
              coszen, ssolar,                                                  &
              alvb,alvd,alnb,alnd,                                             &
              t3d, t8w, tsk,                                                   &
!              prsl_pa,prsi_pa, cldsw, rd_, g_,                                 &
              prsl_pa,prsi_pa, gcldsw, rd_, g_,                                & !jhc (gcldsw)  
              gtaucl, slmsk2, snow,                                            & !jhc (add gtaucl)
              grh2o, gqq(1,kcloud), ncloud,                                    &
              f_qv, f_qc, f_qr, f_qi, f_qs, f_qg,                              &
              f_qndrop,                                                        & !czhao
#ifdef CHEM
              lat,                                                             &
#endif
              gdfvb,gdfvd,gdfnb,gdfnd,                                         &
              ids,ide, jds,jde, kds,kde,                                       &
              ims,ime, jms,jme, kms,kme,                                       &
              its,ite, jts,jte, kts,kte)
#endif
#else
     call rad_sw_gfdl_driver(ipts,ssolar,iswsrc,ibnd,prssi,                    &
              temp,rh2o,qo3,rco2,coszen,taucl,                                 &
              cldsw,cfacsw,icfc,icwp,                                          &
#ifdef ICECLOUD
              cwp,cip,rew,rei,fice,                                            &
#endif
              alvb,alvd,alnb,alnd,paera,                                       &
              hsw,tswuc,tswdc,sswuc,sswdc,                                     &
              tswu0,sswu0,sswd0,                                               &
#ifdef CHEM
              lat,ibeg,                                                        &
#endif
#ifdef NIM_DIAG
              gdfvb,gdfvd,gdfnb,gdfnd,fnet)
#else
              gdfvb,gdfvd,gdfnb,gdfnd)
#endif /* NIM_DIAG end */
#endif
!
!===> ... save toa and sfc fluxes in w/m**2
!
     do i = 1,ipts
       ir = i + ibeg - 1
       tswup (ir) = tswuc(i)
       tswdn (ir) = tswdc(i)
       sswup (ir) = sswuc(i)
       sswdn (ir) = sswdc(i)
#ifdef CLR
       tswup0(ir) = tswu0(i)
       sswdn0(ir) = sswd0(i)
       sswup0(ir) = sswu0(i)
#endif
!
!===> ... diffused downward sfc fluxes (vis,nir)
!
       sswfvd(ir) = gdfvd(i)
       sswfnd(ir) = gdfnd(i)
!
!===> ... direct beam downward sfc fluxes (vis,nir)
!
       sswfvb(ir) = gdfvb(i)
       sswfnb(ir) = gdfnb(i)
     enddo
!
!===> ... convert heating rates to deg/sec
!
     do k = 1,l
       do i = 1,ipts
#ifdef RRTMGSW
         htsw(i+ibeg-1,k) = hsw(i,k)
#else
         htsw(i+ibeg-1,lp1-k) = hsw(i,k) * daysec
#endif
       enddo
     enddo
!
 300 if(itimlw .eq. 0) go to 400
!
#if defined(LWRMDC) || defined (RRTMGLW)
#ifdef LWRMDC
!
!-------------------------------------------------------------
!-------------------------------------------------------------
!     calling uiuc_chou long wave radiation routine
!
#ifdef DBG
#ifdef NIM
     latd=1
#else
     latd=10
#endif
     if(lat.eq.latd)then
       write(6,*)'before rad_lw_nasa_driver  hlw=',hlw
     endif
#endif
     call rad_lw_nasa_driver(ipts,hlw,slwnc,tlwuc,                             &
#ifdef CLR
                 slwn0,tlwu0,                                                  &
#endif
                 pl2,temp,rh2o,qo3,taulw,cld0)
#ifdef DBG
     if(lat.eq.latd)then
       write(6,*)'after rad_lw_nasa_driver  hlw=',hlw
     endif
#endif
#endif
#ifdef RRTMGLW
!	  
#ifdef DBG
     if(lat.eq.latd)then
       write(6,*)'before rad_lw_rrtmg_driver  hlw=',hlw
     endif
#endif
     call rad_lw_rrtmg_driver(                                                 &
              hlw,                                                             &
              slwnc,tlwuc,                                                     &
#ifdef CLR
              slwn0,tlwu0,                                                     &
#endif
              t3d, t8w, tsk,                                                   &
!              prsl_pa,prsi_pa, cldsw, rd_, g_,                                 &
              prsl_pa,prsi_pa, gcldsw, rd_, g_,                                & !jhc (gcldsw)
              gtaucl, slmsk2, snow, vgtype,                                    & !jhc (add gtaucl)
              grh2o, gqq(1,kcloud), ncloud,                                    &
              f_qv, f_qc, f_qr, f_qi, f_qs, f_qg,                              &
              f_qndrop,                                                        & 
              ids,ide, jds,jde, kds,kde,                                       &
              ims,ime, jms,jme, kms,kme,                                       &
              its,ite, jts,jte, kts,kte                                        &
                                                                               )
#ifdef DBG
     if(lat.eq.latd)then
       write(6,*)'after rad_lw_rrtmg_driver  hlw=',hlw
     endif
#endif
!	  
#endif
#else
!************************************************!
!        calling gfdl long wave radiation routine
!************************************************!
!==> ...  get cld factor for lw calculations
!
!  tbef = second()
!
     call rad_lw_gfdl_cloud(ipts,cfaclw,cldlw,nclds,kbtm,ktop)
!
     call rad_lw_gfdl_driver(ipts,hlw,slwnc,tlwuc,                             &
#ifdef CLR
              slwn0,tlwu0,                                                     &
#endif
#ifdef NIM
              prssi,                                                           &
#endif
              prssl,temp,rh2o,qo3,cfaclw,cldlw,nclds,ktop,kbtm)

#endif
!-------------------------------------------------------------
!-------------------------------------------------------------
!
!===> ... save toa and sfc fluxes in w/m**2
!
     do i = 1,ipts
       ir = i + ibeg - 1
#ifdef RRTMGLW
       tlwup (ir) = tlwuc(i)
       slwdn (ir) = slwnc(i)
#ifdef CLR
       tlwup0(ir) = tlwu0(i)
       slwdn0(ir) = slwn0(i)
#endif
#else
       tlwup (ir) = 1.0e-3*tlwuc(i)
       slwup (ir) = 1.0e-3*hsigma*temp(i,lp1)**4
       slwdn (ir) = slwup(ir) - 1.0e-3*slwnc(i)
#ifdef CLR
       tlwup0(ir) = 1.0e-3*tlwu0(i)
       slwdn0(ir) = slwup(ir) - 1.0e-3*slwn0(i)
#endif
#endif
     enddo
!
!===> ... convert heating rates to deg/sec
!
     do k = 1,l
       do i = 1,ipts
#if defined(RRTMGLW)
         htlw(i+ibeg-1,k) = hlw(i,k)
#else
         htlw(i+ibeg-1,lp1-k) = hlw(i,k) * daysec
#endif
       enddo
     enddo
!
 400 continue
!
   enddo
!
   return
   end subroutine rad_main_solver
!
!------------------------------------------------------------------------------
