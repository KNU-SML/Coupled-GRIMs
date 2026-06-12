#include <define.h>
   subroutine rad_cloud_prepare(lons2,                                         &
#ifndef NONHYD
#ifdef SMP
                     dgrs,                                                     &
#else
                     dgrs,plamgr,pphigr,                                       &
#endif
                     ugrs,vgrs,                                                &
#endif /* NONHYD */
                     tgrs,qgrs,pgr,                                            &
#ifdef NONHYD 
                     wgrs,                                                     &
#endif /* NONHYD */
                     albedr,slmskr,rlon,rlat,                                  &
                     tsear,shelgr,tgr,cvr,cvtr,cvbr,rhcl,                      &
#ifndef SWRMDC
                     ozonea,albdoa,cldary,cldtot,cldcnv,                       &
#else
                     ozonea,albdoa,cldary,cldtot,cldcnv,rhrh,                  &
#endif
                     cldsa,mtopa,mbota,                                        &
                     rrs2,lat,latco,latrue,istrat,                             &
#ifdef NIM
                     ips,ipe,si1d,sl1d,prsi_nim,prsl_nim,                      &
#endif
                     kalb,jo3,slag,rsin1,rcos1,rcos2,                          &
                     fjd,dlt,jsno,workr,lworkr                                 &
#ifdef DG
                     ,cldt,cldc                                                &
#endif
#ifdef EXPLICIT_CLOUDINESS
                     ,qcicps, qrscps                                           &
                     ,qcis_cps                                                 &
#endif
                     ,icwp,qcis)
!-------------------------------------------------------------------------------
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!     ************************************************************
!     *  added accumulation of clds and convective cloud in dg3  *
!     *                                        k.a.c sept 1994   *
!     *  changed  h,m,l calculation                              *
!     *     in rad_cloudiness (removed facv)                     *
!     *     and added proper total cloud calculation             *
!     *  changed avecld calc in cldiag (used cldary) and         *
!     *     used total cloud calculated in rad_cloudiness        *
!     *                                         k.a.c. nov94     *
!     *  interpolate o3 profile to each gridpoint, ie use        *
!     *     proper surface pressure                              *
!     *                                         k.a.c. dec94     *
!     *  fix pl1 for operations, where dgz is on and dg3 is off, *
!     *     ....note dg is on if either dgz or dg3 is on         *
!     *                                         k.a.c. jan94     *
!     ************************************************************
!
!        updates made to add oceanic stratus and to fix conv cloud..
!                     to rad_driver - ivv(2),ibl are set=1....
!                     to rad_driver - set min q to 1.e-10,rather than 1.e-6
!                                 to anticipate avoiding cld creation
!                                 in extremely dry,cold (winter) regions
!                                 where 1.e-6 could imply hi valu of rh
!                     to rad_cloudiness - multitude of changes
!        updates made to add grid point diagnostics ..k.a.c...sep 91
!                     to rad_driver -
!        updates made to fix sw approx              ..k.a.c...nov 91
!                     to rad_cos_zenith
!        updates made to pass and receive sib data  ..k.a.c...mar 92
!                     to rad_driver -
!        updates made to fix sw rad diagnostics     ..k.a.c...jun 92
!                     proper diurnal weighting
!                     to rad_driver and rad_cos_zenith
!        updates made to calculate clear-sky "on-the-fly" kac aug 92
!                     to rad_driver,radfs,fst,spa,lwr,swr
!                     ...for cloud forcing....
!        updates made for the completely new cloud routine
!                     (rad_cloudiness),use
!                     flag ivva to control vertical velocity adj.
!                     for low cld (=0: without,  =1: with)
!                     use flag invr to control lapse rate inversion
!                     type of cld (=0: without,  =1: with)
!                     to rad_driver and radfs ...y.h.           ...dec92
!        updates made to call cld optical property routine (rad_cloud_comp),
!                     to give cld emissivity, optical depth, layer
!                     reflectance and transmitance
!                     to rad_driver and radfs ...y.h.           ...feb93
!
!       updates made to allow tuned cld-rh data to be used..ctune
!                   to rad_driver and rad_cloudiness     ..k.a.c...dec 92
!                   spatial interpolation of tables ........may93
!                   use only 1 set of tuning tables for all fcst hrs,
!                    the tuning of the 24hr fcst       .....jan94
!                    old code used 6 tables..see ckc94 .....feb94
!                   since tuning done for h,m,l cld, vertically
!                    blend the relations at old hml bdries..jan94
!        updates made to change definition of h,m,l domains..
!          to rad_driver,rad_cloudiness_init,rad_cloudiness, rad_cloud_comp   
!                                                    ..k.a.c...dec92 + aug93
!          wmo aerosols distributions, and b.p.breiglebs surface
!          albedo scheme.            ..............y.h...sep95
!             kalb is the control flag for surface albedo
!             kalb=0 use the old matthews data, =1 use the new scheme.
!             jo3 is the control flag for ozone climatology
!             jo3=0 use the old gfdl data, =1 use the new nasa data
#if defined (SWRMDC) || defined (RRTMGSW)
!        recoded the default part of cloud scheme patterned after
!          ccm3 (kiehl et al. 1998,j.clim; and 1994,jgr)... y.h. aug98
!        updates made to opac aerosol algorithm (1998)..... y.h. sep99
!        updates made to cloud overlapping and moorthi's liquid
!          water cloud fraction scheme ............y.h...oct00
#endif
!
!  ::: structure :::
!
!    [rad_cloud_prepare]
!      |
!      |--- A)   [rad_cloudiness] - [rad_cloudiness-adjust] *
!      |--- B)   [rad_cloudiness_nasa] *
!      |
!      |--- [rad_ozone_gfdl] *
!      |--- [rad_ozone_nasa] *
!      |--- [rad_albedo8snow] *
!
!-------------------------------------------------------------------------------
   use paramodel, only : LATG2S,LONF2S,levh_,ntotal_,nwmass_,ngases_,levs_,levp1_
   use constant, only  : pi_,cp_,akapa_,rv_,rd_,g_,fv_
#ifdef RMP
#ifdef MP
   use paramodel, only : igrd12p_,jgrd12p_
#else
   use paramodel, only : igrd12_
#endif
   use rscomver
#endif
#ifdef DFS
   use dfsvar, only    : si=>sigmafull,sl=>sigma,del=>delsig
#endif
#if !defined (SWRMDC) && !defined (RRTMGSW)
#ifndef ICECLOUD
   use rdparm
#else
   use rdparm8
#endif
#else
   use rdparm99
#endif
#ifndef RMP
   use comfver       
#endif
#ifdef MP
   use commpi
#endif
#ifdef NIM_DBGM
   use comio,     only : iope,ncpu,ndbg1,kdbg1
#endif
#ifndef RMP
#define LONLENS lons2
#else
#ifdef MP
#ifdef RMPVECTORIZE
#define LONLENS igrd12p_*jgrd12p_
#else
#define LONLENS lonlen(mype)*2
#endif
#else
#undef RMPVECTORIZE
#define LONLENS igrd12_
#endif
#endif
!-------------------------------------------------------------------------------
!
   integer,parameter    ::  mcld=3,nseal=2,nbin=100,nlon=2,nlat=4
!
#ifndef NONHYD
   real                 ::  dgrs(LONF2S,levs_)
#ifndef SMP
   real                 ::  plamgr(LONF2S),pphigr(LONF2S)
#endif
   real                 ::  ugrs(LONF2S,levs_)
   real                 ::  vgrs(LONF2S,levs_)
#endif /* NONHYD */
   real                 ::  tgrs(LONF2S,levs_)
   real                 ::  qgrs(LONF2S,levh_)
   real                 ::  pgr(LONF2S)
#ifdef NONHYD
   real                 ::  wgrs(LONF2S,levp1_)
#else
#if !defined(RMP) && !defined(NIM)
   real                 ::  rrs2(LATG2S)
#endif
#endif /* NONHYD */
   real                 ::  workr(LONF2S,lworkr)
#ifdef NIM
   real                 ::  si1d(levs_+1), sl1d(levs_)
   real                 ::  prsi_nim(LONF2S,levs_+1), prsl_nim(LONF2S,levs_)
#endif
   real                 ::  rlon(LONF2S),rlat(LONF2S)
   real                 ::  albedr(LONF2S),slmskr(LONF2S)
   real                 ::  vvel(LONF2S,levs_)
   real                 ::  ssnow(LONF2S)
!
   real                 ::  slk(levs_)
   real                 ::  sik(levp1_)
   real                 ::  prsi(LONF2S,levs_+1)
   real                 ::  prsik(LONF2S,levs_+1)
   real                 ::  prsl(LONF2S,levs_)
   real                 ::  prslk(LONF2S,levs_)
   real                 ::  phii(LONF2S,levs_+1)
   real                 ::  phil(LONF2S,levs_)
   real                 ::  zl(LONF2S,levs_)
   real                 ::  zi(LONF2S,levs_+1)
   real                 ::  delprsi(LONF2S,levs_)
#ifdef NONHYD
   real                 ::  delz(LONF2S,levs_)
#endif /* NONHYD */
!
!   cldary contains multi layers of cloud
!
#if !defined (SWRMDC) && !defined (RRTMGSW)
   real                 ::  cldary(LONF2S,levs_)
#else
   real                 ::  cldary(LONF2S,levs_)
   real                 ::  rhrh (LONF2S,levs_)
#endif
   real                 ::  cldcnv(LONF2S,levs_)
   real                 ::  cldtot(LONF2S,levs_)
#ifdef DG
   integer              ::  icnv(LONF2S)
   real                 ::  cldc(LONF2S,levs_)
   real                 ::  cldt(LONF2S,levs_)
#endif
   real                 ::  qcis(LONF2S,levs_)
   real                 ::  ozonea(LONF2S,levs_)
   real                 ::  albdoa(LONF2S)
   real                 ::  cldsa(LONF2S,4)
   integer              ::  mtopa(LONF2S,3)
   integer              ::  mbota(LONF2S,3)
   real                 ::  wrkema(LONF2S)
   real                 ::  tsear(LONF2S),shelgr(LONF2S)
   real                 ::  cvr(LONF2S),cvtr(LONF2S)
   real                 ::  cvbr(LONF2S),tgr(LONF2S)
#ifdef EXPLICIT_CLOUDINESS
   real                 ::  qcicps(LONF2S,levs_), qrscps(LONF2S,levs_)
   real                 ::  qcis_cps(LONF2S,levs_)
#endif
#ifdef HYBRID
   real                 ::  sihyb(LONF2S,levp1_),cihyb(LONF2S,levp1_)
   real                 ::  slhyb(LONF2S,levs_),clhyb(LONF2S,levs_)
   real                 ::  delhyb(LONF2S,levs_),slkhyb(LONF2S,levs_)
#endif
#ifdef NIM
   integer              ::  ips,ipe
#endif
   integer,save         ::  invr,ivva,ko3cli
   real,save            ::  rhmax,xlim
!
!  array added for rh-cl calculation
!  indices for lon,lat,cld type(l,m,h), land/sea respectively
!  nlon=1-2, for eastern and western hemispheres
!  nlat=1-4, for 60n-30n,30n-equ,equ-30s,30s-60s
!  land/sea=1-2 for land(and seaice),sea
!
   real                 ::  rhcl (nbin,nlon,nlat,mcld,nseal)
   real                 ::  rhcla(nbin,nlon,mcld)
   real                 ::  rhcld(LONF2S,nbin,mcld)
   real                 ::  temp1, temp2
   real,save            ::  xlabdy(3),xlobdy(3)
!
!  xlabdy = lat bndry between tuning regions,+/- xlim for transition
!  xlobdy = lon bndry between tuning regions
!
   data xlabdy / 30.e0 , 0.e0 , -30.e0 /
   data xlobdy / 0.e0 , 180.e0 , 360.e0 /
   data xlim / 5.e0 /
!
!  initial rh crit. set 1 for ocean, set 2 for land.
!  invr=0 no lapse rate inversion type cld, =1 wiht it
!  ivva=0 no vertical velocity adj. for low cld, =1 with adj.
!
   data rhmax/1.00e0/, invr/1/, ivva/1/
!
! initialize local variables
!
   vvel=0. ; ssnow=0. ; slk=0.  ; sik=0.
   prsi=0. ; prsik=0. ; prsl=0. ; prslk=0.
   phii=0. ; phil=0.
   zl=0.   ; zi=0.    ; delprsi=0.
#ifdef NONHYD
   delz=0.
#endif /* NONHYD */
#ifdef DG
   icnv=0
#endif
   wrkema=0.
#ifdef HYBRID
   sihyb=0. ; cihyb=0. ; slhyb=0. ; clhyb=0. ; delhyb=0. ; slkhyb=0.
#endif
   rhcla=0. ; rhcld=0.
!
!   ko3cli=0,1 for gfdl,nasa climo
!     only for ko3=1(meaning o3 is input to radfs)..see rad_main_driver
!
   ko3cli=jo3
#ifdef RMP
   do j = 1,LONLENS
     pgr(j) = exp(pgr(j))
   enddo
#endif
#ifdef DG
   do j = 1,LONLENS
     icnv(j) = 0.e0
   enddo
#endif
!
#ifdef HYBRID
!
! ak5 : Pa -> cb in dyn_hybrid_setup.F90
!
   do i = 1,LONLENS
     work4 = (pgr(i)/100.)**(akapa_)
     prsik(i,1)     = work4
     do k = 1,levp1_
       sihyb(i,k)   = ak5(levp1_-k+1)/pgr(i)+bk5(levp1_-k+1)
       cihyb(i,k)   = 1.-sihyb(i,k)
     enddo
     do k = 1,levs_
       slhyb(i,k)   = 0.5*(sihyb(i,k)+sihyb(i,k+1))
       slkhyb(i,k)  = slhyb(i,k)**(akapa_)
       delhyb(i,k)  = sihyb(i,k)-sihyb(i,k+1)
       clhyb(i,k)   = 1.-slhyb(i,k)
!
       prsi(i,k)    = sihyb(i,k)*pgr(i)
       prsik(i,k+1) = work4 * (sihyb(i,k+1)**(akapa_))
       prsl(i,k)    = slhyb(i,k)*pgr(i)
       prslk(i,k)   = work4 * slkhyb(i,k)
     enddo
       prsi(i,levs_+1) = pgr(i) * sihyb(i,levs_+1)
   enddo
#else
#ifdef NIM
   do k = 1,levs_
     sl(k)=sl1d(k)
     slk(k)=sl(k)**(akapa_)
   enddo
   do k = 1,levp1_
     si(k)=si1d(k)
     sik(k)=si(k)**(akapa_)
   enddo
!
   do i = 1,LONLENS
     work4=(pgr(i)/100.)**(akapa_)
     do k = 1,levs_
       prsl(i,k)   = prsl_nim(i,k)
       prsi(i,k)   = prsi_nim(i,k)
     enddo
     prsi(i,levs_+1) = prsi_nim(i,levs_+1)
   enddo
#else /* not NIM */
   do k = 1,levs_
     slk(k)=sl(k)**(akapa_)
   enddo
   do k = 1,levp1_
     sik(k)=si(k)**(akapa_)
   enddo
!
   do i = 1,LONLENS
     work4=(pgr(i)/100.)**(akapa_)
     prsik(i,1)=work4
     do k = 1,levs_
       prsl(i,k)   = pgr(i)*sl(k)
       prsi(i,k)   = pgr(i)*si(k)
       prslk(i,k)  = work4   *slk(k)
       prsik(i,k+1)= work4   *sik(k+1)
     enddo
     prsi(i,levs_+1) = pgr(i)*si(levs_+1)
   enddo
#endif   /* NIM */
#endif   /* Hybrid */
!
   call dyn_get_pressure(lons2,LONF2S,levs_,akapa_,cp_,fv_,tgrs,qgrs,pgr,si,   &
                prsi,prsik,prsl,prslk,phii,phil,delprsi)
!
#ifndef NIM
   do k = 1,levs_
     do i = 1,lons2
       zl(i,k)=phil(i,k)/g_
       zi(i,k)=phii(i,k)/g_
     enddo
   enddo
   do i = 1,lons2
     zi(i,levs_+1)=phii(i,levs_+1)/g_
   enddo
#endif
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
#ifdef MP
   if (mype.eq.master) then
#endif
     call print_maxmin_seven(prsi,lons2,LONF2S,levs_+1,1,levs_+1,              &
                            'prsi in rad_cloud_prepare')
     call print_maxmin_seven(delprsi,lons2,LONF2S,levs_,1,levs_,               &
                            'delprsi in rad_cloud_prepare')
     call print_maxmin_seven(prsl,lons2,LONF2S,levs_,1,levs_,                  &
                            'prsl in rad_cloud_prepare')
     call print_maxmin_seven(prsik,lons2,LONF2S,levs_+1,1,levs_+1,             &
                            'prsik in rad_cloud_prepare')
     call print_maxmin_seven(prslk,lons2,LONF2S,levs_,1,levs_,                 &
                            'prslk in rad_cloud_prepare')
     call print_maxmin_seven(phii,lons2,LONF2S,levs_+1,1,levs_+1,              &
                            'phii in rad_cloud_prepare')
     call print_maxmin_seven(phil,lons2,LONF2S,levs_,1,levs_,                  &
                            'phil in rad_cloud_prepare')
#ifdef MP
   endif
#endif
#ifdef OPENMP
!$omp end master
#endif
#endif /* DBG end */
#ifdef NIM_DBGM
   print*, 'rad_cloud_prepare prsi_nim',ips,ndbg1,(prsi_nim(ndbg1,k),k=1,5)
   print*, 'rad_cloud_prepare prsl_nim',ips,ndbg1,(prsl_nim(ndbg1,k),k=1,5)
#endif
!
!   generalized the computation and suitable for all grid ---
!                                             by h.-m.juang
   do i = 1,LONLENS
     isla = 1
     if (slmskr(i).lt.1.e0) isla = 2
     xlatpt = rlat(i) * 180.e0 / pi_
!
!  get rh-cld relation for this lat
!
     kla = 4
     do k = 1,3
       if (xlatpt.gt.xlabdy(k)) then
         kla = k
!
         exit
!
       end if
     enddo
!
     klap=0
     do k = 1,3
       xlnn = xlabdy(k)+xlim
       xlss = xlabdy(k)-xlim
       if (xlatpt.lt.xlnn.and.xlatpt.gt.xlss) then
         kla  = k
         klap = k+1
!
         exit
!
       endif
     enddo
!
     if( klap .eq. 0 ) then
       do kc = 1,mcld
         do lo = 1,nlon
           do nbi = 1,nbin
             rhcla(nbi,lo,kc) = rhcl(nbi,lo,kla,kc,isla)
           enddo
         enddo
       enddo
     else
!
!  linear transition between latitudinal regions...
!
       temp1=(xlatpt-xlss)/(xlnn-xlss)
       do kc = 1,mcld
         do lo = 1,nlon
           do nbi = 1,nbin
             temp2=(rhcl(nbi,lo,kla,kc,isla)-rhcl(nbi,lo,klap,kc,isla))
             rhcla(nbi,lo,kc) = temp1 * temp2 + rhcl(nbi,lo,klap,kc,isla)
           enddo
         enddo
       enddo
     endif
!
!  get rh-cld relation for this lon
!
     xlonpt = rlon(i) * 180.e0 / pi_
     lo=1
     if( xlonpt .gt. 180.0 ) lo=2
     ikn = 0
     do k = 1,3
       diflo = abs(xlonpt-xlobdy(k))
       if (diflo.lt.xlim) then
         ikn = k
         ilft = lo
         irgt = ilft + 1
         if (irgt.gt.nlon) irgt = 1
         xlft = xlobdy(ikn) - xlim
         xrgt = xlobdy(ikn) + xlim
!
         exit
!
       endif
     enddo
     if( ikn.eq.0 ) then
       do k = 1,mcld
         do nbi=1,nbin
           rhcld(i,nbi,k) = rhcla(nbi,lo,k)
         enddo
       enddo
     else
       do k = 1,mcld
         do nbi = 1,nbin
           rhcld(i,nbi,k) =                                                    &
           (rhcla(nbi,ilft,k)-rhcla(nbi,irgt,k))                               &
                  * (xlonpt-xrgt)/(xlft-xrgt)+rhcla(nbi,irgt,k)
         enddo
       enddo
     endif
!
   enddo
!
!   get mean zenith angle for this dtswav-both nh and sh
!   start radfs set-up for both hemispheres
!   get vertical motion (cb/sec) in vvel
!
#ifdef RMPVECTORIZE
#define OMEGAS omegasv
#else
#define OMEGAS dyn_get_omega
#endif
!
#ifdef NONHYD
   do k = 1,levs_
     do i = 1,LONLENS
       wmean=-0.5*(wgrs(i,k)+wgrs(i,k+1))
       vvel(i,k)=wmean*prsl(i,k)*g_/(rd_*tgrs(i,k))
     enddo
   enddo
#else /* else of NONHYD */
#ifndef SMP
#ifndef RMPVECTORIZE
   call dyn_get_omega(LONLENS,LONF2S,levs_,                                    &
               pphigr(1),plamgr(1),workr,ugrs(1,1),vgrs(1,1),                  &
#ifndef RMP
               dgrs(1,1),delprsi,rrs2(lat),vvel,prsi(1,1),prsl)
#else
               dgrs(1,1),delprsi,rrs2     ,vvel,prsi(1,1),prsl)
#endif
#else
   call OMEGAS(LONLENS,LONF2S,levs_,                                           &
               pphigr(1),plamgr(1),workr,ugrs(1,1),vgrs(1,1),                  &
#ifndef RMP
               dgrs(1,1),del,rrs2(lat),vvel,pgr(1),sl)
#else
               dgrs(1,1),del,rrs2     ,vvel,pgr(1),sl)
#endif
#endif
#endif /* endif SMP */
#endif /* endif NONHYD */
!
!  get model diagnosed clds
!
#if !defined (SWRMDC) && !defined (RRTMGSW)
   call rad_cloudiness(lons2,LONF2S,levs_,nbin,mcld,                           &
               qgrs(1,1),tgrs(1,1),vvel,                                       &
               cvr(1),cvtr(1),cvbr(1),                                         &
               prsi,prsl,slmskr(1),                                            &
               cldsa(1,1),mtopa(1,1),mbota(1,1),                               &
               cldary(1,1),ivva,invr,rhmax,                                    &
#ifdef DBG
               rlat(1),rhcld,istrat,ntotal_,lat)
#else
               rlat(1),rhcld,istrat,ntotal_)
#endif
!
!dg3     unpack cldamt and conv cldamt from cldary..(stratus is +2)
!dg3      radiation sees stratiform or convective...not merged.....
!dg3      so cldtot reflects this below
!dg3  note:cldary 2-4 digits to left of decimal=cv cloud
!dg3       cldary 1 digit to left of decimal+fractional part=strat cld
!dg3        get stratiform cld into cldtot, conv cld in cldcnv
!dg3        ..anvil ci is stratiform, so place the 1 lyr in cnv also
!dg3        ..icnv = cv cloud top layer
!
   do k = 1,levs_
     do i = 1,LONLENS
#ifdef INTERACTIVE_STRATUS
       cldtot(i,k) = amod (cldary(i,k),2.e0)
       cldcnv(i,k) = 0.
#else
       cldtot(i,k) = amod (cldary(i,k),2.e0)
       cldcnv(i,k) = float(int(cldary(i,k))/10)*1.e-3
#endif
     enddo
   enddo
#else
   call rad_cloudiness_nasa(lons2,LONF2S,levs_,nbin,mcld,                      &
                   qgrs(1,1),tgrs(1,1),vvel,                                   &
                   cvr(1),cvtr(1),cvbr(1),                                     &
                   prsi,prsl,slmskr(1),                                        &
                   cldsa(1,1),mtopa(1,1),mbota(1,1),                           &
                   cldtot(1,1),cldcnv(1,1),ivva,invr,icwp,qcis,rhrh,           &
#ifdef EXPLICIT_CLOUDINESS
                   qcicps(1,1),qrscps(1,1),                                    &
                   qcis_cps,                                                   &
#endif
#ifdef DBG
                   rlat(1),rhcld,istrat,ntotal_,lat)
#else
                   rlat(1),rhcld,istrat,ntotal_)
#endif
#endif
!
!dg3   anvil ci is stratiform, so place the 1 lyr in cnv also
!dg3   icnv = cv cloud top layer
!
#ifdef DG
   do k = 1,levs_
     do i = 1,LONLENS
       cldt(i,k) = cldtot(i,k)
       cldc(i,k) = cldcnv(i,k)
       if(cldc(i,k).gt.0.e0) then
         cldt(i,k) = cldc(i,k)
         icnv(i) = k
       end if
     end do
   end do
!
!dg3   anvil ci placed into layer above conv cld top
!
   do i = 1,LONLENS
     k = icnv(i)
     if (k.gt.1.and.k.lt.levs_) then
       cldc(i,k+1) = cldt(i,k+1)
     end if
   enddo
#endif
!
!  compute mean cloud diagnostics + h,m,l,total cloud
!     call cvdiag(lons2,avecv(1,lat  ),avecv(1,latco),
!    1                         cvr(1),cvtr(1),cvbr(1))
!     call cldiag(lons2,avecld(1,lat  ),cldl(1,lat  ),
!    1                  avecld(1,latco),cldl(1,latco),
!    2                         cldsa(1,1),cldary(1,1))
!
   if(thour.eq.0.or.ngases_.le.0) then
     if (ko3cli.eq.0)                                                          &
       call rad_ozone_gfdl(LONLENS,LONF2S,levs_,                               &
                     prsl,ozonea(1,1),                                         &
                     rlat(1),rsin1,rcos1,rcos2)
     if (ko3cli.eq.1)                                                          &
       call rad_ozone_nasa(LONLENS,LONF2S,levs_,                               &
                   prsl,ozonea(1,1),                                           &
                   rlat(1),lat)
   endif
!
   do i = 1,LONLENS
     if(slmskr(i).eq.2.0e0) then
       tsear(i) = min(tsear(i),271.2e0)
     else if(slmskr(i).eq.0.0e0) then
       tsear(i) = max(tsear(i),271.21e0)
     else if(slmskr(i).eq.1.0e0 .and. shelgr(i).gt.0.0e0)  then
       tsear(i) = min(tsear(i),273.16e0)
     endif
   enddo
   if (kalb .ge. 1) then
     do i = 1,LONLENS
       albdoa(i)=0.
     enddo
     return
   endif
!
!   the following determines surface albedo (albdoa),where snow exists.
!
   do i = 1,LONLENS
     ssnow(i) = shelgr(i) * 0.1e0
   enddo
#ifndef RMP
   call rad_albedo8snow(LONLENS,latrue,jsno,                                   &
#else
   call rad_albedo8snow(LONLENS,lat,jsno,                                      &
#endif
               albdoa(1),rlat(1),albedr(1),slmskr(1),                          &
               ssnow, tgr(1), tgrs(1,1) )
!
   return
   end subroutine rad_cloud_prepare
!
!-------------------------------------------------------------------------------
   subroutine rad_cloudiness(imx2,imx22,kmx,nbin,mcld,                         &
              q,t,vvel,cv,cvt,cvb,prsi,prsl,slmsk,                             &
              cld,mtop,mbot,cldary,ivva,invr,rhmax,                            &
#ifdef DBG
              xlatrd,rhcld,istrat,ntotal,lat)
#else
              xlatrd,rhcld,istrat,ntotal)
#endif
!-------------------------------------------------------------------------------
!
!   from yh.rad.mdl93(cldnew28).......
!     later updated from yh.rad.mdl94(cldmul28)...22jan94
!     later updated from yh.rad.mdl94(cldml28a)... 1feb94
!     later updated from yh.rad.mdl94(cldml28b)... 5feb94
!               subr rad_cloud_comp replaced
!               added vertical interp of cld-rh relations(istrat gt 1)
!     later updated from yh.rad.mdl94(cldml28e)... 11mar94
!               subr rad_cloud_comp replaced,gcl adjusted
!     later updated from phys_mps_wsm6................... 24mar94
!              subr rad_cloud_comp , low enhanced to old value..0.14..
!              subr gclnew , llyr calculation adj to old valu(kl-1)
!                             llyrl was ok.. ive removed it and
!                             replaced it by its equivalent, klowb
!     later updated from phys_mps_wsm6................... 30mar94
!               subr rad_cloud_comp , low and middle (not cv) enhanced=0.10
!---------------------------------------------------------------------
!     nov., 1992 - y.h., k.a.c., and a.k.
!        cloud parameterization patterned after slingo and slingos
!        work (jgr, 1991).
!     stratiform clouds are allowed in any layer except the surface
!        and upper stratosphere.  the relative humidity criterion may
!        vary in different model layers.
!
!     output cloud amounts are in cldary(i,k), k=1 is the lowest
!        model layer, stratiform (str) and convective (cnv) types of
!        cloud are compressed into one word: camt = str + 1.0e4*cnv
!        low marine stratus amts are flaged by adding 2.
!
!   for istrat = 0, there is rh-cld relation for each layer..
!                      crit rh computed within..
!   for istrat = 1, rh-cld relation from tables created using
!                     mitchell-hahn tuning technique (a.f. rtneph obs)
!                  ...stratus computed similar to old opnl rad_cloudiness.....
!                      except no cloud below layer=klowb..approx 955mb
!
!     convective clouds are from model convective scheme and are
!        no longer broken into .75,.25,.25..rather cc itself is used..
!        convective still takes precedence over stratiform in radfs
!         but here cv+st merge exits in cldary...(in radiation use of
!        cc gives improvement to tropical middle cld (as did st+cv))
!
!     clouds are also divided into 3 atmospheric domains (l,m,h) for
!        diagnostic purposes.  they are computed from random overlap
!        assumption for separated cloud layers and maximum overlap
!        for adjacent cloud layers.  a total cloud fraction is also
!        computed.
!
!     h,m,l domain pressure tops 'ptop1(k)' vary linearly from
!        'ptopc(k,1)' at 45deg to 'ptopc(k,2)' at the pole
!
!     input variables:
!        ps (cb)       - surface pressure
!        q  (kg/kg)    - specific humidity
!        t  (deg k)    - absolute temperature
!        vvel(cb/sec)  - vertical velocity
!        cv,cvt,cvb    - conv cld fraction, top, bottom layer from
!        si,sl         - mdl sigma interface and layer mean
!        slmsk         - sea/land mask array(sea:0.,land:1.,snow:2.)
!        ivva          - flag to control vertical velocity adj.
!                        =1: with, =0: without
!        invr          - flag to control lapse rate inversion cld
!                        =1: with, =0: without
!        rhmax         - upper limit of relative humidity to
!                        form overcast cloud (cld fractn = 1.)
!
!  modify to as an array (h.-m. h. juang)
!       xlatrd        - current latitude in radians (1st data pt)
!                         for models with diff lat at each pt, need to
!                         use the lat of all points....careful.....
!        rhcld         - cloud-rh relations from mitchell+hahn,
!                        using a.f. rtneph analyses
!        istrat        - 0 or 1:for default or 'rhcld' tables
!                        in the stratiform cloud calculation
!
!    output variables:
!       cldary         - vertical column array of cloud fraction
!                        profile
!       cld            - cld fraction in 3 types of domains (l,m,h)
!                          and total in 4th layer
!       mtop,mbot      - top, bottom layers of clouds (l,m,h)
!
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_,icloud_,nwmass_
   use constant, only  : pi=>pi_,rd=>rd_,rv=>rv_,akapa_
   use comcd1
!-------------------------------------------------------------------------------
#ifndef ICE
   integer,parameter    ::  nxx=7501
   real                 ::  tbpvs0(nxx)
   common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
#endif
   real,parameter       ::  eps=rd/rv, epsm1=rd/rv-1.0
   real                 ::  cv(imx22),  cvt(imx22),  cvb(imx22)
   real                 ::  slmsk(imx22), prsi(imx22,kmx+1), prsl (imx22,kmx)
   real                 ::  t(imx22,kmx), vvel(imx22,kmx), q(imx22,kmx*ntotal)
   real                 ::  cld(imx22,4)
   real                 ::  cldary(imx22,kmx),xlatrd(imx22)
   real                 ::  spnt(imx22),cr1se1(imx22),cr1sl1(imx22),cr1(imx22)
   real                 ::  xcrh1(imx22),cr2(imx22),xcrh2(imx22)
   integer              ::  mtop(imx22,3),      mbot(imx22,3)
!
!   rh-cld relationships for each point
!
   dimension rhcld(imx22,nbin,mcld)
!
!   ptopc(k,l): top presure of each cld domain (k=1-4 are sfc,l,m,h;
!       l=1,2 are low-lat (<45 degree) and pole regions)
!
!  workspace
!
#ifdef INTERACTIVE_STRATUS
   real                   :: clwt,clwmin,clwm, tem1
   logical                :: inversn(imx22) 
   integer, dimension(imx22) :: kinver
#endif
   logical              ::  bitx(ILOTS),                                       &
           bity(ILOTS),                                                        &
           bitz(ILOTS),                                                        &
           bit1, bit2,                                                         &
           bitm(ILOTS)
   real                 ::  rhrh (ILOTS,levs_)
   real                 ::  prsly(ILOTS,levs_)
   real                 ::  dthdp(ILOTS,levs_)
   real                 ::  theta(ILOTS,levs_)
   real                 ::  cl1  (ILOTS)
   real                 ::  cl2  (ILOTS)
   real                 ::  omeg (ILOTS)
   real                 ::  ptop1(ILOTS,4)
   integer              ::  kcut (ILOTS),kbase(ILOTS)
   integer              ::  kbt1 (ILOTS),kth1 (ILOTS)
   integer              ::  kbt2 (ILOTS),kth2 (ILOTS)
   integer              ::  kcvb (ILOTS),kcvt (ILOTS)
   integer              ::  ksave(ILOTS)
!-------------------------------------------------------------------------------
!
! initialize local variables
!
   rhrh =0.;  prsly=0.;  dthdp=0.;  theta=0.;  cl1  =0.;  cl2  =0.
   omeg =0.;  ptop1=0.;  kcut =0.;  kbase=0.;  kbt1 =0.;  kth1 =0.
   kbt2 =0.;  kth2 =0.;  kcvb =0.;  kcvt =0.;  ksave=0.
!
!   begin here 
!
   kdim=kmx
   kdimp=kmx+1
   levm1=kmx-1
   levm2=kmx-2
!
!  find top pressure for each cloud domain
!
   do k = 1,4
     do i = 1,imx2
       fac = max (0.0e0, 4.0e0*abs(xlatrd(i))/pi-1.0e0)
       ptop1(i,k) = ptopc(k,1) + (ptopc(k,2)-ptopc(k,1)) * fac
     enddo
   enddo
!
!  low cloud top sigma level, computed for each lat cause
!       domain definition changes with latitude...
   klow=kdim
   do k = kdim,1,-1
     do i = 1,imx2
       if (prsi(i,k)/prsi(i,1) .lt. ptop1(i,2)*1.0e-3) klow = min(klow,k)
     enddo
   enddo
!
!  potential temp and layer relative humidity
!
   do k = 1,kdim
     do i = 1,imx2
       cldary(i,k) = 0.0e0
       prsly(i,k) = prsl(i,k) * 10.0e0
       exnr = (prsly(i,k)*0.001e0) ** (-akapa_)
       theta(i,k) = exnr * t(i,k)
#ifdef ICE
       es = fpvs(t(i,k))
#else
       es = fpvs0(t(i,k))
#endif
       qs = eps * es / (prsl(i,k) + epsm1*es)
       rhrh(i,k) = max (0.0e0, min (1.0e0, q(i,k)/qs))
     enddo
   enddo
!
!   potential temp lapse rate
!
   do k = 1,levm1
     do i = 1,imx2
       dthdp(i,k) = (theta(i,k+1) - theta(i,k)) /(prsly(i,k+1) - prsly(i,k))
     enddo
   enddo
! 
!     find the stratosphere cut off layer for high cloud. it
!      is assumed to be above the layer with dthdp less than
!      -0.25 in the high cloud domain (from looking at 1 case).
!
   do i = 1,imx2
     kcut(i) = levm2
   enddo
   do k = klow+1,levm2
     bit1 = .false.
     do i = 1,imx2
       if (kcut(i).eq.levm2 .and. prsly(i,k).le.ptop1(i,3) .and.               &
               dthdp(i,k).lt.-0.25e0) then
         kcut(i) = k
       end if
       bit1    = bit1 .or. kcut(i).eq.levm2
     enddo
     if (.not. bit1) exit
   enddo
!
! ------------------------------------------------------------------
!
   if (istrat.le.0) then
!
! ------------------------------------------------------------------
!     default scheme ....tuned for 28 lyrs by y-t hou.
!     calculate stratiform cloud and put into array 'cldary'
!     the relative humidity criteria are preset for each model
!     sigma level, (1) for ocean points, and (2) for land points.
! ------------------------------------------------------------------
!
     do k = klowb,levm2
       bit1 = .false.
       do i = 1,imx2
         bitx(i) = k.le.kcut(i)
         bit1 = bit1 .or. bitx(i)
       enddo
!
       if (.not. bit1) cycle
!
       do i = 1,imx2
         spnt(i) = max (0.6e0, min (0.85e0, 0.96e0-0.6e0*prsl(i,k)/prsi(i,1)))
         cr1se1(i) = (0.41e0*prsl(i,k)/prsi(i,1) - 0.71e0)**2 + 0.52e0
         cr1sl1(i) = 0.8e0 - 0.167e0*prsl(i,k)/prsi(i,1)
       enddo
       do i = 1,imx2
         if (slmsk(i).eq.1.0e0) then
           cr1(i) = cr1se1(i)
         else
           cr1(i) = cr1sl1(i)
         end if
         xcrh1(i) = spnt(i) * (rhmax - cr1(i))
         cr2(i) = cr1(i) + xcrh1(i)
         xcrh2(i) = rhmax - cr2(i)
         cl1(i) = max (0.e0, (rhrh(i,k)-cr1(i))/xcrh1(i)) ** 3
         if (cl1(i).gt.1.0e0)                                                  &
           cl1(i) = 1.e0 + sqrt((rhrh(i,k)-cr2(i))/xcrh2(i))
       enddo
       do i = 1,imx2
         if (bitx(i)) then
           cldary(i,k) = min (1.0e0, 0.5e0*cl1(i))
         end if
       enddo
     enddo
!
!  special treatment on low clouds
!
     dvvcld = vvcld(1) - vvcld(2)
     rclap = 1.0e0 / (0.8e0 - crhrh)
     do i = 1,imx2
       kbase(i) = 0
     enddo
!
     do k = klowb,klowt
!
       do i = 1,imx2
         omeg(i) = 10.0e0 * vvel(i,k)
         cl1 (i) = 0.0e0
       enddo
!
       if (ivva .le. 0) go to 250
!
!  vertical velocity adjustment on low clouds
!
       bit1 = .false.
       do i = 1,imx2
         bitx(i) = prsly(i,k).ge.ptop1(i,2) .and. cldary(i,k).gt.0.0e0
         bit1 = bit1 .or. bitx(i)
       enddo
!
       if (.not. bit1) go to 250
!
       do i = 1,imx2
         if (bitx(i)) then
           if (omeg(i).ge.vvcld(1)) then
             cldary(i,k) = 0.0e0
             else if(omeg(i).gt.vvcld(2)) then
             cr1(i) = (vvcld(1) - omeg(i)) / dvvcld
!            cldary(i,k) = cldary(i,k) * cr1(i)
             cldary(i,k) = cldary(i,k) * sqrt(cr1(i))
           endif
         endif
       enddo
!
!  t inversion related stratus clouds
!
250    continue
!
       if (invr .lt. 1) cycle
       bit1 = .false.
       do i = 1,imx2
         bitx(i) = prsly(i,k).ge.pstrt .and. slmsk(i).le.0.0                   &
                .and. dthdp(i,k).le.clapse
         bit1 = bit1 .or. bitx(i)
       enddo
!
       if (.not. bit1) cycle
!
       do i = 1,imx2
         if ( kbase(i).eq.0  .and. rhrh(i,k).gt.crhrh .and. bitx(i) )          &
               kbase(i) = k
       enddo
       do i = 1,imx2
         if (kbase(i).gt.0 .and. bitx(i) .and. cldary(i,k+1).le.0.1e-1         &
                  .and. cldary(i,k+2).le.0.1e-1) then
           cr1(i) = min (1.0e0,                                                &
                max (0.0e0,  16.67e0*(clapse-dthdp(i,k)) ))
           if(rhrh(i,kbase(i)).lt.0.8e0) then
             cr1(i) = cr1(i) * (rhrh(i,kbase(i))-crhrh) * rclap
           endif
!
!  for t inversion type cloud, add flag value of 2.0
!
             cldary(i,k) = max (cldary(i,k), cr1(i)) + 2.0e0
         endif
       enddo
     enddo
!
   end if
!
   if (istrat.gt.0) then
#ifdef INTERACTIVE_STRATUS
     do i = 1, imx2
       inversn(i) = .false.
       kinver (i) = 1
     enddo
     do k = kdim-1, 2, -1
       do i = 1, imx2
         if (prsly(i,k) > 600.0 .and. (.not.inversn(i))) then
           tem1 = t(i,k-1) - t(i,k)
           if (tem1 > 0.1 .and. t(i,k) > 278.0) then
             inversn(i) = .true.
             kinver(i)  = k
           endif
         endif
       enddo
     enddo
     clwmin = 0.0
!
     do k = kdim, 1, -1
       do i = 1, imx2
         clwt = 1.0e-6 * (prsly(i,k)*0.001)
         clwt = 2.0e-6 * (prsly(i,k)*0.001)
         qsum = 0.
!
         do ic = icloud_,nwmass_ 
           kk = k + levs_ * (ic-1)
           qsum = q(i,kk) + qsum
         enddo
!
         if (qsum > clwt .or.                                                  &
           (inversn(i) .and. k >= kinver(i)) ) then
           onemrh= max( 1.e-10, 1.0-rhrh(i,k) )
           clwm  = clwmin / max( 0.01, prsly(i,k)*0.001 )
!
#ifdef ICE              
           es = fpvs(t(i,k))
#else
           es = fpvs0(t(i,k))
#endif
           qs = eps * es / (prsl(i,k) + epsm1*es)
           tem1  = min(max(sqrt((onemrh*qs)),0.0001),1.0)
           tem1  = 100.0 / tem1
           value = max( min( tem1*(qsum-clwm), 50.0 ), 0.0 )
           tem2  = sqrt( sqrt(rhrh(i,k)) )
           cldary(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
         endif
       enddo
     enddo
#else
!
! ------------------------------------------------------------------
!     calculate stratiform cloud and put into array 'cldary' using
!       the cloud-rel.humidity relationship from table look-up..where
!       tables obtained using k.mitchell frequency distribution tuning
!        (observations are daily means from us af rtneph).....k.a.c.
!       tables created without lowest 10 percent of atmos.....k.a.c.
! ------------------------------------------------------------------
!  this loop to retrieve cloud from rh rewritten 950113 -mi
!
     do klev = klowb,levm2
       do i = 1,imx2
         kbase(i)=0
         bitx(i)=.false.
       enddo
       do kc = mcld,1,-1
         do i = 1,imx2
           if(prsly(i,klev).ge.ptop1(i,kc+1)) kbase(i)=kc
         enddo
       enddo
       nx=0
       nhalf=(nbin+1)/2
       do i = 1,imx2
         if(kbase(i).le.0.or.klev.gt.kcut(i)) then
           cldary(i,klev)=0.
         elseif(rhrh(i,klev).le.rhcld(i,1,kbase(i))) then
           cldary(i,klev)=0.
         elseif(rhrh(i,klev).ge.rhcld(i,nbin,kbase(i))) then
           cldary(i,klev)=1.
         else
           bitx(i)=.true.
           ksave(i)=nhalf
           nx=nx+1
         endif
       enddo
       do while(nx.gt.0)
         nhalf=(nhalf+1)/2
         do i = 1,imx2
           if(bitx(i)) then
             crk=rhrh(i,klev)
             cr1(i)=rhcld(i,ksave(i),kbase(i))
             cr2(i)=rhcld(i,ksave(i)+1,kbase(i))
             if(crk.le.cr1(i)) then
               ksave(i)=max(ksave(i)-nhalf,1)
             elseif(crk.gt.cr2(i)) then
               ksave(i)=min(ksave(i)+nhalf,nbin-1)
             else
               cldary(i,klev)=0.01*(ksave(i)+(crk-cr1(i))/(cr2(i)-cr1(i)))
               bitx(i)=.false.
               nx=nx-1
             endif
           endif
         enddo
       enddo
     enddo
!
!     clean out not-suspected marine stratus regions...
!      cause tuning procedure not carried out down to lyr3 and we
!      get too much lo cloud if we don t clean it out..
!
     do i = 1,imx2
       bitm(i) = .true.
     enddo
     do k = klowb,llyr
       do i = 1,imx2
         if(bitm(i)) then
           bitm(i) = prsly(i,k).lt.pstrt                                       &
              .or. slmsk(i).gt.0.0 .or. dthdp(i,k).gt.clapkc.or.               &
               rhrh (i,k+1).gt.0.60e0.or.rhrh (i,k+2).gt.0.60e0
           kbase(i) = k
         endif
       enddo
     enddo
     do k = 1,llyr
       do i = 1,imx2
         if(bitm(i)) cldary(i,k) = 0.0e0
       enddo
     enddo
!
!  special treatment on low clouds
!
     dvvcld = vvcld(1) - vvcld(2)
!
     do k = klowb,klow
!
       do i = 1,imx2
         omeg(i) = 10.0e0 * vvel(i,k)
         cl1 (i) = 0.0e0
       enddo
!
       if (ivva .le. 0) go to 920
!
!  vertical velocity adjustment on low clouds
!
       bit1 = .false.
       do i = 1,imx2
         bitx(i) = prsly(i,k).ge.ptop1(i,2) .and. cldary(i,k).gt.0.0e0
         bit1 = bit1 .or. bitx(i)
       enddo
!
       if (.not. bit1) go to 920
!
       if(k.gt.llyr) then
         do i = 1,imx2
           if (bitx(i)) then
             if(omeg(i).ge.vvcld(1)) then
               cldary(i,k) = 0.0e0
             else if(omeg(i).gt.vvcld(2)) then
               cr1(i) = (vvcld(1) - omeg(i)) / dvvcld
!              cldary(i,k) = cldary(i,k) * cr1(i)
               cldary(i,k) = cldary(i,k) * sqrt(cr1(i))
             endif
           endif
         enddo
       else
!
!  no vvel filter for marine stratus region
!
         do i = 1,imx2
           if (bitm(i)) then
             if (bitx(i)) then
               if(omeg(i).ge.vvcld(1)) then
                 cldary(i,k) = 0.0e0
               else if(omeg(i).gt.vvcld(2)) then
                 cr1(i) = (vvcld(1) - omeg(i)) / dvvcld
                 cldary(i,k) = cldary(i,k) * sqrt(cr1(i))
               endif
             endif
           endif
         enddo
       endif
!
! t inversion related stratus clouds
!
920    continue
!
       if (invr .lt. 1) cycle
       if (k.gt.llyr) cycle
!
       bit1 = .true.
       do i = 1,imx2
       bit1 = bit1 .and. bitm(i)
     enddo
     if (bit1) cycle
       do i = 1,imx2
         if (.not.bitm(i)) then
           if (dthdp(i,kbase(i)).gt.clpse) then
!
!   smooth transition for cloud when dthdp between
!       clapse and clapse+dclps  (-0.05 and -0.06)
!
             cfiltr = 1.0e0 - ((clpse - dthdp(i,kbase(i))) / dclps)
             cldary(i,k) = cldary(i,k)*cfiltr
           end if
!
!  for t inversion type cloud, add flag value of 2.0
!
           cldary(i,k) = cldary(i,k)+2.0e0
         end if
       enddo
     enddo
! 
#endif
   end if
!
! ------------------------------------------------------------------
!     add convective cloud into 'cldary', no merge at this point..
!     two types of clouds are separated by a factor of 1.0e+4
! ------------------------------------------------------------------
!
   bit1 = .false.
   do i = 1,imx2
     bitx(i) = cv(i).gt.0.0e0 .and. cvt(i).ge.cvb(i)
     bit1 = bit1 .or. bitx(i)
   enddo
!
   if (.not. bit1) go to 550
!
   do i = 1,imx2
     if (bitx(i)) then
       kcvb(i) = nint(cvb(i))
       kcvt(i) = min(levm2, nint(cvt(i)))
     else
       kcvb(i) = 1
       kcvt(i) = 1
     end if
   enddo
   do k = klowb,levm2
     bit2 = .false.
     do i = 1,imx2
       bity(i) = bitx(i) .and. kcvb(i).le.k .and. kcvt(i).ge.k
       bit2 = bit2 .or. bity(i)
     enddo
     if (.not. bit2) cycle
     do i = 1,imx2
       if (bity(i)) cldary(i,k) = cldary(i,k)                                  &
                  + 10.0e0 * aint(1.0e3 * cv(i))
     enddo
   enddo
!
!     if mean cvt layer higher than 400mb add anvil cirrus
!
   bit2 = .false.
   do i = 1,imx2
      bitz(i) = bitx(i) .and. prsly(i,kcvt(i)).le.cvtop
      bit2 = bit2 .or. bitz(i)
   enddo
!
   if (.not. bit2) go to 500
!
   do i = 1,imx2
     if (bitz(i)) then
       kk = kcvt(i)
       cr1(i) = max (0.0e0, min (1.0e0, 2.0e0*(cv(i)-0.3e0)))
!
! get stratus back before doing anvil calculation
!
       cr2(i) =  mod(cldary(i,kk),10.e0)
       cldary(i,kk) = cr2(i) + 10.0e0*aint(1.0e3*cr1(i))
     end if
   enddo
!
500 continue
550 continue
!
!  introduce adjustment of clouds (mk)
!
#ifdef CLDADJ
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
   if(lat.eq.10) then
     print *,'before rad_cloudiness_adjust lat=',lat
     print *,(cldary(20,k),k=1,levs_)
   endif
   if(lat.eq.30) then
     print *,'before rad_cloudiness_adjust lat=',lat
     print *,(cldary(20,k),k=1,levs_)
   endif
#ifdef OPENMP
!$omp end master
#endif
#endif
   call rad_cloudiness_adjust(cldary,imx22,imx2,kmx,ptop1,ILOTS,prsi,xlatrd)
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
   if(lat.eq.10) then
     print *,'after rad_cloudiness_adjust lat=',lat
     print *,(cldary(20,k),k=1,levs_)
   endif
   if(lat.eq.30) then
     print *,'after rad_cloudiness_adjust lat=',lat
     print *,(cldary(20,k),k=1,levs_)
   endif
#ifdef OPENMP
!$omp end master
#endif
#endif /* DBG end */
#endif
!
! -------------------------------------------------------------------
!     separate clouds into 3 pressure domains (l,m,h).  within each
!     of the domains, assume separated cloud layers are randomly
!     overlapped and adjacent cloud layers are maximum overlapped.
!     vertical location of each type of cloud is determined by
!     the thickest continuing cloud layers in the domain.
! -------------------------------------------------------------------
!  loop over 3 cloud domains (l,m,h)
!
   do  l =1,3
!
     do i = 1,imx2
       cld (i,l) = 0.0e0
       mtop(i,l) = 1
       mbot(i,l) = 1
       cl1 (i) = 0.0e0
       cl2 (i) = 0.0e0
       kbt1(i) = 1
       kbt2(i) = 1
       kth1(i) = 0
       kth2(i) = 0
     enddo
!
     do k = 2,levm2
       bit1 = .false.
       do i = 1,imx2
         bitx(i) = (prsly(i,k).ge.ptop1(i,l+1)) .and.                          &
                     (prsly(i,k).lt.ptop1(i,l)) .and. (cldary(i,k).gt.0.0e0)
         bit1 = bit1 .or. bitx(i)
       enddo
!
       if (.not. bit1) cycle
!
       do i = 1,imx2
         cr1(i)  =  mod(cldary(i,k), 2.0e0)
         cr2(i)  = float(int(cldary(i,k)) / 10) * 1.0e-3
         if (bitx(i)) then
!
! kth2 le 0 : 1st cld layer.
!
           if(kth2(i).le.0) then
             kbt2(i) = k
             kth2(i) = 1
!
!  kth2 gt 0 : consecutive cld layer.
!
           else
             kth2(i) = kth2(i) + 1
           endif
!
!   physical cloud as seen by radiation..conv takes precedence
!   except anvil cirrus not random overlapped with cv tower as
!   in radiation code(so hi may be slight underestimate)....
!
           if (cr2(i).gt.0.0e0) then
             cl2 (i) = max (cl2(i), cr2(i))
           else
             cl2 (i) = max (cl2(i), cr1(i))
           end if
         endif
       enddo
       bit2 = .false.
!
!  bity=true if next lyr=clear or we change cloud domains..
!
       do i = 1,imx2
         bity(i) = bitx(i) .and. (cldary(i,k+1).le.0.0e0                       &
                            .or.  prsly(i,k+1).lt.ptop1(i,l+1) )
         bit2 = bit2 .or. bity(i)
       enddo
!
       if (.not. bit2) cycle
!
!     at the domain boundary or separated cld lyrs, random overlap.
!     choose the thickest or the largest fraction amt as the cld
!     layer in that domain
!
       do i = 1,imx2
         if (bity(i)) then
           if (cl1(i).gt.0.0e0) then
             kbt1(i) = int( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i))                  &
                          / (cl1(i) + cl2(i)) )
             kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i))                 &
                           / (cl1(i) + cl2(i)) ) + 1
             cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)
           else
             kbt1(i) = kbt2(i)
             kth1(i) = kth2(i)
             cl1 (i) = cl2 (i)
           endif
           kbt2(i) = 1
           kth2(i) = 0
           cl2 (i) = 0.0e0
         endif
       enddo
     enddo
!
!  finish one domain, save effective clouds
!
     do i = 1,imx2
       cld(i,l) =  cl1(i)
       mtop(i,l) = max(kbt1(i), kbt1(i)+kth1(i)-1)
       mbot(i,l) = kbt1(i)
     enddo
   enddo
!
!      calculate total cloud from the multi-lyr cloud array
!      in a manner as seen by the radiation code.
!      where, max overlap is used for vertically adjacent cloud layers
!      a clear layer separates two contiguously layered cloud types.
!      where, for convection any anvil is considered a separate
!         randomly overlapped cloud..
!      ilow=0,1 if no,yes preceeding model layer was cloudy..
!      clow contains the cloudiness of preceeding separate layered cld
!
   do i = 1,imx2
     cld(i,4) = 0.e0
     icvec = 0
     ilow = 0
     clow = 0.e0
     do k = 1,kdim
       ccldy =  mod(cldary(i,k), 2.0e0)
       ccvec = float(int(cldary(i,k)) / 10) * 1.0e-3
       if (ccvec.gt.0.e0) then
         ccldy = ccvec
         icvec = 1
       end if
       if (ccldy.gt.0.e0) then
         if (ilow.eq.0) then
           clow = ccldy
           ilow = 1
         else
           if (icvec.gt.0) then
!
!  if convective and an adjacent lyr=stratiform (ie ccldy changes),
!  then random overlap the preceeding cloud tower...
!
             if (ccldy.ne.clow) then
               cld(i,4) = cld(i,4) + (1.-cld(i,4))*clow
               clow = ccldy
             end if
!
! max overlap for non convective adjacent cld layers...
!
           else
             clow = max(ccldy,clow)
           end if
         end if
       else
         if (ilow.eq.1) then
!
!  if this is first clear layer in a gap betwixt cldlyrs, then
!   random overlap the preceeding clouds with the ones below..
!
           cld(i,4) = cld(i,4) + (1.-cld(i,4))*clow
           ilow = 0
         end if
       end if
     enddo
   enddo
!
   return
   end subroutine rad_cloudiness
!
!-------------------------------------------------------------------------------
#undef ADJCNV
#include <define.h>
   subroutine rad_cloudiness_adjust(cldary,imx22,imx2,kmx,                     &
                                    ptop1,ilots,prsi,xlatrd)
!-------------------------------------------------------------------------------
!
!  adjust full level of model cloudiness using latitude dependent ratio table
!
!-------------------------------------------------------------------------------
   use constant, only    :  pi_
!
   integer,parameter    ::  idmr=19
   real                 ::  cldary(imx22,kmx)
   real                 ::  ptop1(ilots,4)
   real                 ::  prsi(imx22,kmx+1)
   real                 ::  xlatrd(imx2)
   real                 ::  ratio(idmr,3)
!
   logical lstr
!
!     data ratio/1.2, 1.4, 1.4, 1.2, 1.1, 1.2, 1.8, 2.1, 1.8, 1.3,
!    1           1.5, 1.7, 1.5, 1.3, 1.1, 1.0, 1.0, 1.1, 1.3,
!    2           0.8, 0.7, 0.8, 0.8, 0.8, 0.7, 0.6, 0.5, 0.5, 0.6,
!    3           0.5, 0.5, 0.7, 0.8, 0.5, 0.5, 0.5, 0.5, 0.5,
!    4           0.4, 0.2, 0.4, 0.4, 0.5, 0.7, 0.6, 0.5, 0.7, 0.7,
!    5           0.5, 0.4, 0.6, 0.8, 0.6, 0.2, 0.2, 0.2, 0.2/
!
! 1st iteration
!     data ratio/1.0, 0.2, 0.6, 0.6, 0.6, 1.0, 1.5, 1.4, 1.1, 0.6,
!    1           0.8, 1.3, 1.4, 1.0, 0.7, 0.6, 0.8, 0.5, 0.3,
!
!    2           0.1, 0.5, 0.9, 0.8, 0.8, 0.7, 0.6, 0.5, 0.4, 0.4,
!    3           0.4, 0.3, 0.6, 1.0, 1.1, 1.1, 1.1, 1.3, 1.3,
!
!    4           1.2, 0.9, 0.3, 0.4, 0.6, 0.6, 0.6, 0.5, 0.6, 0.5,
!    5           0.5, 0.3, 0.4, 0.7, 0.6, 0.3, 0.2, 0.2, 0.1/
! final iteration
!     data ratio/1.0, 0.2, 0.6, 0.6, 0.7, 1.5, 1.5, 1.3, 1.1, 0.5,
!    1           1.2, 1.3, 1.4, 1.0, 0.6, 0.5, 0.8, 0.4, 0.2,
!
!    2           0.1, 0.7, 1.3, 0.9, 0.7, 0.6, 0.4, 0.2, 0.2, 0.2,
!    3           0.3, 0.2, 0.5, 1.0, 1.2, 1.1, 1.2, 1.4, 1.1,
!
!    4           2.0, 1.6, 0.1, 0.2, 0.5, 0.6, 0.5, 0.3, 0.5, 0.3,
!    5           0.5, 0.2, 0.6, 0.6, 0.4, 0.1, 0.1, 0.1, 0.1/
! yet final iteration
!     data ratio/1.5, 1.5, 1.5, 0.6, 0.7, 1.8, 1.3, 1.4, 1.0, 0.5,
!    1           1.2, 1.2, 1.3, 1.0, 0.5, 0.6, 0.8, 0.8, 0.8,
!
!    2           1.6, 1.6, 1.6, 0.8, 0.7, 0.5, 0.3, 0.1, 0.1, 0.1,
!    3           0.3, 0.1, 0.4, 0.8, 0.9, 1.0, 1.0, 1.2, 1.2,
!
!    4           0.05, 0.05, 0.05, 0.1, 0.6, 0.7, 0.6, 0.3, 0.6, 0.4,
!    5           0.7, 0.2, 1.0, 0.6, 0.3, 0.1, 0.05, 0.05, 0.05/
!
   data ratio/0.4, 0.4, 0.6, 0.6, 0.7, 1.5, 1.5, 1.4, 1.0, 1.0,                &
              1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.8, 0.4, 0.2,                     &
              0.5, 0.5, 0.5, 0.7, 0.7, 0.5, 0.3, 0.02, 0.02, 0.01,             &
              0.01, 0.02, 0.4, 0.8, 0.9, 1.0, 1.0, 1.2, 1.2,                   &
              0.02, 0.02, 0.02, 0.3, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,             &
              0.6, 0.6, 0.6, 0.6, 0.3, 0.05, 0.02, 0.02, 0.02/
!-------------------------------------------------------------------------------
   do k = 1,kmx 
     do i = 1,imx2
!
!  find index of ratio from latitude and level
!
       kc=0
       do kk = 1,3
         if(prsi(i,k)/prsi(i,1).le.ptop1(i,kk).and.                            &
            prsi(i,k)/prsi(i,1).gt.ptop1(i,kk+1)) kc=kk
       enddo
       if(kc.gt.0.and.cldary(i,k).gt.0.) then
!
         rlat=float(idmr-1)/pi_*xlatrd(i)+float(idmr+1)/2.
         ilat=int(rlat)
         ilat=min(max(ilat,1),idmr)
!
! interpolate ratio to the model grid
!
         ilatp=min(ilat+1,idmr)
         ratiox=ratio(ilat,kc)+                                                &
           (ratio(ilatp,kc)-ratio(ilat,kc))*(rlat-float(ilat))
!
! identify and extract convective clouds
!
         if(cldary(i,k).gt.1.e2) then
           ccnv=float(int(cldary(i,k))/10)*1.0e-3
           cstr=cldary(i,k)-10.e0*aint(1.e3*ccnv)
           if(ccnv.lt.0..or.ccnv.gt.1.0) then
             write(6,*)'convective cloud amount incorrect. ccnv=',ccnv
#ifdef MP
#ifdef RMP
             call rmpabort
#else
             call mpabort
#endif
#else
             call abort
#endif
           endif
         else
           ccnv=0.
           cstr=cldary(i,k)
         endif
!
! identify temp inversion clouds
!
         lstr=.false.
         if(cstr.gt.2) then
           lstr=.true.
         endif
         cstr=mod(cstr,2.0e0)
!
         if(cstr.lt.0..or.cstr.gt.1.0) then
           write(6,*) 'stable cloud amount incorrect. cstr=',cstr
#ifdef MP
#ifdef RMP
           call rmpabort
#else
           call mpabort
#endif
#else
           call abort
#endif
         endif
!
!  adjust cloud amount
!
#ifdef ADJCNV
#ifdef ICECLOUD
         ccnv=min(max(ccnv*ratiox,0.),1.0)
#else
         ccnv=ccnv*ratiox
#endif
#else
         ccnv=ccnv
#endif
#ifdef ICECLOUD
         cstr=min(max(cstr*ratiox,0.),1.0)
#else
         cstr=cstr*ratiox
#endif
!
         if(lstr) then
           cstr=cstr+2.0e0
         endif
!
!  reconstruct cldary
!
         cldary(i,k)=cstr+10.0e0*aint(1.0e3*ccnv)
!
       endif
     enddo
   enddo
!
   return
   end subroutine rad_cloudiness_adjust
!
!-------------------------------------------------------------------------------
   subroutine rad_ozone_gfdl(ims2,imx2,kmx,                                    &
                     prsl,qo3,xlat,rsin1,rcos1,rcos2)
!-------------------------------------------------------------------------------
!  subroutine: rad_ozone_gfdl
!
!  abstract:
!    compute model lyr o3 profile from the original gfdl data
!
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS
   use constant, only  : pi_
   use comfcst, only   : nl, prgfdl                                           ,&
                         dduo3n=>xduo3n, ddo3n2=>xdo3n2                       ,&
                         ddo3n3=>xdo3n3, ddo3n4=>xdo3n4
!-------------------------------------------------------------------------------
   real                 ::  prsl(imx2,kmx)
   real                 ::  qo3(imx2,kmx),xlat(imx2)
!
   integer              ::  jjrow(ILOTS)
   real                 ::  tthan(ILOTS)
   real                 ::  qo3o3(ILOTS,nl),psm(ILOTS)
!-------------------------------------------------------------------------------
!
!...       begin here .....
!
   rndg = 180./pi_
   do i = 1,ims2
     th2= 0.2e0*xlat(i)*rndg
     jjrow(i)=19.001e0-th2
     tthan(i)=(19-jjrow(i))-th2
   enddo
10 continue
!
!....   seasonal and spatial interpolation done below.
!
   do k = 1,nl
     do i = 1,ims2
       do3v = dduo3n(jjrow(i),k) + rsin1*ddo3n2(jjrow(i),k)                    &
                +rcos1*ddo3n3(jjrow(i),k)                                      &
                +rcos2*ddo3n4(jjrow(i),k)
       do3vp = dduo3n(jjrow(i)+1,k) + rsin1*ddo3n2(jjrow(i)+1,k)               &
                +rcos1*ddo3n3(jjrow(i)+1,k)                                    &
                +rcos2*ddo3n4(jjrow(i)+1,k)
!
!...   now latitudinal interpolation, and
!          convert o3 into mass mixing ratio(original data mpy by 1.e4)
!      flip vertical coordinate too...
!
       qo3o3(i,nl+1-k) = 1.e-4 * (do3v+tthan(i)*(do3vp-do3v))
     enddo
   enddo
!
!...    vertical (linear in ln p) interpolate for each gridpoint
!
   numitr = 0
   ilog = nl
!
21 continue
!
   ilog = (ilog+1)/2
   if(ilog.eq.1) go to 22
   numitr = numitr + 1
!
   go to 21
!
22 continue
!
   do k = 1,kmx
     nhalf=(nl+1)/2
     do i = 1,ims2
       jjrow(i) = nhalf
       psm(i) = prsl(i,k)
     enddo
     do it=1,numitr
       nhalf=(nhalf+1)/2
       do i = 1,ims2
         if(psm(i).lt.prgfdl(jjrow(i))) then
           jjrow(i) = jjrow(i) + nhalf
         else if(psm(i).ge.prgfdl(jjrow(i)-1)) then
           jjrow(i) = jjrow(i) - nhalf
         endif
         jjrow(i) = min(jjrow(i),nl)
         jjrow(i) = max(jjrow(i),2)
       enddo
     enddo
     do i = 1,ims2
       if(psm(i).gt.prgfdl(1)) then
         qo3(i,k) = qo3o3(i,1)
       else if(psm(i).lt.prgfdl(nl)) then
         qo3(i,k) = qo3o3(i,nl)
       else
         aplo = alog(prgfdl(jjrow(i)-1))
         apup = alog(prgfdl(jjrow(i)))
         qo3(i,k) = qo3o3(i,jjrow(i)) + (alog(psm(i))-apup) /                  &
                    (aplo-apup) *                                              &
                    (qo3o3(i,jjrow(i)-1)-qo3o3(i,jjrow(i)))
       endif
     enddo
   enddo
!
   return
   end subroutine rad_ozone_gfdl
!
!-------------------------------------------------------------------------------
   subroutine rad_ozone_nasa(ims2,imx2,kmx,prsl,qo3,xlat,lat)
!-------------------------------------------------------------------------------
!
! abstract:
!    compute model lyr o3 profile from the original nasa data
! 
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS
   use constant, only  : pi_
#ifdef DFS
   use dfsvar, only    : iope
#endif
   use comfcst, only   : loz, psnasa, o3nasa
!-------------------------------------------------------------------------------
   real                 ::  prsl(imx2,kmx)
   real                 ::  qo3(imx2,kmx),xlat(imx2)
!
   real                 ::  wgts(ILOTS)
   real                 ::  psm(ILOTS)
   real                 ::  qo3o3(ILOTS,loz)
   integer              ::  jlats(ILOTS)
   integer              ::  jjrow(ILOTS)
!-------------------------------------------------------------------------------
!
!...       begin here .....
!
   rndg = 180./pi_
   do 10 i=1,ims2
     rlat=xlat(i)*rndg
     do j = 1,37
       lato3=90-(j-1)*5
       xlato3=lato3
       if (xlato3.lt.rlat) then
!
! add hoon
!
         if (j.eq.1) then
           jlats(i)=j+1
         else
           jlats(i)=j
         endif
         wgts(i) = (rlat-xlato3)/5.
!          go to 1013
!
         go to 10
!
       end if
     enddo
10 continue
!
!....  latitudinal (linear) interpolation done below.
!       flip vertical coordinate too...
!
   do k = 1,loz
     do i = 1,ims2
       qo3o3(i,loz+1-k)=wgts(i)                                                &
                      * (o3nasa(jlats(i)-1,k)-o3nasa(jlats(i),k))              &
                      +  o3nasa(jlats(i),k)
     enddo
   enddo
!
!...    vertical (linear in ln p) interpolate for each gridpoint
!
   numitr = 0
   ilog = loz
!
21 continue
!
   ilog = (ilog+1)/2
!
   if(ilog.eq.1) go to 22
!
   numitr = numitr + 1
!
   go to 21
!
22 continue
!
   do k = 1,kmx
     nhalf=(loz+1)/2
     do i = 1,ims2
       jjrow(i) = nhalf
       psm(i) = prsl(i,k)
     enddo
     do it = 1,numitr
       nhalf=(nhalf+1)/2
       do i = 1,ims2
         if(psm(i).lt.psnasa(jjrow(i))) then
           jjrow(i) = jjrow(i) + nhalf
         else if(psm(i).ge.psnasa(jjrow(i)-1)) then
           jjrow(i) = jjrow(i) - nhalf
         endif
!           jjrow(i) = min(jjrow(i),nl)
         jjrow(i) = min(jjrow(i),loz)
         jjrow(i) = max(jjrow(i),2)
       enddo
     enddo
     do i = 1,ims2
       if(psm(i).gt.psnasa(1)) then
         qo3(i,k) = qo3o3(i,1)
       else if(psm(i).lt.psnasa(loz)) then
         qo3(i,k) = qo3o3(i,loz)
       else
         aplo = alog(psnasa(jjrow(i)-1))
         apup = alog(psnasa(jjrow(i)))
         qo3(i,k) = qo3o3(i,jjrow(i)) + (alog(psm(i))-apup) / (aplo-apup) *    &
                 (qo3o3(i,jjrow(i)-1)-qo3o3(i,jjrow(i)))
       endif
     enddo
   enddo
   !
   return
   end subroutine rad_ozone_nasa
!
!-------------------------------------------------------------------------------
   subroutine rad_albedo8snow(imx2,latrue,jsno,                                &
                              albdoa,rlat,albedr,slmskr,ssnow,tgr,tar)
!-------------------------------------------------------------------------------
   use paramodel
   use constant, only    : pi_
!-------------------------------------------------------------------------------
   real                 ::  albdoa(imx2),albedr(imx2),slmskr(imx2)
   real                 ::  ssnow(imx2),rlat(imx2)
   real                 ::  snodeg, snochk
!
! added by bob grumbine for sea ice albedo algorithm
!
   real                 ::  tgr(imx2), tar(imx2)
!-------------------------------------------------------------------------------
!
! modified by hmh juang for simplicity  and for regional model use
!
   snodeg = 70.*pi_/180.
!
! the following determines surface albedo (albdoa),where snow exists.
!
   snochk = pi_
   if(jsno.eq.0) snochk=snodeg  ! jsno=0 for rmp
!
   do i = 1, imx2
!
! limit background albedo (in case snow leaves greenland)
!
     albdoa(i)=min(albedr(i),0.6e0)
     if(slmskr(i).eq.1.0e0) then
!
! check lat<jsno for global and abs(rlat)>snodeg for regional
!
       if(latrue.lt.jsno .or. abs(rlat(i)).gt.snochk) then
         if(ssnow(i).gt.0.e0) albdoa(i)=0.75e0
       else
         if(ssnow(i).ge.1.e0) then
           albdoa(i)=0.6e0
         else if(ssnow(i).gt.0.0e0) then
           albdoa(i)=albdoa(i)+sqrt(ssnow(i))*(0.6e0-albdoa(i))
         endif
       endif
     else if(slmskr(i).eq.2.0e0) then
       if (ssnow(i) .gt. 0.0) then
         if (tgr(i) .lt. 273.16 - 5.) then
           albdoa(i) = 0.8
         else if (tgr(i) .le. 273.16) then
           albdoa(i) = 0.65 + 0.03*(273.16 - tgr(i))
         else
           albdoa(i) = 0.65
         endif
       else
         if (tgr(i) .lt. 271.2 .or. tar(i) .lt. 273.16) then
           albdoa(i) = 0.65
         else if (tar(i) .lt. 273.16+5.) then
           albdoa(i) = 0.65 - 0.04*(tar(i) -273.16)
         else
           albdoa(i) = 0.45
         endif
       endif
     endif
     albdoa(i) = max(albdoa(i),.06e0)
   enddo
!
   return
   end subroutine rad_albedo8snow
!
!-------------------------------------------------------------------------------
   subroutine rad_cloudiness_nasa(imx2,imx22,kmx,nbin,mcld,                    &
              q,t,vvel,cv,cvt,cvb,prsi,prsl,slmsk,                             &
              cld,mtop,mbot,                                                   &
              cldtot,cldcnv,ivva,invr,icwp,qci,rhrh,                           &
#ifdef EXPLICIT_CLOUDINESS
              qcicps, qrscps,                                                  &
              qci_cps,                                                         &
#endif
#ifdef DBG
              xlatrd,rhcld,istrat,ntotal,lat)
#else
              xlatrd,rhcld,istrat,ntotal)
#endif
!-------------------------------------------------------------------------------
!
!     nov., 1992 - y.h., k.a.c., and a.k.                            
!        cloud parameterization patterned after slingo and slingos  
!        work (jgr, 1991).                                         
!     stratiform clouds are allowed in any layer except the surface 
!        and upper stratosphere.  the relative humidity criterion may 
!        vary in different model layers.                             
!     output cloud amounts are in cldary(i,k), k=1 is the lowest   
!        model layer, stratiform (str) and convective (cnv) types of 
!        cloud are compressed into one word: camt = str + 1.0e4*cnv 
!        low marine stratus amts are flaged by adding 2.           
!
!..   for istrat = 0, there is rh-cld relation for each layer.. 
!                      crit rh computed within..               
!..   for istrat = 1, rh-cld relation from tables created using  
!                     mitchell-hahn tuning technique (a.f. rtneph obs)
!          ...stratus computed similar to old opnl rad_cloudiness.....
!                      except no cloud below layer=klowb..approx 955mb
!
!     convective clouds are from model convective scheme and are    
!        no longer broken into .75,.25,.25..rather cc itself is used..
!        convective still takes precedence over stratiform in radfs  
!         but here cv+st merge exits in cldary...(in radiation use of
!        cc gives improvement to tropical middle cld (as did st+cv))
!                                                                  
!     clouds are also divided into 3 atmospheric domains (l,m,h) for 
!        diagnostic purposes.  they are computed from random overlap
!        assumption for separated cloud layers and maximum overlap 
!        for adjacent cloud layers.  a total cloud fraction is also 
!        computed.                                                 
!                                                                 
!     h,m,l domain pressure tops 'ptop1(k)' vary linearly from   
!        'ptopc(k,1)' at 45deg to 'ptopc(k,2)' at the pole      
!                                                              
!     aug 1998   y-t hou    recoded the default part of cloud scheme
!         patterned after ccm3 (kiehl et al. 1998,j.clim; and 1994,jgr).
!         output rh and conv clouds use seperated arrays cldtot and
!         cldcnv to avoid packing-unpacking process.
!
!     sep 1999   y-t hou 
!         output rhrh variable for opac type aerosol optical
!         properties calc. later in sw radiation.
!
!     aug. 2000 --- y.h.
!        convective clouds are adjusted by rh. the overlapping of
!        thick clouds for lw is determined by the differences between
!        two contiguous cloud layers.  the h,m,l, and total cloud
!        amounts are calculated from lw radiation point of view. the
!        total cloud amount from sw point view is given in the
!        subroutine 'rad_cloud_comp' (e.g. 1-cfac at the surface).
!
!-------------------------------------------------------------------------------
!     input variables:                                               
!        ps (cb)       - surface pressure                           
!        q  (kg/kg)    - specific humidity                         
!        t  (deg k)    - absolute temperature                     
!        vvel(cb/sec)  - vertical velocity                       
!        cv,cvt,cvb    - conv cld fraction, top, bottom layer from  
!        si,sl         - mdl sigma interface and layer mean       
!        icwp          - =0 rh based cloud scheme
!                        =1 liquid water cloud scheme
!        qci           - cloud liqui water
!        slmsk         - sea/land mask array(sea:0.,land:1.,snow:2.) 
!        ivva          - flag to control vertical velocity adj.     
!                        =1: with, =0: without                     
!        invr          - flag to control lapse rate inversion cld 
!                        =1: with, =0: without                   
!********xlatrd        - current latitude in radians (1st data pt)  
!********                 for models with diff lat at each pt, need to
!********                 use the lat of all points....careful.....  
!        rhcld         - cloud-rh relations from mitchell+hahn,     
!                        using a.f. rtneph analyses                
!        istrat        - 0 or 1:for default or 'rhcld' tables     
!                        in the stratiform cloud calculation     
!tune                                                           
!    output variables:                                         
!       cldtot         - vertical column array of stratiform cloud
!       cldcnv         - vertical column array of convective cloud
!       cld            - cld fraction in 3 types of domains (l,m,h)
!                        and total in 4th layer                  
!       mtop,mbot      - top, bottom layers of clouds (l,m,h)     
!       rhrh           - relative humidity
!                                                                
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_,ngases_,icloud_,nwmass_
   use constant, only  : g_,pi_,rd_,rv_,akapa_
   use comcd1
!-------------------------------------------------------------------------------
!                                                                   
   real, parameter      ::  rd=rd_, rv=rv_, eps=rd/rv, epsm1=rd/rv-1.0, pi=pi_ 
!
! --- input variable arrays
!
   real                 ::  cv(imx22),cvt(imx22),cvb(imx22),                   &
                            slmsk(imx22),prsi(imx22,kmx+1),prsl(imx22,kmx),    &
                            t(imx22,kmx),vvel(imx22,kmx),q(imx22,kmx*ntotal),  &
                            xlatrd(imx22),qci(imx22,kmx),                      &
                            delts(imx22),clwt(imx22)
#ifdef EXPLICIT_CLOUDINESS
   real                 ::  qcicps(imx22,kmx), qrscps(imx22,kmx)
   real                 ::  qci_cps(imx22,kmx),tmp
#endif
!...    rh-cld relationships for each point              
   real                 ::  rhcld(imx22,nbin,mcld)                  
!
! --- output variable arrays
!
   real                 ::  cld(imx22,4),cldtot(imx22,kmx),cldcnv(imx22,kmx),  &
                            rhrh(imx22,kmx)
   integer              ::  mtop(imx22,3), mbot(imx22,3) 
!                                                        
! --- ptopc(k,l): top presure of each cld domain (k=1-4 are sfc,l,m,h;
!       l=1,2 are low-lat (<45 degree) and pole regions)             
!
   logical              ::  bitx(ILOTS), bity(ILOTS), bitm(ILOTS), bit1,       &
                            bit2, inversn(ILOTS)
   real                 ::  adcr (ILOTS,levs_),prsly(ILOTS,levs_),             &
                            dthdp(ILOTS,levs_),theta(ILOTS,levs_),             &
                            cl1  (ILOTS),      cl2  (ILOTS),                   &
                            clwp (ILOTS,levs_),omeg (ILOTS),                   &
                            ptop1(ILOTS,4),    qs   (ILOTS,levs_)
   integer              ::  kcut (ILOTS),      kbase(ILOTS),                   &
                            kbt1 (ILOTS),      kth1 (ILOTS),                   &
                            kbt2 (ILOTS),      kth2 (ILOTS),                   &
                            ksave(ILOTS),      ityp1(ILOTS),                   &
                            ityp2(ILOTS),      kinver(ILOTS)
!-------------------------------------------------------------------------------
!===>    begin here ................................................ 
!
   kdim  = kmx                                                      
   kdimp = kmx+1                                                  
   levm1 = kmx-1                                                 
   levm2 = kmx-2                                                
   rhh  = 0.95e0           ! critical value for hi  cld
   rhm  = 0.70e0           ! critical value for mid cld
   rhl  = 0.70e0           ! critical value for low cld
!
!...  find top pressure for each cloud domain                   
!
   do k = 1,4                                                 
     dptop = ptopc(k,2)-ptopc(k,1)
     do i = 1,imx2                                           
       fac = max(0.0e0, 4.0e0*abs(xlatrd(i))/pi-1.0e0)    
       ptop1(i,k) = ptopc(k,1) + dptop * fac 
     enddo
   enddo                                                
!
! --- low cloud top sigma level, computed for each lat cause 
!       domain definition changes with latitude...          
   klow = kdim                                            
   do k = kdim,1,-1
     do i = 1,imx2 
       if (prsi(i,k)/prsi(i,1) .lt. ptop1(i,2)*1.0e-3) klow = min(klow,k)
     enddo
   enddo
!
! --- potential temp and layer relative humidity               
!
   do k = 1,kdim
     do i = 1,imx2
       cldtot(i,k) = 0.0e0
       cldcnv(i,k) = 0.0e0
       prsly(i,k) = prsl(i,k) * 10.0e0
       exnr = (prsly(i,k)*0.001e0) ** (-akapa_)
       theta(i,k) = exnr * t(i,k)
#ifdef ICE
       es = fpvs(t(i,k))
#else
       es = fpvs0(t(i,k))
#endif
       qs(i,k) = eps * es / (prsl(i,k) + epsm1*es)
       rhrh(i,k) = max (0.0e0, min (1.0e0, q(i,k)/qs(i,k)))
     enddo
   enddo
!
! --- potential temp lapse rate                               
!
   do k = 1,levm1
     do i = 1,imx2
       dthdp(i,k) = (theta(i,k+1) - theta(i,k)) /                              &
                    (prsly(i,k+1) - prsly(i,k))
     enddo
   enddo
!
! --- cloud liquid water path
!
   if(icwp .eq. 1) then
     do k = 1,kdim
       do i = 1,imx2
         delts(i) = (prsi(i,k) - prsi(i,k+1))*1.0e6 / prsi(i,1) / g_
       enddo      
       do i = 1,imx2
#ifdef EXPLICIT_CLOUDINESS
#ifdef CPS_QSUM
         tmp = qci(i,k) + qci_cps(i,k)
         clwp(i,k) = max(tmp, 0.0) * delts(i)       ! g/m**2
#else
         clwp(i,k) = max(qci(i,k), 0.0) * delts(i)  ! g/m**2
#endif
#else
         clwp(i,k) = max(qci(i,k), 0.0) * delts(i)  ! g/m**2
#endif
       enddo
     enddo
   endif
! ------------------------------------------------------------------
!     find the stratosphere cut off layer for high cloud. it       
!      is assumed to be above the layer with dthdp less than      
!      -0.25 in the high cloud domain (from looking at 1 case).  
! ------------------------------------------------------------------ 
   do i = 1,imx2
     kcut(i) = levm2                                             
   end do
   do k = klow+1,levm2
     bit1 = .false.
     do i = 1,imx2
       if (kcut(i).eq.levm2 .and. prsly(i,k).le.ptop1(i,3) .and.               &
           dthdp(i,k).lt.-0.25e0) then
         kcut(i) = k
       end if
       bit1    = bit1 .or. kcut(i).eq.levm2
     end do
     if (.not. bit1) go to 85
   enddo
85 continue                                                     
!
! --- get convective cloud top and base level, set logical array bitx
!
   bit1 = .false.
   do i = 1,imx2
     bitx(i) = cv(i).gt.0.0 .and. cvt(i).ge.cvb(i)
     bit1 = bit1 .or. bitx(i)
   enddo
   if (bit1 .and. icwp.eq.0) then
     do i = 1,imx2
       kbt1(i) = 1
       kth1(i) = 1
       if (bitx(i)) then
         kbt1(i) = nint(cvb(i))
         kth1(i) = min( levm2, nint(cvt(i)) )
       endif
     enddo
!
! --- put convective cloud into 'cldcnv', no merge at this point..
!
     do k = 2,levm2
       do i = 1,imx2
         if (bitx(i) .and. kbt1(i).le.k .and. kth1(i).ge.k) then
           cv1 = cv(i)*(0.25+3.75*max(0.0,min(0.2,rhrh(i,k)-crhrh)))
           if (cv1 .ge. 0.001) cldcnv(i,k) = cv1
         endif
       enddo
     enddo
!
! --- if mean cvt layer higher than 400mb add anvil cirrus
!
     do i = 1,imx2
       if (bitx(i) .and. prsly(i,kth1(i)).le.cvtop) then
         kk = kth1(i)
         cldcnv(i,kk) = max(0.0, min(1.0, 2.0*(cv(i)-0.3) ))
       endif
     enddo
   endif
! ------------------------------------------------------------------ 
   if (istrat.le.0) then                                         
     if (icwp .eq. 1) then
       do k = klowb,levm2
         do i = 1,imx2
           onemrh = max(1.0e-10, 1.0-rhrh(i,k))
           value = max(0.0, min(50.0, 1000.0*clwp(i,k)/onemrh ))
           cldtot(i,k) = max(0.0, rhrh(i,k)*(1.0-exp(value)) )
         enddo
       enddo
     else
! ------------------------------------------------------------------
!      patterned after ccm-3 scheme     
!      calculate stratiform cloud and put into array 'cldtot'      
! ------------------------------------------------------------------ 
       xx = -(g_ * g_) / (0.00035 * rd)
       do k = klowb,levm2
         do i = 1,imx2
           if (prsly(i,k) .ge. ptop1(1,2)) then
!
! --- low cloud domain
!
             rhcr = rhl
!             if (nint(slmsk(i)) .eq. 1) rhcr = rhl - 0.10
           else
!
! --- mid-high cloud domain
!
             if (prsly(i,k) .ge. ptop1(1,3)) then
               crh = rhm
               if (nint(slmsk(i)) .eq. 1) crh = rhm - 0.10
             else
               crh = rhh
               if (nint(slmsk(i)) .eq. 1) crh = rhh - 0.10
             endif
             aa = xx*prsly(i,k)*dthdp(i,k) / (t(i,k)-theta(i,k))
             rhcr = 0.999 - (1.0-crh)*(1.0 - min(1.0, max(0.0, aa)))
           endif
           bb = max(0.0, (rhrh(i,k)-rhcr)/(1.0-rhcr) )
           cldtot(i,k) = min(1.0, bb**2)
         enddo
       enddo
! ------------------------------------------------------------------ 
!     special treatment on low clouds                               
! ------------------------------------------------------------------
       do i = 1,imx2 
         bitx(i) = .false.
         kbt1(i) = klowb - 1
       end do
!
       do k = klowb-1,llyr
         do i = 1,imx2
           if (.not.bitx(i)) then
             bitx(i) = nint(slmsk(i)).eq.0                                     &
               .and. prsly(i,k).ge.pstrt .and. dthdp(i,k).le.clapkc            &
               .and. rhrh(i,k+1).le.0.6  .and. rhrh(i,k+2).le.0.6
             kbt1(i) = k
           endif
         enddo
       enddo
!
! --- vertical velocity adjustment on low clouds                   
!
       dvvcld = vvcld(1) - vvcld(2)
       if (ivva.gt.0 .and. dvvcld.gt.1.0e-6) then
         do k = klowb,klow                                       
           do i = 1,imx2                                          
             if (.not.bitx(i)) then
               cr1 = min(1.0, max(0.0,                                         &
                     (vvcld(1) - 10.0*vvel(i,k)) / dvvcld )) 
               cldtot(i,k) = cldtot(i,k) * sqrt(cr1)
             endif 
           enddo
         enddo
       endif
!
! --- t inversion related stratus clouds                        
!
       if (invr .gt. 0) then
         bit1 = .false.
         do i = 1,imx2
           bit1 = bit1 .or. bitx(i) 
         end do
         if (bit1) then
!
! --- smooth transition for cloud when dthdp between clapse
!     and clapse+dclps (-0.05 and -0.06)
!
           do k = klowb,klow+2
             do i = 1,imx2
               if (bitx(i)) then
                 cfiltr = max(0.3,                                             &
                        1.0-((clapse-dthdp(i,kbt1(i))) / dclps))
!
! --- for t inversion type cloud, add flag value of 2.0        
!
                 cldtot(i,k) = min(0.95, cldtot(i,k)*cfiltr)
               endif
             enddo
           enddo
         endif
       endif
!
     endif                  ! end icwp block
   endif                    ! end istrat=0 block
! ------------------------------------------------------------------ 
   if (istrat.gt.0) then                                         
!
     if (icwp .eq. 1) then    ! liquid water based clouds
#ifdef Moorthi
       inversn = .false.
       kinver  = kdim
       do k = 2,kdim
         do i = 1,imx2
           if (prsly(i,k).gt.600.0 .and. (.not.inversn(i))) then
             tem = t(i,k+1) - t(i,k)
             if (tem .gt. 0.1 .and. t(i,k) .gt. 278.0) then
               inversn(i) = .true.
               kinver(i)  = k
             endif
           endif
         enddo
       enddo
       clwmin = 0.0e-6
#endif
       do k = 1,kdim
#ifdef Moorthi
         do i = 1,imx2
!           clwt(i) = 1.0e-6 * prsl(i,k)/prsi(i,1)
           clwt(i) = 2.0e-6 * prsl(i,k)/prsi(i,1)
         enddo
#endif
         do i = 1,imx2
#ifdef Moorthi
           if (clwp(i,k) .gt. clwt(i)) then
#endif
             onemrh = max(1.e-10, 1.-rhrh(i,k))
#ifdef Moorthi
!
! --- Moorthi
!
             tem = 1.0 / max(prsly(i,k)*0.001, 0.01)
             clwm= clwmin * tem
!
             tem1 = min(max(sqrt(sqrt(onemrh*qs(i,k))),0.0001),1.0)
             tem1 = 2000.0 / tem
             if (inversn(i) .and. k.le.kinver(i)) tem = tem * 5.0
             value = max(min(tem*(clwp(i,k)-clwm), 50.),0.)
#endif
!
! --- Xu and Randall 1996
!
             qsum = 0.0
             do ic = icloud_,nwmass_
               kk = k + levs_ * (ic-1)
               qsum = q(i,kk) + qsum
             enddo
#ifdef EXPLICIT_CLOUDINESS
#ifdef CPS_QSUM
             qsum = qsum + qcicps(i,k) + qrscps(i,k) 
#endif
#endif
             tem1 = min( max( sqrt( onemrh * qs(i,k)), 0.0001), 1.0)
             tem1 = 100.0 / tem1
             value = max( min( tem1 * (qsum), 50.0), 0.0)
!
! -- consider to pressure
!
!             value = max( min( tem1 * (qsum-clwm), 50.0), 0.0)
!
             tem2   = sqrt( sqrt( rhrh(i,k)) )
             cldtot(i,k) = max( tem2 * (1.0 - exp(-value)), 0.0)
#ifdef Moorthi
           endif
#endif
         enddo
       enddo
     else                           !   Old Campana Clouds!
!
! ------------------------------------------------------------------
!     calculate stratiform cloud and put into array 'cldtot' using 
!       the cloud-rel.humidity relationship from table look-up..where
!       tables obtained using k.mitchell frequency distribution tuning
!       (observations are daily means from us af rtneph).....k.a.c.
!       tables created without lowest 10 percent of atmos.....k.a.c.
!      (observations are synoptic using -6,+3 window from rtneph)
!       tables are created with lowest 10-percent-of-atmos, and are
!            now used..  25 october 1995 ... kac.
! ------------------------------------------------------------------ 
!  this loop to retrieve cloud from rh rewritten 950113 -mi         
!
!--- an updated table with mcld=4 instead of old version of mcld=3
!
       do klev = 2,levm2                                              
         do i = 1,imx2                                               
           kbase(i) = 0                                             
           bitx(i) = .false.                                       
         enddo                                                  
         do kc = mcld,1,-1
           do i = 1,imx2                                         
             if(prsly(i,klev).ge.ptop1(i,kc+1)) kbase(i) = kc     
           enddo  
         enddo   
!
         nx = 0                                                
         nhalf = (nbin+1)/2                                   
         do i = 1,imx2                                       
           if(kbase(i).le.0.or.klev.gt.kcut(i)) then      
             cldtot(i,klev) = 0.                           
           elseif(rhrh(i,klev).le.rhcld(i,1,kbase(i))) then      
             cldtot(i,klev) = 0.                                  
           elseif(rhrh(i,klev).ge.rhcld(i,nbin,kbase(i))) then 
             cldtot(i,klev) = 1.                                
           else                                              
             bitx(i) = .true.                                 
             ksave(i) = nhalf                                
             nx = nx+1                                      
           endif                                         
         enddo                                          
         do while(nx.gt.0)                              
           nhalf = (nhalf+1)/2                          
           do i = 1,imx2                               
             if(bitx(i)) then                       
               crk = rhrh(i,klev)                    
               cr1 = rhcld(i,ksave(i),kbase(i))     
               cr2 = rhcld(i,ksave(i)+1,kbase(i))  
               if(crk.le.cr1) then              
                 ksave(i) = max(ksave(i)-nhalf,1)          
               elseif(crk.gt.cr2) then                  
                 ksave(i) = min(ksave(i)+nhalf,nbin-1)   
               else
#ifndef INTERACTIVE_STRATUS
                 cldtot(i,klev) = 0.01*(ksave(i)+(crk-cr1)/(cr2-cr1))
#else
                 qsum = 0.
                 do ic = icloud_,nwmass_
                   kk = klev + levs_ * (ic-1)
                   qsum = q(i,kk) + qsum
                 enddo
                 onemrh = max(1.e-10,1.-rhrh(i,klev))
                 value = max(min(1000.*qsum/onemrh,50.),0.)
                 cldtot(i,klev) = max(rhrh(i,klev)*(1.-exp(-value)),0.)
#endif
                 bitx(i) = .false.                                    
                 nx = nx-1                                           
               endif                                              
             endif                                               
           enddo                                                
         enddo                                                 
       enddo                                                  
!
!     clean out not-suspected marine stratus regions...
!      cause tuning procedure not carried out down to lyr3 and we
!      get too much lo cloud if we don t clean it out..
!
      do i = 1,imx2
        bitm(i) = .true.
      enddo
      do k = klowb,llyr
        do i = 1,imx2
          if(bitm(i)) then
            bitm(i) = prsly(i,k).lt.pstrt                                      &
               .or. slmsk(i).gt.0.0 .or. dthdp(i,k).gt.clapkc.or.              &
                rhrh (i,k+1).gt.0.60e0.or.rhrh (i,k+2).gt.0.60e0
            kbase(i) = k
          endif
        enddo
      enddo
      do k = 1,llyr
        do i = 1,imx2
          if(bitm(i)) cldtot(i,k) = 0.0e0
        enddo
      enddo
!
! ------------------------------------------------------------------ 
!     special treatment on low clouds                               
! ------------------------------------------------------------------
       dvvcld = vvcld(1) - vvcld(2)                                 
!                                                                 
       do 950 k = 2,klow                                            
!                                                                 
         do i = 1,imx2                                        
           omeg(i) = 10.0e0 * vvel(i,k)                        
           cl1 (i) = 0.0e0                                    
         enddo                               
         if (ivva .le. 0) go to 920                        
! --- vertical velocity adjustment on low clouds         
         bit1 = .false.                                  
         do i = 1,imx2                                
           bitx(i) = prsly(i,k).ge.ptop1(i,2) .and. cldtot(i,k).gt.0.0e0
           bit1 = bit1 .or. bitx(i)                                    
         enddo                                                    
         if (.not. bit1) go to 920                                  
         if(k.gt.llyr) then                                        
           do i = 1,imx2                                        
             if (bitx(i)) then                                   
               if(omeg(i).ge.vvcld(1)) then                     
                 cldtot(i,k) = 0.0e0                           
               else if(omeg(i).gt.vvcld(2)) then              
                 cr1 = (vvcld(1) - omeg(i)) / dvvcld         
                 cldtot(i,k) = cldtot(i,k) * sqrt(cr1)     
               endif                                      
             endif                                       
           enddo                                     
         endif                                         
!
! --- t inversion related stratus clouds             
!
920      if (invr .lt. 1) go to 950                  
         if (k.gt.llyr) go to 950                   
         bit1 = .true.                             
         bitm(1:ILOTS)=.false.
         do i = 1,imx2                          
           bit1 = bit1 .and. bitm(i)             
         enddo                               
         if (bit1) go to 950                   
         do i = 1,imx2                      
           if (.not.bitm(i)) then             
             if (dthdp(i,kbase(i)).gt.clpse) then   
!---   smooth transition for cloud when dthdp between   
!           clapse and clapse+dclps  (-0.05 and -0.06) 
               cfiltr = 1.0e0 - ((clpse - dthdp(i,kbase(i))) / dclps) 
               cldtot(i,k) = cldtot(i,k)*cfiltr                      
             endif                                                
           endif                                               
         enddo                                             
950    continue                                              
! ------------------------------------------------------------------
     endif                        ! if for liquid water or campana
   endif                        ! istrat gt 0 
! --- reset cloud amount to 0 if less than 0.001
   do k = 1,kdim
     do i = 1,imx2
       if (cldtot(i,k) .lt. 0.001) cldtot(i,k) = 0.0
     end do
   enddo
!
! -------------------------------------------------------------------
!     separate clouds into 3 pressure domains (l,m,h).  within each  
!     of the domains, assume separated cloud layers are randomly    
!     overlapped and adjacent cloud layers are max/ran overlapped. 
!     vertical location of each type of cloud is determined by    
!     the thickest continuing cloud layers in the domain.        
! -------------------------------------------------------------------
! ..... vertically contiguous clouds are treated with ran/max overlap
! ... crtcl- criterion of applying random or max overlaping (0.2-0.5)
! ... adcr - adjustment over land/sea: 0 (more ran) - 0.5 (more max)
   crtcl = 0.3
   do k = 1,kdim
     do i = 1,imx2
       adcr(i,k) = 0.0
       if (nint(slmsk(i)) .ne. 1) then            ! sea/ice points
         if (prsly(i,k) .ge. ptop1(i,2)) then     ! low cloud domain
           adcr(i,k) = 0.2
         endif
       endif
     enddo
   enddo
!
! --- loop over 3 cloud domains (l,m,h)                             
!
   do l = 1,3                                                 
     do i = 1,imx2                                          
       cld (i,l) = 0.0e0                                     
       mtop(i,l) = 1                                        
       mbot(i,l) = 1                                       
       cl1 (i) = 0.0e0                                    
       cl2 (i) = 0.0e0                                   
       kbt1(i) = 1                                      
       kbt2(i) = 1                                     
       kth1(i) = 0                                    
       kth2(i) = 0                                   
       bitm(i) = cldtot(i,2).gt.0.0 .or. cldcnv(i,2).gt.0.0
       ityp1(i)= 0
     enddo                                       
!                                                     
     do 700 k = 2,levm2                             
       bit1 = .false.                            
       do i = 1,imx2                              
!
! --- bitx logical array for layer k
!
         bitx(i) = (prsly(i,k).ge.ptop1(i,l+1)) .and.                          &
                   (prsly(i,k).lt.ptop1(i,l)) .and. bitm(i)
         bit1 = bit1 .or. bitx(i)                     
!
! --- bitm logical array updated for layer k+1
!
         bitm(i) = cldtot(i,k+1).gt.0.0 .or. cldcnv(i,k+1).gt.0.0
       enddo
       if (.not. bit1) go to 700                     
       do i = 1,imx2                                  
         if (bitx(i)) then
!
! ---  physical cloud as seen by radiation..conv takes precedence
! ---  if both are presence, use max value
!
           cr1  = cldtot(i,k)
           cr2  = cldcnv(i,k)
           cr1  = max (cr1, cr2)
           ityp2(2) = 1                  
           if (cr2 .gt. 0.0) ityp2(i) = 2
           if(kth2(i).le.0) then        
!
! --- kth2 le 0 : 1st cld layer.              
!
             kbt2(i) = k                  
             kth2(i) = 1                 
             cl2 (i) = cr1
             ityp1(i) = ityp2(i)
           else                         
!
! --- kth2 gt 0 : consecutive cld layer.  
!
             sdif = 2.0*abs(cl2(i)-cr1)/(cl2(i)+cr1)-adcr(i,k)
             if (ityp1(i).ne.ityp2(i) .or. sdif.ge.crtcl) then
               kbt1(i) = int( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i))                &
                            / (cl1(i) + cl2(i)) )
               kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i))               &
                             / (cl1(i) + cl2(i)) ) + 1
               cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)
               kbt2(i) = k                                     
               kth2(i) = 1                                    
               cl2 (i) = cr1
               ityp1(i) = ityp2(i)
             else
               kth2(i) = kth2(i) + 1
               cl2 (i) = max(cl2(i), cr1)
             endif
           endif                                             
         endif                                              
       enddo
       bit2 = .false.                                      
!....  bity=true if next lyr=clear or we change cloud domains.. 
       do i = 1,imx2                                          
         bity(i) = bitx(i) .and. (.not. bitm(i)                                &
                            .or.  prsly(i,k+1).lt.ptop1(i,l+1) ) 
         bit2 = bit2 .or. bity(i)                               
       enddo
       if (.not. bit2) go to 700                               
! --- at the domain boundary or separated cld lyrs, random overlap. 
!     choose the thickest or the largest fraction amt as the cld   
!     layer in that domain                                        
       do i = 1,imx2                                            
         if (bity(i)) then                                   
           kbt1(i) = int( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i))                    &
                        / (cl1(i) + cl2(i)) )              
           kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i))                   &
                         / (cl1(i) + cl2(i)) ) + 1          
           cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)       
           kbt2(i) = 1                                    
           kth2(i) = 0                                   
           cl2 (i) = 0.0e0                              
         endif                                         
       enddo
!
700  continue                                         
! --- finish one domain, save effective clouds          
     do i = 1,imx2                                    
       cld(i,l) =  cl1(i)                          
       mtop(i,l) = max(kbt1(i), kbt1(i)+kth1(i)-1)  
       mbot(i,l) = kbt1(i)                         
     enddo
   enddo                               
!
!....  calculate total cloud from the multi-lyr cloud array  
!  .......in a manner as seen by the radiation code........ 
!      where, max/ran overlap is used for vertically adjacent cloud
!      layers. for convection any anvil is considered a separate  
!      randomly overlapped cloud..                               
!
!      cl1, ityp1 store previous layer cloud amount and type
!      cl2, ityp2 store current layer cloud amount and type
!      type = 0,1,2 for clear layer, stratiform and convective cloud.
!
   do i = 1,imx2                                                 
     cld(i,4) = 0.e0                                           
     ityp1(i) = 0                                             
     cl1(i) = 0.0e0                                          
   enddo
!
   do k = 1,kdim
     do i = 1,imx2
       cl2(i) = max(cldcnv(i,k), cldtot(i,k))
       ityp2(i) = 0
       if (cldtot(i,k) .gt. 0.0) ityp2(i) =1
       if (cldcnv(i,k) .gt. 0.0) ityp2(i) =2
       bitx(i) = ityp2(i) .gt. 0
       if (k .lt. kdim) then
         bity(i) = bitx(i) .and. cldtot(i,k+1).le.0.0                          &
                           .and. cldcnv(i,k+1).le.0.0
       else
         bity(i) = bitx(i)
       endif
     enddo
!
     do i = 1,imx2
       if (bitx(i)) then
         if (ityp1(i) .eq. 0) then
           cl1(i) = cl2(i)
           ityp1(i) = ityp2(i)
         else
           sdif = 2.0*abs(cl1(i)-cl2(i))/(cl1(i)+cl2(i))-adcr(i,k)
           if (ityp1(i).ne.ityp2(i) .or. sdif.ge.crtcl) then
             cld(i,4) = cld(i,4) + cl1(i) - cld(i,4)*cl1(i)
             cl1(i) = cl2(i)
             ityp1(i) = ityp2(i)
           else
             cl1(i) = max(cl1(i), cl2(i))
           endif
         endif
       endif
!
       if (bity(i)) then
         cld(i,4) = cld(i,4) + cl1(i) - cld(i,4)*cl1(i)
         cl1(i) = 0.0
         ityp1(i) = 0
       endif
     enddo
   enddo
!
   return                                                  
   end subroutine rad_cloudiness_nasa                                                   
!
!-------------------------------------------------------------------------------
