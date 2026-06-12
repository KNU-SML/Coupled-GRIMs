#include <define.h>
   subroutine rad_main_driver(lons2,                                           &
                     tgrs,qgrs,qcis,pgr,                                       &
#ifndef SWRMDC
                     paerr,ozonea,albdoa,                                      &
#else
#ifndef RMP
                     ozonea,albdoa,                                            &
#else
!    &                  tgrs,qgrs,qcis,pgr,paerr,ozonea,albdoa,
                     tgrs,qgrs,qcis,pgr,ozonea,albdoa,                         &
#endif
#endif
                     slmskr,coszro,coszdg,rlat,tsear,                          &
#ifdef NIM
                     si1d,sl1d,prsi_nim,prsl_nim,                              &
#ifdef NIM_DIAG
                     fnet,                                                     &
#endif
#endif
#if defined (RRTMGSW) || defined (RRTMGLW)
                     snoweq,vtype,                                             &
#endif
                     alvbr,alnbr,alvdr,alndr,                                  &
                     cldary,cldtot,cldcnv,cldsa,mtopa,mbota,                   &
#ifdef SWRMDC
                     rhrh,kprf,idxc,cmix,denn,                                 &
#endif
#ifndef RMP
                     lat,latco,sdec,solc,rsin1,rcos1,rcos2,                    &
#else
                     lat,sdec,solc,rsin1,rcos1,rcos2,                          &
#endif
                     raddt,dtlw,itimsw,itimlw,kalb,                            &
                     iswsrc,ibnd,ko3,icwp,icfc,                                &
                     swhr,hlwr,                                                &
#ifdef EXPLICIT_CLOUDINESS
                     qcis_cps,                                                 &
                     taucld,cldwp,cldip,                                       &
#endif
                     sfnswr,sfdlwr,tsflwr,                                     &
#ifdef VIC
                     sfdswr,sfuswr,                                            &
#endif
                     gdfvdr,gdfndr,gdfvbr,gdfnbr)
!-------------------------------------------------------------------------------
!
! subroutine: rad_main_driver
!
! program history log:
!   1988-05-06  kenneth campana        initial development
!   1998-05-06  yu-tai hou             operational implementation
!   1999-09-01  song-you hong          add ice effects and prognostic clouds
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!     ************************************************************
!     *  added accumulation of clds and convective cloud in dg3  *
!     *                                        k.a.c sept 1994   *
!     *  f3d added for clouds..mi new code=f94/source2/diagnew   *
!     *                               b katz + k.a.c oct  1994   *
!     *  changed to store total cloud and all lyrs        *
!     *     of cloud......                                       *
!     *                                         k.a.c. nov94     *
!     *  fix pl1 for operations, where dgz is on and dg3 is off, *
!     *     ....note dg is on if either dgz or dg3 is on         *
!     *                                         k.a.c. jan94     *
!     ************************************************************
!
!        updates made to fix the h2d,h3d files...kac aug 90...
!        updates made to rad_driver - call wrth2d before wrtrad (so ctop ok)
!                     to rad_driver - send work array to wrth3d
!                     to wrth3d - to write proper layers of heat..
!                                 (in wrtrad)
!        updates made to add grid point diagnostics ..k.a.c...sep 91
!                     to rad_driver -
!        updates made to pass and receive sib data  ..k.a.c...mar 92
!                     to rad_driver -
!        updates made to fix sw rad diagnostics     ..k.a.c...jun 92
!                     proper diurnal weighting
!                     to rad_driver and rad_cos_zenith
!        updates made to calculate clear-sky "on-the-fly" kac aug 92
!                     to rad_driver,radfs,fst,spa,lwr,swr
!                     ...for cloud forcing....
!        updates made to call cld optical property routine (rad_cloud_comp),
!                     to give cld emissivity, optical depth, layer
!                     reflectance and transmitance
!                     to rad_driver and radfs ...y.h.           ...feb93
!
!        updates made to change definition of h,m,l domains..
!          to rad_cloud_comp   k.a.c...dec92 + aug93
!
!          wmo aerosols distributions, and b.p.breiglebs surface
!          albedo scheme.            ..............y.h...sep95
!             to rad_main_driver - call rad_main_solver 
!                                    ... to replace the old radfs
!                                        use gfdls lw and chous sw
!                         kalb=0 use old sfc albedo scheme (matthews)
!                             =1 use new sfc albedo scheme (breigleb)
!                         iswsrc flags for selections of sw absorbers
!                                1:aerosol, 2:o2, 3:co2, 4:w.v. 5:o3
!                                 =1: with the absorber, =0: without
#ifdef SWRMDC
!        updates made to opac aerosol algorithm (1998)...y.h..sep99
#endif
!
!-------------------------------------------------------------------------------
#if defined (RRTMGSW) || defined (RRTMGLW)
   use paramodel, only : LONF2S,LATG2S,slvark_,levs_,levh_,levp1_
#else
   use paramodel, only : LONF2S,LATG2S,slvark_,levs_,levp1_
#endif
   use constant, only  : cal_,akapa_
#ifdef DFS
   use dfsvar, only    : si=>sigmafull,sl=>sigma,iope
#else
#ifdef NIM
   use comio, only     : iope,ncpu
#else
   use comio, only     : iope
#endif
#endif
#ifndef RMP
   use comfver
   use radiag
   use comgpd
#else
#ifdef MP
   use paramodel, only : igrd12p_,jgrd12p_
   use commpi
#else
   use paramodel, only : igrd12_
#endif
   use rscomver
   use rscomrad
   use rscomgpd
#endif
#ifndef SWRMDC
#ifndef ICECLOUD
   use rdparm
#else
   use rdparm8
#endif
#else
   use aerparm
   use rdparm99
#endif
#if defined (DG) || defined (DG3)
   use comfcst, only : dfvbr, dfnbr, dfvdr, dfndr
#endif
#ifdef SMP_NUDGING_RAD
   use comfcst, only : psfc,tsfc,qsfc,prec,scmlwup,scmlwdn,scmswup,scmswdn
#endif
#ifdef CHEM
   use dao_mod, only : c_albd=>albd, c_cldfrc=>cldfrc, c_radswg=>radswg
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
#define LONLENS igrd12_
#endif
#endif
!soojin_couple
#ifdef AOMG
  use comfphys, only : cpl_sswdn, cpl_sswup
#endif /*AOMG*/
!-------------------------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------------------------
   real,parameter       ::  cnwatt=-cal_*1.e4/60.,cnprog=1./cnwatt
   real                 ::  tgrs(LONF2S,levs_)
#if defined (RRTMGSW) || defined (RRTMGLW)
   real                 ::  qgrs(LONF2S,levh_)
#else
   real                 ::  qgrs(LONF2S,levs_)
#endif
   real                 ::  pgr(LONF2S)
   real                 ::  qcis(LONF2S,levs_)
   real                 ::  coszer(LONF2S),coszdg(LONF2S)
   real                 ::  rlat(LONF2S),slmskr(LONF2S)
#if defined (RRTMGSW) || defined (RRTMGLW)
   real                 ::  snoweq(LONF2S),vtype(LONF2S)
   real                 ::  fnp(levs_), fnm(levs_), t8w(LONF2S,levp1_), slm(levs_)
#endif
#ifndef SWRMDC
   real                 ::  cldary(LONF2S,levs_)
   real                 ::  cldtot(LONF2S,levs_)
   real                 ::  cldcnv(LONF2S,levs_)
#else
   real                 ::  cldary(LONF2S,levs_)
   real                 ::  rhrh (LONF2S,levs_)
   real                 ::  cldtot(LONF2S,levs_)
   real                 ::  cldcnv(LONF2S,levs_)
   integer              ::  idxc(nxc,LONF2S)
   real                 ::  cmix(nxc,LONF2S)
   real                 ::  denn(ndn,LONF2S)
   integer              ::  kprf(LONF2S)
#endif
   real                 ::  flwup(LONF2S),fswup(LONF2S)
   real                 ::  fswdn(LONF2S)
   real                 ::  sswup(LONF2S),sswdn(LONF2S)
   real                 ::  slwup(LONF2S),slwdn(LONF2S)
!
   real                 ::  prsi(LONF2S,levs_+1)
   real                 ::  prsl(LONF2S,levs_)
   real                 ::  prsik(LONF2S,levs_+1)
   real                 ::  prslk(LONF2S,levs_)
#ifdef CLR
   real                 ::  flwup0(LONF2S),fswup0(LONF2S)
   real                 ::  sswup0(LONF2S),sswdn0(LONF2S)
   real                 ::  slwdn0(LONF2S)
#endif
#ifndef SWRMDC
   real                 ::  paerr(LONF2S,5)
#endif
   real                 ::  ozonea(LONF2S,levs_)
   real                 ::  albdoa(LONF2S)
   real                 ::  cldsa(LONF2S,4)
   integer              ::  mtopa(LONF2S,3)
   integer              ::  mbota(LONF2S,3)
   real                 ::  coszro(LONF2S)
   real                 ::  tsear(LONF2S)
   real                 ::  swhr(LONF2S,levs_)
   real                 ::  hlwr(LONF2S,levs_)
   real                 ::  sfnswr(LONF2S),sfdlwr(LONF2S)
   real                 ::  tsflwr(LONF2S)
#ifdef VIC
   real                 ::  sfdswr(LONF2S),sfuswr(LONF2S)
#endif
   real                 ::  alvbr(LONF2S),alnbr(LONF2S)
   real                 ::  alvdr(LONF2S),alndr(LONF2S)
   real                 ::  gdfvbr(LONF2S),gdfnbr(LONF2S)
   real                 ::  gdfvdr(LONF2S),gdfndr(LONF2S)
#ifdef EXPLICIT_CLOUDINESS
   real                 ::  qcis_cps(LONF2S,levs_)
   real                 ::  taucld(LONF2S,levs_)
   real                 ::  cldwp(LONF2S,levs_),cldip(LONF2S,levs_)
#endif
#ifdef HYBRID
   real                 ::  sihyb(LONF2S,levp1_),cihyb(LONF2S,levp1_)
   real                 ::  slhyb(LONF2S,levs_),clhyb(LONF2S,levs_)
   real                 ::  delhyb(LONF2S,levs_),slkhyb(LONF2S,levs_)
#endif
#ifdef NIM
   real                 ::  si1d(levp1_),sl1d(levs_)
   real                 ::  prsi_nim(LONF2S,levp1_), prsl_nim(LONF2S,levs_)
#ifdef NIM_DIAG
   real                 ::  fnet(LONF2S)
#endif
#endif
   real                 ::  dtlw,raddt,rcos2,rcos1,rsin1,solc,sdec,work4
#if defined (RRTMGSW) || defined (RRTMGLW)
   real                 ::  w1,w2
#endif
   integer              ::  icfc,icwp,ibnd,itimlw,itimsw,i,k,itop,ibtc
   integer              ::  ko3,kalb,latco,lat,lons2
!
!   downward sw fluxes from sw sib rad
!   saved for h2d file
!
   integer              ::  iswsrc(nsrc)
#ifdef DBG
   integer              ::  latd
#endif
!-------------------------------------------------------------------------------
!
! initialize local vairalbes
!
   flwup=0. ; fswup=0. ; fswdn=0. ; sswup=0. ; sswdn=0. ; slwup=0. ; slwdn=0.
   prsi=0.  ; prsl=0.  ; prsik=0. ; prslk=0.
#ifdef CLR
   flwup0=0.; fswup0=0.; sswup0=0.; sswdn0=0.; slwdn0=0.
#endif
#ifdef HYBRID
   sihyb=0. ; cihyb=0. ; slhyb=0. ; clhyb=0. ; delhyb=0.; slkhyb=0.
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
   enddo
   do k = 1,levp1_
     si(k)=si1d(k)
   enddo
   do i = 1,LONLENS
     do k = 1,levs_
       prsl(i,k)   = prsl_nim(i,k)
       prsi(i,k)   = prsi_nim(i,k)
     enddo
     prsi(i,levs_+1) = prsi_nim(i,levs_+1)
   enddo
#else
   do i = 1,LONLENS
     do k = 1,levs_
       prsl(i,k)   = pgr(i)*sl(k)
       prsi(i,k)   = pgr(i)*si(k)
     enddo
     prsi(i,levs_+1) = pgr(i)*si(levs_+1)
   enddo
#endif /* NIM */
#endif /* HYBRID */
!
#if defined (RRTMGSW) || defined (RRTMGLW)
! calculate interface level temperature
!
   do k = 2,levs_
     slm(k)=0.5*(sl(k)+sl(k-1))
     fnp(k)=0.5*sl(k)/slm(k)
     fnm(k)=0.5*sl(k-1)/slm(k)
   enddo
!
   do i = 1,LONLENS
     do k = 2,levs_
       t8w(i,k)=fnm(k)*tgrs(i,k)+fnp(k)*tgrs(i,k-1)
     enddo
   enddo
   do i = 1,LONLENS
     w1=(prsi(i,1)-prsl(i,2))/(prsl(i,1)-prsl(i,2))
     w2= 1. -w1
     t8w(i,1)=w1*tgrs(i,1)+w2*tgrs(i,2)
   enddo
   do i = 1,LONLENS
     w1=(prsi(i,levs_)-prsl(i,levs_-1))/(prsl(i,levs_)-prsl(i,levs_-1))
     w2= 1. -w1
     t8w(i,levp1_)=w1*tgrs(i,levs_)+w2*tgrs(i,levs_-1)
!     t8w(i,levp1_)=tgrs(i,levs_)
   enddo
#endif
#ifdef DBG
#ifdef NIM
   latd=lat
   if(iope) then
#else
   latd=10
   if(lat.eq.latd .and. iope) then
#endif
     write(6,*) 'before rad_main_solver '
     write(6,*) si,sl,lat,rlat(1),solc,rsin1,rcos1,rcos2,slmskr(latd),         &
              qgrs(latd,10),tgrs(latd,10),ozonea(latd,10),tsear(latd),         &
              coszro(latd),cldary(latd,10),albdoa(latd),                       &
              alvbr(latd),alnbr(latd),alvdr(latd),alndr(latd)
     write(6,*) 'itime ',itimsw,itimlw,ko3,kalb,iswsrc,                        &
              hlwr(latd,10),slwup(latd),slwdn(latd),flwup(latd),               &
              swhr(latd,10),sswup(latd),sswdn(latd),fswup(latd),fswdn(latd)
   endif
#endif
!
   call rad_main_solver(LONLENS,LONF2S,prsi,prsl,                              &
               lat,rlat(1),solc,rsin1,rcos1,rcos2,slmskr(1),                   &
               qgrs(1,1),tgrs(1,1),ozonea(1,1),tsear(1),                       &
               coszro(1),cldary(1,1),cldtot(1,1),cldcnv(1,1),                  &
#if defined (RRTMGSW) || defined (RRTMGLW)
               snoweq(1),vtype(1),t8w(1,1),                                    &
#endif
               albdoa(1),                                                      &
#ifndef SWRMDC
               alvbr(1),alnbr(1),alvdr(1),alndr(1),paerr(1,1),                 &
#else
               alvbr(1),alnbr(1),alvdr(1),alndr(1),                            &
               kprf(1),idxc(1,1),cmix(1,1),denn(1,1),rhrh(1,1),                &
#endif
               itimsw,itimlw,ko3,kalb,iswsrc,                                  &
               ibnd,icwp,qcis,icfc,                                            &
               hlwr(1,1),slwup,slwdn,flwup,                                    &
               swhr(1,1),sswup,sswdn,fswup,fswdn,                              &
#ifdef NIM_DIAG
               fnet,                                                           &
#endif
#ifdef CLR
               flwup0,fswup0,sswdn0,sswup0,slwdn0,                             &
#endif
#ifdef EXPLICIT_CLOUDINESS
               qcis_cps,                                                       &
               taucld,cldwp,cldip,                                             &
#endif
               gdfvbr(1),gdfnbr(1),gdfvdr(1),gdfndr(1))
#ifdef DBG
#ifdef NIM
   if(iope) then
#else
   if(lat.eq.latd .and. iope) then
#endif
     write(6,*) 'after rad_main_solver '
     write(6,*) si,sl,lat,rlat(1),solc,rsin1,rcos1,rcos2,slmskr(latd),         &
              qgrs(latd,10),tgrs(latd,10),ozonea(latd,10),tsear(latd),         &
              coszro(latd),cldary(latd,10),albdoa(latd),                       &
              alvbr(latd),alnbr(latd),alvdr(latd),alndr(latd) 
     write(6,*) 'itime ',itimsw,itimlw,ko3,kalb,iswsrc,                        &
              hlwr(latd,10),slwup(latd),slwdn(latd),flwup(latd),               &
              swhr(latd,10),sswup(latd),sswdn(latd),fswup(latd),fswdn(latd)
   endif
#endif
!
!  cnprog is conversion from w/m**2 to phys_lsm_osu1 units
!
   if(itimsw.eq.1) then
     do i = 1,LONLENS
#ifdef SMP_NUDGING_RAD
       print*,'rad_main_driver scmswup scmswdn',scmswup,scmswdn
       print*,'rad_main_driver scmlwup scmlwdn',scmlwup,scmlwdn
       sswdn(i) = scmswdn(i)
#endif
       sfnswr(i) = (sswdn(i)-sswup(i))*cnprog
#ifdef VIC
       sfdswr(i) = sswdn(i)*cnprog
       sfuswr(i) = sswup(i)*cnprog
#endif
     enddo
   endif
!
   if(itimlw.eq.1) then
     do i = 1,LONLENS
       sfdlwr(i) = slwdn(i)*cnprog
       tsflwr(i) = tgrs(i,1)
     enddo
   endif
!
!   save 4 components of downward sw flux
!
#ifdef DG
   if(itimsw.eq.1) then
!
!  accumulate for h2d file..
!
     do i=1,LONLENS
       dfvbr(i,lat) = dfvbr(i,lat) + raddt * gdfvbr(i)
       dfnbr(i,lat) = dfnbr(i,lat) + raddt * gdfnbr(i)
       dfvdr(i,lat) = dfvdr(i,lat) + raddt * gdfvdr(i)
       dfndr(i,lat) = dfndr(i,lat) + raddt * gdfndr(i)
     enddo
   endif
#endif
!
!  grid point monitor-data on radiation grid
!
#ifdef DGP
   if(isave.ne.0.and.npoint.gt.0) then
     do igpt = 1,npoint
       if (lat.eq.jgrd(igpt)) then
         do id = 1,LONLENS
           if (id.eq.igrd(igpt)) then
             svdata( 25,igpt,itnum)= cldsa(id,3)
             svdata( 26,igpt,itnum)= cldsa(id,2)
             svdata( 27,igpt,itnum)= cldsa(id,1)
             svdata( 38,igpt,itnum)= cldsa(id,4)
             if(isshrt.lt.1.and.itimsw.eq.1) then
               svdata( 41,igpt,itnum)= id
               svdata( 42,igpt,itnum)= lat
               svdata( 43,igpt,itnum)= slmskr(id)
               svdata( 44,igpt,itnum)= tsear (id)
               svdata( 45,igpt,itnum)= sswdn(id)
               svdata( 46,igpt,itnum)= sswup(id)
               svdata( 48,igpt,itnum)= fswdn(id)
               svdata( 49,igpt,itnum)= fswup(id)
               svdata( 51,igpt,itnum)= mtopa(id,3)
               svdata( 52,igpt,itnum)= mtopa(id,2)
               svdata( 53,igpt,itnum)= mtopa(id,1)
               svdata( 54,igpt,itnum)= mbota(id,3)
               svdata( 55,igpt,itnum)= mbota(id,2)
               svdata( 56,igpt,itnum)= mbota(id,1)
               svdata( 57,igpt,itnum)= coszro(id)
               svdata( 58,igpt,itnum)= asin(sdec)*180.e0/3.14159265e0
#ifdef CLR
               svdata( 70,igpt,itnum)= flwup0(id)
               svdata( 71,igpt,itnum)= slwdn0(id)
               if (coszro(id).gt.0.) then
                 svdata( 72,igpt,itnum)= fswup0(id)*coszdg(id)/coszro(id)
                 svdata( 73,igpt,itnum)= sswdn0(id)*coszdg(id)/coszro(id)
                 svdata( 74,igpt,itnum)= sswup0(id)*coszdg(id)/coszro(id)
               else
                 svdata( 72,igpt,itnum)=0.
                 svdata( 73,igpt,itnum)=0.
                 svdata( 74,igpt,itnum)=0.
               endif
#endif
             endif
             if(isshrt.lt.1.and.itimlw.eq.1) then
                svdata( 47,igpt,itnum)= slwdn(id)
                svdata( 50,igpt,itnum)= flwup(id)
             endif
             if (ilshrt.lt.1.and.itimsw.eq.1) then
               do kc = 1,levs_
#ifndef SWRMDC
                 cvcl = float(int(cldary(id,kc))/10)*1.e-3
                 if (cvcl.gt.0.0e0) then
                   svdata(kc+slvark_+7*levs_,igpt,itnum) = cvcl
                 else
                   svdata(kc+slvark_+7*levs_,igpt,itnum) =                     &
                   amod(cldary(id,kc),2.e0)
                 end if
#else
                 svdata(kc+slvark_+7*levs_,igpt,itnum) =                       &
                 max(cldcnv(id,kc), cldtot(id,kc))
#endif
               enddo
             endif
           endif
         enddo
       endif
     enddo
   endif
#endif /* DGP end */
!
! save gridded radiative fluxes and cld data
!
   if(itimsw.eq.1) then
     do i = 1,LONLENS
       if(kalb.eq.0) then
         fluxr(i,lat,17) = fluxr(i,lat,17) + raddt * albdoa(i)
       else
         fluxr(i,lat,17) = fluxr(i,lat,17) + raddt *                           &
              (0.466e0 * alvdr(i) + 0.529e0 * alndr(i))
       endif
       fluxr(i,lat,26) = fluxr(i,lat,26) + raddt * cldsa(i,4)
     enddo
   endif
   if(itimlw.eq.1) then
     do i = 1, LONLENS
       fluxr(i,lat,1 ) = fluxr(i,lat,1 ) + dtlw * flwup(i)
       fluxr(i,lat,19) = fluxr(i,lat,19) + dtlw * slwdn(i)
       fluxr(i,lat,20) = fluxr(i,lat,20) + dtlw * slwup(i)
#ifdef CLR
       fluxr(i,lat,21) = fluxr(i,lat,21) + dtlw * flwup0(i)
       fluxr(i,lat,25) = fluxr(i,lat,25) + dtlw * slwdn0(i)
#endif
     enddo
   endif
!
!   proper diurnal sw wgt..coszro=mean cosz over daylight, while
!                           coszdg= mean cosz over entire interval
!
   do i = 1,LONLENS
#ifndef NIM
     if (coszro(i).gt.0.) then
       fluxr(i,lat,2 ) = fluxr(i,lat,2 ) + raddt * fswup(i)                    &
                                     * coszdg(i)/coszro(i)
       fluxr(i,lat,3 ) = fluxr(i,lat,3 ) + raddt * sswup(i)                    &
                                     * coszdg(i)/coszro(i)
       fluxr(i,lat,4 ) = fluxr(i,lat,4 ) + raddt * sswdn(i)                    &
                                     * coszdg(i)/coszro(i)
       fluxr(i,lat,18) = fluxr(i,lat,18) + raddt * fswdn(i)                    &
                                     * coszdg(i)/coszro(i)
#ifdef CLR
       fluxr(i,lat,22) = fluxr(i,lat,22) + raddt * fswup0(i)                   &
                                     * coszdg(i)/coszro(i)
       fluxr(i,lat,23) = fluxr(i,lat,23) + raddt * sswdn0(i)                   &
                                     * coszdg(i)/coszro(i)
       fluxr(i,lat,24) = fluxr(i,lat,24) + raddt * sswup0(i)                   &
                                     * coszdg(i)/coszro(i)
#endif
     end if
#else /* NIM */
     fluxr(i,lat,2 ) = fluxr(i,lat,2 ) + raddt * fswup(i)
     fluxr(i,lat,3 ) = fluxr(i,lat,3 ) + raddt * sswup(i)
     fluxr(i,lat,4 ) = fluxr(i,lat,4 ) + raddt * sswdn(i)
     fluxr(i,lat,18) = fluxr(i,lat,18) + raddt * fswdn(i)
#ifdef CLR
     fluxr(i,lat,22) = fluxr(i,lat,22) + raddt * fswup0(i)
     fluxr(i,lat,23) = fluxr(i,lat,23) + raddt * sswdn0(i)
     fluxr(i,lat,24) = fluxr(i,lat,24) + raddt * sswup0(i)
#endif
#endif
   enddo

#ifdef AOMG
!soojin_coupled
!radiation exchage variable
#endif /*AOMG*/

!
!   save cld frac,toplyr,botlyr and top temp
!   note that order of high, middle and low clouds is
!   reversed for proper output to sflux and h2d files.
!
   do k = 1,3
     do i = 1,LONLENS
       fluxr(i,lat,8-k) = fluxr(i,lat,8-k) + raddt * cldsa(i,k)
!
!  save interface pressure (cb) of top/bot,
!
       itop = mtopa(i,k)
       ibtc = mbota(i,k)
       fluxr(i,lat,11-k) = fluxr(i,lat,11-k) + raddt *                         &
                          prsi(i,itop+1) * cldsa(i,k)
       fluxr(i,lat,14-k) = fluxr(i,lat,14-k) + raddt *                         &
                          prsi(i,ibtc) * cldsa(i,k)
       fluxr(i,lat,17-k) = fluxr(i,lat,17-k) + raddt *                         &
                          tgrs(i,itop) * cldsa(i,k)
     enddo
   enddo
!
#ifdef CHEM
   ! save 2-D variables for chemistry
   do i = 1,LONLENS
     if(kalb.eq.0) then
       c_albd(i,lat) = albdoa(i)
     else
       c_albd(i,lat) = 0.466e0 * alvdr(i) + 0.529e0 * alndr(i)
     endif
     c_cldfrc(i,lat) = cldsa(i,4)
     c_radswg(i,lat) = sswdn(i) * ( 1d0 - c_albd(i,lat) )
   enddo
#endif
!
   return
   end subroutine rad_main_driver
!
!-------------------------------------------------------------------------------
