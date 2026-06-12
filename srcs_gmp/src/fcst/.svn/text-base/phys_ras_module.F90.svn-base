#include <define.h>
   module phys_ras_module
#ifdef RAS
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [phys_cps_ras]
!      |
!      |---[phys_ras_module]
!               |
!               |--- [ras_driver] * --- [ras_solver] * --- [ras_cloud] *
!                                                              |
!                                                       [ras_cwork_crit] *
!    [dyn_sph_driver] or [physics_main_driver]                 |
!      |                                               [ras_cloud_factor] *
!      |---[phys_ras_module]
!               |
!               |--- [ras_setup] * --- [ras_setup_sub] *
!                                        |- [ras_setup_param] *
!
! program history log:
!   1995-01-01  moorthi                development
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ras_driver(len, lenc, km, lm, nstrp, dt,                         &
               prsi, prsl,                                                     &
               krmin, krmax, ncrnd, afac, rannum, ufac,                        &
               rgas, cp, grav, alhl,                                           &
               tin, qin, uin, vin,                                             &
               evapc, kbot, ktop, icps, lat, cd,                               &
               sig, prj, prh, fpk, hpk, sgb, ods, rasal, prns)
!-------------------------------------------------------------------------------
!
! subroutine: ras_driver
!
! program history log:
!   1991-01-01  moorthi                initial mrf
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real,parameter       ::  frac=5.0, crtmsf=0.0
!
!--  input
!
   real                 ::  dt, rgas, cp, grav, alhl, afac, ufac, rkap 
   integer              ::  i, k, l, ll, len, lenc, km, lm, nstrp, lat,        &
                            krmin, krmax, ncrnd, kmp1, ir1, lmm1, lmp1, ia
   real                 ::  tin(lenc,km), qin(lenc,km), uin(lenc,km)
   real                 ::  vin(lenc,km), cd(len)
   real                 ::  evapc(len)
   real                 ::  prsi(lenc,km+1),prsl(lenc,km)
   integer              ::  kbot(len), ktop(len), icps(len)
   real                 ::  rannum(ncrnd*3)
!
   real                 ::  tem, tem1, tem2
   real                 ::  sig(lm+1),  prj(lm+1), prh(lm),  fpk(lm), hpk(lm)
   real                 ::  sgb(lm),  ods(lm), rasal(lm), prns(nstrp)
!
!--  locals
!
   real, allocatable:: pib(:), pik(:), pki(:), pcu(:)                          &
            ,          tns(:), qns(:), uos(:), vos(:)                          &
            ,          gam(:,:)                                                &
            , uoc(:,:), voc(:,:), tn0(:,:), qn0(:,:)                           &
            , tcu(:,:), qcu(:,:), clt(:,:)
!
   allocate ( pib(len),    pik(len),    pki(len),    pcu(len)                  &
            , tns(len),    qns(len),    uos(len),    vos(len)                  &
            , gam(len,lm)                                                      &
            , uoc(len,lm), voc(len,lm), tn0(len,lm), qn0(len,lm)               &
            , tcu(len,lm), qcu(len,lm), clt(len,lm))
!-------------------------------------------------------------------------------
   kmp1 = km + 1
   ir1  = 5 * lm + 1
   lmm1 = lm - 1
   lmp1 = lm + 1
   rkap = rgas / cp
   ia = 1
!
   do i = 1,len
      pib(i) = prsi(i,1) * 10.0
      pik(i) = pib(i) ** rkap
      pki(i) = 1.0 / pik(i)
      uos(i) = 0.0
      vos(i) = 0.0
      tns(i) = 0.0
      qns(i) = 0.0
   enddo
   do l = 1,lmm1
      ll = km - l + 1
      do i = 1,len
         tn0(i,l) = tin(i,ll) * pki(i) * (1.0/prh(l))
         qn0(i,l) = qin(i,ll)
         uoc(i,l) = uin(i,ll)
         voc(i,l) = vin(i,ll)
      enddo
   enddo
   tem2 = 1.0 / (sig(lm+1) - sig(lm))
   do l = 1,nstrp
      do i = 1,len
         tem  = (prsi(i,l) - prsi(i,l+1))/prsi(i,1) * tem2
         tem1 = tem * (1.0/prns(l))
         tns(i) = tns(i) + tin(i,l) * tem1
         qns(i) = qns(i) + qin(i,l) * tem
         uos(i) = uos(i) + uin(i,l) * tem
         vos(i) = vos(i) + vin(i,l) * tem
      enddo
   enddo
   do i = 1,len
      tns(i)    = tns(i) * pki(i)
      tn0(i,lm) = tns(i)
      qn0(i,lm) = qns(i)
      uoc(i,lm) = uos(i)
      voc(i,lm) = vos(i)
      pcu(i)    = 0.0
   enddo
   do l = 1,lm
      do i = 1,len
         tcu(i,l) = 0.0
         qcu(i,l) = 0.0
         clt(i,l) = 0.0
      enddo
   enddo
!
   call ras_solver (len, len, lm,  dt                                          &
             ,ncrnd, krmin, krmax, frac, rasal, .true.                         &
             ,cp,  alhl, grav, crtmsf, rannum                                  &
             ,sgb,    ods, prh, prj, fpk, hpk, sig(2)                          &
             ,pib,    pik, pki, gam                                            &
             ,uoc,    voc, tn0, qn0                                            &
             ,pcu,    clt, tcu, qcu                                            &
             ,cd)
!
   tem = dt / 86400000.0
   do i = 1,len
     evapc(i) = pcu(i) * tem
   enddo
!
   do l = 1,lmm1
      ll = km - l + 1
      do i = 1,len
         tin(i,ll) = tn0(i,l) * pik(i) * prh(l)
         qin(i,ll) = qn0(i,l)
      enddo
   enddo
   do i = 1,len
      icps(i) = 0
      kbot(i) = kmp1
      ktop(i) = 0
      tns(i) = (tn0(i,lm) - tns(i)) * pik(i) * prh(lm)
      qns(i) = qn0(i,lm) - qns(i)
      uos(i) = uoc(i,lm) - uos(i)
      vos(i) = voc(i,lm) - vos(i)
   enddo

   do l = lmm1,1,-1
      do i = 1,len
         if (tcu(i,l) .ne. 0.0) then
            icps(i) = 1
         endif
      enddo
   enddo
!
!  new test for convective clouds ! added ib 07/25/96
!
   do l = lmm1,1,-1
      do i = 1,len
         if (clt(i,l) .gt. 0) ktop(i) = l
      enddo
   enddo
   do l = 1,lmm1
      do i = 1,len
         if (clt(i,l) .gt. 0) kbot(i) = l
      enddo
   enddo
   do i = 1,len
      ktop(i) = kmp1 - ktop(i)
      kbot(i) = kmp1 - kbot(i)
   enddo
!
   do l = 1,nstrp
      do i = 1,len
         tin(i,l) = tin(i,l) + tns(i)
         qin(i,l) = qin(i,l) + qns(i)
      enddo
   enddo
!
!  calling ras is over
!
   deallocate ( pib,    pik,    pki,    pcu                                    &
             ,tns,    qns,    uos,    vos                                      &
             ,gam                                                              &
             ,uoc, voc, tn0, qn0                                               &
             ,tcu, qcu, clt)
!
   return
   end subroutine ras_driver
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine ras_solver (len, lenc,   k,   dt                                 &
                   ,ncrnd, krmin, krmax, frac, rasal, botop                    &
                   ,cp,  alhl, grav, crmsfx, rannum                            &
                   ,sgo, onbdsg, prh, prj, fpk, hpk, sge                       &
                   ,pii, pik,    pki,  gam                                     &
                   ,uoi, voi,    poi, qoi                                      &
                   ,pcu, cln,    q1,  q2                                       &
                   ,cd)
!-------------------------------------------------------------------------------
!
! subroutine: ras_solver
!
!********************************************************************!
!********************************************************************!
!******************** relaxed arakawa-schubert **********************!
!************************ parameterization **************************!
!********************* plug compatible version  *********************!
!*********************  for sigma=p/ps models  **********************!
!************************* 04 august 1992 ***************************!
!********************************************************************!
!************************** developed by ****************************!
!**************************              ****************************!
!************************ shrinivas moorthi *************************!
!************************       and         *************************!
!************************  max j. suarez ****************************!
!************************                ****************************!
!******************** laboratory for atmospheres ********************!
!****************** nasa/gsfc, greenbelt, md 20771 ******************!
!********************************************************************!
!********************************************************************!
!
! program history log:
!   1991-07-01  moorthi          mrf implementation
!
! references : moorthi and suarez (1992, mwr)
!
!------------------------------------------------------------------------------e
   real,parameter       ::  one=1.0, daylen=86400.0, rknob=10.0
   integer,parameter    ::  ibmax=1,icm=100
!
   real                 ::  pii(len),   pik(len),   pki(len), rasal(k-1)
   real                 ::  uoi(len,k), voi(len,k), poi(len,k), qoi(len,k)
   real                 ::  cd(len)
!
   real                 ::  q1(len,k), q2(len,k),  cln(len,k), pcu(len)
!
   real                 ::  gam(len,k)
   real                 ::  rannum(ncrnd*3)
!
   real                 ::  sgo(k),    prj(k),    prh(k),    fpk(k),    hpk(k)
   real                 ::  sge(k),    eps(icm),  del(icm),  sdu(icm),  sdd(icm)
   real                 ::  epss(icm)
   real                 ::  onbdsg(k)
   integer              ::  ic(icm)
   real                 ::  rhfacl(ibmax), rhfacs(ibmax)
   logical              ::  botop
!
   km1   = k  - 1
   albycp = alhl / cp
!
!-clm rhfacs(1) = 0.85
!
   rhfacs(1) = 0.70
!     rhfacs(2) = 0.70
!     rhfacs(3) = 0.70
!
!-clm rhfacl(1) = 0.85
!
   rhfacl(1) = 0.70
!     rhfacl(2) = 0.70
!     rhfacl(3) = 0.70
!
   do l = 1,km1
     eps(l)   = (prj(l+1)-prh(l  ))/prh(l  )
     del(l+1) = (prh(l+1)-prj(l+1))/prh(l+1)
!
     tem = cp / (prh(l+1)-prh(l))
     sdu(l+1) = prh(l+1) * (prj(l+1)-prh(l  )) * tem
     sdd(l  ) = prh(l  ) * (prh(l+1)-prj(l+1)) * tem
   enddo
!
   do ibl = 1,ibmax
     kpbl1 = k - ibl
     kpbl  = kpbl1 + 1
     kcr   = min(kpbl1,krmax)
     kfx   = kpbl1 - kcr
     ncmx  = kfx + ncrnd
!
     if (kfx .gt. 0) then
       if (botop) then
         do nc = 1,kfx
           ic(nc) = kpbl - nc
         enddo
       else
         do nc = kfx,1,-1
           ic(nc) = kpbl - nc
         enddo
       endif
     endif
!
!   warning!! if ibl> 3, this wont work!!!
!   ----------------------------------------
!   rannum dimension not sufficient!
!
     if (ncrnd .gt. 0) then
       ii = (ibl-1)*ncrnd
       do i = 1,ncrnd
         irnd = (rannum(i+ii)-0.0005)*(kcr-krmin+1)
         ic(kfx+i) = irnd + krmin
       enddo
     endif
!
     fracs  = (sge(kpbl) - sge(kpbl1)) * frac
!
     do nc=1,ncmx
!
       ib = ic(nc)
!
       do i = ib,k
         epss(i) = eps(i)
       enddo
       epss(ib) = eps(ib) * 0.5
!
       call ras_cloud(len, lenc, kpbl-ib+1, k-ib+1                             &
         ,dt, rasal(ib), fracs, rhfacl(ibl), rhfacs(ibl)                       &
         ,cp,  alhl, grav, crmsfx, prj(ib)                                     &
         ,cd                                                                   &
         ,sgo(ib),   onbdsg(ib), prh(ib),   fpk(ib), hpk(ib), sge(ib)          &
         ,epss(ib),  del(ib),    sdu(ib),   sdd(ib)                            &
         ,pii,       pik,     pki, gam(1,ib)                                   &
         ,uoi(1,ib), voi(1,ib),  poi(1,ib), qoi(1,ib)                          &
         ,pcu,       cln(1,ib),  q1(1,ib),  q2(1,ib))
!
     enddo
   enddo
!
   return
   end subroutine ras_solver
!
!-------------------------------------------------------------------------------
   subroutine ras_cloud (len, lenc, k, km, dt, rasal, frac                     &
                    ,rhfacl, rhfacs                                            &
                    ,cp,  alhl, grav, crmsfx, prj                              &
                    ,cd                                                        &
                    ,sgo, onbdsg, prh,  fpk, hpk, sge                          &
                    ,eps, del,    sdu,  sdd                                    &
                    ,pii, pik, pki, gam                                        &
                    ,uoi, voi, poi, qoi                                        &
                    ,pcu, cln, tcu, qcu)
!-------------------------------------------------------------------------------
   real, allocatable :: hol(:,:), qol(:,:), eta(:,:), hst(:,:)                 &
            ,gmh(:,:)
   real, allocatable :: tx1(:), tx2(:), tx3(:), tx4(:), tx5(:)                 &
            ,tx6(:), alm(:), wlq(:), wfn(:), akm(:)                            &
            ,txc(:), txe(:), clf(:), uht(:), vht(:)                            &
            ,qs1(:)
   integer, allocatable :: ia(:),  i1(:),  i2(:)
   real, allocatable :: clp(:)
!
!   subroutine ras_cloud
!     16 march   1988
!
   real,parameter       ::  zero=0.0, one=1.0, half=0.5, cmb2pa=100.0
   real,parameter       ::  daylen=86400.0,    epsln=0.001, qmin=1.0e-8
   real,parameter       ::  rhram=0.35
   integer,parameter    ::  nx=7501
   integer,parameter    ::  ic=1
#ifndef ICE
   real                 ::  tbpvs0(nx)
   common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
#endif
!
! from mrf model
!
   real,parameter       ::  rd= 2.8705e+2 , rv= 4.6150e+2
   real,parameter       ::  epsa=rd/rv,     epsam1=rd/rv-1.
!
   real,parameter       ::  zm1p04 = -1.04e-4, two89= 2.89e-5, zp44= 0.44
   real,parameter       ::  zp01= 0.01, zp1 = 0.1, zp001= 0.001
   real,parameter       ::  zp578 = 0.578, z3600 = 3600.
   real,parameter       ::  thousand = 1000.
!
   real                 ::  gam(len,k)
   real                 ::  pii(len),    pik(len),    pki(len)
   real                 ::  uoi(len,k),  voi(len,k),  poi(len,km),  qoi(len,km)
   real                 ::  cd(len)
!
   real                 ::  tcu(len,km),  qcu(len,km),  cln(len,k), pcu(len)
!
   real                 ::  onbdsg(km), fpk(km), hpk(km),  prh(km)
   real                 ::  prj(km),    sgo(km), sge(km)
   real                 ::  eps(km),    del(km), sdu(km),  sdd(km),  bta(km)
!
!
! from ras2_liquid_water
!
   real,parameter       ::  h2omw  = 18.01, airmw  = 28.97
   real,parameter       ::  esfac = h2omw/airmw       
   real,parameter       ::  erfac = (1.0-esfac)/esfac 
!
   real                 ::  tbsvp(nx), tbdsvp(nx)
   common/comsvp/ c1xsvp,c2xsvp,tbsvp, tbdsvp
!-------------------------------------------------------------------------------
   allocate (hol(len ,k), qol(len,k), eta(len,k), hst(len,k)                   &
            ,gmh(len,k))
   allocate (tx1(lenc), tx2(lenc), tx3(lenc), tx4(lenc), tx5(lenc)             &
            ,tx6(lenc), alm(lenc), wlq(lenc), wfn(lenc), akm(lenc)             &
            ,txc(lenc), txe(lenc), clf(lenc), uht(lenc), vht(lenc)             &
            ,qs1(lenc))
   allocate (ia(lenc),  i1(lenc),  i2(lenc))
   allocate (clp(lenc))
!
   km1 = k  - 1
   ic1 = ic + 1
   ip1 = 1
!
   onebg  = 1.0  / grav
   twobal = 2.0 / alhl
   albycp = alhl / cp
   el2orc = alhl*alhl / (rv*cp)
!
   do i = 1,k
     bta(i) = fpk(i) * onebg
   enddo
   bta(ic)  = hpk(ic)  * onebg
!
   do l = ic,k
     do i = 1,lenc
       tl = poi(i,l) * pik(i) * prh(l)
       pl = pii(i) * sgo(l)
!
       xj  = min(max(c1xsvp+c2xsvp*tl,1.),float(nx))
       jx  = min(xj,nx-1.)
       qx  = tbsvp(jx)  + (xj-jx)*(tbsvp(jx+1)-tbsvp(jx))
       dqx = tbdsvp(jx) + (xj-jx)*(tbdsvp(jx+1)-tbdsvp(jx))
       d = (pl-erfac*qx)
       if (d .gt. 0.) then
         d    = 1.0 / d
         gmh(i,l) = amin1(qx * d,1.0)
         dqs = (1.0 + erfac*gmh(i,l)) * d * dqx
       else
         gmh(i,l) = 1.0
         dqs = 0.
       endif
!
       gam(i,l) = 1.0 / (1.0 + dqs * albycp)
      enddo
   enddo
!
   do i = 1,lenc
     tem      = pik(i) * poi(i,k)
     qs1(i)   = gmh(i,k)
     qol(i,k) = amax1(amin1(qs1(i),qoi(i,k)), qmin)
     hol(i,k) = tem * (hpk(k)+cp*prh(k)) + qol(i,k)*alhl
     tx2(i)   = tem * fpk(k)
     eta(i,k) = zero
!
     clp(i)    = zero
     tem = rhfacs
!
     tem1   = qol(i,k) / qs1(i) - tem
     tx3(i) = tem1 + rhram
     if (tx3(i) .gt. 0.)                                                       &
         clp(i) = max(zero, min(one, exp(20.0*tem1) )) ! 08/08/96
!
   enddo
!
   do l = km1,ic,-1
     do i = 1,lenc
       tem      = pik(i) * poi(i,l)
       qs1(i)   = gmh(i,l)
       qol(i,l) = amax1(amin1(qs1(i),qoi(i,l)), qmin)
       tem1     = tx2(i) + tem * (hpk(l)+cp*prh(l))
       hol(i,l) = tem1   + qol(i,l) * alhl
       hst(i,l) = tem1   + qs1(i)   * alhl
       eta(i,l) = eta(i,l+1) + tem * bta(l)
       tx2(i)   = tx2(i)     + tem * fpk(l)
     enddo
   enddo
!
   do i = 1,lenc
     tx2(i) = hol(i,k)  - hst(i,ic)
     tx1(i) = zero
   enddo
!
!  entrainment parameter alm
!
   do l = ic,km1
     do i = 1,lenc
       tx1(i) = tx1(i) + (hst(i,ic)-hol(i,l)) * (eta(i,l)-eta(i,l+1))
     enddo
   enddo
!
   len1 = 0
   len2 = 0
   isav = 0
   do i = 1,lenc
     if (tx1(i) .gt. epsln .and. (tx2(i)  .gt. zero)                           &
                           .and. (tx3(i)  .gt. zero) ) then
        len1      = len1 + 1
        ia(len1)  = i
        alm(len1) = tx2(i) / tx1(i)
     endif
   enddo
   len2 = len1
!
!  for non-entraining cloud lamda=0
!
   if (k .gt. 2) then
     do i = 1,lenc
       if (tx2(i) .le. 0.0 .and. (hol(i,k) .gt. hst(i,ic1))                    &
                           .and. (tx3(i)   .gt. zero) ) then
          len2      = len2 + 1
          ia(len2)  = i
          alm(len2) = 0.0
       endif
     enddo
   endif
!
   if (len2 .eq. 0) then
     deallocate (hol, qol, eta, hst,gmh)
     deallocate (tx1, tx2, tx3, tx4, tx5                                       &
             , tx6, alm, wlq, wfn, akm                                         &
             , txc, txe, clf, uht, vht                                         &
             , qs1)
     deallocate (ia,  i1,  i2)
     deallocate (clp)
     return
   endif
   len11 = len1 + 1
!
!  normalized massflux
!
   do i = 1,len1
     eta(i,k) = 1.0
     tx2(i) = pii(ia(i)) * sgo(ic)
     txc(i) = 1.0
   enddo
!
   do i = len11,len2
     eta(i,k) = 1.0
     wfn(i)   = 0.0
     ii       = ia(i)
     if (hst(ii,ic1) .lt. hst(ii,ic)) then
       txc(i) = (hst(ii,ic1)-hol(ii,k))/(hst(ii,ic1)-hst(ii,ic))
     else
       txc(i) = 0.0
     endif
     tx2(i) = pii(ii) * (sgo(ic1) - txc(i) * (sgo(ic1)-sgo(ic)))
   enddo
!
   call ras_cwork_crit(len2, tx2, tx3)
!
   do l = km1,ic,-1
     do i = 1,len2
       tx5(i) = 1.0 + alm(i) * eta(ia(i),l)
     enddo
     do i = 1,len2
       eta(i,l) = tx5(i)
     enddo
   enddo
!
   do i = 1,len2
     tx5(i) = 1.0
   enddo
!
!  cloud workfunction
!
   if (len1 .gt. 0) then
     do i = 1,len1
       wfn(i) = - gam(ia(i),ic) * eps(ic) * hst(ia(i),ic) * eta(i,ic1)
     enddo
   endif
!
   do i = 1,len2
     tx1(i) = hol(ia(i),k)
     tem    = sgo(ic1) - txc(i) * (sgo(ic1)-sgo(ic))
     tx3(i) = tx3(i) * pii(ia(i)) * (sge(km1)-tem)
     akm(i) = 0.0
   enddo
!
   if (ic1 .le. km1) then
     do l = km1,ic1,-1
       do i = 1,len2
         tem    = tx1(i) + (eta(i,l) - eta(i,l+1)) * hol(ia(i),l)
!
         wfnt   = gam(ia(i),l) * (tx1(i)*eps(l) + tem*del(l)                   &
               - (eta(i,l+1)*eps(l) + eta(i,l)*del(l)) * hst(ia(i),l))
         wfn(i) = wfn(i) + wfnt
         akm(i) = akm(i) - min(wfnt, 0.0)
!
         tx1(i) = tem
       enddo
     enddo
   endif
!
   if (len1 .gt. 0) then
     do i = 1,len1
       wfn(i) = wfn(i) + tx1(i) * gam(ia(i),ic) * eps(ic)
     enddo
   endif
   do i = 1,len2
     if (wfn(i) .gt. 0.0) then
       akm(i) = akm(i) / wfn(i)
       tem    =  min(cd(ia(i))*100.00,0.50)
       if (akm(i) .gt. tem) wfn(i) = 0.0
     endif
     wfn(i) = wfn(i) - tx3(i)
   enddo
!
   lena = 0
   if (len1 .gt. 0) then
     do i = 1,len1
       if (wfn(i) .gt. 0.0 .and. tx5(i) .ge. 0.0) then
         lena = lena + 1
         i1(lena) = ia(i)
         i2(lena) = i
         tx1(lena) = wfn(i)
         tx2(lena) = qs1(ia(i))
         txc(lena) = 1.0
       endif
     enddo
   endif
!
   lenb = lena
   do i = len11,len2
     if (wfn(i) .gt. 0.0 .and. txc(i) .gt. 0.0                                 &
                         .and. tx5(i) .ge. 0.0) then
       lenb = lenb + 1
       i1(lenb)  = ia(i)
       i2(lenb)  = i
       tx1(lenb) = wfn(i)
       tx2(lenb) = qs1(ia(i))
       tx4(lenb) = txc(i)
     endif
   enddo
!
   if (lenb .eq. 0) then
     deallocate (hol, qol, eta, hst,gmh)
     deallocate (tx1, tx2, tx3, tx4, tx5                                       &
              ,tx6, alm, wlq, wfn, akm                                         &
              ,txc, txe, clf, uht, vht                                         &
              ,qs1)
     deallocate (ia,  i1,  i2)
     deallocate (clp)
     return
   endif
!
   do i = 1,lenb
     wfn(i) = tx1(i)
     qs1(i) = tx2(i)
   enddo
!
   do l = ic,k
#ifdef CRAY_THREAD
!dir$ ivdep
#endif
     do i = 1,lenb
       eta(i,l) = eta(i2(i),l)
     enddo
   enddo
!
   lena1 = lena + 1
!
   do i = 1,lena
     ii = i1(i)
     tx6(i) = hst(ii,ic) - hol(ii,ic)
     uht(i) = uoi(ii,ic)
     vht(i) = voi(ii,ic)
   enddo
!
   do i = lena1,lenb
     ii = i1(i)
     txc(i) = tx4(i)
     tem    = txc(i) * (hol(ii,ic)-hol(ii,ic1)) + hol(ii,ic1)
     tx6(i) = hol(ii,k) - tem
!
     tem1   = txc(i) * (qol(ii,ic)-qol(ii,ic1)) + qol(ii,ic1)
     tx5(i) = tem    - tem1 * alhl
     qs1(i) = tem1   + tx6(i)*(one/alhl)
     tx3(i) = hol(ii,ic)
!
     uht(i) = txc(i) * (uoi(ii,ic)-uoi(ii,ic1)) + uoi(ii,ic1)
     vht(i) = txc(i) * (voi(ii,ic)-voi(ii,ic1)) + voi(ii,ic1)
   enddo
!
   do i = 1,lenb
     ii = i1(i)
     wlq(i) = qol(ii,k) - qs1(i) * eta(i,ic)
     uht(i) = uoi(ii,k) - uht(i) * eta(i,ic)
     vht(i) = voi(ii,k) - vht(i) * eta(i,ic)
     txe(i) = hol(ii,k)
   enddo
!
   do l = km1,ic,-1
     do i = 1,lenb
       ii = i1(i)
       tem    = eta(i,l) - eta(i,l+1)
       wlq(i) = wlq(i) + tem * qol(ii,l)
       uht(i) = uht(i) + tem * uoi(ii,l)
       vht(i) = vht(i) + tem * voi(ii,l)
     enddo
   enddo
!
!     calculate gs and part of akm (that requires eta)
!
   do i = 1,lenb
     ii = i1(i)
     tem        = (poi(ii,km1) - poi(ii,k)) * pik(ii)
     hol(i,k)   = tem    * (sdu(k)*onbdsg(k))
     hol(i,km1) = tem    * (sdd(km1)*onbdsg(km1))

     akm(i)     = zero
     tx2(i)     = pii(ii)*sgo(ic)
   enddo
!
   if (ic1 .le. km1) then
      do l = km1,ic1,-1
         do i = 1,lenb
            ii = i1(i)
            tem        = (poi(ii,l-1) - poi(ii,l)) * pik(ii) * eta(i,l)
!
            hol(i,l)   = tem * (sdu(l  )*onbdsg(l))  + hol(i,l)
            hol(i,l-1) = tem * (sdd(l-1)*onbdsg(l-1))
!
            akm(i) = akm(i) - hol(i,l)*(eta(i,l)*del(l)+eta(i,l+1)*eps(l))
         enddo
      enddo
   endif
!
   call ras_cloud_factor(lenb, tx2, tx1, clf)
!
   do i = 1,lenb
     tx2(i) = (one - tx1(i)) * wlq(i)
     wlq(i) = tx1(i) * wlq(i)
     tx1(i) = hol(i,ic)
   enddo
!
   do i = lena1,lenb
     ii = i1(i)
     tx1(i) = tx1(i) + (tx5(i)-tx3(i)+qol(ii,ic)*alhl) * onbdsg(ic)
   enddo
!
   do i = 1,lenb
     hol(i,ic) = tx1(i) - tx2(i) * (alhl*onbdsg(ic))
   enddo
!
   if (lena .gt. 0) then
     do i = 1,lena
       akm(i) = akm(i) - eta(i,ic1) * eps(ic) * tx1(i)
     enddo
   endif
!
   temu = (prh(k) - prj(k)) / (prh(k) - prh(k-1))
   temd = 1.0 - temu
!
   do i = 1,lenb
     ii = i1(i)
     tx3(i)   =  temu*qol(ii,km1) + temd*qol(ii,k)
     gmh(i,k) = hol(i,k) + (tx3(i)-qol(ii,k)) * (alhl*onbdsg(k))
     akm(i) = akm(i) + gam(ii,km1)*eps(km1) * gmh(i,k)
   enddo
!
   if (ic1 .le. km1) then
     do l = km1,ic1,-1
       temu = (prh(l) - prj(l)) / (prh(l) - prh(l-1))
       temd = 1.0 - temu
!
       do i = 1,lenb
         ii = i1(i)
         tx2(i) = tx3(i)
         tx3(i)   =  temu*qol(ii,l-1) + temd*qol(ii,l)
!
         gmh(i,l) = hol(i,l) + ((qol(ii,l) - tx2(i)) * eta(i,l+1)              &
                          +  (tx3(i)    - qol(ii,l)) * eta(i,l))               &
                          *  (alhl*onbdsg(l))
       enddo
     enddo
   endif
!
   do i = lena1,lenb
     tx3(i) = (txe(i)-tx6(i)-tx5(i)) - tx3(i) * alhl
   enddo
   do i = 1,lena
     tx3(i) = (qol(i1(i),ic) - tx3(i)) * alhl
   enddo
   do i = 1,lenb
     gmh(i,ic) = tx1(i) + onbdsg(ic)                                           &
                         * (eta(i,ic1)*tx3(i) + eta(i,ic)*tx6(i))
   enddo
!
!     calculate hc part of akm
!
   if (ic1 .le. km1) then
     do i = 1,lenb
       tx1(i) = gmh(i,k)
     enddo
     do l = km1,ic1,-1
       do i = 1,lenb
         ii = i1(i)
         tx1(i) = tx1(i) + (eta(i,l) - eta(i,l+1)) * gmh(i,l)
         tx2(i) = gam(ii,l-1) * eps(l-1)
       enddo
!
       if (l .eq. ic1) then
         do i = lena1,lenb
           tx2(i) = zero
         enddo
       endif
!
       do i = 1,lenb
         akm(i) = akm(i) + (tx2(i) + gam(i1(i),l)*del(l)) * tx1(i)
       enddo
     enddo
   endif
!
   do i = lena1,lenb
     tx1(i) = sgo(ic) + (sgo(ic1) - sgo(ic)) * (one-txc(i))
     if ((tx1(i) .ge. sge(ic)) .and. (tx1(i) .le. sgo(ic1))) then
       txc(i)     = one - (tx1(i) - sge(ic)) / (sgo(ic1) - sge(ic))
!
       hol(i,ic1) = hol(i,ic1) + hol(i,ic) * (onbdsg(ic1)/onbdsg(ic))
       hol(i,ic)  = zero
!
       gmh(i,ic1) = gmh(i,ic1) + gmh(i,ic)*(onbdsg(ic1)/onbdsg(ic))
       gmh(i,ic)  = zero
     elseif (tx1(i) .lt. sge(ic)) then
       txc(i) = 1.0
     else
       txc(i) = 0.0
     endif
   enddo
!
   tem0 = pcu(i1(ip1))
   do i = 1,lenb
     ii = i1(i)
     if (akm(i) .lt. zero .and. wlq(i) .ge. 0.0) then
       wfn(i) = - clp(ii) * txc(i) * wfn(i) * rasal / akm(i)
     else
        wfn(i) = zero
     endif
     wfn(i) = amin1(wfn(i), frac)
!
!     compute cloud amount
!
     if (txc(i) .gt. 0.9999) then
        tx1(i) = cln(ii,ic)
     else
        tx1(i) = cln(ii,ic1)
     endif
     if (wfn(i) .gt. crmsfx)  tx1(i) = tx1(i) + clf(i)
     if (tx1(i) .gt. one)  tx1(i) = one
!
!     precipitation
!
     wlq(i) = (cmb2pa/grav) * wlq(i) * wfn(i) * pii(ii)
!
   enddo
!
   do i = 1,lenb
     ii = i1(i)
     if (txc(i) .gt. 0.9999) then
        cln(ii,ic) = tx1(i)
     else
        cln(ii,ic1) = tx1(i)
     endif
   enddo
!
!     theta and q change due to cloud type ic
!

   do l = ic,k
#ifdef CRAY_THREAD
!dir$ ivdep
#endif
     do i = 1,lenb
       ii = i1(i)
       tem      = (gmh(i,l) - hol(i,l)) * wfn(i)
       tem1     =  hol(i,l) * wfn(i)
       poi(ii,l) = poi(ii,l) + tem1 * pki(ii) * (one/(cp*prh(l)))
       qoi(ii,l) = qoi(ii,l) + tem * (one/alhl)
!
       tcu(ii,l) = tcu(ii,l) + tem1 * (daylen/(cp*dt))
       qcu(ii,l) = qcu(ii,l) + tem  * (daylen/(alhl*dt))
     enddo
!
   enddo
   do i = 1,lenb
     ii = i1(i)
!
!       cumulus precipitation
!
!       cumulus friction at the bottom level
!
     tem      = wfn(i) * (0.5*onbdsg(k))
     tx1(i)   = uoi(ii,km1) - uoi(ii,k)
     tx2(i)   = voi(ii,km1) - voi(ii,k)
!
     uoi(ii,k) = uoi(ii,k) + tx1(i) * tem
     voi(ii,k) = voi(ii,k) + tx2(i) * tem
   enddo
!
!     cumulus friction at all other levels
!
   do l = km1,ic1,-1
#ifdef CRAY_THREAD
!dir$ ivdep
#endif
     do i = 1,lenb
       ii = i1(i)
       tem = wfn(i) * (0.5*onbdsg(l))
       tem1   = tx1(i)
       tem2   = tx2(i)
       tx1(i) = (uoi(ii,l-1) - uoi(ii,l)) * eta(i,l)
       tx2(i) = (voi(ii,l-1) - voi(ii,l)) * eta(i,l)
!
       uoi(ii,l) = uoi(ii,l) + (tx1(i) + tem1) * tem
       voi(ii,l) = voi(ii,l) + (tx2(i) + tem2) * tem
     enddo
   enddo
!
#ifdef CRAY_THREAD
!dir$ ivdep
#endif
   do i = 1,lenb
     ii = i1(i)
     if (txc(i) .ge. 1.0) then
       tem    = wfn(i) * (onbdsg(ic)*0.5)
       uoi(ii,ic) = uoi(ii,ic) + (tx1(i) + uht(i) + uht(i)) * tem
       voi(ii,ic) = voi(ii,ic) + (tx2(i) + vht(i) + vht(i)) * tem
     elseif (txc(i) .gt. 0.0) then
       tem    = wfn(i) * (onbdsg(ic1)*0.5)
       uoi(ii,ic1) = uoi(ii,ic1) - (tx1(i) - uht(i) + uht(i)) * tem
       voi(ii,ic1) = voi(ii,ic1) - (tx2(i) - vht(i) + vht(i)) * tem
     endif
   enddo
!
   ipa = 1
   if (ipa .gt. lenb) ipa = 1
   ipo = i1(ipa)
!
   rphf = z3600/dt
   tem = onbdsg(ic)
   do i = 1,lenb
     tx5(i) = tem * wfn(i) * eta(i,ic)
     if (tx5(i) .gt. 1.0) tx5(i) = 1.0
   enddo
!
    do  i = 1,lenb
      ii = i1(i)
      tx6(i) = wlq(i)
      tx4(i) = sqrt(pii(ii))
      txe(i) = (grav/cmb2pa) / pii(ii)
    enddo
!
!      do loop for moisture evaporation ability and convec evaporation.
!
    do 5000 l=ic,km
      tem    = rphf * sqrt(sgo(l)*zp001)
      dsgbg  = (onebg*cmb2pa) / onbdsg(l)
      rknobt = 5.0
!
      albyp = albycp / prh(l)
#ifdef CRAY_THREAD
!dir$ ivdep
#endif
    do i = 1,lenb
       ii = i1(i)
       tx1(i) = poi(ii,l) * pik(ii) * prh(l)
       tx3(i) = qoi(ii,l)
       txc(i) = sgo(l)*pii(ii)
!
!   first iteration - increment halved
!
#ifdef ICE
       es   = 10.0*fpvs(tx1(i))
#else
       xj1=min(max(c1xpvs0+c2xpvs0*tx1(i),1.),float(nx))
       jx1=min(xj1,nx-1.)
       fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
       es   = 10.0*fpvs01
!         es   = 10.0*fpvs0(tx1(i))
#endif
       qs   = epsa * es / (txc(i) + epsam1*es)
       tsq  = tx1(i) * tx1(i)
       delq = 0.5 * (qs - tx3(i)) * tsq / (tsq + el2orc * qs)
!
       tx2(i) = delq
       tx1(i) = tx1(i) - delq * albycp
       tx3(i) = tx3(i) + delq
!
!   second iteration
!
#ifdef ICE
       es   = 10.0*fpvs(tx1(i))
#else
       xj1=min(max(c1xpvs0+c2xpvs0*tx1(i),1.),float(nx))
       jx1=min(xj1,nx-1.)
       fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
       es   = 10.0*fpvs01
!         es   = 10.0*fpvs0(tx1(i))
#endif
       qs   = epsa * es / (txc(i) + epsam1*es)
       tsq  = tx1(i) * tx1(i)
       delq = (qs - tx3(i)) * tsq / (tsq + el2orc * qs)
!
       tx2(i) = tx2(i) + delq
       tx1(i) = tx1(i) - delq * albycp
       tx3(i) = tx3(i) + delq
!
!   third iteration
!
#ifdef ICE
       es   = 10.0*fpvs(tx1(i))
#else
       xj1=min(max(c1xpvs0+c2xpvs0*tx1(i),1.),float(nx))
       jx1=min(xj1,nx-1.)
       fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
       es   = 10.0*fpvs01
!         es   = 10.0*fpvs0(tx1(i))
#endif
       qs   = epsa * es / (txc(i) + epsam1*es)
       tsq  = tx1(i) * tx1(i)
       delq = (qs - tx3(i)) * tsq / (tsq + el2orc * qs)
!
       tx2(i) = tx2(i) + delq
!
!   "evaporation efficiency" using  kesslers parameterization         !
!

#ifdef CRAY_THREAD
!dir$ ivdep
#endif
      ii = i1(i)
!
      tx2(i) = max(tx2(i) * pii(ii) * dsgbg, 0.0)
      exparg  = zm1p04 * dt * (tx6(i)*tx4(i)*tem) ** zp578
!
      tx3(i) = tx2(i) * (1.0 - exp(exparg)) * min(tx5(i)*rknobt,1.0)
!
      tx3(i) = min(tx6(i), tx3(i))
!
      tx6(i) = tx6(i) - tx3(i)
!
      tx3(i)    = tx3(i) * txe(i) * onbdsg(l)
!
      qoi(ii,l) = qoi(ii,l) + tx3(i)
      poi(ii,l) = poi(ii,l) - tx3(i) * pki(ii) * albyp
!
      tcu(ii,l) = tcu(ii,l) - tx3(i) * (albycp*daylen/dt)
      qcu(ii,l) = qcu(ii,l) + tx3(i) * (daylen/dt)
!
   enddo
!
 5000  continue
!
   do i = 1,lenb
     ii = i1(i)
     pcu(ii) = pcu(ii) + tx6(i) * (daylen/dt)
   enddo
!
   deallocate (hol, qol, eta, hst,gmh)
   deallocate (tx1, tx2, tx3, tx4, tx5                                         &
              ,tx6, alm, wlq, wfn, akm                                         &
              ,txc, txe, clf, uht, vht                                         &
              ,qs1)
   deallocate (ia,  i1,  i2)
   deallocate (clp)
!
   return
   end subroutine ras_cloud
!
!-------------------------------------------------------------------------------
   subroutine ras_cwork_crit(len, pl, acr)
!-------------------------------------------------------------------------------
!
! subroutine:  ras_cwork_crit
!
!   1991-02-01
!
!-------------------------------------------------------------------------------
   real                 ::  pl(len), acr(len)
   integer, allocatable ::  iwk(:)
   real                 ::  a(15), ac(15), ad(15)
   common /rasacr/ a, ac, ad, actop
!
   allocate (iwk(len))
!
   do i = 1,len
     iwk(i) = pl(i) * 0.02 - 0.999999999
     if (iwk(i) .gt. 1) then
       if (iwk(i) .le. 15) then
            acr(i) = ac(iwk(i)) + pl(i) * ad(iwk(i))
       else
            acr(i) = a(15)
       endif
     else
         acr(i) = actop
     endif
   enddo
!
   deallocate (iwk)
!
   return
   end subroutine ras_cwork_crit
!
!-------------------------------------------------------------------------------
   subroutine ras_cloud_factor(len, pl, rno, clf)                                        
!-------------------------------------------------------------------------------
!     parameter (p5=50.0,  p8=900.0, pt8=0.8, pt2=0.2)
!     parameter (p5=200.0,  p8=800.0, pt8=0.4, pt2=0.6)
!     parameter (p5=400.0,  p8=800.0, pt8=0.5, pt2=0.5)
!     parameter (p5=200.0,  p8=800.0, pt8=0.8, pt2=0.2)
   real,parameter       ::  p5=400.0,  p8=800.0, pt8=0.8, pt2=0.2
!     parameter (p5=500.0,  p8=800.0, pt8=0.5, pt2=0.5)
   real,parameter       ::  pfac=pt2/(p8-p5)
   real,parameter       ::  p4=400.0,    p6=950.0
   real,parameter       ::  crtm=5.0e-7, crtmi=1.0/crtm
   real,parameter       ::  cfac=1.0/(p6-p4)
   real                 ::  pl(len),  rno(len), clf(len)                                    
!-------------------------------------------------------------------------------
   do i = 1,len                                                             
     rno(i) = 1.0                                                              
     clf(i) = 1.0                                                              
!                                                                               
     if (pl(i) .ge. p5 .and. pl(i) .le. p8) then                               
       rno(i) = (p8-pl(i))*pfac + pt8                                        
     elseif (pl(i) .gt. p8 ) then                                              
       rno(i) = pt8                                                          
     endif                                                                     
     if (pl(i) .ge. p4 .and. pl(i) .le. p6) then                               
       clf(i) = (p6-pl(i))*cfac                                               
     elseif (pl(i) .gt. p6 ) then                                              
         clf(i) = 0.0                                                           
     endif                                                                     
   enddo
!                                                                               
   return                                                                    
   end subroutine ras_cloud_factor                                                                      
!
!-------------------------------------------------------------------------------
#include "abort.h"

   subroutine ras_setup(lev, si, sl, del, cp, rgas, dt, nsphys, fh             &
                    ,sig, sgb, prh, prj, hpk, fpk, ods, prns                   &
                    ,rasal, lm, krmin, krmax, nstrp                            &
                    ,ncrnd, rannum, afac, ufac)
!-------------------------------------------------------------------------------
   real,parameter       ::  actp=1.7,   facm=1.50
!
   real                 ::  si(lev+1), sl(lev), del(lev)
   real                 ::  sig(lev+1), prj(lev+1), prh(lev),fpk(lev), hpk(lev)
   real                 ::  sgb(lev),   ods(lev),   rasal(lev), prns(lev/2)
   real                 ::  rannum(200,nsphys)
   real                 ::  ph(15), a(15), ac(15), ad(15), au(15)
#ifdef MAC
   integer, allocatable ::  nrnd(:)
#endif
!
#ifdef SX6
   real,parameter       ::  seed=1.0
#endif
#ifdef NEC
   real,parameter       ::  seed=1.0
#endif
#ifdef ES
   real,parameter       ::  seed=1.0
#endif
!
   common /rasacr/ a, ac, ad, actop
!
   real                 ::  adata(15)
!
   data ph/150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0              &
          ,550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/
!
   data adata/ 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677                  &
          ,0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664                      &
          ,0.0553, 0.0445, 0.0633/
!
   logical,save         ::  first
   integer,save         ::  iseed
   data iseed/0/, first/.true./, fh0/0./
!
#ifdef NEC
   real*16 rseed
#endif
#ifdef ES
   real*16 rseed
#endif
#ifdef SX6
   real*16 rseed
#endif
!-------------------------------------------------------------------------------
!
!     set critical workfunction arrays
!
   if (first) then
     do i = 1,15
       a(i)=adata(i)
     enddo
     first   = .false.
     actop   = actp*facm
     do l = 1,15
       a(l) = a(l)*facm
     enddo
     do l = 2,15
       tem   = ph(l) - ph(l-1)
       au(l) = a(l-1) / tem
       ad(l) = a(l)   / tem
       ac(l) = ph(l)*au(l) - ph(l-1)*ad(l)
       ad(l) = ad(l) - au(l)
     enddo
     call ras_setup_sub
   endif
!
!      set other parameters
!
   if (fh .eq. fh0) then
     iseed = iseed + 1
   else
     iseed = fh*3600.0 + 0.0001
     iseed = max(iseed, 1)
     fh0   = fh
   endif
!
   rkap = rgas / cp
!
   nstrp = 1
   do l = 2, lev
     if (sl(l) .gt. 0.95) nstrp = nstrp + 1
   enddo
!  nstrp = 1
   do l = 1,nstrp+1
     prj(l) = (si(l)/1000.0) ** rkap
   enddo
!
   afac = 0.0
   do l = 1,nstrp
     prns(l) =  (si(l)*prj(l)-si(l+1)*prj(l+1)) /                              &
                   ((si(l) - si(l+1)) * (rkap+1.0) )
     afac    = afac + prns(l) * del(l)
   enddo
   afac = afac / (si(1) - si(nstrp+1))
!
   lm   = lev - nstrp + 1
   lmm1 = lm - 1
   do l = 1,lm
     sig(l) = si(lev-l+2)
     prj(l) = (sig(l)/1000.0) ** rkap
   enddo
   sig(lm+1) = si(1)
   prj(lm+1) = (sig(lm+1)/1000.0) ** rkap
!
   krmin = 1
   krmax = lmm1
   do l = 1,lmm1
     sgb(l) = sl(lev-l+1)
     if (sgb(l) .le. 0.067) krmin = l
     if (sgb(l) .le. 0.76) krmax = l
   enddo
   sgb(lm) = 0.5 * (sig(lm) + sig(lm+1))
!
   ncrnd   = 48 * (dt/3600) + 0.50001
   ncrnd   = max(ncrnd, 1)
   rasalf  = 0.30
!
   if (ncrnd*3 .gt. 100) then
     write(6,*)' dimension of rannum too small -- job terminated'
     call MPABORT
   endif
#ifdef LINUX_INTEL
   call random_seed(iseed)
   call random_number(rannum)
#else

#define DEFAULT
#ifdef RANF
#undef DEFAULT
   call ranset(iseed)
#endif
#ifdef DEFAULT
#ifdef NEC
!  call xsrandom(iseed)
#else
#ifdef ES
!  call xsrandom(iseed)
#else
#ifdef SX6
!  call xsrandom(iseed)
#else
#ifndef MAC
   call srand(iseed)
#endif
#endif
#endif
#endif
#endif

#ifndef MAC
   do nn = 1,nsphys
     do i = 1,ncrnd*3
#define DEFAULT
#ifdef RANF
#undef DEFAULT
       rannum(i,nn) = ranf()
#endif
#ifdef DEFAULT
#ifdef NEC
       rseed=iseed        
       rannum(i,nn) = frand(rseed)
#else
#ifdef ES
       rseed=iseed
       rannum(i,nn) = frand(rseed)
#else
#ifdef SX6
       rseed=iseed
       rannum(i,nn) = frand(rseed)
#else
       rannum(i,nn) = rand()
#endif
#endif
#endif
#endif
     enddo
   enddo
#endif
#endif

#ifdef MAC
   call random_seed(size=krsize)
   allocate ( nrnd(krsize) )
   nrnd=iseed
   call random_seed(put=nrnd)
   call random_number(rannum)
   deallocate ( nrnd )
#endif
!
   critsig=0.95
!
   do l = 1,lm
     prh(l  ) = (sig(l+1)*prj(l+1)-sig(l)*prj(l)) /                            &
                  ( (sig(l+1)-sig(l))*(rkap+1.0) )
     fpk(l  ) = cp * (prj(l+1)-prj(l))
     hpk(l  ) = cp * (prj(l+1)-prh(l))
     ods(l) = sig(l+1) - sig(l)
     rasal(l) = rasalf
!
     if(sig(l).gt.critsig) rasal(l)=0.
!
   enddo
!
   do l = 1,lm
     ods(l) = 1.0 / ods(l)
   enddo
   afac = 1.0 / (afac - prh(lmm1))
   ufac = 0.0
!
   return
   end subroutine ras_setup
!
!-------------------------------------------------------------------------------
   subroutine ras_setup_sub     
!-------------------------------------------------------------------------------
   integer,parameter    ::  nx=7501
   real                 ::  tbsvp(nx), tbdsvp(nx)                                           
   common/comsvp/ c1xsvp,c2xsvp,tbsvp, tbdsvp                                
!-------------------------------------------------------------------------------
   xmin=180.0                                                                
   xmax=330.0                                                                
   xinc=(xmax-xmin)/(nx-1)                                                   
   c1xsvp=1.-xmin/xinc                                                       
   c2xsvp=1./xinc                                                            
!
   do jx = 1,nx                                                                
     x=xmin+(jx-1)*xinc                                                      
     t=x                                                                     
     call ras_setup_param(t, tbsvp(jx), tbdsvp(jx))  
   enddo                                                                     
!
   return                                                                    
   end subroutine ras_setup_sub                                                                      
!
!-------------------------------------------------------------------------------
   subroutine ras_setup_param(tt, es, desdt)       
!-------------------------------------------------------------------------------
   implicit none                                                             
   real                 ::  tt, es, desdt                                                        
 !                                                                            
   real                 ::  airmw, h2omw, lice, runiv, rgas                               
   parameter ( airmw  = 28.97      )                                         
   parameter ( h2omw  = 18.01      )                                         
   parameter ( lice   = 2.834e6    )                                         
   parameter ( runiv  = 8314.3     )                                         
   parameter ( rgas   = runiv/airmw)                                         
 !                                                                            
   real                 ::  tmin, tmax, zeroc, tmix                                              
   real                 ::  tstarr1, tstarr2, tstarr3, tstarr4                                   
   parameter ( zeroc   = 273.16 )                                            
   parameter ( tmin    = -95. )                                              
   parameter ( tstarr1 = -75. )                                              
   parameter ( tstarr2 = -65. )                                              
   parameter ( tstarr3 = -50. )                                              
   parameter ( tstarr4 = -40. )                                              
   parameter ( tmix    = -20. )                                              
   parameter ( tmax    = +60. )                                              
 !                                                                            
   real           ::  esfac, erfac                                                  
   parameter ( esfac = h2omw/airmw       )                                   
   parameter ( erfac = (1.0-esfac)/esfac )                                   
 !                                                                            
   real                 ::  b0, b1, b2, b3, b4, b5, b6                                     
   parameter (b6= 6.136820929e-11*esfac)                                     
   parameter (b5= 2.034080948e-8 *esfac)                                     
   parameter (b4= 3.031240396e-6 *esfac)                                     
   parameter (b3= 2.650648471e-4 *esfac)                                     
   parameter (b2= 1.428945805e-2 *esfac)                                     
   parameter (b1= 4.436518521e-1 *esfac)                                     
   parameter (b0= 6.107799961e+0 *esfac)                                     
!                                                                             
   real                 ::  bi0, bi1, bi2, bi3, bi4, bi5, bi6                                    
   parameter (bi6= 1.838826904e-10*esfac)                                    
   parameter (bi5= 4.838803174e-8 *esfac)                                    
   parameter (bi4= 5.824720280e-6 *esfac)                                    
   parameter (bi3= 4.176223716e-4 *esfac)                                    
   parameter (bi2= 1.886013408e-2 *esfac)                                    
   parameter (bi1= 5.034698970e-1 *esfac)                                    
   parameter (bi0= 6.109177956e+0 *esfac)                                    
!                                                                             
   real                 ::  s10, s11, s12, s13, s14, s15, s16                                    
   parameter (s16= 0.516000335e-11*esfac)                                    
   parameter (s15= 0.276961083e-8 *esfac)                                    
   parameter (s14= 0.623439266e-6 *esfac)                                    
   parameter (s13= 0.754129933e-4 *esfac)                                    
   parameter (s12= 0.517609116e-2 *esfac)                                    
   parameter (s11= 0.191372282e+0 *esfac)                                    
   parameter (s10= 0.298152339e+1 *esfac)                                    
!                                                                             
   real                 ::  s20, s21, s22, s23, s24, s25, s26                                    
   parameter (s26= 0.314296723e-10*esfac)                                    
   parameter (s25= 0.132243858e-7 *esfac)                                    
   parameter (s24= 0.236279781e-5 *esfac)                                    
   parameter (s23= 0.230325039e-3 *esfac)                                    
   parameter (s22= 0.129690326e-1 *esfac)                                    
   parameter (s21= 0.401390832e+0 *esfac)                                    
   parameter (s20= 0.535098336e+1 *esfac)                                    
!                                                                             
   real                 ::  c1, c2, c3, c4, c5, c6                                         
   parameter (c1= b1   )                                                     
   parameter (c2= b2*2.)                                                     
   parameter (c3= b3*3.)                                                     
   parameter (c4= b4*4.)                                                     
   parameter (c5= b5*5.)                                                     
   parameter (c6= b6*6.)                                                     
!                                                                             
   real                 ::  ci1, ci2, ci3, ci4, ci5, ci6                                   
   parameter (ci1= bi1   )                                                   
   parameter (ci2= bi2*2.)                                                   
   parameter (ci3= bi3*3.)                                                   
   parameter (ci4= bi4*4.)                                                   
   parameter (ci5= bi5*5.)                                                   
   parameter (ci6= bi6*6.)                                                   
!                                                                             
   real                 ::  d11, d12, d13, d14, d15, d16                                  
   parameter (d11= s11   )                                                   
   parameter (d12= s12*2.)                                                   
   parameter (d13= s13*3.)                                                   
   parameter (d14= s14*4.)                                                   
   parameter (d15= s15*5.)                                                   
   parameter (d16= s16*6.)                                                   
!                                                                            
   real                 ::  d21, d22, d23, d24, d25, d26
   parameter (d21= s21   )                                                   
   parameter (d22= s22*2.)                                                   
   parameter (d23= s23*3.)                                                   
   parameter (d24= s24*4.)                                                   
   parameter (d25= s25*5.)                                                   
   parameter (d26= s26*6.)                                                   
!                                                                             
   real                 ::  t, w, qx, dqx                                                        
!-------------------------------------------------------------------------------
   t = amax1(amin1(tt-zeroc,tmax),tmin)                                      
   dqx = 0.                                                                  
   qx = 0.                                                                   
!
   if(t.gt.0.) then                                                          
      qx = (t*(t*(t*(t*(t*(t*b6+b5)+b4)+b3)+b2)+b1)+b0)                        
      dqx = (t*(t*(t*(t*(t*c6+c5)+c4)+c3)+c2)+c1)                              
   elseif(t.lt.tstarr1) then                                                 
      qx = (t*(t*(t*(t*(t*(t*s16+s15)+s14)+s13)+s12)+s11)+s10)                 
      dqx = (t*(t*(t*(t*(t*d16+d15)+d14)+d13)+d12)+d11)                        
   elseif(t.lt.tstarr2) then                                                 
      w = (tstarr2 - t)/(tstarr2-tstarr1)                                      
      qx =     w *(t*(t*(t*(t*(t*(t*s16+s15)+s14)+s13)+s12)+s11)+s10)          &
       + (1.-w)*(t*(t*(t*(t*(t*(t*s26+s25)+s24)+s23)+s22)+s21)+s20)          
    dqx =     w *(t*(t*(t*(t*(t*d16+d15)+d14)+d13)+d12)+d11)                   &
        + (1.-w)*(t*(t*(t*(t*(t*d26+d25)+d24)+d23)+d22)+d21)                 
   elseif(t.lt.tstarr3) then                                                 
    qx = (t*(t*(t*(t*(t*(t*s26+s25)+s24)+s23)+s22)+s21)+s20)                 
    dqx = (t*(t*(t*(t*(t*d26+d25)+d24)+d23)+d22)+d21)                        
   elseif(t.lt.tstarr4) then                                                 
    w = (tstarr4 - t)/(tstarr4-tstarr3)                                      
    qx =     w *(t*(t*(t*(t*(t*(t*s26+s25)+s24)+s23)+s22)+s21)+s20)            &
       + (1.-w)*(t*(t*(t*(t*(t*(t*bi6+bi5)+bi4)+bi3)+bi2)+bi1)+bi0)          
    dqx =     w *(t*(t*(t*(t*(t*d26+d25)+d24)+d23)+d22)+d21)                   &
        + (1.-w)*(t*(t*(t*(t*(t*ci6+ci5)+ci4)+ci3)+ci2)+ci1)                 
   elseif(t.lt.tmix) then                                                    
    qx = (t*(t*(t*(t*(t*(t*bi6+bi5)+bi4)+bi3)+bi2)+bi1)+bi0)                 
    dqx = (t*(t*(t*(t*(t*ci6+ci5)+ci4)+ci3)+ci2)+ci1)                        
   else                                                                      
    w = (tmix - t)/tmix                                                      
    qx =     w *(t*(t*(t*(t*(t*(t*b6+b5)+b4)+b3)+b2)+b1)+b0)                   &
       + (1.-w)*(t*(t*(t*(t*(t*(t*bi6+bi5)+bi4)+bi3)+bi2)+bi1)+bi0)          
    dqx =     w *(t*(t*(t*(t*(t*c6+c5)+c4)+c3)+c2)+c1)                         &
         + (1.-w)*(t*(t*(t*(t*(t*ci6+ci5)+ci4)+ci3)+ci2)+ci1)                
   endif                                                                     
!
   es    = qx                                                                
   desdt = dqx                                                               
!                                                                             
   return                                                                    
   end subroutine ras_setup_param
!-------------------------------------------------------------------------------
#endif /* RAS end */
   end module phys_ras_module
