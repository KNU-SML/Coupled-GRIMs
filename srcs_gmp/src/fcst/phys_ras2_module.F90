#include <define.h>
   module phys_ras2_module
#ifdef RASV2
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [phys_cps_ras2]
!      |
!      |---[phys_ras2_module]
!               |
!               |--- [ras2_solver] --- [ras2_cloud_main] *
!                                           |--- [ras2_liquid_water] *
!                                           |--- [ras2_cwork_crit] *
!                                           |--- [ras2_down_draft] *
!                                                   |- [ras2_down_step1] *
!                                                   |- [ras2_down_step2] *
!                                                   |- function : [qrpf] *
!                                                   |- function : [vtpf] *
!    [dyn_sph_driver] or [phys_main_driver]
!      |
!      |---[phys_ras2_module]
!               |
!               |--- [ras2_setup] --- [ras2_setup_sub1] *
!               |                  |- [ras2_setup_sub2] *
!               |
!               |--- [ras2_setup_hybrid] *
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
!-------------------------------------------------------------------------------
   subroutine ras2_solver(k,   dt, pdd                                         &
                   ,ncrnd, krmin, krmax, kfmax, frac, rasal, botop             &
                   ,revap, max_neg_bouy, alfint, cumfrc                        &
!    *,                revap, max_neg_bouy, alfint, alfinq
                   ,cp,  alhl, alhf, grav, rkap, rannum                        &
                   ,mct, kctop                                                 &
!    *,                mct, kctop, sgh, sig
!
                   ,poi, qoi, uvi, kbl                                         &
!    *,                qli, qii
                   ,qli, qii, clw                                              &
                   ,rain, cln,    q1,  q2, pcu, flx, rcu                       &
                   ,cdrag, dsfc                                                &
                   ,prs,   psj, wrkfun, calkbl, crtfun, updret)
!-------------------------------------------------------------------------------
!
! subroutine: ras2_solver
!
!********************************************************************!
!********************************************************************!
!******************** relaxed arakawa-schubert **********************!
!************************ parameterization **************************!
!********************* plug compatible driver  *********************!
!************************* 27 february 1998 **************************!
!********************************************************************!
!************************** developed by ****************************!
!**************************              ****************************!
!************************ shrinivas moorthi *************************!
!************************                ****************************!
!************************      emc/ncep  ****************************!
!********************************************************************!
!********************************************************************!
!
! program history log:
!   1999-09-01  moorthi         initial implementation
!
! references : moorthi and suarez (1992, mwr)
!
!-------------------------------------------------------------------------------
   real,parameter       ::  one=1.0, daylen=86400.0, pfac=1.0/450.0
   integer,parameter    ::  icm=100
!
   integer              ::  k, mct, ncrnd
   real                 ::  rasal(k-1), uvi(k,2), poi(k), qoi(k)
!    *,    qli(k),     qii(k)
   real                 ::  qli(k),     qii(k), clw(k), cln(k)
   real                 ::  cdrag,      wfnc,   rcu(k,2), dsfc
   real                 ::  dt, pdd,    frac,   max_neg_bouy
   real                 ::  cp, alhl,   alhf,   grav, rkap
   integer              ::  kbl, kctop(mct+1), krmin, krmax, kfmax
!
   real                 ::  q1(k), q2(k), flx(k), pcu(k), rain
   real                 ::  prs(k+1), psj(k+1),  alfint(k), alfinq(k)
   real                 ::  rannum(ncrnd)
!
   integer              ::  ic(icm), iidiag(100)
   logical              ::  botop, revap, lprnt, wrkfun, calkbl, crtfun
   logical              ::  dndrft, updret, cumfrc
!
   real                 ::  cfac, savt, savq, savw(2), st, sq
   real                 ::  sw(2), sprs, sprj, wr1
   real                 ::  rhfacs, rhfacl, tem
!
   integer              ::  kcr, kfx, ncmx, nc, ktem, i, l, lct, lcb, ntrc
!fpp$ noconcur r
   rhfacs = 0.70
   rhfacl = 0.70
!
!     lend = max(1, len / 200)
!     if (lend .le. 100) then
!        do i = 1,lend
!           iidiag(i) = (i-1)*100 + 1
!        enddo
!     endif
!
   kcr   = min(mct,krmax)
   ktem  = min(mct,kfmax)
   kfx   = ktem - kcr
   ncmx  = kfx + ncrnd
!
   if (kfx .gt. 0) then
     if (botop) then 
       do nc = 1,kfx
         ic(nc) = ktem + 1 - nc
       enddo
     else   
       do nc = kfx,1,-1
         ic(nc) = ktem + 1 - nc
       enddo
     endif  
   endif
!
   if (ncrnd .gt. 0) then
     do i = 1,ncrnd
       irnd = (rannum(i)-0.0005)*(kcr-krmin+1)
       ic(kfx+i) = irnd + krmin
     enddo
   endif
!
   ntrc = 0
   if (cumfrc) ntrc = 2
   do nc = 1,ncmx
!
     ib = ic(nc)
!
     lct  = kctop(ib)
     lcb  = kctop(ib+1) - 1
!
!        if (me .eq. 0) then
!        write(6,*) ' lct=',lct,' lcb=',lcb,' ib=',ib
!    *,' kctop=',kctop(ib),kctop(ib+1),' prs=',prs(1,lcb),prs(1,lcb+1)
!        write(6,*) ' ib=',ib,' sgh=',sgh(ib),sgh(ib+1)
!    *,' kctop=',kctop(ib),kctop(ib+1),' sig=',sig(lcb),sig(lcb+1)
!        endif
!
     if (lct .lt. lcb) then
!          tem1    = 1.0 / (sgh(ib+1) - sgh(ib))
!          tem     = (sig(lcb+1) - sig(lcb)) * tem1
!          tem2    = tem
       wr1     = 1.0 / (prs(lcb+1) - prs(lct))
       tem     = (prs(lcb+1) - prs(lcb)) * wr1
!
       savt    = poi(lcb)
       savq    = qoi(lcb)
       savw(1) = uvi(lcb,1)
       savw(2) = uvi(lcb,2)
       sprs    = prs(lcb)
       sprj    = psj(lcb)
!
       poi(lcb)   = poi(lcb)   * tem
       qoi(lcb)   = qoi(lcb)   * tem
       uvi(lcb,1) = uvi(lcb,1) * tem
       uvi(lcb,2) = uvi(lcb,2) * tem
!      prs(lcb)   = prs(lct)
!      psj(lcb)   = psj(lct)
       do l = lct,lcb-1
!        tem  = (sig(l+1) - sig(l)) * tem1 
!        tem2 = tem2 + tem
         tem        = (prs(l+1) - prs(l))   * wr1
         poi(lcb)   = poi(lcb)   + poi(l)   * tem
         qoi(lcb)   = qoi(lcb)   + qoi(l)   * tem
         uvi(lcb,1) = uvi(lcb,1) + uvi(l,1) * tem
         uvi(lcb,2) = uvi(lcb,2) + uvi(l,2) * tem
!
!        if (me .eq. 0) then
!           print *,' tem=',tem,' poi=',poi(i,lcb),poi(i,l)
!    *,' l=',l,' lcb=',lcb,' lct=',lct
!        endif
       enddo
!          write(0,*) ' tem2=',tem2
       st       = poi(lcb)
       sq       = qoi(lcb)
       sw(1)    = uvi(lcb,1)
       sw(2)    = uvi(lcb,2)
       prs(lcb) = prs(lct)
       psj(lcb) = psj(lct)
     endif
!
     dndrft = pdd .gt. 0.0
!
!     write(0,*)' ib=',ib,' lct=',lct,' lcb=',lcb,' kdd=',kdd,' dndrft='
!    *,       dndrft
!
!     if (me .eq. 0) then
!     print *,' calling cloud type ib=',ib,' kbl=',kbl
!    *,' kpbl=',kpbl,' alfint=',alfint,' frac=',frac
!    *,' ntrc=',ntrc
!
!     if (lprnt) then
!     ix = 1
!     print *,' poi=',(poi(ix,l),l=1,k)
!     print *,' qoi=',(qoi(ix,l),l=1,k),' kbl=',kbl
!     endif
!
     do i = 1,k-lcb+1
!      q1(ib)   = 0.0
!      q2(ib)   = 0.0
!      pcu(ib)  = 0.0
       flx(lcb) = 0.0
     enddo
     wfnc = 0.0
!
     do i = lcb,k-1
       alfinq(i) = alfint(i)
     enddo
     alfinq(k) = 0.0
!    alfinq = 1.0     ! 96/09/20
!
!    if (me .eq. 0) then
!    print *,' calling cloud type ib= ', ib,' dt=',dt
!    print *,' poi=',(poi(1,l),l=1,k)
!    print *,' qoi=',(qoi(1,l),l=1,k)
!    endif
!    print *,' alft=',alfint
!
     tla = -10.0
!
     call ras2_cloud_main(k, lcb, ntrc                                         &
                 ,rasal(lcb), frac,  max_neg_bouy                              &
                 ,alfint, alfinq, rhfacl, rhfacs                               &
                 ,cp,  alhl, alhf, rkap, grav                                  &
!
                 ,poi, qoi, uvi, prs, psj                                      &
                 ,qli, qii, kbl, dsfc                                          &
                 ,cdrag                                                        &
                 ,idiag                                                        &
                 ,q1, q2, rcu, pcu, flx                                        &
                 ,rain, revap, dt                                              &
                 ,wfnc, wrkfun, calkbl, crtfun, tla, dndrft, pdd)
!    *,              wfnc, wrkfun, calkbl, crtfun, tla, dndrft, updret)
!     if (me .eq. 0) then
!     print *,' after calling cloud type ib= ', ib
!     print *,' poi=',(poi(l),l=1,k)
!     print *,' qoi=',(qoi(l),l=1,k)
!     endif
! 
!     compute cloud amounts for the goddard radiation
!
     if (flx(kbl) .gt. 0.0) then
       pl = 0.5 * (prs(lct) + prs(lcb+1))
       cfac = min(1.0, max(0.0, (850.0-pl)*pfac))
     else
       cfac = 0.0
     endif
!
     if (lct .lt. lcb) then
       st        = poi(lcb)   - st
       sq        = qoi(lcb)   - sq
       sw(1)     = uvi(lcb,1) - sw(1)
       sw(2)     = uvi(lcb,2) - sw(2)
!
       poi(lcb)   = savt
       qoi(lcb)   = savq
       uvi(lcb,1) = savw(1)
       uvi(lcb,2) = savw(2)
       prs(lcb)   = sprs
       psj(lcb)   = sprj
     
       do l = lct,lcb
         poi(l)   = poi(l)     + st
         qoi(l)   = qoi(l)     + sq
         uvi(l,1) = uvi(l,1)   + sw(1)
         uvi(l,2) = uvi(l,2)   + sw(2)
         cln(l)   = min(cln(l) + cfac, 1.0)
       enddo
       do l = lct,lcb-1
         qli(l) = qli(lcb)
         qii(l) = qii(lcb)
       enddo
     else
       cln(lcb) = min(cln(lcb) + cfac, 1.0)
     endif
!
#ifdef DBG
!     print *,' rain=',rain(1),' ib=',ib
!     tem1 = 0.0
!     tem2 = 0.0
!     tem3 = 0.0
!     tem4 = 0.0
!     do l = 1,k
!        cf1a = (100.0*cp) / (alhl*grav)
!        cf2a = 100.0 / grav
!        print *,' cp=',cp,' alhl=',alhl,' l=',l
!    *,' cf1a=',cf1a,' cf2a=',cf2a
!        print *,' tcu=',(tcu(i,l),i=1,len)
!        print *,' qcu=',(qcu(i,l),i=1,len)
!        do i = 1,len
!           tem = prs(i,l+1)-prs(i,l)
!           cf1 = cf1a * tem
!           cf2 = cf2a * tem
!           tem1 = tem1 + q1(i,l)*cf1
!           tem2 = tem2 + (q2(i,l) - (alhf/alhl)*qii(i,l))*cf2
!           tem4 = tem4 + (qli(i,l)+(1.0+alhf/alhl)*qii(i,l))*cf2

!           tem1 = tem1 + (q1(i,l)-(alhf/cp)*qii(i,l))*cf1
!           tem2 = tem2 + q2(i,l)*cf2
!           tem4 = tem4 + (qli(i,l)+qii(i,l))*cf2
!        enddo
!     enddo
!     print *,' tem1=',tem1,' tem2=',tem2,' tem3=',tem3
!     do i = 1,len
!        tem3 = tem3 + rain(i)
!     enddo
!     tem4 = tem4 + tem3
!     print *,' after cloud: tem1=',tem1,' tem2=',tem2,' tem3=',tem3
!    *,                    ' tem4=',tem4
!     print *,' rain=',(rain(i),i=1,len)
!
!     if (lprnt) then
!     print *,' pcu=',pcu(ix)
!     print *,' poi=',(poi(ix,l),l=1,k)
!     print *,' qoi=',(qoi(ix,l),l=1,k)
!     endif
!
!     print *,' rain=',rain(ix)
#endif

!   warining!!!!
!   by doing the following, cloud does not contain environmental
!   condensate!
!
     do l = 1,k
       clw(l ) = clw(l) + qli(l) + qii(l)
       qli(l)  = 0.0
       qii(l)  = 0.0
!
!      q1(l) = 0.0
!      q2(l) = 0.0
     enddo
!    rain = 0.0
!
   enddo
!
!        cln(ib) = min(cln(ib), 1.0)
!
   return
   end subroutine ras2_solver
!
!-------------------------------------------------------------------------------
   subroutine ras2_cloud_main(                                                 &
                     k, kd, m                                                  &
                    ,rasalf, fracbl, max_neg_bouy                              &
                    ,alfint, alfinq, rhfacl, rhfacs                            &
                    ,cp, alhl, alhf, rkap, grav                                &
                    ,toi, qoi, roi, prs, prj                                   &
                    ,qli, qii, kpbl, dsfc                                      &
                    ,cd                                                        &
                    ,idiag                                                     &
                    ,tcu, qcu, rcu, pcu, flx                                   &
!    u,                 tcd, qcd
                    ,cup, revap, dt                                            &
                    ,wfnc, wrkfun, calkbl, crtfun, tla, dndrft, pdd)
!-------------------------------------------------------------------------------
!
!**********************************************************************!
!******************** relaxed  arakawa-schubert ***********************!
!****************** plug compatible scalar version ********************!
!************************ subroutine cloud  ***************************!
!************************  23 june 1999     ***************************!
!**************************  version 2.0  *****************************!
!*******************  s. moorthi and m.j. suarez **********************!
!**********************************************************************!
!   reference:
!     noaa technical report nws/ncep 99-01:
!     documentation of version 2 of relaxed-arakawa-schubert
!     cumulus parameterization with convective downdrafts, june 1999.
!
!**********************************************************************!
!
!===>    updates cloud tendencies due to a single cloud
!===>    detraining at level kd.
!
!**********************************************************************!
!************* shrinivas.moorthi@noaa.gov (301) 763 8000(x7233) *******!
!***************  max.suarez@gsfc.nasa.gov (301) 286 7373 *************!
!**********************************************************************!
!**********************************************************************!
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!===>  toi(k)     inout   temperature            kelvin
!===>  qoi(k)     inout   specific humidity      non-dimensional
!===>  roi(k,m)   inout   tracer                 arbitrary
!===>  qli(k)     inout   liquid water           non-dimensional
!===>  qii(k)     inout   ice                    non-dimensional

!===>  prs(k+1)   input   pressure @ edges       mb
!===>  prj(k+1)   input   (p/p0)^kappa  @ edges  non-dimensional

!===>  k      input   the rise & the index of the subcloud layer
!===>  kd     input   detrainment level ( 1<= kd < k )          
!===>  m      input   number of gases. may be zero.
!===>  dndrft input   logical .true. or .false.
!===>  pdd    input   pressure level above which downdrfat can exist hpa
!
!===>  tcu(k  )   update  temperature tendency       deg
!===>  qcu(k  )   update  water vapor tendency       (g/g)
!===>  rcu(k,m)   update  tracer tendencies          nd
!===>  pcu(k-1)   update  precip @ base of layer     kg/m^2
!===>  flx(k  )   update  mass flux @ top of layer   kg/m^2
!===>  cup        update  precipitation at the surface kg/m^2
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   real,parameter       ::  almin1=0.00e-6, almin2=0.00e-5, almax=1.0e-1
!     parameter (almin1=0.10e-4, almin2=0.15e-4, almax=1.0e-1)
!
!
!  input arguments
!
   integer              ::  kpbl,   kbl,     kb1
   integer              ::  k, kd, m
   integer              ::  idiag
   logical              ::  revap, dndrft, wrkfun, calkbl, crtfun, calcup
   real                 ::  toi(k  ),  qoi(k ),  prs(k+1)
   real                 ::  qli(k),    qii(k)
   real                 ::  prj(k+1),  roi(k,m)
   real                 ::  cd,        ufn,     dsfc
   real                 ::  rasalf, fracbl, max_neg_bouy, alfint(k)
   real                 ::  rhfacl, rhfacs
   real                 ::  cp, alhl, alhf, rkap, grav, alfinq(k),   pdd
!
!  update arguments
!
   real                 ::  tcu(k), qcu(k), rcu(k,m)
   real                 ::  tcd(k), qcd(k)
   real                 ::  pcu(k), flx(k)
   real                 ::  cup
!
!  temporary work space
!
   real                 ::  hol(kd:k),  qol(kd:k),   gaf(kd:k+1)
   real                 ::  hst(kd:k),  qst(kd:k),   tol(kd:k)
   real                 ::  gmh(kd:k),  gms(kd:k+1), gam(kd:k+1)
   real                 ::  akt(kd:k),  akc(kd:k),   bkc(kd:k)
   real                 ::  ltl(kd:k),  rnn(kd:k),   fco(kd:k)
   real                 ::  prh(kd:k),  pri(kd:k)
   real                 ::  qil(kd:k),  qll(kd:k)
   real                 ::  dlb(kd:k+1),dlt(kd:k+1), eta(kd:k+1)
   real                 ::  prl(kd:k+1)
!
   real                 ::  alm,   det,    hcc,  clp
   real                 ::  hsu,   hsd,    qtl,  qtv
   real                 ::  akm,   wfn,    hos,  qos
   real                 ::  amb,   tx1,    tx2,  tx3
   real                 ::  tx4,   tx5,    qis,  qls
   real                 ::  hbl,   qbl,    rbl(m)
   real                 ::  qlb,   qib,    pris
   real                 ::  wfnc,  tx6,    acr
   real                 ::  tx7,   tx8,    tx9
!
   integer              ::  ia,  i1,  i2, id1, id2
   integer              ::  ib,  i3
!
   logical              ::  unsat
   logical              ::  lowest, skpdd
!
   real                 ::  tl, pl, ql, qs, dqs, st1, sgn
   real                 ::  tau, qtvp, hb, qb, tb
   real                 ::  hccp, ds, dh, ambmax, x00, epp, qtlp, qqq
   real                 ::  dpi, dphib, dphit, del_eta, detp
   real                 ::  qudfac, testmb
   real                 ::  tem, tem1, tem2, tem3, tem4, onebg, rkpp1i, st2
   real                 ::  errh, errw, erre
   real                 ::  tem5, tem6, hbd, qbd
   integer              ::  i, l,  lend1, lend2, n,  kd1, ii
   integer              ::  kp1, it, km1, ktem, kk, kk1, lm1, ll, lp1
!
!     reevaporation
!
   real                 :: clfrac, dt
   real                 :: actevap,arearat,deltaq,mass,massinv,potevap
   real                 :: teq,qsteq,dqdt,qeq
   real                 :: elocp,gravcon,afc, rknob, elfocp
!
   real,parameter       ::  two=2.0,  half=0.5,  one=1.0, zero=0.0
   real,parameter       ::  cmb2pa=100.0   !  conversion mb to pascals
   real,parameter       ::  rhmax=1.0      !  max relative humidity
   real,parameter       ::  quad_lam=1.0   !  mask for quadratic lambda
   real,parameter       ::  rhram=0.15     !  pbl relative humidity ramp
   real,parameter       ::  hcrit=3000.0   !  critical moist static energy
!
   real,parameter       ::  airmw  = 28.97
   real,parameter       ::  h2omw  = 18.01
   real,parameter       ::  esfac = airmw/h2omw-one
!
   real,parameter       ::  c0i=2.0e-3,c0=2.0e-3
   real,parameter       ::  errmin=0.0001, errmi2=0.1*errmin
   integer,parameter    ::  kblmx=10
!     parameter (c0=0.5e-3, kblmx=7, errmin=0.00001, errmi2=0.1*errmin)
   real,parameter       ::  tf=230.16, tcr=260.16, tcrf=1.0/(tcr-tf)
!     parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))
!
!     temporary workspace and parameters needed for downdraft
!
   real                 ::  tla, gmf
!
!     parameter (face=0.0, delx=10000.0, angf=0.5)
!
   real                 ::  buy(kd:k+1), qrb(kd:k),   qrt(kd:k)
   real                 ::  etd(kd:k+1), hod(kd:k+1), qod(kd:k+1)
   real                 ::  ghd(kd:k),   gsd(kd:k),   evp(kd:k)
   real                 ::  etz(kd:k)
   real                 ::  train, dof, cldfrd
   integer              ::  idh
   real                 ::  fac, rsum1, rsum2, rsum3
   logical              ::  ddft, updret
!-----------------------------------------------------------------------

!fpp$ expand (ras2_liquid_water, ras2_cwork_crit)
!fpp$ noconcur r
   tcd = 0.0
   qcd = 0.0
!
   qudfac = quad_lam*half
   onebg  = 1.0 / grav
   rkpp1i = one / (one+rkap)
   kd1    = kd + 1
   kp1    = k  + 1
   km1    = k  - 1
!
   do l = kd, k
     rnn(l) = 0.0
   enddo
!
   elocp  = alhl / cp
   elfocp = (alhl+alhf) / cp
!
   dlt(kp1) = prj(kp1)
   cldfrd   = 0.0
   dof      = 0.0
   prl(kp1) = prs(kp1)
   do l = kd,k
     dlt(l) = prj(l)
     tol(l) = toi(l)
     qol(l) = qoi(l)
     prl(l) = prs(l)
     buy(l) = 0.0
   enddo
!
   do l = kd, k
     dpi = one / (prl(l+1) - prl(l))
     pri(l) = (grav/cmb2pa) * dpi
     prh(l) = (dlt(l+1)*prl(l+1) - dlt(l)*prl(l))*dpi * rkpp1i
     pl     = (prl(l+1) + prl(l))*half
!    tl     = pol(l)*prh(l)
     tl     = tol(l)
     akt(l) = (prl(l+1) - pl) * dpi
!    akt(l) = half
!
     call ras2_liquid_water(tl, pl, qs, dqs)
! 
     qst(l) = qs
     gam(l) = dqs * elocp
     gaf(l) = (one/alhl)*(gam(l)/(one+gam(l)))
     ql     = amax1(amin1(qs*rhmax,qol(l)), 1.0e-10)
     qol(l) = ql
     tem    = cp * tl
     ltl(l) = tem*(one+gam(l))/(one+esfac*(qst(l)+tl*dqs))
     eta(l) = one / (ltl(l) * (1.0 + esfac*ql))
     hol(l) = tem + ql * alhl
     hst(l) = tem + qs * alhl
   enddo
!
   eta(k+1) = zero
   gms(k)   = zero
!
   akt(kd)  = half
   gms(kd)  = zero
!
   clp      = zero
!
   gam(k+1) = gam(k)
   gaf(k+1) = gaf(k)
!
   do l = k,kd1,-1
     tem    = cp * tol(l) / prh(l)
     tem1   = (one + esfac*qol(l)) * tem
! 
     dphib  = (prj(l+1) - prh(l)) * tem1
     dphit  = (prh(l  ) - prj(l)) * tem1
!
     tem    = tem * eta(l)
     dlb(l) = (prj(l+1) - prh(l)) * tem
     dlt(l) = (prh(l  ) - prj(l)) * tem
!
     qrb(l) = dphib
     qrt(l) = dphit
!
     eta(l) = eta(l+1) + dphib
!
     hol(l) = hol(l) + eta(l)
     hst(l) = hst(l) + eta(l)
!
     eta(l) = eta(l) + dphit
   enddo
!
!     for the cloud top layer
!
   l = kd
   tem    = cp * tol(l) * (prj(l+1) - prh(l)) / prh(l)
   dphib  = tem * (one + esfac*qol(l))
!
   dlb(l) = tem * eta(l)
!
   qrb(l) = dphib
   qrt(l) = dphib
!
   eta(l) = eta(l+1) + dphib
!
   hol(l) = hol(l) + eta(l)
   hst(l) = hst(l) + eta(l)
!
!     modification of the environmental moist static energy if ice present
!
   do l = kd,k
     hol(l) = hol(l) - alhf * qii(l)
   enddo
!
!     to determine kbl internally -- if kbl is defined externally
!     the following two loop should be skipped
!
   if (calkbl) then
     kbl   = k
     tx1   = zero
     unsat = .false.
     ktem = max(kd, k-kblmx-2)
     do l = km1,ktem,-1
       tem = hol(k) - hol(l)
       tx3 = (hol(l) - hol(l+1)) / (prl(l+2) - prl(l))
!
       if (tx3 .lt. tx1 .and. tem .lt. hcrit) then
            tx1   = tx3
            kbl   = l+1
            unsat = .true.
       elseif (unsat .and.                                                     &
              ( ((kbl .lt. k-1) .and. tx3 .gt. 0.5*tx1)                        &
                 .or. tem .gt. hcrit) ) then
            tx1 = -1.0e20
       endif
     enddo
!
     do l = ktem,k
       if (hol(k) .gt. hst(l)) kbl = l
     enddo
     kbl = min(k, max(kbl+1, kd-1))
!
     kbl  = max(kbl,k-kblmx)
!
     kpbl = kbl
   else
     kbl  = kpbl
   endif
!
   kbl      = max(kbl,kd)
   kb1      = kbl - 1
   tx1      = eta(kbl)
!
   pris     = one / (prl(k+1)-prl(kbl))
!
   gms(kbl) = 0.0
   gmh(kbl) = 0.0
   rnn(kbl) = 0.0
!
   do l = k,kd,-1
     if (l .ge. kbl) then
       eta(l) = (prl(k+1)-prl(l)) * pris
     else
       rnn(l) = (eta(l) - tx1) * onebg
       eta(l) =  rnn(l) - rnn(l+1)
       gmh(l) =  rnn(l) * rnn(l) * qudfac
       gms(l) =  gmh(l) - gmh(l+1)
     endif
   enddo
!
   hbl = hol(k) * eta(k)
   qbl = qol(k) * eta(k)
   qlb = qli(k) * eta(k)
   qib = qii(k) * eta(k)
   tx1 = qst(k) * eta(k)
!
   do l = km1,kbl,-1
     tem = eta(l) - eta(l+1)
     hbl = hbl + hol(l) * tem
     qbl = qbl + qol(l) * tem
     qlb = qlb + qli(l) * tem
     qib = qib + qii(l) * tem
     tx1 = tx1 + qst(l) * tem
   enddo
!                                   find min value of hol in tx2
   tx2 = hol(kd)
   idh = kd + 1
   do l = kd1,kb1
     if (hol(l) .lt. tx2) then
       tx2 = hol(l)
       idh = l             ! level of minimum moist static energy!
     endif
   enddo
!  idh = 1
   idh = max(kd+1, idh)
!
   tem1 = hbl - hol(kd)
   tem  = hbl - hst(kd1)                                                       &
                - ltl(kd1) *( esfac *(qol(kd1)-qst(kd1)))
   lowest = kd .eq. kb1
!
   tx1   = qbl / tx1
   unsat = (tem .gt. zero .or. (lowest .and. tem1 .ge. zero))                  &
            .and. (tx1 .gt. rhfacs-rhram)                                      &
            .and. (kbl .gt. kd)
!
!===>  if no sounding meets first condition, return
!
   if (.not. unsat) return
!
   do n = 1,m
     rbl(n) = roi(k,n) * eta(k)
     do l = km1,kbl,-1
       rbl(n) = rbl(n) + roi(l,n)*(eta(l)-eta(l+1))
     enddo
   enddo
!
!                   linear case!
!     tem = 1.0 / rhram
!     qtl    = max(zero, min(1.0, (tx1-rhfacs+rhram)*tem))
!
!  commented on 09/23/98 to test linear case!
!
   tem1   = tx1 - rhfacs
   qtl    = max(zero, min(one, exp(20.0*tem1) ))
!
   tx4    = 0.0
   tx5    = 0.0
!
   tx3      = qst(kbl) - gaf(kbl) * hst(kbl)
   qil(kbl) = max(0.0, min(1.0, (tcr-tol(kbl))*tcrf))
!
   do l = kb1,kd1,-1
     tem      = qst(l) - gaf(l) * hst(l)
     tem1     = (tx3 + tem) * 0.5
     st2      = (gaf(l)+gaf(l+1)) * 0.5
!
     fco(l+1) =            tem1 + st2 * hbl
     rnn(l+1) = rnn(l+1) * tem1 + st2 * tx4
     gmh(l+1) = gmh(l+1) * tem1 + st2 * tx5
!
     tx3      = tem
     tx4      = tx4 + eta(l) * hol(l)
     tx5      = tx5 + gms(l) * hol(l)
!
     qil(l)   = max(0.0, min(1.0, (tcr-tol(l))*tcrf))
     qll(l+1) = (0.5*alhf) * st2 * (qil(l)+qil(l+1)) + one
   enddo
!
!     for the cloud top -- l=kd
!
   l = kd
!
   tem      = qst(l) - gaf(l) * hst(l)
   tem1     = (tx3 + tem) * 0.5
   st2      = (gaf(l)+gaf(l+1)) * 0.5
!
   fco(l+1) =            tem1 + st2 * hbl
   rnn(l+1) = rnn(l+1) * tem1 + st2 * tx4
   gmh(l+1) = gmh(l+1) * tem1 + st2 * tx5
!
   fco(l)   = tem + gaf(l) * hbl
   rnn(l)   = tem * rnn(l) + (tx4 + eta(l)*hol(l)) * gaf(l)
   gmh(l)   = tem * gmh(l) + (tx5 + gms(l)*hol(l)) * gaf(l)
!
!   replace fco for the bottom
!
   fco(kbl) = qbl
   qis      = qii(kd)
   qls      = qli(kd)
   qil(kd)  =  max(0.0, min(1.0, (tcr-tol(kd))*tcrf))
   qll(kd1) = (0.5*alhf) * st2 * (qil(kd) + qil(kd1)) + one
   qll(kd ) = alhf * gaf(kd) * qil(kd) + one
!
   do l = kd,kb1
     fco(l) = fco(l+1) - fco(l)
     rnn(l) = rnn(l+1) - rnn(l) + eta(l)*qol(l)
     gmh(l) = gmh(l+1) - gmh(l) + gms(l)*qol(l)
!
!    tem    = c0 * eta(l)
     tem    = (c0*(1.0-qil(l)) + c0i*qil(l)) * eta(l)
     bkc(l) = qll(l+1) - tem * (1.0-akt(l))
     akt(l) = qll(l)   + tem * akt(l)
     akc(l) = 1.0 / akt(l)
   enddo
!
   rnn(kbl) = 0.0
   tx3      = 0.0
   tx4      = 0.0
   tx5      = 0.0
   do l = kb1,kd1,-1
     tem    = bkc(l-1)       * akc(l)
     tx3    = (tx3 + fco(l)) * tem
     tx4    = (tx4 + rnn(l)) * tem
     tx5    = (tx5 + gmh(l)) * tem
   enddo
   if (kd .lt. kb1) then
     st1   = tx3 / bkc(kd)
     hsd   = hst(kd1)                                                          &
            + ltl(kd1) * (esfac *(qol(kd1)-qst(kd1))                           &
                     +  st1 * (1.0 + esfac*qol(kd1)) )                         &
            - (0.5*alhf) * st1 * (qil(kd1)+qil(kd))
   else
      hsd   = hbl
   endif
!
   tx3 = (tx3 + fco(kd)) * akc(kd)
   tx4 = (tx4 + rnn(kd)) * akc(kd)
   tx5 = (tx5 + gmh(kd)) * akc(kd)
   alm = alhf*qil(kd) - ltl(kd) * (1.0 + esfac*qol(kd))
!
   hsu = hst(kd) + ltl(kd) * esfac * (qol(kd)-qst(kd))
!
!===> vertical integrals needed to compute the entrainment parameter
!
   tx1 = alm * tx4
   tx2 = alm * tx5
!
   do l = kd,kb1
     tau = hol(l) - hsu
     tx1 = tx1 + tau * eta(l)
     tx2 = tx2 + tau * gms(l)
   enddo
!
!     modify hsu to include cloud liquid water and ice terms
!
   hsu   = hsu - alm * tx3
!
   clp   = zero
   alm   = -100.0
   hos   = hol(kd)
   qos   = qol(kd)
   unsat = hbl .gt. hsu .and. abs(tx1) .gt. 1.0e-4
!**********************************************************************!
   st1    = half*(hsu + hsd)
   if (unsat) then
!
!  standard case:
!   cloud can be neutrally bouyant at middle of level kd w/ +ve lambda.
!   epp < .25 is required to have real roots.
!
     st1 = one / tx1
     x00 = -(hbl-hsu) * st1
     epp = -x00 * tx2 * st1
     akm = epp
     epp = epp + epp
!
!        if((tx1.lt.zero .and. epp.lt..5) .or.
!    *      (tx1.gt.zero .and. epp.lt.0.0)      ) then
!          if (epp.eq.zero) epp = 1.e-8
!          sgn = sign(1.,tx1)
!          alm = x00*(one/epp)*(one+sgn*sqrt(one-epp-epp))
!          clp = one
!        endif
!
     if (epp .lt. 0.5) then
       if (epp.eq.zero) epp = 1.e-8
       tem1   = x00 / epp
       st2    = tem1 * sqrt(one - epp - epp)
       tem2   = tem1 - st2
       tem1   = tem1 + st2
       if (tem1 .gt. almax) tem1 = -100.0
       if (tem2 .gt. almax) tem2 = -100.0
       alm    = max(tem1, tem2)
       clp    = one
     endif
!
!
!  clip case:
!   non-entrainig cloud detrains in lower half of top layer.
!   no clouds are allowed to detrain below the top layer.
!
   elseif ( (hbl .le. hsu) .and.(hbl .gt. st1   )     ) then
     alm = zero
     clp = (hbl-st1) / (hsu-st1)
   endif
!
   unsat = .true.
   if (almin1 .gt. 0.0) then
     if (alm .ge. almin1) unsat = .false.
   else
     lowest   = kd .eq. kb1
     if ( (alm .gt. zero) .or.                                                 &
         (.not. lowest .and. alm .eq. zero) ) unsat = .false.
   endif
!
!===>  if no sounding meets second condition, return
!
   if (unsat) return 
!
   if(clp.gt.zero .and. clp.lt.one) then
     st1     = half*(one+clp)
     st2     = one - st1
     hst(kd) = hst(kd)*st1 + hst(kd1)*st2
     hos     = hol(kd)*st1 + hol(kd1)*st2
     qst(kd) = qst(kd)*st1 + qst(kd1)*st2
     qos     = qol(kd)*st1 + qol(kd1)*st2
!
     dlb(kd) = dlb(kd)*st1 + dlb(kd1)*st2
     ltl(kd) = ltl(kd)*st1 + ltl(kd1)*st2
     eta(kd) = eta(kd)*clp
     gms(kd) = gms(kd)*clp
   endif
!
   if (almin2 .ne. 0.0) then
     st1 = 0.0
     if (almin1 .ne. almin2) st1 = 1.0 / max(1.0e-10,(almin2-almin1))
       if (alm .lt. almin2) then
          clp = clp * (alm - almin1) * st1
       endif
   endif
!
!**********************************************************************!
!
!    critical workfunction is included in this version
!
   if (crtfun) then
     acr = 0.0
     tem = prl(kd1) - (prl(kd1)-prl(kd)) * clp * half
     call ras2_cwork_crit(tem, st1)
     acr = (prl(k) - tem) * st1
     clp = clp * qtl
   else
     acr = 0.0
     clp = clp * qtl
   endif
!
!===>  normalized massflux
!
!  eta is the thickness coming in and the mass flux going out.
!  gms is the thickness of the square; it is later reused for gamma_s
!
!     eta(k) = one
   do l = kb1,kd,-1
     eta(l) = eta(l+1) + alm * (eta(l) + alm * gms(l))
   enddo
!
!===>  cloud workfunction
!
   wfn   = zero
   akm   = zero
   det   = zero
   hcc   = hbl
   unsat = .false.
   qtl   = qst(kb1) - gaf(kb1)*hst(kb1)
   qtv   = qbl
   tx1   = hbl
!
   tem   = qst(kbl) - gaf(kbl)*hst(kbl)
   qtv   = 0.5 * ((tem+qtl) + (gaf(kbl)+gaf(kb1))*hbl)
   det   = max(0.0, qbl-qtv)
   qtv   = qbl - det
!
   do l = kb1,kd1,-1
     del_eta = eta(l) - eta(l+1)
     hccp = hcc + del_eta*hol(l)
!
     qtlp = qst(l-1) - gaf(l-1)*hst(l-1)
     qtvp = 0.5 * ((qtlp+qtl)*eta(l)                                           &
                + (gaf(l)+gaf(l-1))*hccp)
     detp = (bkc(l)*det - (qtvp-qtv) + del_eta*qol(l)) * akc(l)
!
     tem1   = akt(l)   - qll(l)
     tem2   = qll(l+1) - bkc(l)
     rnn(l) = tem1*detp  + tem2*det
!
     if (detp .le. zero) unsat = .true.
!
     st1  = hst(l) - ltl(l)*esfac*(qst(l)-qol(l))
     tem2 = hccp   + detp   * (qil(l-1)+qil(l)) * (alhf*0.5)
!
     st2  = ltl(l) * (1.0 + esfac*qol(l))
     tem5 = qli(l) + qii(l)
     tem3 = (tx1  - eta(l+1)*st1 - st2*(det-tem5))  * dlb(l)
     tem4 = (tem2 - eta(l  )*st1 - st2*(detp-tem5)) * dlt(l)

     st1  = tem3 + tem4

     wfn = wfn + st1       
     akm = akm - min(st1,zero)
!
!    buy(l) = eta(l+1)*tem3 + eta(l)*tem4
     buy(l) = 0.5 * (eta(l+1) + eta(l)) * st1
!
     hcc = hccp
     det = detp
     qtl = qtlp
     qtv = qtvp
     tx1 = tem2
   enddo

   del_eta = eta(kd) - eta(kd1)
   hccp    = hcc + del_eta*hos
!
   qtlp = qst(kd) - gaf(kd)*hst(kd)
   qtvp = qtlp*eta(kd) + gaf(kd)*hccp
   detp = (bkc(kd)*det - (qtvp-qtv) + del_eta*qos) * akc(kd)
!
   tem1    = akt(kd)  - qll(kd)
   tem2    = qll(kd1) - bkc(kd)
   rnn(kd) = tem1*detp  + tem2*det
!
   if (detp.le.zero) unsat = .true.
!
!#ifdef DBG
   st1 = half * (hst(kd)  - ltl(kd)*esfac*(qst(kd)-qos)                        &
       + hst(kd1) - ltl(kd1)*esfac*(qst(kd1)-qol(kd1)))
   st2 = half * (ltl(kd)  * (1.0+esfac*qol(kd))                                &
              +  ltl(kd1) * (1.0+esfac*qol(kd1)) )
!#endif
!     st1 = hst(kd)  - ltl(kd)*esfac*(qst(kd)-qos)
!     st2 = ltl(kd)  * (1.0+esfac*qol(kd))
!
!     st1 = half*(tx1-eta(kd1)*st1-st2*det)*dlb(kd)
   tem5 = half * (qli(kd1)+qii(kd1)+qls+qis)
   st1  = half * (tx1-eta(kd1)*st1-st2*(det-tem5))*dlb(kd)
!
   wfn = wfn + st1
   akm = akm - min(st1,zero)
!
!     buy(kd) = eta(kd1) * st1
   buy(kd) = 0.5 * (eta(kd1) + eta(kd)) * st1
!
   det = detp
   hcc = hccp
   akm = akm / wfn

!**********************************************************************!
!
!     if only to calculate workfunction save it and return
!
   if (wrkfun) then
     if (wfn .ge. 0.0) wfnc = wfn
     return
   elseif (.not. crtfun) then
     acr = wfnc
   endif
!
!===>  third check based on cloud workfunction
!
   calcup = .false.
!
   tem  =  min(cd*50.0, max_neg_bouy)
   if (wfn .gt. acr .and.  (.not. unsat) .and. akm .le. tem) then
     calcup = .true.
   endif
!
!===>  if no sounding meets third condition, return
!
   if (.not. calcup) return
!
   do l = kbl,k 
     rnn(l) = 0.0 
   enddo
!
!     if downdraft is to be invoked, do preliminary check to see
!     if enough rain is available and then call ras2_down_draft.
!
   ddft = .false.
   if (dndrft) then
     train = 0.0
     if (clp .gt. 0.0) then
       do l = kd,kb1
         train = train + rnn(l)
       enddo
     endif

     pl = (prl(kd1) + prl(kd))*half
     if (train .gt. 1.0e-4 .and. pl .le. pdd) ddft  = .true.
     do l = kd,km1
       if (l .lt. kbl-2 .and. ddft) then
         if (buy(l) .lt. 0.1) ddft = .false.
       endif
     enddo
   endif
!
   if (ddft) then
!
!     call downdraft scheme based on (cheng and arakawa, 1997)
!
     call ras2_down_draft(                                                     &
                 k, kd                                                         &
                ,cp, alhl, alhf, rkap, grav, tla, alfint                       &
                ,tol, qol, hol, prl, qst, hst, gam, gaf, hbl, qbl              &
                ,qrb, qrt, buy, kbl, idh, eta, rnn                             &
                ,alm, wfn, train, ddft                                         &
                ,etd, hod, qod, evp, dof, cldfrd, etz                          &
                ,gms, gsd, ghd)
!    *,             tx1, tx2, tx3, tx4, tx5, tx6, tx7, tx8, tx9)

   endif
!
!  no downdraft case (including case with no downdraft soln)
!
   if (.not. ddft) then
     do l = kd,k+1
       etd(l) = 0.0
       hod(l) = 0.0
       qod(l) = 0.0
     enddo
     do l = kd,k
       evp(l) = 0.0
       etz(l) = 0.0
     enddo
   endif
!
!===> calculate gammas  i.e. tendencies per unit cloud base massflux
!           includes downdraft terms!
!
   tx1 = qib
   tx2 = qlb
   do l = kb1,kd,-1
     del_eta = eta(l) - eta(l+1)
     tx1     = tx1 + del_eta * qii(l)
     tx2     = tx2 + del_eta * qli(l)
   enddo
!
   tem = det * qil(kd)
   st1 = (hcc+alhf*tem-eta(kd)*hst(kd)) / (1.0+gam(kd))
   ds  = eta(kd1) * (hos- hol(kd)                                              &
       - alhl*(qos - qol(kd)) + alhf*(qis-qii(kd)) )
   dh  = eta(kd1) * (hos- hol(kd))
   gms(kd) = (ds + st1) * pri(kd)
   gmh(kd) = pri(kd) * (hcc-eta(kd)*hos + dh)
!
!      tendency for suspended environmental ice and/or liquid water
!
   qil(kd) =     (tem + eta(kd1)*(qis-qii(kd))                                 &
                      + tx1 - eta(kd)*qis ) * pri(kd)
   qll(kd) = (det-tem + eta(kd1)*(qls-qli(kd))                                 &
                      + tx2 - eta(kd)*qls ) * pri(kd)
!
   ghd(kd) = 0.0
   gsd(kd) = 0.0
!
   do l = kd1,k
     st1 = one - alfint(l)
     st2 = one - alfinq(l)
     if (l .lt. kbl) then
       hb       = alfint(l)*hol(l-1) + st1*hol(l)
       qb       = alfint(l)*qol(l-1) + st1*qol(l)
!
       tem      = alfinq(l)*qii(l-1) + st2*qii(l)
       tem2     = alfinq(l)*qli(l-1) + st2*qli(l)
! 
       tem1     = eta(l) * (tem - qii(l))
       tem3     = eta(l) * (tem2 - qli(l))
!
       tem5     =  etd(l) * (hod(l) - hb)
       tem6     =  etd(l) * (qod(l) - qb)
!
       dh       = eta(l) * (hb - hol(l)) + tem5
       ds       = dh - alhl * (eta(l) * (qb - qol(l)) + tem6)
!
       gmh(l)   = dh * pri(l)
       gms(l)   = ds * pri(l)
!
       ghd(l)   = tem5 * pri(l)
       gsd(l)   = (tem5 - alhl * tem6) * pri(l)
!
       qil(l)   = tem1 * pri(l)
       qll(l)   = tem3 * pri(l)
!
       tem1     = eta(l) * (qii(l-1) - tem)
       tem3     = eta(l) * (qli(l-1) - tem2)
!
       dh       = eta(l) * (hol(l-1) - hb) - tem5
       ds       = dh - alhl * eta(l) * (qol(l-1) - qb)                         &
                    + alhl * (tem6 - evp(l-1))
!
       gmh(l-1) = gmh(l-1) + dh * pri(l-1)
       gms(l-1) = gms(l-1) + ds * pri(l-1)
!
       ghd(l-1) = ghd(l-1) - tem5 * pri(l-1)
       gsd(l-1) = gsd(l-1) - (tem5-alhl*(tem6-evp(l-1))) * pri(l-1)
!
       qil(l-1) = qil(l-1) + tem1 * pri(l-1)
       qll(l-1) = qll(l-1) + tem3 * pri(l-1)
     elseif (l .eq. kbl) then
        hb       = alfint(l)*hol(l-1) + st1*hbl
        qb       = alfint(l)*qol(l-1) + st1*qbl
!
!       hb       = hbl
!       qb       = qbl
        hbd      = alfint(l)*hol(l-1) + st1*hol(l)
        qbd      = alfint(l)*qol(l-1) + st1*qol(l)
!
        tem      = alfinq(l)*qii(l-1) + st2*qib
        tem2     = alfinq(l)*qli(l-1) + st2*qlb
!
        tem1     = eta(l) * (tem - qib)
        tem3     = eta(l) * (tem2 - qlb)
!
        tem5     =  etd(l) * (hod(l) - hbd)
        tem6     =  etd(l) * (qod(l) - qbd)
!
        tem4     = (grav/cmb2pa) * pris
        tx1      = eta(l) * (hb - hbl) * tem4
        tx2      = tx1 - alhl * eta(l) * (qb - qbl) * tem4
        dh       = tem5
!
        ds       =  dh - alhl * (tem6 + evp(l))
!
        gmh(l)   = tx1 + dh * pri(l)
        gms(l)   = tx2 + ds * pri(l)
!
        ghd(l)   = tem5 * pri(l)
        gsd(l)   = (tem5 - alhl * (tem6+evp(l))) * pri(l)
!
        qil(l)   = tem1 * tem4
        qll(l)   = tem3 * tem4
!
        tem1     = eta(l) * (qii(l-1) - tem)
        tem3     = eta(l) * (qli(l-1) - tem2)
!
        dh       = eta(l) * (hol(l-1) - hb) - tem5
        ds       = dh - alhl * eta(l) * (qol(l-1) - qb)                         &
                    + alhl * (tem6 - evp(l-1))
!
        gmh(l-1) = gmh(l-1) + dh * pri(l-1)
        gms(l-1) = gms(l-1) + ds * pri(l-1)
!
        ghd(l-1) = ghd(l-1) - tem5 * pri(l-1)
        gsd(l-1) = gsd(l-1) - (tem5-alhl*(tem6-evp(l-1)))                       &
                                     * pri(l-1)
!
        qil(l-1) = qil(l-1) + tem1 * pri(l-1)
        qll(l-1) = qll(l-1) + tem3 * pri(l-1)
     else
       hbd      = alfint(l)*hol(l-1) + st1*hol(l)
       qbd      = alfint(l)*qol(l-1) + st1*qol(l)
       tem5     =  etd(l) * (hod(l) - hbd)
       tem6     =  etd(l) * (qod(l) - qbd)
       dh       =  tem5
       ds       =  dh - alhl * (tem6 + evp(l))
!
       gmh(l)   = tx1 + dh * pri(l)
       gms(l)   = tx2 + ds * pri(l)
       ghd(l)   = dh * pri(l)
       gsd(l)   = ds * pri(l)
!
       dh       = - tem5
       ds       = dh  + alhl * tem6
       gmh(l-1) = gmh(l-1) + dh * pri(l-1)
       gms(l-1) = gms(l-1) + ds * pri(l-1)
!
       ghd(l-1) = ghd(l-1) + dh * pri(l-1)
       gsd(l-1) = gsd(l-1) + ds * pri(l-1)
!
       qil(l)   = qil(l-1)
       qll(l)   = qll(l-1)
     endif
   enddo
!
   hbd  = hol(k)
   qbd  = qol(k)
   tem5 =  etd(k) * (hod(k+1) - hbd)
   tem6 =  etd(k) * (qod(k+1) - qbd)
   dh   = - tem5
   ds   = dh  + alhl * tem6
   tem1 = dh * pri(k)
   tem2 = ds * pri(k)
   gmh(k) = gmh(k) + tem1
   gms(k) = gms(k) + tem2
   ghd(k) = ghd(k) + tem1
   gsd(k) = gsd(k) + tem2
!
   tem4   = - (grav/cmb2pa) * pris
   tx1    = dh * tem4
   tx2    = ds * tem4
!
   do l = kbl,k
     gmh(l) = gmh(l) + tx1
     gms(l) = gms(l) + tx2
     ghd(l) = ghd(l) + tx1
     gsd(l) = gsd(l) + tx2
   enddo
!
!**********************************************************************!
!
!
!===>  kernel (akm) calculation begins
!
!===>  modify sounding with unit mass flux
!
   testmb = 0.1
!
   do l = kd,k
     tem1   = gmh(l) - ghd(l)
     tem2   = gms(l) - gsd(l)
     hol(l) = hol(l) +  tem1*testmb
!        qol(l) = qol(l) + (tem1-tem2) * (testmb/alhl)
     qol(l) = qol(l) + (tem1-tem2                                              &
                          + alhf*qil(l))  * (testmb/alhl)
     hst(l) = hst(l) +  tem2*(one+gam(l))*testmb
     qst(l) = qst(l) +  tem2*gam(l)*(testmb/alhl)
   enddo
!
   hos = hos + (gmh(kd)-ghd(kd))  * testmb
   qos = qos + (gmh(kd)-gms(kd)                                                &
!    *          -  ghd(kd)+gsd(kd) ) * (testmb/alhl)
             -  ghd(kd)+gsd(kd)+ alhf*qil(kd)) * (testmb/alhl)
!
   tem = prl(k+1) - prl(k)
   hbl = hol(k) * tem
   qbl = qol(k) * tem
   do l = km1,kbl,-1
     tem = prl(l+1) - prl(l)
     hbl = hbl + hol(l) * tem
     qbl = qbl + qol(l) * tem
   enddo
   hbl = hbl * pris
   qbl = qbl * pris
!
!**********************************************************************!
!===>  cloud workfunction for modified sounding, then kernel (akm)
!
   akm = zero
   tx1 = zero
   qtl = qst(kb1) - gaf(kb1)*hst(kb1)
   qtv = qbl
   hcc = hbl
   tx2 = hcc
   tx4 = (alhf*0.5)*max(0.0,min(1.0,(tcr-tol(kb1))*tcrf))
!
   tem   = qst(kbl) - gaf(kbl)*hst(kbl)
   qtv   = 0.5 * ((tem+qtl) + (gaf(kbl)+gaf(kb1))*hbl)
   tx1   = max(0.0, qbl-qtv)
   qtv   = qbl - tx1
!

   do l = kb1,kd1,-1
     del_eta = eta(l) - eta(l+1)
     hccp = hcc + del_eta*hol(l)
!
     qtlp = qst(l-1) - gaf(l-1)*hst(l-1)
     qtvp = 0.5 * ((qtlp+qtl)*eta(l)                                           &
                   +(gaf(l)+gaf(l-1))*hccp)
     detp = (bkc(l)*tx1 - (qtvp-qtv) + del_eta*qol(l)) * akc(l)
!
     st1 = hst(l) - ltl(l)*esfac*(qst(l)-qol(l))
!
     tem2 = (alhf*0.5)*max(0.0,min(1.0,(tcr-tol(l-1))*tcrf))
     tem1 = hccp + detp * (tem2+tx4)
!
     st2  = ltl(l) * (1.0 + esfac * qol(l))
     tem5   = qli(l) + qii(l)
     akm  = akm +(  (tx2  -eta(l+1)*st1-st2*(tx1-tem5))  * dlb(l)              &
         + (tem1 -eta(l  )*st1-st2*(detp-tem5)) * dlt(l) )
!
      hcc  = hccp
      tx1  = detp
      tx2  = tem1
      qtl  = qtlp
      qtv  = qtvp
      tx4  = tem2
   enddo
!
!  eventhough we ignore the change in lambda, we still assume
!  that the cloud-top contribution is zero; as though we still
!  had non-bouyancy there.
!
!#ifdef DBG
   st1 = half * (hst(kd)  - ltl(kd)*esfac*(qst(kd)-qos)                        &
       +  hst(kd1) - ltl(kd1)*esfac*(qst(kd1)-qol(kd1))) 
   st2 = half * (ltl(kd)  * (1.0+esfac*qol(kd))                                &
              +  ltl(kd1) * (1.0+esfac*qol(kd1)) )
!#endif
!
!     st1 = hst(kd)  - ltl(kd)*esfac*(qst(kd)-qos)
!     st2 = ltl(kd)  * (1.0+esfac*qol(kd))
! 
   tem5 = half * (qli(kd1)+qii(kd1)+qls+qis)
   akm  = akm + half * (tx2-eta(kd1)*st1-st2*(tx1-tem5))*dlb(kd)
!
   akm = (akm - wfn) * (one/testmb)
!
!**********************************************************************!
!
!===>   mass flux
!
   amb = - (wfn-acr) / akm
!
!===>   relaxation and clipping factors
!     tem    = max(min((prl(kd1)-prl(kd))*(rasalf*0.02), 0.5), 0.1)
!     amb(i) = amb * clp * tem
!
   amb = amb * clp * rasalf
!    
!===>   sub-cloud layer depth limit on mass flux
!
   ambmax = (prl(kp1)-prl(kbl))*(cmb2pa*fracbl/grav)
   amb    = amax1(amin1(amb, ambmax),zero)
!
!**********************************************************************!
!*************************results**************************************!
!**********************************************************************!
!===>  precipitation and clw detrainment
!
   dsfc = dsfc + amb * etd(k) * (1.0/dt)
!
   do l = kbl,kd,-1
     pcu(l) = pcu(l) + amb*rnn(l)      !  (a40)
   enddo
!
!===> temparature and q change and cloud mass flux due to cloud type kd
!
   tx1 = amb * (one/cp)
   tx2 = amb * (one/alhl)
   do l = kd,k
     st1    = gms(l)*tx1
     toi(l) = toi(l) + st1
     tcu(l) = tcu(l) + st1
     tcd(l) = tcd(l) + gsd(l) * tx1
!
     flx(l) = flx(l) + eta(l)*amb
!
     qii(l) = qii(l) + qil(l) * amb
     tem    = 0.0
!
     qli(l) = qli(l) + qll(l) * amb + tem
!
     st1          = (gmh(l)-gms(l)+alhf*qil(l)) * tx2
!
     qoi(l) = qoi(l) + st1
     qcu(l) = qcu(l) + st1
     qcd(l) = qcd(l) + (ghd(l)-gsd(l)) * tx2
   enddo
!
   tx1 = 0.0
   tx2 = 0.0
!
!     reevaporation of falling convective rain
!
   if (revap) then
     gravcon = cmb2pa/grav
     afc     = -(1.04e-4*dt)*(3600./dt)**0.578
     rknob   = 5.0
!
     tx3    = amb*eta(kd)*pri(kd)
!
     cldfrd = min(amb*cldfrd, 1.0)
!
     do l = kd,k
!                                                 for l=kd,k
       if (l .ge. idh .and. ddft) then
         tx2 = tx2 + amb * rnn(l)
       else
         tx1 = tx1 + amb * rnn(l)
       endif
!
       clfrac = min(tx3*rknob, 1.0)
!
       if (tx1 .gt. 0. .or. tx2 .gt. 0.0) then
         teq     = toi(l)
         qeq     = qoi(l)
         pl      = 0.5 * (prl(l+1)+prl(l))
!
         st1     = max(0.0, min(1.0, (tcr-teq)*tcrf))
         st2     = st1*elfocp + (1.0-st1)*elocp
!
         call ras2_liquid_water ( teq,pl,qsteq,dqdt )
         deltaq = 0.5 * (qsteq-qeq) / (1.+st2*dqdt)
!
         qeq    = qeq + deltaq
         teq    = teq - deltaq*st2
!
         tem1   = max(0.0, min(1.0, (tcr-teq)*tcrf))
         tem2   = tem1*elfocp + (1.0-tem1)*elocp
!
         call ras2_liquid_water ( teq,pl,qsteq,dqdt )
         deltaq = (qsteq-qeq) / (1.+tem2*dqdt)
!
         qeq    = qeq + deltaq
         teq    = teq - deltaq*tem2
!
         if (qeq .gt. qoi(l)) then
           potevap = (qeq-qoi(l))*(prl(l+1)-prl(l))*gravcon
!
           tem3    = sqrt(pl*0.001)
           tem4    = potevap * (1. - exp( afc*sqrt(tx1*tem3) ) )
           actevap = min(tx1, tem4*clfrac)
!
           tem4    = potevap * (1. - exp( afc*sqrt(tx2*tem3) ) )
           tem4    = min(tx2, tem4*cldfrd)
!
           tx1     = tx1 - actevap
           tx2     = tx2 - tem4
           st1      = (actevap+tem4) * pri(l)
           qoi(l)  = qoi(l) + st1
           qcu(l)  = qcu(l) + st1
!

           st1     = st1 * elocp
           toi(l)  = toi(l) - st1 
           tcu(l)  = tcu(l) - st1
         endif
       endif
     enddo
!
   else
     do l = kd,k
       tx1 = tx1 + amb * rnn(l)
     enddo
   endif
!
   cup = cup + tx1 + tx2 + dof * amb
!
!    mixing of passive gases
!
   do n = 1,m
!
     do l = kd,k
       hol(l) = roi(l,n)
     enddo
!
     hcc     = rbl(n)
     hod(kd) = hol(kd)
!
!      compute downdraft properties for the tracer
!
     do l = kd1,k
       st1 = one - alfint(l)
       hb  = alfint(l)  * hol(l-1) + st1 * hol(l)
       if (etz(l-1) .ne. 0.0) then
         del_eta = etd(l) - etd(l-1)
         tem     = 1.0 / etz(l-1)
         if (del_eta .gt. 0.0) then
           hod(l) = (etd(l-1)*(hod(l-1)-hol(l-1))                              &
                   +  etd(l)  *(hol(l-1)-hb)+  etz(l-1)*hb) * tem
         else
           hod(l) = (etd(l-1)*(hod(l-1)-hb) + etz(l-1)*hb) * tem
         endif
       else
         hod(l) = hb
       endif
     enddo
!          
     do l = kb1,kd,-1
       hcc = hcc + (eta(l)-eta(l+1))*hol(l)
     enddo
!
     gmh(kd) = pri(kd) * (hcc-eta(kd)*hol(kd))
     do l = kd1,k
       st1 = one - alfint(l)
       if (l .lt. kbl) then
         hb       = alfint(l) * hol(l-1) + st1 * hol(l)
         tem5     = etd(l)    * (hod(l) - hb)
         dh       = eta(l)    * (hb - hol(l)) + tem5
         gmh(l  ) = dh * pri(l)
         dh       = eta(l)    * (hol(l-1) - hb) - tem5
         gmh(l-1) = gmh(l-1)  + dh * pri(l-1)
       elseif (l .eq. kbl) then
         hb       = alfint(l) * hol(l-1) + st1 * rbl(n)
         hbd      = alfint(l) * hol(l-1) + st1 * hol(l)
         dh       = etd(l)    * (hod(l) - hbd)
         tem4     = (grav/cmb2pa) * pris
         tx1      = eta(l)    * (hb - rbl(n)) * tem4
         gmh(l)   = tx1       + dh * pri(l)
         dh       = eta(l)    * (hol(l-1) - hb) - dh
         gmh(l-1) = gmh(l-1)  + dh * pri(l-1)
       else
         hbd      = alfint(l) * hol(l-1) + st1 * hol(l)
         dh       = etd(l)    * (hod(l) - hbd)
         gmh(l)   = tx1       + dh * pri(l)
         gmh(l-1) = gmh(l-1)  - dh * pri(l-1)
       endif
     enddo
!
     do l = kd,k
       st1      = gmh(l)*amb
       roi(l,n) = hol(l)    + st1
       rcu(l,n) = rcu(l,n) + st1
     enddo
   enddo                             ! tracer loop m
!
   return
   end subroutine ras2_cloud_main
!
!-------------------------------------------------------------------------------
   subroutine ras2_liquid_water(tt,p,q,dqdt)
!-------------------------------------------------------------------------------
!     implicit none
   real                 ::  tt, p, q, dqdt

   real,parameter       ::   airmw  = 28.97       
   real,parameter       ::   h2omw  = 18.01       
   real,parameter       ::   lice   = 2.834e6     
   real,parameter       ::   runiv  = 8314.3      
   real,parameter       ::   rgas   = runiv/airmw 
!

   real,parameter       ::   esfac = h2omw/airmw        
   real,parameter       ::   erfac = (1.0-esfac)/esfac  

   integer,parameter    ::  nx=7501
   real                 ::  tbsvp(nx), tbdsvp(nx)
   common/comsvp/ c1xsvp,c2xsvp,tbsvp, tbdsvp
!fpp$ noconcur r
   real d, qx, dqx
!-------------------------------------------------------------------------------
   xj  = min(max(c1xsvp+c2xsvp*tt,1.),float(nx))
   jx  = min(xj,nx-1.)
   qx  = tbsvp(jx)  + (xj-jx)*(tbsvp(jx+1)-tbsvp(jx))
   dqx = tbdsvp(jx) + (xj-jx)*(tbdsvp(jx+1)-tbdsvp(jx))
!
   d = (p-erfac*qx)
   if (d .gt. 0.) then
     d    = 1.0 / d
     q    = amin1(qx * d,1.0)
     dqdt = (1.0 + erfac*q) * d * dqx
   else
     q    = 1.0
     dqdt = 0.
   endif
!
   return
   end subroutine ras2_liquid_water
!
!-------------------------------------------------------------------------------
   subroutine ras2_cwork_crit(pl, acr)
!-------------------------------------------------------------------------------
   real                 ::  pl, acr
   integer              ::  iwk
   real                 ::  ac(16), ad(16)
   common /rasacr/ ac, ad
!-------------------------------------------------------------------------------
   iwk = pl * 0.02 - 0.999999999
   iwk = max(1, min(iwk,16))
   acr = ac(iwk) + pl * ad(iwk)
!
   return
   end subroutine ras2_cwork_crit
!
!-------------------------------------------------------------------------------   
   subroutine ras2_down_draft(                                                 &
                     k, kd                                                     &
                    ,cp, alhl, alhf, rkap, grav, tla, alfint                   &
                    ,tol, qol, hol, prl, qst, hst, gam, gaf, hbl, qbl          &
                    ,qrb, qrt, buy, kbl, idh, eta, rnn                         &
                    ,alm, wfn, train, ddft                                     &
                    ,etd, hod, qod, evp, dof, cldfrd, wcb                      &
                    ,gms, gsd, ghd)
!    *,                 tx1, tx2, tx3, tx4, tx5, tx6, tx7, tx8, tx9)
!-------------------------------------------------------------------------------   
!**********************************************************************!
!******************** cumulus downdraft subroutine ********************!
!****************** based on cheng and arakawa (1997)  ****** *********!
!********************* subroutine ras2_down_draft  ********************!
!************************  30 april 1998  *****************************!
!*******************  s. moorthi and m.j. suarez **********************!
!**********************************************************************!
!**********************************************************************!
!************* shrinivas.moorthi@noaa.gov (301) 763 8000(x7233) *******!
!***************  max.suarez@gsfc.nasa.gov (301) 286 7373 *************!
!**********************************************************************!
!**********************************************************************!
!
!===>  tol(k)     input   temperature            kelvin
!===>  qol(k)     input   specific humidity      non-dimensional
!===>  prl(k+1)   input   pressure @ edges       mb
!===>  k     input   the rise & the index of the subcloud layer
!===>  kd    input   detrainment level ( 1<= kd < k )          
!-------------------------------------------------------------------------------   
   implicit none
!-------------------------------------------------------------------------------   
!
!  input arguments
!
   integer              ::  k, kd
   real                 ::  cp, alhl, alhf, rkap, grav, alfint(k)
   integer              ::  kbl, kb1
!
   logical              ::  skpdd, skpup
!
   real                 ::  hol(kd:k),   qol(kd:k),   gaf(kd:k+1)
   real                 ::  hst(kd:k),   qst(kd:k),   tol(kd:k)
   real                 ::  buy(kd:k+1), qrb(kd:k),   qrt(kd:k)
   real                 ::  gam(kd:k+1), rnn(kd:k)
   real                 ::  eta(kd:k+1), prl(kd:k+1)
   real                 ::  hbl,     qbl,        pris
   real                 ::  train,   wfn,        alm
!
!  temporary work space
!
   real                 ::  gms(kd:k+1)
   real                 ::  tx1,    tx2,  tx3, tx4
   real                 ::  tx5,    tx6,  tx7, tx8, tx9
   logical              ::  unsat
!
   real                 ::  tl, pl, ql, qs, dqs, st1,  hb, qb, tb
   real                 ::  qqq, picon, piinv, del_eta
   real                 ::  tem, tem1, tem2, tem3, tem4, onebg, st2
   real                 ::  errh, errw, erre, tem5, tem6, hbd, qbd
   integer              ::  i, l,  n, ix, kd1, ii
   integer              ::  kp1, it, km1, ktem, kk, kk1, lm1, ll, lp1
   integer              ::  ip1
!
   real,parameter       ::  two=2.0,  half=0.5,  one=1.0, zero=0.0
   real,parameter       ::  cmb2pa=100.0  !  conversion mb to pascals
!
   real,parameter       ::  airmw  = 28.97
   real,parameter       ::  h2omw  = 18.01
   real,parameter       ::  esfac = airmw/h2omw-one
   real,parameter       ::  errmin=0.0001, errmi2=0.1*errmin
!     parameter (errmin=0.00001, errmi2=0.1*errmin)
!
   real                 ::  tla, stla,  ctl2, ctl3
   real                 ::  vtp, rgas, ctla, vtrm, vtpexp
   real                 ::  wcmin, wcbase, f2
   real                 ::  qraf, qrbf, cmpor
!
   real,parameter       ::  face=0.0, delx=10000.0, angf=0.5
   real,parameter       ::  onpg=1.0+0.5, gmf=1.0/onpg, rpart=0.0
!     parameter (aa1=1.0, bb1=1.5, cc1=1.1, dd1=0.85, f3=cc1, f5=2.5)
   real,parameter       ::  aa1=2.0, bb1=1.5, cc1=1.1, dd1=0.85, f3=cc1, f5=2.5
!     parameter (aa1=1.0, bb1=1.0, cc1=1.0, dd1=1.0, f3=cc1,  f5=1.0)
   real,parameter       ::  qrmin=1.0e-6,wc2min=0.01,gmf1=gmf/aa1,gmf5=gmf/f5
!     parameter (qrmin=1.0e-6, wc2min=1.00, gmf1=gmf/aa1, gmf5=gmf/f5)
!
   real,parameter       ::  pi=3.1415926535897931
   integer,parameter    ::  itrmu=14, itrmd=18, itrmin=7
!     parameter (itrmu=10, itrmd=10, itrmin=5)
   real                 ::  qrp(kd:k+1), wvl(kd:k+1), al2
!
   real                 ::  rnf(kd:k),   etd(kd:k+1), wcb(kd:k)
   real                 ::  hod(kd:k+1), qod(kd:k+1), evp(kd:k)
   real                 ::  ror(kd:k+1), stlt(kd:k)
   real                 ::  ghd(kd:k),   gsd(kd:k)
   real                 ::  rnt,        rnb
   real                 ::  errq,       cldfrd,     rntp
   integer              ::  idn, idh,itr
   real                 ::  elm(k)
!     real em(k*k), elm(k)
   real                 ::  edz, ddz, ce, qhs, fac, facg, asin, rsum1, rsum2, rsum3
   logical              ::  ddft, updret, ddlgk
!
   real                 ::  aa(kd:k,kd:k+1), qw(kd:k,kd:k)
   real                 ::  bud(kd:k), vt(2), vrw(2), trw(2)
   real                 ::  gqw(kd:k)
   real                 ::  qa(3),     wa(3),    dof, dofw
!    *,    gqw(kd:k), wcb(kd:k)
!   real                 ::  qrpf, vtpf
!-------------------------------------------------------------------------------   

!fpp$ expand (qrpf, vtpf)
!fpp$ noconcur r
   onebg  = 1.0 / grav
   kd1    = kd + 1
   kp1    = k  + 1
   km1    = k  - 1
   kb1    = kbl - 1
   rgas   = rkap * cp
   cmpor  = cmb2pa / rgas
   vtp    = 36.34*sqrt(1.2)* (0.001)**0.1364
   vtpexp = -0.3636
   piinv  = 1.0 / pi
   picon  = pi * onebg * 0.5
!
!     compute rain water budget of the updraft (cheng and arakawa, 1997)
!
   cldfrd = 0.0
   rntp   = 0.0
   dof    = 0.0
   errq   = 10.0
   rnb    = 0.0
   rnt    = 0.0
   tx2    = prl(kbl)
!
   tx1      = (prl(kd) + prl(kd1)) * 0.5
   ror(kd)  = cmpor*tx1 / (tol(kd)*(1.0+esfac*qol(kd)))
!     gms(kd)  = vtp * ror(kd) ** vtpexp
   gms(kd)  = vtp * vtpf(ror(kd))
!
   qrp(kd)  = qrmin
!
   tem        = tol(k) * (1.0 + esfac * qol(k))
   ror(k+1) = 0.5 * cmpor * (prl(k+1)+prl(k)) / tem
   gms(k+1) = vtp * vtpf(ror(k+1))
   qrp(k+1) = qrmin
   buy(kd)  = max(buy(kd), 0.1)
!
   do l = kd1,k
     tem = 0.5 * (tol(l)+tol(l-1))                                             &
         * (1.0 + 0.5*esfac * (qol(l)+qol(l-1)))
     ror(l) = cmpor * prl(l) / tem
!       gms(l) = vtp * ror(l) ** vtpexp
     gms(l) = vtp * vtpf(ror(l))
     qrp(l) = qrmin
     buy(l) = max(buy(l), 0.1)
   enddo
!
   call ras2_down_step1(tx1, alm, stla, ctl2, al2, pi, tla, tx2, wfn, tx3)
!
!    following ucla approach for rain profile
!
   f2      = 2.0*bb1*onebg/(pi*0.2)
   wcmin   = sqrt(wc2min)
   wcbase  = wcmin
!
   stla = f2     * stla * al2
   ctl2 = dd1    * ctl2
   ctl3 = 0.1364 * ctl2
!
   do l = kd,k
     rnf(l)   = 0.0
     wvl(l)   = 0.0
     stlt(l)  = 0.0
     gqw(l)   = 0.0
     do n = kd,k
       qw(n,l) = 0.0
     enddo
   enddo
!
!-----qw(n,l) = d(w(n)*w(n))/dqr(l)
!
   kk = kbl
   wvl(kk)    = wcbase
   qw(kd,kd)  = -qrb(kd)  * gmf1
   ghd(kd)    = eta(kd)   * eta(kd)
   gqw(kd)    = qw(kd,kd) * ghd(kd)
   gsd(kd)    = 1.0 / ghd(kd)
!
   gqw(kk)    = -  qrb(kk-1) * (gmf1+gmf1)
!
   wvl(kk)    = wcbase
   wcb(kk)    = wcbase * wcbase
   tx1        = wcb(kk)
   gsd(kk)    = 1.0
   ghd(kk)    = 1.0
!
   tem        = gmf1 + gmf1
   do l = kb1,kd1,-1
     ghd(l)  = eta(l) * eta(l)
     gsd(l)  = 1.0 / ghd(l)
     gqw(l)  = - ghd(l) * (qrb(l-1)+qrt(l)) * tem
     qw(l,l) = - qrt(l) * tem
!
     tx1     = tx1 + buy(l) * tem
     wcb(l)  = tx1 * gsd(l)
   enddo
!
   tem1        = (qrb(kd) + qrt(kd1) + qrt(kd1)) * gmf1
   gqw(kd1)    = - ghd(kd1) * tem1
!     qw(l,kd1)   = - qrt(kd1) * tem
   qw(kd1,kd1) = - qrt(kd1) * tem
   wcb(kd)     = (tx1 + buy(kd)*tem) * gsd(kd)
!
   do l = kd1,kbl
     do n = kd,l-1
        qw(n,l) = gqw(l) * gsd(n)
     enddo
   enddo
   qw(kbl,kbl) = 0.0
!
   wvl(kbl)    = wcbase
   stlt(kbl)   = 1.0 / wcbase
!
   do l = kd,k+1
     do n = kd,k
       aa(n,l) = 0.0
     enddo
   enddo
!
   skpup = .false.
!
   do itr = 1,itrmu               ! rain profile iteration starts!
     if (.not. skpup) then
!
!-----calculating the vertical velocity
!
       tx1      = 0.0
       ghd(kbl) = 1.0 / qrp(kbl)
       do l = kb1,kd,-1
         tx1     = tx1    + qrp(l+1) * gqw(l+1)
         st1     = wcb(l) + qw(l,l)  * qrp(l)+ tx1      * gsd(l)
         wvl(l)  = sqrt(max(st1,wc2min))
         stlt(l) = 1.0 / wvl(l)
         ghd(l)  = 1.0 / qrp(l)
       enddo
!
!-----calculating trw, vrw and of
!
!         vt(1)   = gms(kd) * qrp(kd)**0.1364
       vt(1)   = gms(kd) * qrpf(qrp(kd))
       trw(1)  = eta(kd) * qrp(kd) * stlt(kd)
       tx6     = trw(1) * vt(1)
       vrw(1)  = f3*wvl(kd) - ctl2*vt(1)
       bud(kd) = stla * tx6 * qrb(kd) * 0.5
       rnf(kd) = bud(kd)
       dof     = 1.1364 * bud(kd) * ghd(kd)
       dofw    = -bud(kd) * stlt(kd)
!
       rnt     = trw(1) * vrw(1)
       tx2     = 0.0
       tx4     = 0.0
       rnb     = rnt
       tx1     = 0.5
       tx8     = 0.0
!
       if (rnt .ge. 0.0) then
         tx3 = (rnt-ctl3*tx6) * ghd(kd)
         tx5 = ctl2 * tx6 * stlt(kd)
       else
         tx3 = 0.0
         tx5 = 0.0
         rnt = 0.0
         rnb = 0.0
       endif
!
       do l = kd1,kb1
         ktem    = max(l-2, kd)
         ll      = l - 1
!
!           vt(2)   = gms(l) * qrp(l)**0.1364
         vt(2)   = gms(l) * qrpf(qrp(l))
         trw(2)  = eta(l) * qrp(l) * stlt(l)
         vrw(2)  = f3*wvl(l) - ctl2*vt(2)
         qqq     = stla * trw(2) * vt(2)
         st1     = tx1  * qrb(ll)
         bud(l)  = qqq * (st1 + qrt(l))
!
         qa(2)   = dof
         wa(2)   = dofw
         dof     = 1.1364 * bud(l) * ghd(l)
         dofw    = -bud(l) * stlt(l)
!
         rnf(ll) = rnf(ll) + qqq * st1
         rnf(l)  =           qqq * qrt(l)
!
         tem3    = vrw(1) + vrw(2)
         tem4    = trw(1) + trw(2)
!
         tx6     = .25 * tem3 * tem4
         tem4    = tem4 * ctl3
!
!-----by qr above
!
!           tem1    = .25*(trw(1)*tem3 - tem4*vt(1))*tx7
         tem1    = .25*(trw(1)*tem3 - tem4*vt(1))*ghd(ll)
         st1     = .25*(trw(1)*(ctl2*vt(1)-vrw(2))                             &
                     * stlt(ll) + f3*trw(2))
!-----by qr below
         tem2    = .25*(trw(2)*tem3 - tem4*vt(2))*ghd(l)
         st2     = .25*(trw(2)*(ctl2*vt(2)-vrw(1))                             &
                    * stlt(l)  + f3*trw(1))
         !
         !      from top to  the kbl-2 layer
         !
         qa(1)   = tx2
         qa(2)   = qa(2) + tx3 - tem1
         qa(3)   = -tem2
!
         wa(1)   = tx4
         wa(2)   = wa(2) + tx5 - st1
         wa(3)   = -st2
!
         tx2     = tem1
         tx3     = tem2
         tx4     = st1
         tx5     = st2
!
         vt(1)   = vt(2)
         trw(1)  = trw(2)
         vrw(1)  = vrw(2)
!
         if (wvl(ktem) .eq. wcmin) wa(1) = 0.0
         if (wvl(ll)   .eq. wcmin) wa(2) = 0.0
         if (wvl(l)    .eq. wcmin) wa(3) = 0.0
         do n = ktem,kbl
           aa(ll,n) = (wa(1)*qw(ktem,n) * stlt(ktem)                           &
                    +  wa(2)*qw(ll,n)   * stlt(ll)                             &
                    +  wa(3)*qw(l,n)    * stlt(l) ) * 0.5
         enddo
         aa(ll,ktem) = aa(ll,ktem) + qa(1)
         aa(ll,ll)   = aa(ll,ll)   + qa(2)
         aa(ll,l)    = aa(ll,l)    + qa(3)
         bud(ll)     = (tx8 + rnn(ll)) * 0.5                                   &
                       - rnb + tx6 - bud(ll)
         aa(ll,kbl+1) = bud(ll)
         rnb = tx6
         tx1 = 1.0
         tx8 = rnn(ll)
       enddo
       l  = kbl
       ll = l - 1
!         vt(2)   = gms(l) * qrp(l)**0.1364
       vt(2)   = gms(l) * qrpf(qrp(l))
       trw(2)  = eta(l) * qrp(l) * stlt(l)
       vrw(2)  = f3*wvl(l) - ctl2*vt(2)
       st1     = stla * trw(2) * vt(2) * qrb(ll)
       bud(l)  = st1
!
       qa(2)   = dof
       wa(2)   = dofw
       dof     = 1.1364 * bud(l) * ghd(l)
       dofw    = -bud(l) * stlt(l)
!
       rnf(ll) = rnf(ll) + st1
!
       tem3    = vrw(1) + vrw(2)
       tem4    = trw(1) + trw(2)
!
       tx6     = .25 * tem3 * tem4
       tem4    = tem4 * ctl3
!
!-----by qr above
!
       tem1    = .25*(trw(1)*tem3 - tem4*vt(1))*ghd(ll)
       st1     = .25*(trw(1)*(ctl2*vt(1)-vrw(2))                               &
                   * stlt(ll) + f3*trw(2))
!-----by qr below
       tem2    = .25*(trw(2)*tem3 - tem4*vt(2))*ghd(l)
       st2     = .25*(trw(2)*(ctl2*vt(2)-vrw(1))                               &
                    * stlt(l)  + f3*trw(1))
!
!      for the layer next to the top of the boundary layer
!
       qa(1)   = tx2
       qa(2)   = qa(2) + tx3 - tem1
       qa(3)   = -tem2
!
       wa(1)   = tx4
       wa(2)   = wa(2) + tx5 - st1
       wa(3)   = -st2
!
       tx2     = tem1
       tx3     = tem2
       tx4     = st1
       tx5     = st2
!
       idn     = max(l-2, kd)
!
       if (wvl(idn) .eq. wcmin) wa(1) = 0.0
       if (wvl(ll)  .eq. wcmin) wa(2) = 0.0
       if (wvl(l)   .eq. wcmin) wa(3) = 0.0
!
       kk = idn
       do n = kk,l
         aa(ll,n) = (wa(1)*qw(kk,n) * stlt(kk)                                 &
                  +  wa(2)*qw(ll,n) * stlt(ll)                                 &
                  +  wa(3)*qw(l,n)  * stlt(l) ) * 0.5
       enddo
!
       aa(ll,idn) = aa(ll,idn) + qa(1)
       aa(ll,ll)  = aa(ll,ll)  + qa(2)
       aa(ll,l)   = aa(ll,l)   + qa(3)
       bud(ll)    = (tx8+rnn(ll)) * 0.5 - rnb + tx6 - bud(ll)
!
       aa(ll,l+1) = bud(ll)
!
       rnb        = trw(2) * vrw(2)
!
!      for the top of the boundary layer
!
       if (rnb .lt. 0.0) then
          kk    = kbl
          tem   = vt(2) * trw(2)
          qa(2) = (rnb - ctl3*tem) * ghd(kk)
          wa(2) = ctl2 * tem * stlt(kk)
       else
          rnb   = 0.0
          qa(2) = 0.0
          wa(2) = 0.0
       endif
!
       qa(1) = tx2
       qa(2) = dof + tx3 - qa(2)
       qa(3) = 0.0
!
       wa(1) = tx4
       wa(2) = dofw + tx5 - wa(2)
       wa(3) = 0.0
!
       kk = kbl
       if (wvl(kk-1) .eq. wcmin) wa(1) = 0.0
       if (wvl(kk)   .eq. wcmin) wa(2) = 0.0
!
       do ii = 1,2
          n = kk + ii - 2
          aa(kk,n) = (wa(1)*qw(kk-1,n) * stlt(kk-1)                            &
                   +  wa(2)*qw(kk,n)   * stlt(kk)) * 0.5
       enddo
       fac = 0.5
       ll  = kbl
       l   = ll + 1
       lm1 = ll - 1
       aa(ll,lm1)  = aa(ll,lm1) + qa(1)
       aa(ll,ll)   = aa(ll,ll)  + qa(2)
       bud(ll)     = 0.5*rnn(lm1) - tx6 + rnb - bud(ll)
       aa(ll,ll+1) = bud(ll)
!
!-----solving the budget equations for dqr
!
       do l = kd1,kbl
         lm1  = l - 1
         unsat = abs(aa(lm1,lm1)) .lt. abs(aa(l,lm1))
         do  n=lm1,kbl+1
           if (unsat) then
             tx1       = aa(lm1,n)
             aa(lm1,n) = aa(l,n)
             aa(l,n)   = tx1
           endif
         enddo
         tx1 = aa(l,lm1) / aa(lm1,lm1)
         do n = l,kbl+1
           aa(l,n) = aa(l,n) - tx1 * aa(lm1,n)
         enddo
       enddo     
!
!-----back substitution and check if the solution converges
!
       kk = kbl
       kk1 = kk + 1
       aa(kk,kk1) = aa(kk,kk1) / aa(kk,kk)      !   qr correction !
       tx2        = abs(aa(kk,kk1)) * ghd(kk)   !   error measure !
!
       kk = kbl + 1
       do l = kb1,kd,-1
         lp1   = l + 1
         tx1  = 0.0
         do n = lp1,kbl
            tx1  = tx1 + aa(l,n) * aa(n,kk)
         enddo
         aa(l,kk) = (aa(l,kk) - tx1) / aa(l,l)       ! qr correction !
         tx2      = max(tx2, abs(aa(l,kk))*ghd(l))   ! error measure !
       enddo
!
       do l = kd,kbl
         qrp(l) = max(qrp(l)+aa(l,kbl+1), qrmin)
       enddo
!
       if (itr .lt. itrmin) then
         tem = abs(errq-tx2) 
         if (tem .ge. errmi2 .and. tx2 .ge. errmin) then 
           errq  = tx2                              ! further iteration !
         else 
           skpup = .true.                           ! converges      !
           errq  = 0.0                              ! rain profile exists!
         endif 
       else
         tem = errq - tx2
         if (tem .lt. zero .and. errq .gt. 0.5) then
           skpup = .true.                           ! no convergence !
           errq = 10.0                              ! no rain profile!
         elseif (abs(tem).lt.errmi2 .or. tx2.lt.errmin) then
           skpup = .true.                           ! converges      !
           errq = 0.0                               ! rain profile exists!
         else
           errq = tx2                               ! further iteration !
         endif
       endif
!
     endif                                           ! skpup  endif!
!
   enddo                                          ! end of the itr loop!!
!
   if (errq .lt. 0.1) then
     ddft = .true.
     rnb  = - rnb
   else
     ddft = .false.
   endif
!
!     caution !! below is an adjustment to rain flux to maintain
!                conservation of precip!
!
   if (ddft) then
     tx1 = 0.0
     do l = kd,kb1
       tx1 = tx1 + rnf(l)
     enddo
     tx1 = train / (tx1+rnt+rnb)
     if (abs(tx1-1.0) .lt. 0.2) then
        rnt = max(rnt*tx1, 0.0)
        rnb = rnb * tx1
     else
        ddft = .false.
     endif
   endif
!
   dof = 0.0
   if (.not. ddft) return     ! rain profile did not converge!
!
   do l = kd,kb1
     rnf(l) = rnf(l) * tx1
   enddo
!
!     adjustment is over
!
!
!     downdraft
!
   do l = kd,k
     wcb(l) = 0.0
   enddo
!
   skpdd = .not. ddft
!
   errq  = 10.0
   if (.not. skpdd) then
!
!     calculate downdraft properties
!
     kk = max(kb1,kd1)
     do l = kk,k
       stlt(l) = stlt(l-1)
     enddo
     tem1 = 1.0 / bb1
     do l = kd,k
       evp(l) = 0.0
       if (l .lt. kbl) then
         tem     = stla * tem1
         stlt(l) = eta(l) * stlt(l) * tem / ror(l)
       else
         stlt(l) = 0.0
       endif
     enddo
!
     idn      = 99
     etd(kd)  = 0.0
     wvl(kd)  = 0.0
     hod(kd)  = hol(kd)
     qod(kd)  = qol(kd)
     tx1      = min(stlt(kd)*qrb(kd)*1.0, 1.0)    ! sigma at the top
     qrp(kd)  = 0.0
     rntp     = 0.0
     evp(kd)  = 0.0
     buy(kd)  = 0.0
     tx5      = tx1
!
!       here we assume rpart of detrained rain rnt goes to pd
!
     if (rnt .gt. 0.0) then
       rntp    = (1.0 - rpart) * rnt
       qrp(kd) = (rpart*rnt / (ror(kd)*tx1*gms(kd))) ** (1.0/1.1364)
       buy(kd) = - ror(kd) * tx1 * qrp(kd)
     endif
!
!     l-loop for the downdraft iteration from kd1 to k+1 (bottom surface)
!
!     bud(kd) = ror(kd)
     do l = kd1,k+1
!
       if (l .le. k) then
         st1   = 1.0 - alfint(l)
         wa(1) = alfint(l)*hol(l-1) + st1*hol(l)
         wa(2) = alfint(l)*qol(l-1) + st1*qol(l)
         wa(3) = alfint(l)*tol(l-1) + st1*tol(l)
         qa(2) = alfint(l)*hst(l-1) + st1*hst(l)
         qa(3) = alfint(l)*qst(l-1) + st1*qst(l)
       else
         wa(1) = hol(k)
         wa(2) = qol(k)
         wa(3) = tol(k)
         qa(2) = hst(k)
         qa(3) = qst(k)
       endif
!
       fac = 2.0
       if (l .eq. kd1) fac = 1.0
!
       facg    = fac * 0.5 * gmf5     !  12/17/97
!
       ddlgk   =  idn .eq. 99
       bud(kd) = ror(l)
!
       if (ddlgk) then
         tx1    = tx5
         wvl(l) = max(wvl(l-1), 0.1)
         qrp(l) = max(qrp(l-1),qrp(l))
!
!           vt(1)  = gms(l-1) * qrp(l-1) ** 0.1364
         vt(1)  = gms(l-1) * qrpf(qrp(l-1))
         rnt    = ror(l-1) * (wvl(l-1)+vt(1))*qrp(l-1)
!
         tem    = max(alm, 2.5e-4) * max(eta(l), 1.0)
         trw(1) = picon*tem*(qrb(l-1)+qrt(l-1))
         trw(2) = 1.0 / trw(1)
!
         vrw(1) = 0.5 * (gam(l-1) + gam(l))
         vrw(2) = 1.0 / (vrw(1) + vrw(1))
!
         tx4    =  (qrt(l-1)+qrb(l-1))*(onebg*fac*500.00)
!
         dofw   = 1.0 / (wa(3) * (1.0 + esfac*wa(2)))      !  1.0 / tvbar!
!
         etd(l) = etd(l-1)
         hod(l) = hod(l-1)
         qod(l) = qod(l-1)
!
         errq   = 10.0
!
         if (l .le. kbl) then
           tx3 = stlt(l-1) * qrt(l-1) * (0.5*fac)
           tx8 = stlt(l)   * qrb(l-1) * (0.5*fac)
           tx9 = tx8 + tx3
         else
           tx3 = 0.0
           tx8 = 0.0
           tx9 = 0.0
         endif
!
         tem  = wvl(l-1) + vt(1)
         if (tem .gt. 0.0) then
           tem1 = 1.0 / (tem*ror(l-1))
           tx3 = vt(1) * tem1 * ror(l-1) * tx3
           tx6 = tx1 * tem1
         else
           tx6 = 1.0
         endif
       endif
!
       if (l .eq. kd1) then
         if (rnt .gt. 0.0) then
           tem    = max(qrp(l-1),qrp(l))
           wvl(l) = tx1 * tem * qrb(l-1)*(facg*5.0)
         endif
         wvl(l) = max(0.01, wvl(l))
         trw(1) = trw(1) * 0.5
         trw(2) = trw(2) + trw(2)
       else
         if (ddlgk) evp(l-1) = evp(l-2)
       endif
!
!       no downdraft above level idh
!
       if (l .lt. idh) then
!
         etd(l)   = 0.0
         hod(l)   = wa(1)
         qod(l)   = wa(2)
         evp(l-1) = 0.0
         wvl(l)   = 0.0
         qrp(l)   = 0.0
         buy(l)   = 0.0
         tx5      = tx9
         errq     = 0.0
         rntp     = rntp + rnt * tx1
         rnt      = 0.0
         wcb(l-1) = 0.0
       endif
!
!      bud(kd) = ror(l)
!
!       iteration loop for a given level l begins
!
       do itr = 1,itrmd
!
         unsat =  ddlgk .and. (errq .gt. errmin)
         if (unsat) then
!
!             vt(1)  = gms(l) * qrp(l) ** 0.1364
           vt(1)  = gms(l) * qrpf(qrp(l))
           tem    =  wvl(l) + vt(1)
!
           if (tem .gt. 0.0) then
             st1    = ror(l) * tem * qrp(l) + rnt
             if (st1 .ne. 0.0) st1 = 2.0 * evp(l-1) / st1
             tem1   = 1.0 / (tem*ror(l))
             tem2   = vt(1) * tem1 * ror(l) * tx8
           else
             tem1   = 0.0
             tem2   = tx8
             st1    = 0.0
           endif
!
           tem = ror(l)*wvl(l) - ror(l-1)*wvl(l-1)
!
           tx5   = (tx1 -st1 + tem2 + tx3) / (1.0 + tem*tem1)
           tx5   = max(tx5, 0.0)
!
           tem1   = etd(l)
           etd(l) = ror(l) * tx5 * max(wvl(l), 0.0)
!
           del_eta = etd(l) - etd(l-1)
           tem       = del_eta * trw(2)
           tem2      = max(min(tem, 1.0), -1.0)
           if (abs(tem) .gt. 1.0 .and. etd(l) .gt. 0.0 ) then
             del_eta = tem2 * trw(1)
             etd(l)  = etd(l-1) + del_eta
           endif
           if (wvl(l) .gt. 0.0) tx5 = etd(l) / (ror(l)*wvl(l))
           erre  = etd(l) - tem1
!
           tem1 = sqrt(max((trw(1)+del_eta)*(trw(1)-del_eta),0.0))
           edz  = (0.5 + asin(tem2)*piinv)*del_eta + tem1*piinv
!
           ddz   = edz - del_eta
           wcb(l-1) = etd(l) + ddz
!
           tem1  = hod(l)
           if (del_eta .gt. 0.0) then
             qqq    = 1.0 / (etd(l) + ddz)
             hod(l) = (etd(l-1)*hod(l-1) + del_eta*hol(l-1)                    &
                                               + ddz*wa(1)) * qqq
             qod(l) = (etd(l-1)*qod(l-1) + del_eta*qol(l-1)                    &
                                               + ddz*wa(2)) * qqq
           else
             qqq    = 1.0 / (etd(l-1) + edz)
             hod(l) = (etd(l-1)*hod(l-1) + edz*wa(1)) * qqq
             qod(l) = (etd(l-1)*qod(l-1) + edz*wa(2)) * qqq
           endif
           errh  = hod(l) - tem1
           errq  = abs(errh/hod(l))  + abs(erre/max(etd(l),0.00001))
           dof   = ddz
           vt(2) = qqq
!
           ddz  = dof
           tem4 = qod(l)
           tem1 = vrw(1)
!
           qhs  = qa(3) + 0.5 * (gaf(l-1)+gaf(l))* (hod(l)-qa(2))
!
!                                           first iteration       !
!
           st2  = prl(l) * (qhs + tem1 * (qhs-qod(l)))
           tem2 = ror(l) * qrp(l)
           call ras2_down_step2(tem2,qraf,qrbf)
           tem6 = tx5 * (1.6 + 124.9 * qraf) * qrbf * tx4
!
           ce   = tem6 * st2 / ((5.4e5*st2 + 2.55e6)*(etd(l)+ddz))
!
           tem2   = - ((1.0+tem1)*(qhs+ce) + tem1*qod(l))
           tem3   = (1.0 + tem1) * qhs * (qod(l)+ce)
           tem    = max(tem2*tem2 - 4.0*tem1*tem3, 0.0)
           qod(l) = max(tem4, (- tem2 - sqrt(tem)) * vrw(2))
!
!                                            second iteration   !
!
           st2  = prl(l) * (qhs + tem1 * (qhs-qod(l)))
           ce   = tem6 * st2 / ((5.4e5*st2 + 2.55e6)*(etd(l)+ddz))
!
           tem2   = - ((1.0+tem1)*(qhs+ce) + tem1*tem4)
           tem3   = (1.0 + tem1) * qhs * (tem4+ce)
           tem    = max(tem2*tem2 - 4.0*tem1*tem3, 0.0)
           qod(l) = max(tem4, (- tem2 - sqrt(tem)) * vrw(2))
!
!                                              evaporation in layer l-1
!
           evp(l-1) = (qod(l)-tem4) * (etd(l)+ddz)
!
!                                              calculate pd (l+1/2)
!
           qa(1)    = tx1*rnt + rnf(l-1) - evp(l-1)
!
           if (etd(l) .gt. 0.0) then
             tem    = qa(1) / (etd(l)+ror(l)*tx5*vt(1))
             qrp(l) = max(tem, 0.0)
           elseif (tx5 .gt. 0.0) then
             qrp(l) = (max(0.0,qa(1)/(ror(l)*tx5*gms(l))))                     &
                                             ** (1.0/1.1364)
           else
             qrp(l) = 0.0
           endif
!
!                                              compute buoyancy
!
           tem1   = wa(3)+(hod(l)-wa(1)-alhl*(qod(l)-wa(2)))                   &
                                                     * (1.0/cp)
           tem1   = tem1 * (1.0 + esfac*qod(l))
           ror(l) = cmpor * prl(l) / tem1
           tem1   = tem1 * dofw
!!!        tem1   = tem1 * (1.0 + esfac*qod(l)) * dofw
!
           buy(l) = (tem1 - 1.0 - qrp(l)) * ror(l) * tx5
!
!                                              compute w (l+1/2)
!
           tem1   = wvl(l)
           if (etd(l) .gt. 0.0) then
             wvl(l) = vt(2) * (etd(l-1)*wvl(l-1) - facg                        &
                    * (buy(l-1)*qrt(l-1)+buy(l)*qrb(l-1)))
           endif
!
           wvl(l) = 0.5 * (wvl(l) + tem1)
!
           errw   = wvl(l) - tem1
!
           errq   = errq + abs(errw/max(wvl(l), 0.00001))
!
           if (itr .ge. min(itrmin,itrmd/2)) then
             if (etd(l-1) .eq. 0.0 .and. errq .gt. 0.2) then
               ror(l)   = bud(kd)
               etd(l)   = 0.0
               wvl(l)   = 0.0
               errq     = 0.0
               hod(l)   = wa(1)
               qod(l)   = wa(2)
               tx5      = tx1 + tx9
!
               evp(l-1) = 0.0
               tem      = max(tx1*rnt+rnf(l-1), 0.0)
               qa(1)    = tem - evp(l-1)
               qrp(l)   = (qa(1) / (ror(l)*tx5*gms(l)))                        &
                                               ** (1.0/1.1364)
               buy(l)   = - ror(l) * tx5 * qrp(l)
               wcb(l-1) = 0.0
             endif
!
             del_eta = etd(l) - etd(l-1)
             if(del_eta .lt. 0.0 .and. errq .gt. 0.1) then
               ror(l)   = bud(kd)
               etd(l)   = 0.0
               wvl(l)   = 0.0
               tx5      = tx1 + tx9
               cldfrd   = tx5
!
               del_eta  = - etd(l-1)
               edz      = 0.0
               ddz      = -del_eta
               wcb(l-1) = ddz

!
               hod(l)   = hod(l-1)
               qod(l)   = qod(l-1)

!
               tem4     = qod(l)
               tem1     = vrw(1)
!
               qhs      = qa(3) + 0.5 * (gaf(l-1)+gaf(l))                      &
                                      * (hod(l)-qa(2))
!
!                                           first iteration       !
!
               st2  = prl(l) * (qhs + tem1 * (qhs-qod(l)))
               tem2 = ror(l) * qrp(l-1)
!
               call ras2_down_step2(tem2,qraf,qrbf)
!
               tem6 = tx5 * (1.6 + 124.9 * qraf) * qrbf * tx4
!
               ce   = tem6*st2/((5.4e5*st2 + 2.55e6)*(etd(l)+ddz))
!
               tem2   = - ((1.0+tem1)*(qhs+ce) + tem1*qod(l))
               tem3   = (1.0 + tem1) * qhs * (qod(l)+ce)
               tem    = max(tem2*tem2 - 4.0*tem1*tem3, 0.0)
               qod(l) = max(tem4, (- tem2 - sqrt(tem)) * vrw(2))
!
!                                            second iteration   !
!
               st2  = prl(l) * (qhs + tem1 * (qhs-qod(l)))
               ce   = tem6*st2/((5.4e5*st2 + 2.55e6)*(etd(l)+ddz))
!
               tem2   = - ((1.0+tem1)*(qhs+ce) + tem1*tem4)
               tem3   = (1.0 + tem1) * qhs * (tem4+ce)
               tem    = max(tem2*tem2 - 4.0*tem1*tem3, 0.0)
               qod(l) = max(tem4, (- tem2 - sqrt(tem)) * vrw(2))
!
!                                              evaporation in layer l-1
!
               evp(l-1) = (qod(l)-tem4) * (etd(l)+ddz)
!
!                                               calculate pd (l+1/2)
!
               rnn(l-1) = tx1*rnt + rnf(l-1) - evp(l-1)
 
               errq = 0.0
             endif
           endif
         endif
!
       enddo                ! end of the iteration loop  for a given l!
!
       if (l .le. k) then
         if (etd(l-1) .eq. 0.0                                                 &
            .and. errq .gt. errmin*10.0) then
           ror(l)   = bud(kd)
           tx5      = tx1 + tx9
           evp(l-1) = 0.0
           qrp(l)   = 0.0
           qa(1)    = tx1*rnt + rnf(l-1) - evp(l-1)
           if (qa(1) .gt. 0.0) then
             qrp(l) = (qa(1) / (ror(l)*tx5*gms(l)))                            &
                                            ** (1.0/1.1364)
           endif
           etd(l)   = 0.0
           wvl(l)   = 0.0
           st1      = 1.0 - alfint(l)
           hod(l)   = wa(1)
           qod(l)   = wa(2)
!
           errq     = 0.0
           buy(l)   = - ror(l) * tx5 * qrp(l)
           wcb(l-1) = 0.0
         endif
       endif
!
       ll = min(idn, k+1)
       if (errq .lt. 1.0 .and. l .le. ll) then
         if (etd(l-1) .gt. 0.0 .and. etd(l) .eq. 0.0) then
           idn = l
         endif
       endif
! 
!     if downdraft properties are not obtainable, (i.e.solution does
!      not converge) , no downdraft is assumed
!
       if (errq .gt. errmin*10.0 .and. idn .eq. 99)                            &
                             ddft = .false.
!
       dof = 0.0
       if (.not. ddft) return
!
     enddo                      ! end of the l loop of downdraft !
!
     tx1 = 0.0
!
     dof = qa(1)
!
   endif                       ! skpdd endif
!
   rnn(kd) = rntp
   tx1     = evp(kd)
   tx2     = rntp + rnb
!
   ii = idh
!
   if (ii .ge. kd1+1) then
       rnn(kd)   = rnn(kd) + rnf(kd)
       tx2       = tx2 + rnf(kd)
       rnn(ii-1) = 0.0
       tx1       = evp(ii-1)
   endif
!
   do l = kd,k
     ii = idh
     if (l .gt. kd1 .and. l .lt. ii) then
       rnn(l-1) = rnf(l-1)
       tx2      = tx2 + rnn(l-1)
     elseif (l .ge. ii .and. l .lt. idn) then
       rnn(l)   = 0.0
       tx1      = tx1 + evp(l)
     elseif (l .gt. idn) then
       etd(l)   = 0.0
       hod(l)   = 0.0
       qod(l)   = 0.0
       evp(l-1) = 0.0
       rnn(l-1) = rnf(l-1)
       tx2      = tx2 + rnn(l-1)
     endif
   enddo
!
   if (k+1 .gt. idn) then
     etd(k+1) = 0.0
     hod(k+1) = 0.0
     qod(k+1) = 0.0
     evp(k)   = 0.0
     rnn(k)   = rnf(k)
   endif
!
!      for downdraft case the rain is that falls thru the bottom
!
   l = kbl
!
   rnn(l)  = rnb
   cldfrd  = tx5
!
   tx2     = tx2 + dof
!
!     caution !! below is an adjustment to rain flux to maintain
!                conservation of precip!
!
   if (tx1 .gt. 0.0) then
     tx1 = (train - tx2) / tx1
   else
     tx1 = 0.0
   endif
!
   tx5      = evp(kbl)
   evp(kbl) = evp(kbl) * tx1
   tx3      = rnn(kbl) + evp(kbl) + dof
   tx2      = rnn(kbl)
   tx4      = evp(kbl)
!
   do l = kd,kb1
     tx5    = tx5 + evp(l)
     evp(l) = evp(l) * tx1
     tx3    = tx3 + evp(l) + rnn(l)
     tx2    = tx2 + rnn(l)
     tx4    = tx4 + evp(l)
   enddo
!
   return
   end subroutine ras2_down_draft
!
!-------------------------------------------------------------------------------
   subroutine ras2_down_step1( pres, alm, stla, ctl2, al2                      &
                    , pi, tla, prb, wfn, ufn)
!-------------------------------------------------------------------------------
   real                 ::  refp(6), refr(6), tlac(8), plac(8)
   real                 ::  pres, stla, ctl2
   real                 ::  alm, al2, tla, tem, tem1
   real                 ::  prb, acr, wfn, ufn
!
   data plac/100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0/
   data tlac/ 25.0,  22.5,  20.0,  17.5,  15.0,   12.5,   10.0, 10.0/
!
   data refp/500.0, 300.0, 250.0, 200.0, 150.0, 100.0/
   data refr/ 0.25,   0.5,   0.75,  1.0,  1.5,  2.0/
!
!-------------------------------------------------------------------------------
!fpp$ expand (ras2_cwork_crit)
   tem = pi / 180.0
!
   if (tla .lt. 0.0) then
     if (pres .le. plac(1)) then
       tla = tlac(1)
     elseif (pres .le. plac(2)) then
       tla = tlac(2) + (pres-plac(2))*((tlac(1)-tlac(2))/                      &
                                              (plac(1)-plac(2)))
     elseif (pres .le. plac(3)) then
       tla = tlac(3) + (pres-plac(3))*((tlac(2)-tlac(3))/                      &
                                              (plac(2)-plac(3)))
     elseif (pres .le. plac(4)) then
       tla = tlac(4) + (pres-plac(4))*((tlac(3)-tlac(4))/                      &
                                              (plac(3)-plac(4)))
     elseif (pres .le. plac(5)) then
       tla = tlac(5) + (pres-plac(5))*((tlac(4)-tlac(5))/                      &
                                              (plac(4)-plac(5)))
     elseif (pres .le. plac(6)) then
       tla = tlac(6) + (pres-plac(6))*((tlac(5)-tlac(6))/                      &
                                              (plac(5)-plac(6)))
     elseif (pres .le. plac(7)) then
       tla = tlac(7) + (pres-plac(7))*((tlac(6)-tlac(7))/                      &
                                              (plac(6)-plac(7)))
     elseif (pres .le. plac(8)) then
       tla = tlac(8) + (pres-plac(8))*((tlac(7)-tlac(8))/                      &
                                              (plac(7)-plac(8)))
     else
       tla = tlac(8)
     endif
     stla = sin(tla*tem)
     tem1 = cos(tla*tem)
     ctl2 = tem1 * tem1
   else
     stla = sin(tla*tem)
     tem1 = cos(tla*tem)
     ctl2 = tem1 * tem1
   endif
!
   if (pres .ge. refp(1)) then
     tem = refr(1)
   elseif (pres .ge. refp(2)) then
     tem = refr(1) + (pres-refp(1)) *                                          &
                      ((refr(2)-refr(1))/(refp(2)-refp(1)))
   elseif (pres .ge. refp(3)) then
     tem = refr(2) + (pres-refp(2)) *                                          &
                      ((refr(3)-refr(2))/(refp(3)-refp(2)))
   elseif (pres .ge. refp(4)) then
     tem = refr(3) + (pres-refp(3)) *                                          &
                      ((refr(4)-refr(3))/(refp(4)-refp(3)))
   elseif (pres .ge. refp(5)) then
     tem = refr(4) + (pres-refp(4)) *                                          &
                      ((refr(5)-refr(4))/(refp(5)-refp(4)))
   elseif (pres .ge. refp(6)) then
     tem = refr(5) + (pres-refp(5)) *                                          &
                      ((refr(6)-refr(5))/(refp(6)-refp(5)))
   else
     tem = refr(6)
   endif
   al2 = max(alm, 2.0e-4/tem)
!
   return
   end subroutine ras2_down_step1
!
!-------------------------------------------------------------------------------
   subroutine ras2_down_step2(qrp,qraf,qrbf)
!-------------------------------------------------------------------------------
   real                 ::  qrp, qraf, qrbf, xj, c1xqrp, c2xqrp
   integer              ::  jx
   integer,parameter    ::  nqrp=500001
   real                 ::  tbqrp(nqrp), tbqra(nqrp), tbqrb(nqrp)
   common/comqrp/ c1xqrp,c2xqrp,tbqrp,tbqra,tbqrb
!-------------------------------------------------------------------------------
   xj   = min(max(c1xqrp+c2xqrp*qrp,1.),float(nqrp))
   jx   = min(xj,nqrp-1.)
   xj   = xj - jx
   qraf = tbqra(jx)  + xj * (tbqra(jx+1)-tbqra(jx))
   qrbf = tbqrb(jx)  + xj * (tbqrb(jx+1)-tbqrb(jx))
!-------------------------------------------------------------------------------
   return
   end subroutine ras2_down_step2
!
!-------------------------------------------------------------------------------   
   function qrpf(qrp)
!-------------------------------------------------------------------------------
   real                 ::  qrp, qrpf, xj, c1xqrp, c2xqrp
   integer              ::  jx
   integer,parameter    ::  nqrp=500001
   real                 ::  tbqrp(nqrp), tbqra(nqrp), tbqrb(nqrp)
   common/comqrp/ c1xqrp,c2xqrp,tbqrp,tbqra,tbqrb
!-------------------------------------------------------------------------------
   xj   = min(max(c1xqrp+c2xqrp*qrp,1.),float(nqrp))
   jx   = min(xj,nqrp-1.)
   qrpf = tbqrp(jx)  + (xj-jx) * (tbqrp(jx+1)-tbqrp(jx))
! 
   return
   end function qrpf
!
!-------------------------------------------------------------------------------
   function vtpf(ror)
!-------------------------------------------------------------------------------
   real                ::  ror, vtpf, xj, c1xvtp, c2xvtp
   integer             ::  jx
   integer, parameter  ::  nvtp=10001
   real                ::  tbvtp(nvtp)
   common/comvtp/ c1xvtp,c2xvtp,tbvtp
!-------------------------------------------------------------------------------
   xj   = min(max(c1xvtp+c2xvtp*ror,1.),float(nvtp))
   jx   = min(xj,nvtp-1.)
   vtpf = tbvtp(jx)  + (xj-jx) * (tbvtp(jx+1)-tbvtp(jx))
!   
   return
   end function vtpf
!
!-------------------------------------------------------------------------------   
#include "abort.h"
   subroutine ras2_setup(lm, si, sl, cp, rgas, dt, nsphys, fh                  &
                    ,sig, sgb, prj                                             &
                    ,rasal, krmin, krmax, kfmax, ncrnd, rannum                 &
                    ,mct, kctop, dtf, dsfc, len)
!-------------------------------------------------------------------------------
!
!     this version of ras2_setup assumes same number of levels in the vertical
!     but levels are reverse order
!
!-------------------------------------------------------------------------------
   integer              ::  lm,nsphys,krmin, krmax, kfmax, ncrnd,mct
   real                 ::  si(lm+1),  sl(lm),    dsfc(len)
   real                 ::  sig(lm+1), sgb(lm),   prj(lm+1)
   real                 ::  rasal(lm), rannum(200,nsphys)
!
   integer              ::  kctop(mct+1)
   data pctp/70.0/, pcbt/900.0/, pinc/5.0/
!
   logical              ::  first
!-------------------------------------------------------------------------------
!fpp$ noconcur r
   call ras2_setup_sub1(lm, si, sl, cp, rgas, sig, sgb, prj)
   call ras2_setup_sub2(lm, dt, nsphys, fh,  sig                               &
              ,rasal, krmin, krmax, kfmax, ncrnd, rannum                       &
              ,mct, pctp, pcbt, kctop, dtf, dsfc, len)
!
   return
   end subroutine ras2_setup
!
!-------------------------------------------------------------------------------
   subroutine ras2_setup_sub1(lm, si, sl, cp, rgas, sig, sgb, prj)
!-------------------------------------------------------------------------------
!
!     this version of ras_setup assumes same number of levels in the vertical
!     but levels are reverse order
!
!-------------------------------------------------------------------------------
   real  ::  si(lm+1),  sl(lm)
   real  ::  sig(lm+1), sgb(lm),   prj(lm+1)
!
   logical first
   data first/.true./
   save first
!fpp$ noconcur r
!
   if (first) then
     rkap = rgas / cp
     lmm1 = lm - 1
     do l = 1,lm
       sig(l) = si(lm-l+2)
       prj(l) = (0.001*sig(l)) ** rkap
       sgb(l) = sl(lm-l+1)
     enddo
     sig(lm+1) = si(1)
     prj(lm+1) = (0.001*sig(lm+1)) ** rkap
!
   endif
!
   return
   end subroutine ras2_setup_sub1
!
!-------------------------------------------------------------------------------
   subroutine ras2_setup_sub2(lm, dt, nsphys, fh,  prs                         &
                    ,rasal, krmin, krmax, kfmax, ncrnd, rannum                 &
                    ,mct, pctp, pcbt, kctop, dtf, dsfc, len)
!-------------------------------------------------------------------------------
!
!     this version of ras_setup assumes same number of levels in the vertical
!     but levels are reverse order
!
!-------------------------------------------------------------------------------
#ifndef DFS
   use comio, only   : iope
#else
   use dfsvar, only  : iope
#endif
   use phys_ras_module, only : ras_setup_sub
!-------------------------------------------------------------------------------
   real, parameter      ::  actp=1.7,   facm=1.00
   real                 ::  dsfc(len), sig(lm+1), prs(lm+1)
   real                 ::  rasal(lm), rannum(200,nsphys)
   real                 ::  ph(15),    a(15),     ac(16), ad(16)
   integer, allocatable,save  ::  nrnd(:)
!
   common /rasacr/ ac, ad
!
   data ph/150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0              &
          ,550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/
!
   data a/ 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677                      &
          ,0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664                      &
          ,0.0553, 0.0445, 0.0633/
!
   real                 ::  pctop(mct+1), sgh(mct+1)
   integer              ::  kctop(mct+1)
!
   real,save            ::  pinc
   data pinc/5.0/
!  data pctp/70.0/, pcbt/900.0/, pinc/5.0/
!
   real                 ::  seed
   integer,save         ::  krsize,iseed
   logical,save         ::  first
   real,save            ::  fh0
   data iseed/0/, first/.true./, fh0/0.0/
!-------------------------------------------------------------------------------
!fpp$ noconcur r
!
   if (first) then
     first = .false.
!                                   set critical workfunction arrays
     actop = actp*facm
     do l = 1,15
       a(l) = a(l)*facm
     enddo
     do l = 2,15
       tem   = 1.0 / (ph(l) - ph(l-1))
       ac(l) = (ph(l)*a(l-1) - ph(l-1)*a(l)) * tem
       ad(l) = (a(l) - a(l-1)) * tem
     enddo
     ac(1)  = actop
     ac(16) = a(15)
     ad(1)  = 0.0
     ad(16) = 0.0
!
     call ras_setup_sub
     call ras2_setup_param1
     call ras2_setup_param2
!                                   set other parameters
     lmm1 = lm - 1
     tem  = 1.0 / prs(lm+1)
     do l = 1,lm+1
       sig(l) = prs(l) * tem
     enddo
!
     tem  = (pcbt - pctp - (mct-1)*mct*pinc/2) / mct
     pctop(1) = pctp * 0.001
     do i = 1,mct
       pctop(i+1) = 0.001 * (pctop(i)*1000.0 + tem + pinc*(i-1))
     enddo
     write(6,*)' pctop=',pctop
     do l = 1,lm
       if (sig(l) .le. pctop(1)) kctop(1) = l
     enddo
     tem2 = 0.0
     do i = 2,mct+1
       if (kctop(i-1) .lt. lm) then
         ii       = kctop(i-1) + 1
         kctop(i) = ii
         do l = ii,lm
           if (sig(l) .le. pctop(i)) kctop(i) = l
         enddo
         if (pctop(i) .gt. sig(kctop(i))) then
           tem1 = sig(kctop(i)) - sig(kctop(i-1))
           if (tem1 .lt. tem2) kctop(i) = min(lm, kctop(i)+1)
           tem2 = sig(kctop(i)) - sig(kctop(i-1))
         endif
       else
         kctop(i) = lm + 1
       endif
     enddo
     do i = 1,mct+1
       if (kctop(i) .le. lm) mctm = i
     enddo
!
!    do i = 1,mctm
!      sgh(i) = sig(kctop(i))
!    enddo
!    do i = 1,mct
!      sgc(i)  = 0.5 * (sgh(i) + sgh(i+1))
!      dsgh(i) = sgh(i+1) - sgh(i)
!    enddo
     write(6,*)' mct=',mct,' mctm=',mctm
     write(6,*)' kctop=',(kctop(i),i=1,mctm)
!    write(6,*)' sgh=',(sig(kctop(i)),i=1,mctm)
!    write(6,*)' sgc=',sgc
!    write(6,*)' dsgh=',dsgh
!
     do i = 1,mctm-1
       if (kctop(i) .eq. kctop(i+1)) then
         write(6,*)' for cloud-type = ',i, ' kctop same at top',               &
                  ' and bottom - run stopped'
         call MPABORT
       endif
     enddo
!
     krmin = 1
     krmax = mctm - 1
     kfmax = mctm - 1
     do l = 1,mctm-1
       sgc  = 0.5 * (sig(kctop(l)) + sig(kctop(l+1)))
       if (sgc .le. 0.760) krmax = l
       if (sgc .le. 0.930) kfmax = l
     enddo
!
!        write(6,*)' krmin=',krmin,' krmax=',krmax,' kfmax=',kfmax
!
     do i = 1,len
       dsfc(i) = 0.0
     enddo
!
!    rasalf  = 0.25
     rasalf  = 0.30
!
     do l = 1,lm
       rasal(l) = rasalf
     enddo
     call random_seed(size=krsize)
     write(6,*)' krsize=',krsize
     allocate (nrnd(krsize))
!
   endif
!
!     compute ncrnd and the random numbers
!
!        ncrnd   = 42 * (dt/3600) + 0.50001
!
   ncrnd   = (krmax-krmin+1) * (dtf/1200) + 0.50001
   if (dt .gt. dtf) ncrnd = (5*ncrnd) / 4
!        ncrnd   = (krmax-krmin+1) * (dt/1200) + 0.50001
   ncrnd   = max(ncrnd, 1)
!
!        write(6,*)' ncrnd=',ncrnd,' dt=',dt,' dtf=',dtf,  &
!                ' rasalf=',rasalf
!
   if (fh .eq. fh0) then
     iseed = iseed + 1
   else
     iseed = fh*3600.0 + 0.0001
     iseed = max(iseed, 1)
     fh0   = fh
   endif
!
   nrnd = iseed
   call random_seed(put=nrnd)
!
   if (iope) then
     call random_number(rannum)
   endif                
#ifdef MP
#ifndef RMP
   call mpbcastr(rannum,200*nsphys)
#else
   call rmpbcastr(rannum,200*nsphys)
#endif    
#endif
!
   return
   end subroutine ras2_setup_sub2
!
!-------------------------------------------------------------------------------
   subroutine ras2_setup_param1
!-------------------------------------------------------------------------------
   integer,parameter    ::  nqrp=500001
   real                 ::  tbqrp(nqrp), tbqra(nqrp), tbqrb(nqrp)
   common/comqrp/ c1xqrp,c2xqrp,tbqrp,tbqra,tbqrb
!-------------------------------------------------------------------------------
!fpp$ noconcur r
!     xmin=1.0e-6
   xmin=0.0
   xmax=5.0
   xinc=(xmax-xmin)/(nqrp-1)
   c1xqrp=1.-xmin/xinc
   c2xqrp=1./xinc
   tem1 = 0.001 ** 0.2046
   tem2 = 0.001 ** 0.525
   do jx = 1,nqrp
     x         = xmin + (jx-1)*xinc
     tbqrp(jx) =        x ** 0.1364
     tbqra(jx) = tem1 * x ** 0.2046
     tbqrb(jx) = tem2 * x ** 0.525
   enddo
!
   return
   end subroutine ras2_setup_param1
!
!-------------------------------------------------------------------------------
   subroutine ras2_setup_param2
!-------------------------------------------------------------------------------
   parameter(nvtp=10001, vtpexp=-0.3636)
   real tbvtp(nvtp)
   common/comvtp/ c1xvtp,c2xvtp,tbvtp
!-------------------------------------------------------------------------------
!fpp$ noconcur r
   xmin=0.05
   xmax=1.5
   xinc=(xmax-xmin)/(nvtp-1)
   c1xvtp=1.-xmin/xinc
   c2xvtp=1./xinc
   do jx = 1,nvtp
     x         = xmin + (jx-1)*xinc
     tbvtp(jx) =        x ** vtpexp
   enddo
!
   return
   end subroutine ras2_setup_param2
!
!-------------------------------------------------------------------------------
   subroutine ras2_setup_hybrid(lm, prsi, prsl, cp, rgas,                      &
                     sig, sgb, prj,                                            &
                     krmin, krmax, kfmax,                                      &
                     mct, kctop)
!-------------------------------------------------------------------------------
   dimension prsi(lm+1),prsl(lm)
   dimension sig(lm+1), sgb(lm),   prj(lm+1)
!
   dimension pctop(mct+1), kctop(mct+1)
   data pctp/70.0/, pcbt/900.0/, pinc/5.0/ 
!-------------------------------------------------------------------------------
   rkap = rgas / cp
   lmm1 = lm - 1
!
   do l = 1,lm
     sig(l) = prsi(lm-l+2)/prsi(1)
     prj(l) = (0.001*sig(l)) ** rkap
     sgb(l) = prsl(lm-l+1)/prsi(1)
   enddo
!
   sig(lm+1) = prsi(1)/prsi(1)
   prj(lm+1) = (0.001*sig(lm+1)) ** rkap
!
   tem  = (pcbt - pctp - (mct-1)*mct*pinc/2) / mct
   pctop(1) = pctp * 0.001
!
   do i = 1,mct
     pctop(i+1) = 0.001 * (pctop(i)*1000.0 + tem + pinc*(i-1))
   enddo
!
   do l = 1,lm
     if (sig(l) .le. pctop(1)) kctop(1) = l
   enddo
!
   tem2 = 0.0
   do i = 2,mct+1
     if (kctop(i-1) .lt. lm) then
       ii       = kctop(i-1) + 1
       kctop(i) = ii
       do l = ii,lm
         if (sig(l) .le. pctop(i)) kctop(i) = l
       enddo
       if (pctop(i) .gt. sig(kctop(i))) then
         tem1 = sig(kctop(i)) - sig(kctop(i-1))
         if (tem1 .lt. tem2) kctop(i) = min(lm, kctop(i)+1)
         tem2 = sig(kctop(i)) - sig(kctop(i-1))
       endif
     else
       kctop(i) = lm + 1
     endif
   enddo
!
   do i = 1,mct+1
      if (kctop(i) .le. lm) mctm = i
   enddo
!
   krmin = 1
   krmax = mctm - 1
   kfmax = mctm - 1
   do l = 1,mctm-1
      sgc  = 0.5 * (sig(kctop(l)) + sig(kctop(l+1)))
      if (sgc .le. 0.760) krmax = l
      if (sgc .le. 0.930) kfmax = l
   enddo      
!
   return
   end subroutine ras2_setup_hybrid
!
!-------------------------------------------------------------------------------
#endif /* RASV2 end */
   end module phys_ras2_module
