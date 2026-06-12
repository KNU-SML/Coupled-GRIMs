#include <define.h>
   subroutine phys_cps_ras2(lm, dt                                             &
              ,krmin, krmax, kfmax, ncrnd, rannum, pdd                         &
              ,mct, kctop                                                      &
              ,rgas, cp, grav, alhl, alhf                                      &
              ,tin, qin, uin, vin, pin, kpbl                                   &
              ,rainc, kbot, ktop, icps, lat, cdrag                             &
              ,sig, prj, sgb, rasal, clw, clt, dsfc) 
#ifdef RASV2
!-------------------------------------------------------------------------------
!
! program history log:
!   1995-01-01  moorthi                development
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use phys_ras2_module, only : ras2_solver
!-------------------------------------------------------------------------------
   implicit none
!
   real,parameter       ::  frac=0.5, crtmsf=0.0
   real,parameter       ::  tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf)
   real,parameter       ::  max_neg_bouy=0.20
!
   logical,parameter    ::  revap=.true., cumfrc=.false.
   logical,parameter    ::  wrkfun = .false.,  updret = .false.
   logical,parameter    ::  crtfun = .true.,   calkbl = .true.
!
! input
!
   integer              ::  lm, mct, lat, kbot, ktop, icps, kbl, kpbl
   integer              ::  krmin, krmax, kfmax, ncrnd
   real                 ::  tin(lm), qin(lm), uin(lm), vin(lm), pin
   real                 ::  rainc,  clw(lm), cdrag, rannum(ncrnd), clt(lm), dsfc
   real                 ::  dt, pdd, rgas, cp, grav, alhl, alhf
!
   real                 ::  sig(lm+1),  prj(lm+1), alfint(lm)
   real                 ::  sgb(lm), rasal(lm)
   integer              ::  kctop(mct+1)
!
! locals
!
   logical              ::  lprnt
   real                 ::  pib,  pik, rain,  uvc(lm,2), tn0(lm), qn0(lm) 
   real                 ::  tcu(lm),   qcu(lm), pcu(lm)
   real                 ::  flx(lm),   qii(lm), qli(lm),  prs(lm+1), psj(lm+1)
   real                 ::  rcu(lm,2), rkap
   integer              ::  i, ia, l, ll, lmm1, lmp1
!fpp$ noconcur r
!
   lmp1 = lm + 1
   lmm1 = lm - 1
   lmp1 = lm + 1
   rkap = rgas / cp
!
   alfint(1) = 0.5
   do l = 2,lm
!    alfint(l) = (sgb(l) - sig(l)) / (sgb(l) - sgb(l-1))
     alfint(l) = 0.5
   enddo
!
   ia = 1
!
!     assumes initial value of cloud liquid water to be zero
!     i.e. no environmental liquid water used!!!
!
   do l = 1,lm
     clw(l) = 0.0
     clt(l) = 0.0
     qii(l)   = 0.0
     qli(l)   = 0.0
!
!      initialize cloudiness,heating, drying, cloudiness etc.
!
     tcu(l)   = 0.0
     qcu(l)   = 0.0
     pcu(l)   = 0.0
     flx(l)   = 0.0
     rcu(l,1) = 0.0
     rcu(l,2) = 0.0
   enddo
!
   pib = pin * 10.0
   pik = pib ** rkap
!
!     transfer input prognostic data in to local variable
!
   do l = 1,lm
     ll = lm - l + 1
     tn0(l)   = tin(ll)
     qn0(l)   = qin(ll)
     uvc(l,1) = uin(ll)
     uvc(l,2) = vin(ll)
   enddo
!
   kbl  = max(min(lm, lm+1-kpbl), lm/2)
   rain = 0.0
!
   do l = 1,lmp1
     prs(l) = pib * sig(l)
     psj(l) = pik * prj(l)
   enddo
!
   call ras2_solver(lm,  dt, pdd                                               &
             ,ncrnd, krmin, krmax, kfmax, frac, rasal, .true.                  &
             ,revap, max_neg_bouy, alfint, cumfrc                              &
!    *,          revap, max_neg_bouy, alfint, alfinq
             ,cp,  alhl, alhf, grav, rkap, rannum                              &
             ,mct, kctop                                                       &
!    *,          mct, kctop, sgh, sig
!
             ,tn0, qn0, uvc, kbl                                               &
!    *,          qli,    qii
             ,qli,    qii, clw                                                 &
             ,rain,   clt, tcu, qcu, pcu, flx, rcu                             &
             ,cdrag,  dsfc                                                     &
             ,prs,    psj, wrkfun, calkbl, crtfun, updret)
!
   rainc = rain * 0.001
!
!  if (lprnt) then
!    print *,' tn0',(tn0(imax,l),l=1,lm)
!    print *,' qn0',(qn0(imax,l),l=1,lm)
!  endif
!
   do l = 1,lm
     ll = lm - l + 1
     tin(ll) = tn0(l)
     qin(ll) = qn0(l)
     uin(ll) = uvc(l,1)
     vin(ll) = uvc(l,2)
!    clw(l ) = qli(l) + qii(i,l)  ! now done in rasv2!
   enddo
   icps  = 0
!  kbot = lmp1
!  ktop = 0
   ktop = lmp1
   kbot = 0
!
   do l = lmm1,1,-1
     if (tcu(l) .ne. 0.0) then
       icps = 1
     endif
   enddo
!
!  new test for convective clouds ! added in 08/21/96
!
   do l = lmm1,1,-1
     if (clw(l) .gt. 0.0) ktop = l
   enddo
!
   do l = 1,lmm1
     if (clw(l) .gt. 0.0) kbot = l
   enddo
   ktop = lmp1 - ktop
   kbot = lmp1 - kbot
!
   return
#endif /* RASV2 end */
   end subroutine phys_cps_ras2
