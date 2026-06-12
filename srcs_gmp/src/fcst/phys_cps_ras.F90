!
   subroutine phys_cps_ras(len, lenc, km, lm, nstrp, dt,                       &
               prsi , prsl,                                                    &
               krmin, krmax, ncrnd, afac, rannum, ufac,                        &
               rgas, cp, grav, alhl,                                           &
               tin, qin, uin, vin,                                             &
               evapc, kbot, ktop, kuo, lat, cd,                                &
               sig, prj, prh, fpk, hpk, sgb, ods, rasal, prns)
#ifdef RAS
!-------------------------------------------------------------------------------
!
! program history log:
!   1995-01-01  moorthi                development
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only        :  ILOTS, levs_
   use phys_ras_module, only  :  ras_driver
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  len, lenc, km, lm, nstrp, krmin, krmax, ncrnd,lat
   real                 ::  dt, rgas, cp, grav, alhl, afac, ufac
   real                 ::  tin(lenc,km), qin(lenc,km), uin(lenc,km)
   real                 ::  vin(lenc,km), pin(len), cd(len)
   real                 ::  evapc(len)
   real                 ::  rannum(ncrnd*3)
   real                 ::  prsi(lenc,levs_+1), prsl(lenc,levs_)
   integer              ::  kbot(len), ktop(len), kuo(len)
!
   real                 ::  sig(lm+1),  prj(lm+1), prh(lm),  fpk(lm), hpk(lm)
   real                 ::  sgb(lm),  ods(lm), rasal(lm), prns(nstrp)
!-------------------------------------------------------------------------------
   call ras_driver(len, lenc, km, lm, nstrp, dt,                               &
               prsi , prsl,                                                    &
               krmin, krmax, ncrnd, afac, rannum, ufac,                        &
               rgas, cp, grav, alhl,                                           &
               tin, qin, uin, vin,                                             &
               evapc, kbot, ktop, kuo, lat, cd,                                &
               sig, prj, prh, fpk, hpk, sgb, ods, rasal, prns)
!
   return
#endif /* RAS end */
   end subroutine phys_cps_ras
