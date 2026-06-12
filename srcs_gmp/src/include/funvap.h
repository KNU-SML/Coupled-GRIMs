   real :: ftv, ft, fq, fpv, fpvs, fqs
   real :: xt, xq, xtv, xpv, xp
!
!  water vapor functions
!
!  virtual temperature from temperature and specific humidity
!
   ftv(xt,xq)=xt*                                                              &
     (1.e0+(rv_/rd_-1.e0)*xq)
!
!  temperature from virtual temperature and specific humidity
!
   ft(xtv,xq)=xtv/                                                             &
     (1.e0+(rv_/rd_-1.e0)*xq)
!
!  specific humidity from vapor pressure and air pressure
!
   fq(xpv,xp)=(rd_/rv_)*xpv/                                                   &
      (xp+(rd_/rv_-1.e0)*xpv)
!
!  vapor pressure from specific humidity and air pressure
!
   fpv(xq,xp)=xp*xq/                                                           &
      ((rd_/rv_)-                                                              &
       (rd_/rv_-1.e0)*xq)
!
!  saturation vapor pressure (kpa) from temperature
!
   fpvs(xt)=(psat_*1.e-3)*                                                    &
       exp((hvap_/rv_)/ttp_                                                    &
           -(hvap_/rv_)/xt)
!
!  saturation specific humidity from temperature and pressure (kpa)
!
   fqs(xt,xp)=fq(fpvs(xt),xp)
