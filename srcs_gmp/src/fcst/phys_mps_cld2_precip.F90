!
   subroutine phys_mps_cld2_precip (im,ix,km,deltim,delprsi,prsi,prsl,         &
                                        q,cwm,t,rn,lat,rcs)
!-------------------------------------------------------------------------------
!
! subroutine: phys_mps_cld2_precip
!
! abstract:
!    subroutine for precipitation processes 
!    from suspended cloud water/ice
! 
! program history log:
!    1995-01   q. zhao                    originally created
!    1996-01   s.-y. hong                 modified
!    1998-10   s. moorthi and h.-l. pan   modified and rewritten
!    2000-04   s.-y. hong                 modified
!
! references:
!    zhao and carr (1997), monthly weather review (august)       
!    sundqvist et al., (1989) monthly weather review. (august)  
!
!-------------------------------------------------------------------------------
!
!     in this code vertical indexing runs from surface to top of the
!     model
!
!     argument list:
!     --------------
!       im         : inner dimension over which calculation is made
!       ix         : maximum inner dimension
!       km         : number of vertical levels
!       dt         : time step in seconds
!       del(km)    : sigma layer thickness (bottom to top)
!       sl(km)     : sigma values for model layers (bottom to top)
!       ps(im)     : surface pressure (centibars)
!       q(ix,km)   : specific humidity (updated in the code)
!       cwm(ix,km) : condensate mixing ratio (updated in the code)
!       t(ix,km)   : temperature       (updated in the code)
!       rn(im)     : precipitation over one time-step dt (m/dt)
!       sr(im)     : index (=-1 snow, =0 rain/snow, =1 rain)
!       tcw(im)    : vertically integrated liquid water (kg/m**2)
!       cll(ix,km) : cloud cover
!
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
   use constant, only  : cp=>cp_,g=>g_,hsub_,hvap_,psat_,rd_,rv_,ttp_
!-------------------------------------------------------------------------------
   real,parameter       ::  r=rd_
   real,parameter       ::  eps=rd_/rv_,     epsm1=eps-1.0
   real,parameter       ::  h1=1.e0,     h2=2.e0,         h1000=1000.0
   real,parameter       ::  d00=0.e0,    d125=.125e0,     d5=0.5e0
   real,parameter       ::  a1=psat_,    a2=17.2693882e0, a3=273.16e0
   real,parameter       ::  tresh=.95e0, pq0=379.90516e0
   real,parameter       ::  elwv=hvap_,  eliv=hsub_,      row=1.e3
   real,parameter       ::  epsq=2.e-12, dldt=2274.e0,    tm10=ttp_-10.0
   real,parameter       ::  h1000g=h1000/g,a4=35.86e0,eliw=eliv-elwv
   real,parameter       ::  arcp=a2*(a3-a4)/cp, rcp=h1/cp
   real,parameter       ::  pq0c=pq0*tresh, rrog=h1/(row*g)
   real,parameter       ::  rrow=h1/row
!
   real                 ::  q(ix,km), t(ix,km),  cwm(ix,km)
   real                 ::  rn(ix)
!
!  local variables and arrays
!
   real                 ::  cll(ILOTS,levs_)
   real                 ::  wmin(ILOTS,levs_)
   real                 ::  wmini(ILOTS,levs_)
   real                 ::  sr(ILOTS),tcw(ILOTS)
   real                 ::  tt(ILOTS),qq(ILOTS)
   real                 ::  ww(ILOTS)
   real                 ::  err(ILOTS),ers(ILOTS)
   real                 ::  precrl(ILOTS),precsl(ILOTS)
   real                 ::  precrl1(ILOTS),precsl1(ILOTS)
   real                 ::  const(ILOTS), rq(ILOTS)
   real                 ::  condt(ILOTS)
   real                 ::  conde(ILOTS), rconde(ILOTS)
   real                 ::  tmt0(ILOTS)
   real                 ::  wmink(ILOTS),pres(ILOTS)
   real                 ::  u00k(levs_),u00i(levs_),cclim(levs_)
   real                 ::  ke, mi0
   real                 ::  delprsi(ILOTS,km), prsl(ILOTS,km), prsi(ILOTS,ix)
   logical              ::  comput(ILOTS)
   integer              ::  iw(ILOTS,levs_)
   integer              ::  ipr(ILOTS)
   integer              ::  iwl(ILOTS),iwl1(ILOTS)
!
!-----------------------preliminaries ------------------------------------------
!
   do k = 1,km
     do i = 1,im
       cll(i,k) = 0.0
     enddo
   enddo
!
   dt      = 2.0*deltim
   rdt     = h1 / dt
   ke      = 2.0e-6
   us      = h1
   cclimit = 1.0e-3
   climit  = 1.0e-20
   cws     = 0.025
   csm1    = 5.0000e-8
   crs1    = 5.00000e-6
   crs2    = 6.66600e-10
   cr      = 5.0e-4
   aa2     = 1.25e-3
   dtcp    = dt * rcp
   c00 = 3.0e-1 * dt
   cmr = 1.0 / 3.0e-4
   u00b = 0.85
   u00t = 0.85
   tem  = (u00t - u00b) / float(km-1)
!
!--------calculate c0 and cmr using lc at previous step-------------------------
!
   rcsi = 1.0 / rcs
   do k = 1,km
     u00k(k)  = u00b + tem * float(k-1)
     u00k(k)  = min(0.99, u00k(k)+(1.0-u00k(k))*(1.0-rcsi))
     u00i(k)  = 1.0 / u00k(k)
   enddo
!
   do k = 1,km
     do i = 1,im
       iw(i,k)   = 0.0
       wmin(i,k) = 1.0e-10
       wmini(i,k) = 5.0e-6
     enddo
   enddo
!
   do i = 1,im
     iwl1(i)    = 0
     precrl1(i) = d00
     precsl1(i) = d00
     comput(i)  = .false.
     rn(i)      = d00
     sr(i)      = d00
   enddo
!
!------------select columns where rain can be produced--------------------------
!
   do k = 1,km-1
     do i = 1,im
       if (cwm(i,k) .gt. wmin(i,k)) comput(i) = .true.
     enddo
   enddo
!
   ihpr = 0
   do i = 1,im
     if (comput(i)) then
       ihpr      = ihpr + 1
       ipr(ihpr) = i
     endif
   enddo
!
!-------------------------------------------------------------------------------
!
   do n = 1,ihpr
     const(n) = prsi(ipr(n),1) * (h1000*dt/g)
   enddo
!
!*******************************************************************************
!
   do k = km,1,-1
     do n = 1,ihpr
       precrl(n) = precrl1(n)
       precsl(n) = precsl1(n)
       err  (n)  = d00
       ers  (n)  = d00
       iwl  (n)  = 0
       i         = ipr(n)
       tt(n)     = t(i,k)
       qq(n)     = q(i,k)
       ww(n)     = cwm(i,k)
       wmink(n)  = wmin(i,k)
       pres(n)   = prsl(i,k) * h1000
       precrk = max(0.,    precrl1(n))
       precsk = max(0.,    precsl1(n))
       wwn    = max(ww(n), climit)
       if (wwn .gt. wmink(n) .or. (precrk+precsk) .gt. d00) then
         comput(n) = .true.
       else
         comput(n) = .false.
       endif
     enddo
!
     do n = 1,ihpr
       if (comput(n)) then
         wwn    = max(ww(n), climit)
         conde(n)  = const(n) * delprsi(n,k)/prsi(n,1)
         condt(n)  = conde(n) * rdt
         rconde(n) = h1 / conde(n)
         qk        = max(epsq,  qq(n))
         tmt0(n)   = tt(n) - 273.16
!
!  the global qsat computation is done in cb
!
         pres1   = pres(n) * 0.001
!
!           qw      = fpvs(tt(n))
!
         qw      = fpvs0(tt(n))
         qw      = eps * qw / (pres1 + epsm1 * qw)
         qw      = max(qw,epsq)
         tmt15 = min(tmt0(n), -15.)
         ai    = 0.008855
         bi    = 1.0
         if (tmt0(n) .lt. -20.0) then
           ai = 0.007225
           bi = 0.9674
         endif
         qi   = qw * (bi + ai*min(tmt0(n),0.))
         qint = qw * (1.-0.00032*tmt15*(tmt15+15.))
         if (tmt0(n).le.-40.) qint = qi
!
!-------------------ice-water id number iw--------------------------------------
!
         if(tmt0(n).lt.-15.) then
           fi = qk - u00k(k)*qi
           if(fi.gt.d00.or.wwn.gt.climit) then
             iwl(n) = 1
           else
             iwl(n) = 0
           endif
         elseif (tmt0(n).ge.0.) then
           iwl(n) = 0
         else
           iwl(n) = 0
           if(iwl1(n).eq.1.and.wwn.gt.climit) iwl(n)=1
         endif
!
!----------------the satuation specific humidity--------------------------------
!
         fiw   = float(iwl(n))
         qc    = (h1-fiw)*qint + fiw*qi
!
!----------------the relative humidity------------------------------------------
!
         if(qc .le. 1.0e-10) then
           rq(n) = d00
         else
           rq(n) = qk / qc
         endif
       endif
     enddo
!
!---   precipitation production --  auto conversion and accretion
!
     do n = 1,ihpr
       if (comput(n)) then
         wws    = ww(n)
         cwmk   = max(0.0, wws)
         if (iwl(n) .eq. 1) then                 !  ice phase
           amaxcm = max(0.0, cwmk - wmini(ipr(n),k))
           expf      = dt * exp(0.025*tmt0(n))
           psaut     = min(cwmk, 1.0e-3*expf*amaxcm)
           ww(n)     = ww(n) - psaut
           cwmk      = max(0.0, ww(n))
           psaci     = min(cwmk, aa2*expf*precsl1(n)*cwmk)
           ww(n)     = ww(n) - psaci
           precsl(n) = precsl(n) + (wws - ww(n)) * condt(n)
         else                                    !  liquid water
           amaxcm = max(0.0, cwmk - wmink(n))
           tem1      = precsl1(n) + precrl1(n)
           praut     = min(cwmk, c00*amaxcm*amaxcm)
           ww(n)     = ww(n) - praut
           cwmk      = max(0.0, ww(n))
           pracw     = min(cwmk, cr*dt*tem1*cwmk)
           ww(n)     = ww(n) - pracw
           precrl(n) = precrl(n) + (wws - ww(n)) * condt(n)
         endif
       endif
     enddo
!
!-----evaporation of precipitation----------------------------------------------
!**** err & ers positive--->evaporation-- negtive--->condensation
!
     do n = 1,ihpr
       if (comput(n)) then
         qk     = max(epsq,  qq(n))
         tmt0k  = max(-30.0, tmt0(n))
         precrk = max(0.,    precrl(n))
         precsk = max(0.,    precsl(n))
         amaxrq = max(0.,    u00k(k)-rq(n)) * conde(n)
!
!-------------------------------------------------------------------------------
! increase the evaporation for strong/light prec
!-------------------------------------------------------------------------------
!
         ppr    = ke * amaxrq * sqrt(precrk)
         if (tmt0(n) .ge. 0.) then
           pps = 0.
         else
           pps = (crs1+crs2*tmt0k) * amaxrq * precsk * u00i(k)
         end if
!
!---------------correct if over-evapo./cond. occurs-----------------------------
!
         erk=precrk+precsk
         if(rq(n).ge.1.0e-10)  erk = amaxrq * qk * rdt / rq(n)
         if (ppr+pps .gt. abs(erk)) then
           rprs   = erk / (precrk+precsk)
           ppr    = precrk * rprs
           pps    = precsk * rprs
         endif
         ppr       = min(ppr, precrk)
         pps       = min(pps, precsk)
         err(n)    = ppr * rconde(n)
         ers(n)    = pps * rconde(n)
         precrl(n) = precrl(n) - ppr
         precsl(n) = precsl(n) - pps
       endif
     enddo
!
!--------------------melting of the snow----------------------------------------
!
     do n = 1,ihpr
       if (comput(n)) then
         if (tmt0(n) .gt. 0.) then
           amaxps = max(0.,    precsl(n))
           psm1   = csm1 * tmt0(n) * tmt0(n) * amaxps
           psm2   = cws * cr * max(0.0, ww(n)) * amaxps
           ppr    = (psm1 + psm2) * conde(n)
           if (ppr .gt. amaxps) then
             ppr  = amaxps
             psm1 = amaxps * rconde(n)
           endif
           precrl(n) = precrl(n) + ppr
           precsl(n) = precsl(n) - ppr
         else
           psm1 = d00
         endif
!
!---------------update t and q--------------------------------------------------
!
         tt(n) = tt(n) - dtcp * (elwv*err(n)+eliv*ers(n)+eliw*psm1)
         qq(n) = qq(n) + dt * (err(n)+ers(n))
       endif
     enddo
!
     do n = 1,ihpr
       iwl1(n)    = iwl(n)
       precrl1(n) = max(0.0, precrl(n))
       precsl1(n) = max(0.0, precsl(n))
       i          = ipr(n)
       t(i,k)     = tt(n)
       q(i,k)     = qq(n)
       cwm(i,k)   = ww(n)
       iw(i,k)    = iwl(n)
     enddo
!
!  move water from vapor to liquid should the liquid amount be negative
!
     do i = 1,im
       if (cwm(i,k) .lt. 0.) then
         q(i,k)   = q(i,k) + cwm(i,k)
         t(i,k)   = t(i,k) - elwv * rcp * cwm(i,k)
         cwm(i,k) = 0.
       endif
     enddo
   enddo                               ! k loop ends here!
!
!*******************************************************************************
!-----------------------end of precipitation processes--------------------------
!*******************************************************************************
!
   do n = 1,ihpr
     i = ipr(n)
     rn(i) = (precrl1(n)  + precsl1(n)) * rrow  ! precip at surface
!
!----sr=1 if sfc prec is rain ; ----sr=-1 if sfc prec is snow
!----sr=0 for both of them or no sfc prec
!
     rid = 0.
     sid = 0.
     if (precrl1(n) .ge. 1.e-13) rid = 1.
     if (precsl1(n) .ge. 1.e-13) sid = -1.
     sr(i) = rid + sid  ! sr=1 --> rain, sr=-1 -->snow, sr=0 -->both
   enddo
!
   return
   end subroutine phys_mps_cld2_precip
