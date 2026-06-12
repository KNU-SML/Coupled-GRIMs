#include <define.h>
!-------------------------------------------------------------------------------
   subroutine lsm_vic_solver
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [lsm_vic_solver]
!      |
!      |--- [vic_exch] *
!      |--- [vic_rainorsnow] *
!      |--- [vic_runoff_solver] *
!      |--- [vic_snow_solver] *
!      |         |--- [vic_snow_intercept] *
!      |         |--- [vic_snow_melt] *
!      |         |--- [vic_snow_canopy] *
!      |         |--- function : [snow_density] *
!      |--- [vic_soil_moisture_solver] *
!      |         |--- function : [heat_capcty] *
!      |         |--- function : [soil_condvty] *
!      |--- [vic_surface_flux] *
!      |         |--- [vic_surface_temp] *
!      |         |--- function : [snow_albedo] *
!      |--- [vic_stop_run] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   end subroutine lsm_vic_solver
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vic_exch(vtype,   month,   slrough, snrough,   vz0,              &
                        zref,      wind,       u,      ra)
!-------------------------------------------------------------------------------
!
! subprogram: vic_exch
!
! abstract: calculate the vic_exch resistance for each vegetation 
!   layer, and the wind 2m above the layer boundary.  in case of 
!   an overstory, also calculate the wind in the overstory.
!   the values are normalized based on a reference height wind 
!   speed, uref, of 1 m/s.  to get wind speeds and vic_exch 
!   resistances for other values of uref, you need to multiply 
!   the here calculated wind speeds by uref and divide the 
!   here calculated vic_exch resistances by uref
!   
! comments: 
!   please study the paper by wigmosta, vail and lettenmaier (wrr 1994)
!   for understanding the equations used to compute the ra
!
! program history log:
!   2003-06~09  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from 4.0.3 uw vic (calcvic_exch.c)
!
!-------------------------------------------------------------------------------
   use vic_veglib
#include <vartyp.h>
!-------------------------------------------------------------------------------
! ------------------- input variables ----------------------------------
   real     ::  slrough        ! soil roughness (m)
   real     ::  snrough        ! snow roughness (m)
   real     ::  zref           ! reference height for windspeed (m)
   real     ::  wind           ! wind speed at zref (m/s)
   integer  ::  vtype          ! vegetation type
   integer  ::  month          ! month of current time step
! ----------------------------------------------------------------------
!
! ------------------- output variables ---------------------------------
!
! vector of length 3, with u/ra for vegetation layers
! if overstory == true the first value is u/ra in vic vegetation,
! and the second value the u/ra in the overstory. otherwise the first
! value is u/ra in the vic vegetation and the second value
! is not used. the third value is u/ra over snow.
!
   real  ::  u(3)           ! wind for vegetation layers (m/s)
   real  ::  ra(3)          ! vic_exch resistance (s/m)
   real  ::  vz0            ! vegetation roughness (m)
! ----------------------------------------------------------------------
!
! ---------------------- local variables -------------------------------
   integer          ::  ost     ! 1 with overstory  0 no overstory
   real             ::  z0      ! roughness length (m)
   real             ::  d       ! vegetation displacement (m)
   real             ::  height  ! height of the veg layers (top layer first) (m)
   real             ::  trunk   ! multiplier for height that indictaes the tunk
   real             ::  refh, d_lower, d_upper, k2
   real             ::  uh, ut, uw, z0sl, z0sn
   real             ::  z0_lower, z0_upper, zt, zw, tmp_wind, rati
   integer          ::  i
!
!     vnwd          ::  attenuation coeff for wind in the overstory
!     von_k         ::  von constant
!
   real, parameter  ::  vnwd=0.5,von_k = 0.41, huge_resist = 1.e2
! ----------------------------------------------------------------------
!-------------------------------------------------------------------------------
! ----------------------------------------------------------------------
! set parameter values
! ----------------------------------------------------------------------
!
   z0sl  = slrough
   z0sn  = snrough
!
   k2    = von_k * von_k
   z0    = veg_rough(month,vtype)
!
! z0 is very important to determine the vic_exch resistance
! however, in according with vic_veglib.h, the value of veg_rough
! is too big for the some vegetations (6/2004,ji)
!
   z0 = z0/5.0     ! 6/2004
   ost   = veg_ost(vtype)
   trunk = veg_trunk(vtype)
   d     = veg_d(month,vtype)
   height= d/0.67
   vz0 = z0 + d
!
! ----------------------------------------------------------------------
! no overstory, thus maximum one soil layer
! ----------------------------------------------------------------------
!
   if(d.lt.zref) then
     refh = zref
   else
     refh = d + zref + z0
     wind = wind*log(refh)/log(zref)
   end if
!
   do i = 1,3         ! initialize u and ra
     u(i)  = 0
     ra(i) = 0
   end do
!
   if(ost.eq.0) then
     z0_lower = z0
     d_lower  = d
!    
! ----------------------------------------------------------------------
! without snow 
! ----------------------------------------------------------------------
!
     u(1)  = log((2. + z0_lower)/z0_lower)/                                    &
             log((refh - d_lower)/z0_lower)
!
! original vic
!
!         ra(1) = log((2.+(1.0/0.63-1.0)*d_lower)/z0_lower) *                  &
!             log((2.+(1.0/0.63-1.0)*d_lower)/(0.1*z0_lower))/k2
! 
     ra(1) = log((2.+(1.0/0.63-1.0)*d_lower)/z0_lower) *                       &
             log((refh+(1.0/0.63-1.0)*d_lower)/z0_lower)/k2
!
! ----------------------------------------------------------------------
! if snow-covered ground (if currently no snow, we still need compute 
! u(3) and ra(3), because there may have snowfalling onto no snow-ground.)
! ----------------------------------------------------------------------
!
     u(3)  = log((2.+z0sn)/z0sn)/log(refh/z0sn)
     ra(3) = log((2.+z0sn)/z0sn)*log(refh/z0sn)/k2
!
   else
!
! ----------------------------------------------------------------------
! overstory present, one or two vegetation layers possible 
! ----------------------------------------------------------------------
!
     z0_upper = z0
     d_upper  = d
!    
     z0_lower = z0sl
     d_lower  = 0
!    
     zw = 1.5 * height - 0.5 * d_upper
     zt = trunk * height
!
     if (zt .lt. (z0_lower+d_lower)) then
       write(6,*) 'error: trunk space height below lower boundary'
       call vic_stop_run
     end if
!
! ----------------------------------------------------------------------
! resistance for overstory
! ----------------------------------------------------------------------
!
     ra(2) = log((refh-d_upper)/z0_upper)/k2                                   &
             * (height/(vnwd*(zw-d_upper))                                     &
             * (exp(vnwd*(1-(d_upper+z0_upper)/height))-1)                     &
             + (zw-height)/(zw-d_upper)                                        &
             + log((refh-d_upper)/(zw-d_upper)))
!    
! ----------------------------------------------------------------------
! wind at different levels in the profile 
! ----------------------------------------------------------------------
!
     uw=log((zw-d_upper)/z0_upper)/log((refh-d_upper)/z0_upper)
     uh=uw - (1-(height-d_upper)/(zw-d_upper))                                 &
        / log((refh-d_upper)/z0_upper)
     u(2) = uh * exp(vnwd * ((z0_upper+d_upper)/height - 1.))
     ut   = uh * exp(vnwd * (zt/height - 1.))
!    
! ----------------------------------------------------------------------
! resistance at the lower boundary 
! ----------------------------------------------------------------------
!
     u(1)  = log((2.+z0_upper)/z0_upper)/log((refh-d_upper)                    &
             /z0_upper)
!         ra(1) = log((2. + (1.0/0.63 - 1.0) * d_upper) / z0_upper)            &
!                 *log((2.+(1.0/0.63-1.0)*d_upper)/(0.1*z0_upper))/k2
     ra(1) = log((2. + (1.0/0.63 - 1.0) * d_upper) / z0_upper)                 &
             *log((refh+(1.0/0.63-1.0)*d_upper)/z0_upper)/k2
!
! ----------------------------------------------------------------------
! if snow covered ground
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! case 1: the wind profile to a height of 2m above the lower 
!         boundary is entirely logarithmic
! ----------------------------------------------------------------------
!
     if (zt .gt. (2. + z0sn)) then
       u(3) =ut*log((2.+z0sn)/z0sn)/log(zt/z0sn)
       ra(3)=log((2.+z0sn)/z0sn)*log(zt/z0sn)/(k2*ut) 
!            
! ----------------------------------------------------------------------
! case 2: the wind profile to a height of 2m above the lower boundary
!         is part logarithmic and part exponential, but the top of
!         the overstory is more than 2 m above the lower boundary
! ----------------------------------------------------------------------
!
     else if (height .gt. (2. + z0sn)) then
       u(3)  = uh * exp(vnwd * ((2. + z0sn)/height - 1.))
       ra(3) = log(zt/z0sn) * log(zt/z0sn)/(k2*ut) +                           &
               height * log((refh-d_upper)/z0_upper) /                         &
               (vnwd*k2*(zw-d_upper)) *(exp(vnwd*(1-zt/height))                &
               - exp(vnwd*(1-(z0sn+2.)/height)))
!         
! ----------------------------------------------------------------------
! case 3: the top of the overstory is less than 2 m above the lower 
!         boundary.  the wind profile above the lower boundary is 
!         part logarithmic and part exponential, but only extends to 
!         the top of the overstory
! ----------------------------------------------------------------------
!
     else 
       u(3)  = uh
       ra(3) = log(zt/z0sn) * log(zt/z0sn)/(k2*ut) +                           &
               height*log((refh-d_upper)/z0_upper) /                           &
               (vnwd*k2*(zw-d_upper))*(exp(vnwd*(1-zt/height))-1)
     end if
   end if
!
! -- to refrain the vic_exch resistance, ra=ra*log(ref/2.0) 6/2004 ji
!
   tmp_wind = max(1.0,min(10.0,wind))
   rati = log(refh/2.0)
   do i = 1,3
     ra(i) = ra(i)*rati
   enddo
!
   if(wind.gt.0.) then
     u(1) = u(1)*wind
     ra(1)= ra(1)/tmp_wind
     if(ost.ne.0) then
       u(2) = u(2)*wind
       ra(2)= ra(2)/tmp_wind
     end if
     u(3) = u(3)*wind
     ra(3)= ra(3)/tmp_wind
   else
     u(1) = u(1)*wind
     ra(1) = huge_resist
     if(ost.ne.0) then
       u(2) = u(2)*wind
       ra(2)= huge_resist
     end if
     u(3) = u(3)*wind
     ra(3)= huge_resist
   end if
!
#ifdef DBGVIC
   write(6,*)'in vic_exch ra ',ra,' u ',u
   write(6,*)'   vic_exch ',vtype, month, slrough, snrough,                    &
                             zref,    wind
#endif
!
   return
   end subroutine vic_exch
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vic_rainorsnow(msl,   nsl,  msub,    md,                         &
                        nveg,   prc,  tgcm,   fmu,                             &
                       lstsn, rfall, sfall,   cwt,                             &
                         csn,   swq,   rsn,   tsf,                             &
                         tpk,   sfw,   pkw,   smc,                             &
                         sic,   tnd)
!-------------------------------------------------------------------------------
!
! subprogram: vic_rainorsnow
!
! abstract: decides rainfall or snowfall, and redistribute hydrologic
!   properties in the wet and dry districts according to 
!   rainfall intensity
!
! comments: if snow falling, mu = 1.0. if no snow falling, mu varies
!   according to the rainfall intensity, and the related 
!   hydrologic variables, such as ccsnow, cwater, raintf, snowtf, 
!   and other soil variables are recalcuated. if no precipitation,
!   mu doesnt change and its value is decided by the last 
!   precipitation event. 
!
! program history log:
!   2003-06~09  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from uw_vic (calc_rainonly.c,dist_prec.c,...)
!
!-------------------------------------------------------------------------------
   use constant, only : t0c_
!-------------------------------------------------------------------------------
#include <vartyp.h>
! ------------------- input variables ----------------------------------
   integer  ::  msl           ! number of soil layer
   integer  ::  nsl           ! number of soil node
   integer  ::  msub          ! maximum number of subgrids
   integer  ::  md            ! number of precipitation district
   integer  ::  nveg          ! number of subgrids
   real     ::  prc           ! precipitation (rain+snow) (m/time step)
   real     ::  tgcm          ! air temperature (k) 
! ----------------------------------------------------------------------
!
! ------------------- modified variables -------------------------------
   integer  ::  lstsn            ! time steps since last snowfall (time step)
                                 ! lstsn = -1 (no snow), 0 (new snow)
   real     ::  fmu              ! fraction of wet district
   real     ::  cwt(md,msub)     ! canopy intercepted water 1:wet, 2:dry (m)
   real     ::  csn(md,msub)     ! canopy intercepted snow  (m)
   real     ::  swq(md,msub)     ! snow water equivalent (m)
   real     ::  rsn(md,msub)     ! snow density (kg/m3)
   real     ::  tsf(md,msub)     ! snow surface temperature (k)
   real     ::  tpk(md,msub)     ! snow pack temperature (k)
   real     ::  sfw(md,msub)     ! snow surface water (m)
   real     ::  pkw(md,msub)     ! snow pack water (m)
   real     ::  smc(msl,md,msub) ! soil moisture (mm)
   real     ::  sic(msl,md,msub) ! soil ice (mm)
   real     ::  tnd(nsl,md,msub) ! soil temperature (k)
! ----------------------------------------------------------------------
!
! ------------------- output variables ---------------------------------
   real     ::  rfall            ! rainfall (m/time step)
   real     ::  sfall            ! snowfall (m/time step)
! ----------------------------------------------------------------------
!
! ---------------------- local variables -------------------------------
   real, parameter  ::  tmaxsnow=0.5    ! maxi temperature at which snow can fall (c)
   real, parameter  ::  tminrain=-0.5   ! mini temperature at which rain can fall (c)
   real, parameter  ::  prec_expt = 0.6 ! exponential that controls the fraction of a
                                        ! grid cell that receives rain during a storm
                                        ! of given intensity
   integer          ::  nv, n, m, idi ! loop indexes
   real             ::  ta               ! air temperature (c)
   real             ::  oldmu, oldwu, wu
   logical          ::  logicswq
! ----------------------------------------------------------------------
!
   if(prc.gt.0.0) then
     oldmu = fmu
     oldwu = 1.0 - oldmu
     ta = tgcm - t0c_
!
     if(ta.lt.tmaxsnow.and.ta.gt.tminrain)then
       rfall = (ta-tminrain)/(tmaxsnow-tminrain)*prc
     else if(ta .gt.tmaxsnow) then
       rfall = prc
     end if
!
     sfall = prc - rfall
!      
! ----------------------------------------------------------------------
! the fractional coverage of precipitation over an area or grid cell, 
!     mu, is estimated using the equation from fan et. al. (wrr 1996)
!     if snowfalling, mu = 1.0
! ----------------------------------------------------------------------
!
     if(sfall.gt.0.0) then
       lstsn = 0
       fmu = 1.0
     else
       if(md.eq.1) then
         fmu = 1.0
       else
         fmu = 1.0 - exp(-prec_expt*rfall) 
       end if
     end if
!      
! ----------------------------------------------------------------------
! reorganize variables in the wet and dry districts for continous run
! ----------------------------------------------------------------------
!
     if(md.ne.1) then      ! include precipitation distribution
       if(fmu.ne.oldmu) then
         wu = 1.0 - fmu
         do nv = 1,nveg
!
           if(fmu.gt.oldmu) then             ! part of dry -> wet
             if(csn(1,nv).gt.0.or.csn(md,nv).gt.0) then
               csn(1,nv) = (csn(1,nv)*oldmu+csn(md, nv)*(fmu-oldmu))/fmu
             end if
             if(cwt(1,nv).gt.0.or.cwt(md,nv).gt.0) then
               cwt(1,nv) = (cwt(1,nv)*oldmu+cwt(md, nv)*(fmu-oldmu))/fmu
             end if
           else                              ! part of wet -> dry
             if(csn(1,nv).gt.0.or.csn(md,nv).gt.0) then
               csn(md,nv) = (csn(md,nv)*oldwu+csn(1, nv)*(wu-oldwu))/wu
             end if
             if(cwt(1,nv).gt.0.or.cwt(md,nv).gt.0) then
               cwt(md,nv) = (cwt(md,nv)*oldwu+cwt(1, nv)*(wu-oldwu))/wu
             end if
           end if
!
           if(swq(1,nv).gt.0.or.swq(md,nv).gt.0) then
             if(fmu.gt.oldmu) then
               swq(1,nv) = (swq(1,nv)*oldmu+swq(md,nv)*(fmu-oldmu))/fmu
               rsn(1,nv) = (rsn(1,nv)*oldmu+rsn(md,nv)*(fmu-oldmu))/fmu
               tsf(1,nv) = (tsf(1,nv)*oldmu+tsf(md,nv)*(fmu-oldmu))/fmu
               tpk(1,nv) = (tpk(1,nv)*oldmu+tpk(md,nv)*(fmu-oldmu))/fmu
               sfw(1,nv) = (sfw(1,nv)*oldmu+sfw(md,nv)*(fmu-oldmu))/fmu
               pkw(1,nv) = (pkw(1,nv)*oldmu+pkw(md,nv)*(fmu-oldmu))/fmu
             else
               swq(md,nv) = (swq(md,nv)*oldwu+swq(1,nv)*(wu-oldwu))/wu
               rsn(md,nv) = (rsn(md,nv)*oldwu+rsn(1,nv)*(wu-oldwu))/wu
               tsf(md,nv) = (tsf(md,nv)*oldwu+tsf(1,nv)*(wu-oldwu))/wu
               tpk(md,nv) = (tpk(md,nv)*oldwu+tpk(1,nv)*(wu-oldwu))/wu
               sfw(md,nv) = (sfw(md,nv)*oldwu+sfw(1,nv)*(wu-oldwu))/wu
               pkw(md,nv) = (pkw(md,nv)*oldwu+pkw(1,nv)*(wu-oldwu))/wu
             end if
           end if
!
           do m = 1,msl
             if(fmu.gt.oldmu) then             ! part dry -> wet
               smc(m,1,nv) = (smc(m,1,nv)*oldmu+smc(m,md,nv)*(fmu-oldmu))/fmu
               sic(m,1,nv) = (sic(m,1,nv)*oldmu+sic(m,md,nv)*(fmu-oldmu))/fmu
             else                              ! part wet -> dry
               smc(m,md,nv) = (smc(m,md,nv)*oldwu+smc(m,1,nv)*(wu-oldwu))/wu
               sic(m,md,nv) = (sic(m,md,nv)*oldwu+sic(m,md,nv)*(wu-oldwu))/wu
             end if
           end do
!
           do n = 1,nsl
             if(fmu.gt.oldmu) then
               tnd(n,1,nv) = (tnd(n,1,nv)*oldmu+tnd(n,md,nv)*(fmu-oldmu))/fmu
             else
               tnd(n,md,nv) = (tnd(n,md,nv)*oldwu+tnd(n,1,nv)*(wu-oldwu))/wu
             end if
           end do
!
         end do      ! loop of nveg
       end if
     end if        ! end of precip distribution
   end if        ! end of prc>0
!      
   if(sfall.le.0.0) then
     logicswq = .false.
     do nv = 1,nveg         !whether or not snow exists on ground
       if(swq(1,nv).gt.0.or.swq(md,nv).gt.0) then
         logicswq = .true.
       end if
     end do
!
     if(logicswq) then
       if(lstsn.ge.0) then
         lstsn = lstsn + 1
       else
         lstsn = 0
       end if
     else
       lstsn = -1
     end if
   end if
!
   return
   end subroutine vic_rainorsnow
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vic_runoff_solver(msl,  dtime,    ppt,  binf,                    &
                      ws,     ds,    dsm,   cef,                               &
                     smx,    smr,    kst,  expt,                               &
                     dph,    bub,                                              &
                    fsmc,   fsic,  sflow, bflow)
!-------------------------------------------------------------------------------
!
! subroutine: vic_runoff_solver.f
!
! abstract: this subroutine calculates infiltration and runoff from 
!            the surface, gravity driven drainage between all soil 
!            layers, and generates bflow from the bottom layer.
!
! sublayer indecies are always is the current vic model moisture layer
!    1 = thawed sublayer, 2 = frozen sublayer, and 3 = unfrozen sublayer.
!
! program history log:
!           2003-09-01 ji chen modified from 4.0.3 uw vic (runoff.c)
!           2008-08-01 kyeong-hee seol debugged, scm option
!
!-------------------------------------------------------------------------------
#include <vartyp.h>
! ------------------- input variables ----------------------------------
   integer  ::  msl         ! number of soil layer
   real     ::  dtime       ! time step (second)
   real     ::  ppt         ! incoming precip and snow melt (m/time step)
   real     ::  binf        ! vic parameter
   real     ::  ws          ! fraction of max soilm nonlinear baseflow
   real     ::  ds          ! fraction of dsm begin nonlinar bflow
   real     ::  dsm         ! maximum velocity of baseflow (mm/day)
   real     ::  cef         ! exponent used in infiltration curve (n/a)
   real     ::  fsic(msl)   ! soil layer ice (mm)
   real     ::  kst(msl)    ! saturated hydraulic conductivity (mm/day)
   real     ::  expt(msl)   ! parameter describing kst with soil moisture
   real     ::  dph(msl)    ! soil layer thickness (m)
   real     ::  bub(msl)    ! bubbling pressure of soil layer (cm)
   real     ::  smx(msl)    ! soil layer maximum soil moisture (mm)
   real     ::  smr(msl)    ! soil layer residue (mm)
! ----------------------------------------------------------------------
!
! ------------------- modified variables -------------------------------
   real     ::  fsmc(msl)   ! soil layer moisture (mm)
! ----------------------------------------------------------------------
!
! ------------------- output variables ---------------------------------
   real     ::  sflow       ! surface runoff (mm/time step)
   real     ::  bflow       ! base flow (mm/time step)
! ----------------------------------------------------------------------
!
! ---------------------- local variables -------------------------------
   real, parameter  ::  secpday=86400.0
   integer  ::  k, tmp
   real     ::  frcday          ! fractional day of the model time step
   real     ::  inflow, max_infil, ex, a, i_0, basis
   real     ::  moist(msl), q12(msl-1)
   real     ::  top_moist, tmp_moist, tmp_smcn
   real     ::  top_smx, tmp_inflow, tmpws, frac
   real     ::  avgmatri, avgsmr, avgsmx, avgbub, avgb, avgkst, tmpdph
   real     ::  matri(msl), dif(msl), b(msl)
   real     ::  smcon(msl), smrcon(msl), smxcon(msl)
   real     ::  d1, d2, d3, d4
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! ppt = amount of liquid water coming to the surface
! ----------------------------------------------------------------------
!
   inflow = ppt*1000.0
   frcday = dtime / secpday
!
! ----------------------------------------------------------------------
! initialize variables
! ----------------------------------------------------------------------
!
   do k = 1,msl
!
! ----------------------------------------------------------------------
! set layer unfrozen moisture content
! ----------------------------------------------------------------------
!
     moist(k) = fsmc(k) - fsic(k)
     if(moist(k).lt.0) moist(k) = 0.0
   end do
!   
! ----------------------------------------------------------------------
! surf_runoff based on soil moisture level of upper layers
! ----------------------------------------------------------------------
!
   top_moist = 0.
   top_smx   = 0.
!
   do k = 1,msl-1
     top_moist = top_moist + fsmc(k)
     top_smx   = top_smx + smx(k)
   end do
!
   if(top_moist.gt.top_smx) top_moist=top_smx
!   
! ----------------------------------------------------------------------
! calculate surf_runoff from surface--surf_runoff calculations 
!   for top layer only
!   a and i_0 as in wood et al. in jgr 97, d3, 1992 equation (1)
! ----------------------------------------------------------------------
!
   if(binf.le.0) then
     write(6,*) 'error vic_runoff_solver vic binf must be a positive value '
     write(6,*) ' binf = ',binf
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
   max_infil = (1.0 + binf) * top_smx
   ex        = binf / (1.0 + binf)
!
   if(abs(top_moist-top_smx).lt.1.e-6) then
     a   = 1.0
     i_0 = max_infil
   else
     a   = 1.0 - (1.0 - top_moist / top_smx)**ex
     i_0 = max_infil*(1.0-(1.0 - a)**(1.0/binf))
   endif
!
! ----------------------------------------------------------------------
! maximum inflow ---- equation (3a) wood et al.
! ----------------------------------------------------------------------
!
   if(inflow.le.0.0) sflow = 0.0
!
   if(inflow.gt.0) then
     if (max_infil.eq.0.0) then
       sflow = inflow
     else if((i_0 + inflow) .ge. max_infil) then
       sflow = inflow - top_smx + top_moist
     else 
       basis = 1.0 - (i_0 + inflow) / max_infil
       sflow = inflow - top_smx + top_moist + top_smx                          &
                      *(basis**(1.0*(1.0 + binf)))
     endif
   endif
!
   if(sflow.lt.0.) sflow = 0.0
   if(sflow.gt.inflow) sflow = inflow
!
   inflow = inflow - sflow
!
#ifdef DBGVIC
   write(6,*) 'in vic_runoff_solver before computing bflow, sflow=',sflow
#endif
!
! ----------------------------------------------------------------------
! compute flow between soil layers 
! ----------------------------------------------------------------------
!
#ifdef RUNSMDIF
   do k = 1,msl
     if(expt(k).gt.3.0) then
       b(k) = (expt(k)-3.0)/2.0
     else
       b(k) = 1.e-6
     endif
     if(moist(k).gt.smr(k))then
       matri(k)=10*bub(k)*((moist(k)-smr(k))/(smx(k)-smr(k)))**(-b(k))
     else
       matri(k)=1.e6
     endif
     smcon(k) =max(1.e-6, moist(k)/1000.0/dph(k))   ! mm -> unitless
     smrcon(k)=max(1.e-6, smr(k)/1000.0/dph(k))
     smxcon(k)=smx(k)/1000.0/dph(k)
   enddo
#endif
!
! ----------------------------------------------------------------------
! compute drainage between sublayers 
! ----------------------------------------------------------------------
!
   do k = 1,msl-1
!
! ----------------------------------------------------------------------
! brooks & corey relation for hydraulic conductivity
! ----------------------------------------------------------------------
!
#ifdef RUNSMDIF
     if(moist(k).gt.smr(k))then
       tmpdph = dph(k)+dph(k+1)
       avgmatri=10**((dph(k+1)*log10(matri(k))+dph(k)*                         &
                       log10(matri(k+1)))/tmpdph)
       avgsmr=10**((dph(k+1)*log10(smrcon(k))+dph(k)*                          &
                       log10(smrcon(k+1)))/tmpdph)
       avgsmx=10**((dph(k+1)*log10(smxcon(k))+dph(k)*                          &
                       log10(smxcon(k+1)))/tmpdph)
       avgbub=10**((dph(k+1)*log10(bub(k))+dph(k)*                             &
                       log10(bub(k+1)))/tmpdph)
       avgb  =10**((dph(k+1)*log10(b(k))+dph(k)*                               &
                       log10(b(k+1)))/tmpdph)
       avgkst=10**((dph(k+1)*log10(kst(k))+dph(k)*                             &
                       log10(kst(k+1)))/tmpdph)
!
       tmp_smcn = avgsmr + (avgsmx-avgsmr)*                                    &
                     ((avgmatri/avgbub)**(-1.0/avgb))
!
       q12(k)= kst(k) * ((moist(k)-smr(k)) / (smx(k)-smr(k)))**expt(k)
!
       if(abs(smcon(k)-smcon(k+1)).lt.1.e-6) then
         dif(k) = 0.0
       else
         dif(k) = (avgb*avgkst*avgbub/avgsmx)*(((tmp_smcn-                     &
                      avgsmr)/(avgsmx-avgsmr))**(avgb+2))
         dif(k) = dif(k)*(smcon(k)-smcon(k+1))/(tmpdph/2.0)
       endif
     else
       q12(k) = 0.0
       dif(k) = 0.0
     end if
#else
     dif(k) = 0.0
     if(moist(k).gt.smr(k))then
       q12(k) = kst(k)* ((moist(k)-smr(k)) / (smx(k)-smr(k)))**expt(k)
     else
       q12(k) = 0.0
     endif
!    if(fsic(k).lt.1.e-6 .and. q12(k).lt.1.e-6)then
!      if(expt(k).gt.13)then
!        write(6,101)k,q12(k),kst(k),moist(k)/dph(k)/1000.0,
!     &               fsic(k)/dph(k)/1000.0,smr(k)/dph(k)/1000.0,
!     &               smx(k)/dph(k)/1000.0,expt(k),dph(k)
! 101        format('k=',i2,8(1x,e10.4))
!         endif
#endif
!
   end do
!
! -- q12(k) (mm/day) -> (mm/time step)
!   
   do k = 1,msl-1
     q12(k) = q12(k) * frcday
     dif(k) = dif(k) * frcday
   end do
!
! ----------------------------------------------------------------------
! solve for current soil layer moisture, and check versus maximum 
! and minimum moisture contents.  
! ----------------------------------------------------------------------
!
   do k = 1,msl-1
     tmp_inflow = 0.
!
! ----------------------------------------------------------------------
! update soil layer moisture content
! ----------------------------------------------------------------------
!
     moist(k) = moist(k) + inflow - q12(k) - dif(k)
!
! ----------------------------------------------------------------------
! verify that soil layer moisture is less than maximum
! ----------------------------------------------------------------------
!
     if((moist(k)+fsic(k)).gt.smx(k))then
       tmp_inflow=moist(k)+fsic(k)-smx(k)
       moist(k) = smx(k)-fsic(k)
       if(k.eq.1) then
         q12(k) = q12(k) + tmp_inflow
         tmp_inflow = 0
       else
         tmp = k
         do while (tmp_inflow.gt.0)
           tmp=tmp-1
           if(tmp.eq.0) then
!
! ----------------------------------------------------------------------
! if top layer saturated, add to surf_runoff
! ----------------------------------------------------------------------
!
             sflow = sflow + tmp_inflow
             tmp_inflow = 0.0
           else
!
! ----------------------------------------------------------------------
! else add excess soil moisture to next higher layer
! ----------------------------------------------------------------------
!
             moist(tmp) = moist(tmp)+tmp_inflow
             if((moist(tmp)+fsic(tmp)).gt.smx(tmp))then
               tmp_inflow=(moist(tmp)+fsic(tmp))-smx(tmp)
               moist(tmp)=smx(tmp)-fsic(tmp)
             else
               tmp_inflow= 0
             endif
           endif
         enddo    ! end check if excess moisture in top layer
       endif
     endif
!
! ----------------------------------------------------------------------
! verify that current layer moisture is greater than minimum
! ----------------------------------------------------------------------
!
     if ((moist(k)+fsic(k)).lt.smr(k))then
!
! ----------------------------------------------------------------------
! moisture cannot fall below residual moisture content
! ----------------------------------------------------------------------
!
       q12(k)  = q12(k)+moist(k)+fsic(k)-smr(k)
       moist(k)=smr(k)-fsic(k)
       if(moist(k).lt.0) moist(k) = 0.0
     end if
     inflow = q12(k)+dif(k)
   end do                    ! end loop through (msl-1) soil layers 
!
! ----------------------------------------------------------------------
! compute bflow
!    arno model for the bottom soil layer (based on bottom
!    soil layer moisture from previous time step)
! ----------------------------------------------------------------------
!   
   k = msl
#ifdef VICGLB
!
! when using vic global soil data, please check the following paper for
! the detail of computing baseflow algorithm.
! nijssen b., g.m. odonnell, d.p. lettenmaier, d. lohmann, and e.f. wood,
! predicting the discharge of global rivers. j. climate, 3307-3322,2001.
!
! d1=ds, d2=dsm, d3=ws, d4=cef
!
   if(moist(k).ge.smr(k))then
     d1 = ds
     d2 = dsm
     d3 = ws
     d4 = cef

     bflow = d1*moist(k)
     if(moist(k).gt.d3) then
       bflow = bflow + d2*(moist(k)-d3)**d4
     endif
     bflow = bflow * frcday
   else
     bflow = 0.0
   endif
#else
   if(ws.gt.1.0.or.ws.lt.0.0) then
     write(6,*) 'error in runoff_diffu wrong ws: ',ws
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
   if(moist(k).ge.smr(k))then
     tmpws = ws*(smx(k)-smr(k))
     frac  = frcday * ds * dsm /tmpws
     bflow = frac * (moist(k)-smr(k))             ! mm/time step
     if (moist(k).gt.tmpws) then
       frac = (moist(k)-tmpws)/(smx(k)-tmpws)
       bflow= bflow+frcday*dsm*(1.0-ds/ws)*(frac**cef) ! mm/time step
     end if
   else
     bflow = 0.0
   endif
#endif /* VICGLB end */
!
! turn off bflow
!      bflow = 0
!
! ----------------------------------------------------------------------
! extract bflow from the bottom soil layer 
! ----------------------------------------------------------------------
!
   moist(k) = moist(k) + inflow - bflow
!      write(6,100)inflow,bflow,moist(k)
! 100  format('runof ',3f12.7)
!
! ----------------------------------------------------------------------
! check lower sub-layer moistures 
! ----------------------------------------------------------------------
!
   tmp_moist = 0
!   
   if((moist(k)+fsic(k)).lt.smr(k)) then
!
! ----------------------------------------------------------------------
! soil moisture is below minimum 
! ----------------------------------------------------------------------
!
     bflow = bflow+moist(k)+fsic(k)-smr(k)
     moist(k) = smr(k)-fsic(k)
   end if
!
   if(bflow .lt. 0) bflow = 0
!
   if((moist(k)+fsic(k)) .gt. smx(k)) then
!
! ----------------------------------------------------------------------
! soil moisture above maximum 
! ----------------------------------------------------------------------
!
     tmp_moist= moist(k)+fsic(k)-smx(k)
     moist(k) = smx(k) - fsic(k)
     tmp = k
     do while (tmp_moist .gt. 0)
       tmp = tmp - 1
       if(tmp.eq.0) then
!
! ----------------------------------------------------------------------
! if top layer saturated, add to surf_runoff
! ----------------------------------------------------------------------
!
         sflow = sflow + tmp_moist
         tmp_moist = 0
       else 
!
! ----------------------------------------------------------------------
! else if sublayer exists, add excess soil moisture
! ----------------------------------------------------------------------
!
         moist(tmp) = moist(tmp) + tmp_moist 
         if((moist(tmp)+fsic(tmp)).gt.smx(tmp))then
           tmp_moist = moist(tmp) + fsic(tmp) - smx(tmp)
           moist(tmp)= smx(tmp)-fsic(tmp)
         else
           tmp_moist=0
         end if
       endif
     enddo
   end if
!
   do k = 1,msl
     fsmc(k) = moist(k)+fsic(k)
!         fsmc(k) = fsmc(k)
   end do
!
   return
   end subroutine vic_runoff_solver
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vic_snow_solver(msl,     dtime,    month,     rfall,             &
                            sfall,      pgcm,     tair,     longw,             &
                          netshts,   netshtv,   rhoair,     vpair,             &
                              vpd,        ra,        u,       wcr,             &
                              wpw,     vtype,      ost,       frt,             &
                            fflai,     tgrnd,     fsmc,      fsic,             &
                             fswq,      frsn,     ftsf,      ftpk,             &
                             fsfw,      fpkw,     fcwt,      fcsn,             &
                             frtf,    fsnowm,   sncvfr,     laths,             &
                            senhs)
!-------------------------------------------------------------------------------
!
! subroutine: vic_snow_solver
!
! abstract: this routine is to handle the various calls for 
!           solving the various components of the vic snow scheme.
!
! program history log:
!   2003-09-01  ji chen                modified from 4.0.3 uw vic (runoff.c)
!   2008-08-01  kyeong-hee seol        debugged, scm option
!
!-------------------------------------------------------------------------------
   use vic_veglib
#include <vartyp.h>
! ------------------- input variables ----------------------------------
!
! -- model basic parameters
!
   integer  ::  msl              ! number of soil layer
   integer  ::  month            ! current month
   real     ::  dtime            ! time step (s)
!
! -- atmosphere variables
!
   real     ::  rfall            ! rainfall (m/time step)
   real     ::  sfall            ! snowfall (m/time step)
   real     ::  pgcm             ! pressure (pa)
   real     ::  tair             ! air temperature (k)
   real     ::  longw            ! long wave (w/m2)
   real     ::  netshts          ! net short wave radiation on snow ground (w/m2)
   real     ::  netshtv          ! net short wave radiation no snow ground (w/m2)
   real     ::  vpair            ! actual vapor pressure of air (pa)
   real     ::  vpd              ! vapor pressure deficit (pa)
   real     ::  rhoair           ! air density (kg/m3)
   real     ::  u(3)             ! wind speed (m/s)
   real     ::  ra(3)            ! vic_exch resistance (s/m)
!
! -- soil parameters
!
   real     ::  wcr(msl)         ! ~70% of field capacity (mm)
   real     ::  wpw(msl)         ! wilting point soil moisture (mm)
!
! -- vegetation parameters
!
   integer  ::  vtype            ! vegetation class
   integer  ::  ost              ! vegetation overstory (1: with, 0: no)
   real     ::  frt(msl)         ! root content (fraction)
   real     ::  fflai            ! leaf area index
!
! -- land surface variable
!
   real     ::  tgrnd            ! ground surface temperature (k)
! ----------------------------------------------------------------------
!
! ------------------- modified variables -------------------------------
   real     ::  fswq             ! snow water equivalent (m h2o)
   real     ::  frsn             ! snow density (kg/m3)
   real     ::  ftsf             ! snow surface temperature (k)
   real     ::  ftpk             ! snow pack temperature (k)
   real     ::  fsfw             ! snow surface water equivalent (m h2o)
   real     ::  fpkw             ! snow pack water equivalent (m h2o)
   real     ::  fcwt             ! canopy intercepted water (m h2o)
   real     ::  fcsn             ! canopy intercepted snow (m h2o)
   real     ::  fsmc(msl)        ! soil moisture (liquid + ice) (mm)
   real     ::  fsic(msl)        ! soil ice (mm)
! ----------------------------------------------------------------------
!
! ------------------- output variables ---------------------------------
   real     ::  frtf             ! rain through fall (m/time step)
   real     ::  fsnowm           ! snow melt (m/time step)
   real     ::  sncvfr           ! snow cover fraction (fraction)
   real     ::  laths            ! latent heat flux (w/m2)
   real     ::  senhs            ! sensible heat flux (w/m2)
! ----------------------------------------------------------------------
!
! ------------------- local variables ----------------------------------
   real     ::  fstf             ! snow through fall (m/time step)
   real     ::  cnplongw         ! longwave radiation to canopy (w/m2)
   real     ::  longout          ! land emitted longwave radiation (w/m2)
   real     ::  netsht           ! net short wave radiation (w/m2)
   real     ::  vpcnp            ! canopy air saturated vapor pressure
   real     ::  le, ls           ! latent heat of vaporization/sublimation (j/kg)
   real     ::  lathsc           ! latent heat fluxes over canopy (w/m2)
   real     ::  senhsc           ! sensible heat fluxes over canopy (w/m2)
   real     ::  evapc            ! canopy evapotranspiration (w/m2)
   real     ::  snowd            ! snowpack depth (m)
   real     ::  rad              ! land surface net radiation (w/m2)
   real     ::  sfat             ! radiation attentuation factor
   real     ::  old_fswq         ! fswq before computing snow melt (m)
   real     ::  tmptsf4          ! ftsf**4
   real     ::  tmptgr4          ! tgrnd**4
   real     ::  tmp1             ! sncvfr*sigma
   real     ::  mfcwt            ! temporary canopy intercepted water (m h2o)
   real     ::  mfsmc(msl)       ! temporary soil moisture (liquid + ice) (mm)
   integer  ::  m                ! loop index
   real     ::  canopy_evap      ! function name
   real     ::  snow_density     ! function name
!
!          mxcvswq  ::  minimum swq for fully covering a grid (m)
!
   real, parameter  ::  mxcvswq=0.00
!
!     rhoh2o        :: water density (kg/m^3)
!     t0c           :: ice/water mix temperature (k)
!     sigma         ::  stefan-boltzmann constant (w/m2/k-4)
!
   real, parameter  ::  t0c=273.15, sigma=5.67e-8,rhoh2o=1.e3
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! initialize the output and other variables
! ----------------------------------------------------------------------
!
   lathsc = 0.0
   senhsc = 0.0
   evapc  = 0.0
   fstf   = 0.0
!
! ----------------------------------------------------------------------
! check for thin snowpack which only partially covers grid cell
! ----------------------------------------------------------------------
!
   if(fswq .lt. mxcvswq) then
     sncvfr = fswq/mxcvswq
   else 
     sncvfr = 1.0
   endif
!
! ----------------------------------------------------------------------
! if vegetation overstory exists
! ----------------------------------------------------------------------
!
   if(ost.eq.1) then
!
! ----------------------------------------------------------------------
! compute the snow interceptation over whole grid
! ----------------------------------------------------------------------
!
     tmp1 = sncvfr*sigma
     tmptsf4 = ftsf**4
     tmptgr4 = tgrnd**4
!
     cnplongw = longw + tmp1*tmptsf4 + (sigma-tmp1)*tmptgr4
!
     if(fcsn.gt.0.or.sfall.gt.0) then
       call vic_snow_intercept(dtime,    fflai, rhoair,   tair,               &
                               vpair,      vpd,  ra(2),   u(2),               &
                             netshts, cnplongw,  rfall,  sfall,               &
                                pgcm,     fcwt,   fcsn,   frtf,               &
                                fstf,   lathsc, senhsc)
     else                            ! canopy no snow
!
! ----------------------------------------------------------------------
! calculate the net radiation at the canopy surface, using the canopy 
!    temperature. the outgoing longw is subtracted twice, because the 
!    canopy radiates in two directions
! ----------------------------------------------------------------------
!
       longout = sigma * tair**4
       rad   = netshtv + cnplongw - 2*longout
!
       evapc = canopy_evap(msl,  vtype,   dtime,   fflai,                      &
                          fcwt,  rfall,    fsmc,    fsic,                      &
                           wcr,    wpw,     frt,     rad,                      &
                        rhoair,  vpair,     vpd,   ra(3),                      &
                       netshtv,   tair,    pgcm,   mfcwt,                      &
                          frtf,  mfsmc)
!
       fcwt = mfcwt
       do m = 1,msl
         fsmc(m) = mfsmc(m)
       enddo
!
       if(tair.gt.t0c) then
         le = (2.501 - 0.002361 * (tair-t0c)) * 1.0e6
         evapc = le*evapc*rhoh2o
       else
         ls = (677.-0.07*(tair-t0c))*4.1868*1000.0
         evapc = ls*evapc*rhoh2o
       end if
!
     end if
   else
     frtf = rfall
     fstf = sfall
   end if            
!
#ifdef DBGVIC
   print *,'in vic_snow_solver lathsc/evapc=',lathsc,evapc
   print *,'in vic_snow_solver       senhsc=',senhsc
#endif
!
   if((fswq+fstf).gt.0.0)then
!      
! ----------------------------------------------------------------------
! check for thin snowpack which only partially covers grid cell
! the snow melt is modeled over the sncvfr
! ----------------------------------------------------------------------
!
     if((fswq+fstf) .lt. mxcvswq) then
       sncvfr = (fswq+fstf)/mxcvswq
     else
       sncvfr = 1.0
     endif
!
! ----------------------------------------------------------------------
! shortwave radiation attentuation factor because of canopy
! ----------------------------------------------------------------------
!
     if(ost.eq.1)then
#define CORR
#ifdef CORR
       sfat = exp(-0.1*veg_sfat(vtype)*fflai)
#else
       sfat = exp(-veg_sfat(vtype)*fflai)
#endif
     else
       sfat = 1.0
     endif
!
     netsht = netshts*sfat
!
! ----------------------------------------------------------------------
! compute snow melt over the snow_coverage_fraction (sncvfr)
! convert the snow variables to the fraction of sncvfr
! ----------------------------------------------------------------------
!
     fstf = fstf/sncvfr
     fswq = fswq/sncvfr
     fsfw = fsfw/sncvfr
     fpkw = fpkw/sncvfr
!
     old_fswq = fswq
!
! ----------------------------------------------------------------------
! calculate snow depth (h.b.h. 7.2.1) 
! ----------------------------------------------------------------------
!
     snowd = 1000. * fswq / frsn        ! in the units of m
!
! ----------------------------------------------------------------------
! call snow pack accumulation and ablation algorithm
! ----------------------------------------------------------------------
!
     call vic_snow_melt(dtime,   vtype,   month,  rhoair,                      &
                        ra(3),    u(3),   vpair,     vpd,                      &
                         tair,   tgrnd,  netsht,   longw,                      &
                         pgcm,    frtf,    fstf,   snowd,                      &
                         frsn,    fswq,    ftsf,    ftpk,                      &
                         fsfw,    fpkw,  fsnowm,   laths,                      &
                        senhs)
!
     laths = laths*sncvfr + lathsc + evapc
     senhs = senhs*sncvfr + senhsc
!
#ifdef DBGVIC
     print *,'in vic_snow_solver laths=',laths,lathsc,evapc
     print *,'in vic_snow_solver senhs=',senhs,senhsc
#endif
!
! ----------------------------------------------------------------------
! compute snow parameters 
! ----------------------------------------------------------------------
!
     if(fswq .gt. 0.) then

#ifdef DBGVIC
       print *,'in vic_snow_solver befsndens',frsn,fstf,old_fswq,snowd
#endif
!
! ----------------------------------------------------------------------
! calculate snow density
! ----------------------------------------------------------------------
!
       frsn = snow_density(fstf, tair, old_fswq, snowd,                        &
                            dtime, ftsf)

#ifdef DBGVIC
       print *,'in vic_snow_solver after snow_density',frsn
#endif
     end if
!
! ----------------------------------------------------------------------
! compute snow melt over the snow_coverage_fraction (sncvfr)
! convert the snow parameters to sncvfr
! ----------------------------------------------------------------------
!
     fswq  = fswq*sncvfr
     fsfw  = fsfw*sncvfr
     fpkw  = fpkw*sncvfr
     fsnowm= fsnowm*sncvfr
!
   else
!
! ----------------------------------------------------------------------
! ground snow not present
! ----------------------------------------------------------------------
!
     laths = lathsc + evapc
     senhs = senhsc
   end if
!
   if(fswq.le.1.e-6) then
     frsn = 50.0
     fswq =  0.0
     fsfw =  0.0
     fpkw =  0.0
     ftsf = min(tair,t0c)
     ftpk = min(tgrnd,t0c)
   end if
!
   return
   end subroutine vic_snow_solver
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vic_snow_intercept(dtime,     lai, rhoair,    tair,              &
                                 vpair,     vpd,     ra,    wind,              &
                                netsht,   longw,  rainf,   snowf,              &
                                  pgcm,  cwater, ccsnow,  raintf,              &
                                snowtf,  lathsc, senhsc)
!-------------------------------------------------------------------------------
!
! subprogram: snow_intercept
!
! abstract: calculates the interception and subsequent release
!   of snow by the forest canopy using energy balance
!
! program history log:
!   2003-06~09  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from 4.0.3 uw vic (snow_intercept.c)
!
! comments: only the top canopy layer is taken into account for snow
!   interception. snow interception by lower canopy is
!   disregarded. rain water can be intercepted by lower canopy
!   layers. of course:  no vegetation -> no interception
!
!-------------------------------------------------------------------------------
   use constant, only : t0c_, rhoh2o_
!-------------------------------------------------------------------------------
#include <vartyp.h>
! ------------------- input variables ----------------------------------
   real    ::  dtime         ! time step (second)
   real    ::  lai           ! vegetation leaf area index
   real    ::  rhoair        ! air density (kg/m^3)
   real    ::  tair          ! air temperature (k)
   real    ::  vpair         ! actual vapor pressure of air (pa)
   real    ::  vpd           ! air vapor pressure deficit (pa)
   real    ::  ra            ! vic_exch resistance (s/m)
   real    ::  wind          ! wind speed related to ra(2) (m/s)
   real    ::  netsht        ! net short wave (w/m^2)
   real    ::  longw         ! long wave (w/m^2)
   real    ::  rainf         ! rain fall (m/time step)
   real    ::  snowf         ! snow fall (m/time step)
   real    ::  pgcm          ! air pressure (pa)
! ----------------------------------------------------------------------
!
! ------------------- modified variables -------------------------------
   real    ::  cwater        ! canopy intercepted water (m)
   real    ::  ccsnow        ! canopy intercepted snow (m)
! ----------------------------------------------------------------------
!
! ---------------------- output variables ------------------------------
   real    ::  raintf        ! rain throughfall (m/time step)
   real    ::  snowtf        ! snow throughfall (m/time step)
   real    ::  lathsc        ! leaf snow latent heat (w/m2) (+: air to snow)
   real    ::  senhsc        ! leaf snow sensible heat flux (w/m2) (+: to snow)
!-----------------------------------------------------------------------
!
! ---------------------- local variables -------------------------------
   real, parameter  ::  cpair=1004.0        ! specific heat at constant pressure 
                                            !    of air (j/deg/k)
   real, parameter  ::  lai_sm=0.0005       ! lai snow multiplier (m)
   real, parameter  ::  lai_wf=0.0002       ! leaf water factor 
                                            !    for interception storage (m)
   real, parameter  ::  liq_wc=0.035        ! liquid water capacity
   real, parameter  ::  min_is=0.005        ! min snow interception storage (m)
   real, parameter  ::  sigma=5.67e-8       ! stefan-boltzmann constant (w/m^2/k^-4)
   real, parameter  ::  ch_water = 4186.8e3 ! volumetric heat capacity (j/(m3*c))
                                            !    of water
   real, parameter  ::  lf=3.337e5          ! latent heat of freezing (j/kg) at 0c
   real    ::  imax1         ! maxium water intecept regardless of temp (m)
   real    ::  maxint        ! maximum rainf interception storage (m)
   real    ::  blownsnow     ! depth of snow blown of the canopy (m) 
   real    ::  delsnwint     ! change in swe of snow interceped on branches (m)
   real    ::  maxwatint     ! water interception capacity (m)  
   real    ::  maxsnwint     ! snow interception capacity (m) 
   real    ::  intrf         ! fraction of intercepted water which is liquid
   real    ::  intsf         ! fraction of intercepted water which is solid
   real    ::  overload      ! temp variable to calculated structural overloading
   real    ::  ta            ! air temperature (c)
   real    ::  tcanopy       ! canopy temperature (k)
   real    ::  longout       ! emitted longwave radiation (w/m2)
   real    ::  advectede     ! advaected energy (w/m2)
   real    ::  rad           ! net radiation (w/m2)
   real    ::  essnow        ! vapor pressure (pa)
   real    ::  vaporflux     ! vapor flux from canopy (m/s) [+ air -> canpoy]
   real    ::  vaporf        ! vapor from canopy (m) (+ air -> canopy)
   real    ::  ls            ! latent heat of sublimation (j/kg)
   real    ::  refreezee     ! refreeze energy (j)
   real    ::  potsnwmt      ! potential snow melt (m)
   real    ::  exssnwmt      ! excess snow melt (m)
   real    ::  tmpintstr     ! temporary intercepted storage (m)
   real    ::  releasmss     ! released snow mass from leaf (m)
   real    ::  drip          ! leaf drip (m)
   real    ::  inicsnow      ! initial canopy snow storage (m)
   real    ::  tmp_csw       ! temporary canopy water (m)
   real    ::  svp           ! function name
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! initialize variables
! ----------------------------------------------------------------------
!
   raintf = 0.0
   snowtf = 0.0
!
   imax1  = 4.0* lai_sm * lai          ! m
   maxint = lai_wf * lai               ! m
!
! ----------------------------------------------------------------------
! during snow falling, all the snow and water intercepted by vegetation,
! and soil moisture are averaged over the whole grid cell,
! and we only compute the interception over the wet district (mu=1.0);
! but we still need separate canopy snow interception into the wet
! and dry districts in case of rain fall after snowing, resulting in 
! heterogeneity of snow coverage. ji 2003
! ----------------------------------------------------------------------
!      
   if(snowf .gt. 0.) then  ! calculate snow interception
!
! ----------------------------------------------------------------------
! determine the maximum snow interception water equivalent.           
! kobayashi, d., 1987, snow accumulation on a narrow board,           
! cold regions science and technology, (13), pp. 239-245. figure 4.
! ----------------------------------------------------------------------
!
     ta = tair - t0c_
!
     if (ta.lt.-1.0.and.ta.gt.-3.0) then
       maxsnwint = (ta*3.0/2.0) + (11.0/2.0)
     else if (ta.gt.-1.0) then
       maxsnwint = 4.0
     else
       maxsnwint = 1.0
     end if
!      
! ----------------------------------------------------------------------
! therefore lai_ratio decreases as temp decreases
! ----------------------------------------------------------------------
!
     maxsnwint = maxsnwint*lai_sm*lai
!  
! ----------------------------------------------------------------------
! calculate snow interception
! ----------------------------------------------------------------------
!
     delsnwint = (1.0-ccsnow/maxsnwint)*snowf
     if ((delsnwint+ccsnow).gt.maxsnwint) then
       delsnwint = maxsnwint - ccsnow
     end if
     if (delsnwint .lt. 0.0) delsnwint = 0.0
!         
! ----------------------------------------------------------------------
! reduce the amount of intercepted snow if snowing, windy and cold.
! (< -3 to -5 c). schmidt and troendle 1992 western snow conference paper.
! ringyo shikenjo tokyo, #54, 1952. bulletin of the govt. 
! forest exp. station, govt. forest exp. station, meguro, tokyo, japan.
! forstx 634.9072 r475r #54. page 146, figure 10. 
! ----------------------------------------------------------------------
!
     if(ta.lt.-3.0.and.delsnwint.gt.0.0.and. wind.gt.1.0) then
       blownsnow = (0.2 * wind-0.2)*delsnwint
       if(blownsnow.ge.delsnwint) blownsnow = delsnwint
       delsnwint = delsnwint - blownsnow
     end if
!  
! ----------------------------------------------------------------------
! now update snowf and total accumulated intercepted snow amounts
! ----------------------------------------------------------------------
!
     if ((ccsnow+delsnwint).gt.imax1) delsnwint=0.0 
!  
! ----------------------------------------------------------------------
! pixel depth snow through fall (snowtf)
! ----------------------------------------------------------------------
!
     snowtf = snowf - delsnwint
!
! ----------------------------------------------------------------------
! physical depth
! ----------------------------------------------------------------------
!
     ccsnow  = ccsnow + delsnwint
   end if
!
   if(rainf.gt.0) then
!
! ----------------------------------------------------------------------
! calc amount of rain intercepted on branches and stored in intercepted snow
!
! if snow fall and rain fall mixing, newmu = 1.0, just compute rain 
! interception over the wet district; otherwise, we need both the 
! dry and wet districts. however, since there is now rain fall over
! the dry district, we also only compute rain fall interception over
! the wet district. ji 2003
!
! before a new rain fall event occuring, the canopy water and rainthrough 
! fall are averaged over the whole grid cell; during the rain fall event, 
! the canopy water and rainthrough fall are updated according to the 
! fraction of the wet district on where rain fall falls. ji 2003
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! physical depth
! ----------------------------------------------------------------------
!
     maxwatint=liq_wc*ccsnow+maxint
!         
     if ((cwater+rainf).le.maxwatint)then
       cwater = cwater + rainf                ! physical depth
     else
       raintf = cwater + rainf - maxwatint
       cwater = maxwatint                     ! physical depth
     end if
   end if
!
! ----------------------------------------------------------------------
! at this point we have calculated the amount of snow fall intercepted
! and the amount of rain fall intercepted. these values have been 
! appropriately subtracted from snow fall and rain fall to determine 
! snow throughfall and rain throughfall. however, we can end up with the 
! condition that the total intercepted rain plus intercepted snow is 
! greater than the maximum bearing capacity of the tree regardless of air 
! temp (imax1). the following routine will adjust intrain and intsnow 
! by triggering mass release due to overloading. of course since intrain
! and intsnow are mixed, we need to slough them off as fixed fractions  
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! trigger structural unloading
! ----------------------------------------------------------------------
!
   if(lai.gt.0.0) then
     tmp_csw = cwater + ccsnow
     if (tmp_csw.gt.imax1)then 
       overload = tmp_csw - imax1
       intrf = cwater/tmp_csw
       intsf = ccsnow/tmp_csw
       cwater= cwater - overload*intrf
       ccsnow = ccsnow  - overload*intsf
       raintf= raintf + overload*intrf
       snowtf= snowtf + overload*intsf
     end if
   end if
!
   inicsnow = ccsnow
!
! ----------------------------------------------------------------------
! the canopy temperature is assumed to be equal to the air temperature if 
! the air temperature is below t0c, otherwise the canopy temperature is 
! equal to t0c (273.15 k)
! ----------------------------------------------------------------------
!
   if(tair .gt. t0c_) then
     tcanopy = t0c_
   else
     tcanopy = tair
   end if
!
! ----------------------------------------------------------------------
! calculate the net radiation at the canopy surface, using the canopy 
! temperature.  the outgoing longw is subtracted twice, because the 
! canopy radiates in two directions
! ----------------------------------------------------------------------
!
   longout = sigma * tcanopy**4
   rad     = netsht + longw - 2*longout
   if(rad.lt.0.) then        ! in case of setting too high tcanopy
     rad = 0.0
     tcanopy = ((netsht + longw)/2.0/sigma)**(0.25)
     if(tcanopy.lt.(min(t0c_,tair)-5.0)) then
       tcanopy = min(t0c_,tair)-5.0
       longout = sigma * tcanopy**4
       rad     = netsht + longw - 2*longout
     endif
   endif
!
! ----------------------------------------------------------------------
! calculate the vapor mass flux between the canopy and the surrounding 
! air mass
! ----------------------------------------------------------------------
!
   essnow = svp(tcanopy) 
!
! ----------------------------------------------------------------------
! added division by 10 to incorporate change in canopy resistance due
! to smoothing by intercepted snow
! ----------------------------------------------------------------------
!
   vaporflux = rhoair*(0.622/pgcm)*(vpair-essnow)/ra/10.0
   vaporflux = vaporflux/rhoh2o_
!
   if (vpd.eq.0.and.vaporflux.lt.0.0) vaporflux = 0.0
!
! ----------------------------------------------------------------------
! calculate the latent heat flux 
! ----------------------------------------------------------------------
!
   ls = (677.-0.07*(tcanopy-t0c_))*4.1868*1000.0
   lathsc = ls * vaporflux * rhoh2o_
!
! ----------------------------------------------------------------------
! calculate the sensible heat flux
! in original c-vic there is no division of 10 after ra
! however, sometimes senhsc is too big, so the division is added.
! also, the lathsc and senhsc should not too high. so the maximum
! of lathsc and senhsc of 50 w/m2 is assigned.
! for completeness, we should use the energy balance mode for computing
! the canopy energy balance. june 2004/ji
! ----------------------------------------------------------------------
!
   senhsc = rhoair * cpair * (tair-tcanopy)/ra/10.0
!  
! ----------------------------------------------------------------------
! calculate the advected energy
! ----------------------------------------------------------------------
!
   advectede = (ch_water*(tair-tcanopy)*rainf)/dtime
!
! ----------------------------------------------------------------------
! control the quality of computation
! ----------------------------------------------------------------------
!
   lathsc = max(min(lathsc, 50.0),-50.0)
   senhsc = max(min(senhsc, 50.0),-50.0)
   advectede = max(min(advectede, 50.0),-50.0)
!
! ----------------------------------------------------------------------
! calculate the amount of energy available for refreezing
! ----------------------------------------------------------------------
!
   refreezee = senhsc + lathsc + rad + advectede
   refreezee = refreezee * dtime             ! units in (j/m2)
!
! ----------------------------------------------------------------------
! if refreezee is positive it means energy is available to melt
! the intercepted snow in the canopy. if it is negative, it means that 
! intercepted water will be refrozen
! ----------------------------------------------------------------------
!
   vaporf = vaporflux * dtime   ! m/s -> m (+ air -> canopy)
!
   if (refreezee .gt. 0.0) then
!
     if (-(vaporf) .gt. cwater) then
       vaporf =  - cwater
       cwater = 0.0
     else
       cwater = cwater + vaporf
     end if
!      
     potsnwmt = min((refreezee/lf/rhoh2o_), ccsnow)
     maxwatint = liq_wc*ccsnow + maxint
     if ((cwater + potsnwmt) .le. maxwatint) then
       ccsnow  = ccsnow - potsnwmt
       cwater = cwater + potsnwmt
       potsnwmt = 0
     else
       exssnwmt = potsnwmt + cwater - maxwatint
       ccsnow = ccsnow + cwater - maxwatint
       cwater = maxwatint
         
       tmpintstr = 0.0
       releasmss = 0.0
       drip = 0.0
!         
       if (snowtf.gt.0.0.and.inicsnow.le.min_is) then
!
! ----------------------------------------------------------------------
! water in excess of maxwatint has been generated.  if it is 
! snowing and there was little intercepted snow at the beginning 
! of the time step ( <= min_is), then allow the
! snow to melt as it is intercepted
! ----------------------------------------------------------------------
!
         drip  = exssnwmt
         ccsnow = ccsnow - exssnwmt
       else 
!
! ----------------------------------------------------------------------
! else, snowthroughfall = 0.0 or snowthroughfall > 0.0 and there is a 
! substantial amount of intercepted snow at the beginning of the time 
! step ( > min_is). snow melt may generate mass release.
! ----------------------------------------------------------------------
!
         tmpintstr = exssnwmt
       end if
!         
       call vic_snow_canopy(ccsnow,tmpintstr,releasmss,drip)
!
       raintf = raintf + drip
       snowtf = snowtf + releasmss
     end if
!
! ----------------------------------------------------------------------
! if intercepted snow has melted, add the water it held to drip
! ----------------------------------------------------------------------
!
     maxwatint = liq_wc*ccsnow + maxint
     if (cwater .gt. maxwatint) then
       drip = cwater - maxwatint
       cwater = maxwatint
       raintf = raintf + drip
     end if
!  
   else                      ! else (refreezee <= 0.0) 
!      
! ----------------------------------------------------------------------
! refreeze as much surface water as you can
! ----------------------------------------------------------------------
!
     if (refreezee .ge. (- cwater*lf))then
       ccsnow = ccsnow + abs(refreezee)/lf
       cwater = cwater - abs(refreezee)/lf
       refreezee = 0.0
     else 
!
! ----------------------------------------------------------------------
! all of the water in vegetation has been frozen.
! ----------------------------------------------------------------------
!
       ccsnow = ccsnow + cwater
       cwater = 0.0
     end if
!      
     if (-(vaporf) .gt. ccsnow)then
       vaporf = - ccsnow
       ccsnow = 0.0
     else
       ccsnow = ccsnow + vaporf
     end if
   end if
!
   return
   end subroutine vic_snow_intercept
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vic_snow_melt(dtime,    vtype,   month, rhoair,                  &
                           Ra,     wind,   vpair,    vpd,                      &
                         Tair,    Tgrnd,  netsht,  longw,                      &
                         pgcm,     frTF,    fsTF,  snowd,                      &
                         frsn,      swq,     Tsf,    Tpk,                      &
                      sfwater,  pkwater,   snowm,  latHs,                      &
                        senHs)
!-------------------------------------------------------------------------------
!
! subprogram: vic_snow_melt
!
! abstract: calculate snow accumulation and melt using an energy 
!   balance approach for a two layer snow model
!
! program history log:
!   2003-06~09  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from 4.0.3 uw vic (vic_snow_melt.c)
!
!-------------------------------------------------------------------------------
#ifdef VICLSM1
#include <vartyp.h>
! ------------------- input variables ----------------------------------
! 
! -- model basic parameter
!
   real dtime          ! time step (second)
   integer month       ! month of current time step
!
! -- atmosphere data
!
   real fsTF           ! Amount of snow (m/time step)
   real frTF           ! Rain through fall (m/time step)
   real Tair           ! Air temperature (K)
   real rhoair         ! air density  (kg/m3)
   real vpair          ! Actual vapor pressure of air (Pa)
   real vpd            ! Vapor pressure deficit (Pa)
   real Ra             ! vic_exch resistance (s/m)
   real wind           ! Wind speed (m/s)
   real pgcm           ! Air pressure (Pa)

   real netsht         ! net shortwave radiation (w/m2)
   real longw          ! Longwave radiation (w/m2)
!
! -- vegetation parameters
!
   integer vtype       ! vegetation type
!
! -- land surface parameters
!
   real Tgrnd          ! ground surface temperature (K)
   real frsn           ! snow density (kg/m3)
   real snowd          ! snow depth (m)
! ----------------------------------------------------------------------
!
! ------------------- modified variables -------------------------------
   real swq            ! Snow water equivalent at current pixel (m)
   real Tsf            ! Temperature of snow pack surface layer (K)
   real Tpk            ! Temperature of snow pack (K)
   real pkwater        ! Liquid water content of snow pack (m)
   real sfwater        ! Liquid water in snow surface (m)
! ----------------------------------------------------------------------
!
! ------------------- output variables ---------------------------------
   real snowm          ! Amount of snowpack melt (m/time step)
   real latHs          ! Latent heat exchange at surface (W/m2)
   real senHs          ! Sensible heat exchange at surface (W/m2)
   real Refreeze       ! Refreeze energy (W/m2)
! ----------------------------------------------------------------------
!
! ------------------- common block -------------------------------------
#include <vic_snoweb.h>
! ----------------------------------------------------------------------
!
! ---------------------- local variables -------------------------------
   real Ice, Le, Ls
   real oldTsf           ! Snow Surface temp for previous time (K)
   real sfswq            ! Snow water equivalent on snow surface (m)
   real Pkswq            ! Snow pack snow water equivalent (m)
   real SfCC             ! Cold content of surface layer snow pack (J)
   real PkCC             ! Cold content of snow pack (J)
   real SnowTFCC         ! Cold content of new snowfall (J)
   real DeltPkswq        ! Change in snow water equivalent of pack (m)
   real DeltPkCC         ! Change in cold content of the pack
   real vaporf           ! vapor flux (m/time step) (+: to snow surface)

   real MaxLiqWat        ! Maximum liquid water content of pack (m)
   real Qnet             ! Net energy exchange at the surface (W/m2)
   real RefrozenW        ! Amount of refrozen water (m)
   real tmp              ! temporary variable

   character*80 ctrfct   ! determine function

   real func_snoweb      ! function name
   real root_brent       ! function name

   real rhoh2o           ! water density (kg/m^3)
   real Lf               ! Latent heat of freezing (J/kg) at 0!
   real t0c              ! ice/water mix temperature (K)
   real CH_ICE           ! Volumetric heat capacity (J/(m3*C)) of ice
   parameter (Lf=3.337e5, t0c=273.15, CH_ICE=2100.0e3,rhoh2o=1.e3)

   real SNOW_DT          ! snow surface temperature DT (K)
   real MAX_SS           ! MAX SURFACE SWE (first snow layer) (m)
   real LIQ_WC           ! LIQUID WATER CAPACITY (faction)
   parameter (SNOW_DT=1.0, MAX_SS=0.125, LIQ_WC=0.035)
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! convert units ( K -> C)
! ----------------------------------------------------------------------
!
   Tair  = Tair - t0c             ! k -> C
   Tgrnd = Tgrnd - t0c            ! k -> C
   Tsf   = Tsf - t0c              ! k -> C
   Tpk   = Tpk - t0c              ! k -> C
!
   oldTsf= Tsf
!
! ----------------------------------------------------------------------
! Initialize snowpack variables
! ----------------------------------------------------------------------
!
   Ice = swq - pkwater - sfwater
   if(Ice .lt. 0) Ice = 0.0
!
! ----------------------------------------------------------------------
! Reconstruct snow pack  (MAX SURFACE SWE (first snow layer) = 0.125 m)
! ----------------------------------------------------------------------
!
   if (Ice .gt. MAX_SS) then
     Sfswq = MAX_SS
   else
     Sfswq = Ice
   end if
!
   Pkswq = Ice - Sfswq
! 
! ----------------------------------------------------------------------
! Calculate cold contents
! ----------------------------------------------------------------------
!
   SfCC = CH_ICE * Sfswq * Tsf
   PkCC = CH_ICE * Pkswq * Tpk
!
   if(fsTF.gt.0) then            ! snowthrough fall
     if (Tair .gt. 0.0) then
       SnowTFCC = 0.0
     else
       SnowTFCC = CH_ICE * fsTF * Tair
     end if
!
! ----------------------------------------------------------------------
! Distribute fresh snowfall
! ----------------------------------------------------------------------
!
     if (fsTF .gt. (MAX_SS - Sfswq))then
       DeltPkswq = Sfswq + fsTF - MAX_SS
!         
       if (DeltPkswq .gt. Sfswq) then
         DeltPkCC = SfCC + (fsTF - MAX_SS)/                                    &
                    fsTF * SnowTFCC
       else
         DeltPkCC = DeltPkswq/Sfswq * SfCC
       end if
!         
       Sfswq = MAX_SS
       SfCC  = SfCC + SnowTFCC - DeltPkCC
       Pkswq = Pkswq + DeltPkswq
       PkCC  = PkCC + DeltPkCC
     else
       Sfswq = Sfswq + fsTF
       SfCC  = SfCC + SnowTFCC
     end if
   end if
!
   Tsf = SfCC/(CH_ICE * Sfswq)
!
   if (Pkswq .gt. 0.0) then   
     Tpk = PkCC/(CH_ICE * Pkswq)
   else
     Tpk = min(Tgrnd, min(0.0,Tair))
   end if
!
! ----------------------------------------------------------------------
! Adjust ice and snow surf_water
! ----------------------------------------------------------------------
!
   Ice = Ice + fsTF
   sfwater = sfwater + frTF
!
! ----------------------------------------------------------------------
! intialize vicsnoweb common block
! ----------------------------------------------------------------------
!
   Edtime  = dtime
   Emonth  = month
   EfrTF   = frTF
   ETair   = Tair
   Erhoair = rhoair
   Evpair  = vpair
   Evpd    = vpd
   ERa     = Ra
   Ewind   = wind
   Epgcm   = pgcm
   Enetsht = netsht
   Elongw  = longw
   Evtype  = vtype
   EoldTsf = oldTsf
   ETgrnd  = Tgrnd
   Efrsn   = frsn
   Esnowd  = snowd
   Esfswq  = sfswq
   Esfwater= sfwater
!
! ----------------------------------------------------------------------
! Calculate the snow surface energy balance for Tsf = 0.0 
! ----------------------------------------------------------------------
!
   Qnet = func_snoweb(0.0)
!
   Refreeze = Orefreeze
   latHs = OlatHs
   senHs = OsenHs
!
! ----------------------------------------------------------------------
! If Qnet is 0.0, then set the surface temperature to 0.0
! ----------------------------------------------------------------------
!
   if (abs(Qnet).lt.1.e-6)then
     Tsf = 0.0
     if (Refreeze .ge. 0.0) then    ! positive: frozen surface water
       RefrozenW = Refreeze/(Lf * rhoh2o)* dtime
       if (RefrozenW .gt. sfwater) then
         RefrozenW = sfwater
         Refreeze = RefrozenW * Lf * rhoh2o / dtime
       end if
       Sfswq   = Sfswq + RefrozenW
       Ice     = Ice + RefrozenW
       sfwater = sfwater - RefrozenW
       snowm   = 0.0
     else
!
! ----------------------------------------------------------------------
! Calculate snow melt
! ----------------------------------------------------------------------
!
       snowm  = abs(Refreeze)/(Lf *rhoh2o)*dtime      ! m/time step
     end if
!
! ----------------------------------------------------------------------
! Convert vapor mass flux to a depth per timestep and adjust surf_water
! ----------------------------------------------------------------------
!
     Le = 2.501e6
     vaporf = latHs/(Le * rhoh2o)
!
! ----------------------------------------------------------------------
! Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990)
! ----------------------------------------------------------------------
!
     vaporf = vaporf * dtime   ! m/s -> m/time step (+: to snow surface)
!     
     if (sfwater .lt. (-vaporf)) then
       vaporf  = -sfwater
       sfwater = 0.0
     else
       sfwater = sfwater + vaporf
     end if
!
! ----------------------------------------------------------------------
! If snowm < Ice, there was incomplete melting of the pack
! assume: melting snow from bottom to snow surface
! ----------------------------------------------------------------------
!
     if (snowm .lt. Ice) then
       if (snowm .le. Pkswq) then
         sfwater = sfwater + snowm
         Pkswq   = Pkswq - snowm
         Ice     = Ice   - snowm
       else 
         sfwater = sfwater + snowm + pkwater
         pkwater = 0.0
         Pkswq   = 0.0
         Ice     = Ice - snowm
         Sfswq   = Ice
       end if
!    
! ----------------------------------------------------------------------
! Else, Snowm > Ice and there was complete melting of the pack
! ----------------------------------------------------------------------
!
     else
       snowm  = Ice
       sfwater= sfwater + Ice + pkwater
       Ice    = 0.0
       Sfswq  = 0.0
       Pkswq  = 0.0
       pkwater= 0.0
       Tsf    = min(0.0,Tair)
       Tpk    = min(0.0,Tair)
     end if
!
! ----------------------------------------------------------------------
! Else, SnowPackEB(T=0.0) <= 0.0 
! ----------------------------------------------------------------------
!
   else  
!
! ----------------------------------------------------------------------
! Calculate surface layer temperature using "Brent method"
! ----------------------------------------------------------------------
!
     ctrfct = 'func_snoweb'
     Tsf = root_brent(0.0, Tsf-SNOW_DT, ctrfct)
!
     Refreeze = Orefreeze
     latHs = OlatHs
     senHs = OsenHs
!
! ----------------------------------------------------------------------
! since we iterated, the surface layer is below freezing and no snowm
! ----------------------------------------------------------------------
!
     snowm  = 0.0
!    
! ----------------------------------------------------------------------
! Since updated snow_temp < 0.0, all of the liquid water in the surface
! layer has been frozen   
! ----------------------------------------------------------------------
!
     SfSwq   = SfSwq + sfwater
     Ice     = Ice   + sfwater
     sfwater = 0.0
!    
! ----------------------------------------------------------------------
! Convert vapor mass flux to a depth per timestep and adjust surf_water
! ----------------------------------------------------------------------
!
     Ls = (677. - 0.07 * Tsf) * 4.1868e3
!
! ----------------------------------------------------------------------
! Calculate latent heat flux
! Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990)
! ----------------------------------------------------------------------
!
     vaporf = latHs/(Ls * rhoh2o)
     vaporf = vaporf * dtime      ! m/time step (+: to snow surface)
!
     if (Sfswq .lt. -vaporf) then ! in one time the maxi vaporf=Sfswq
       vaporf = -Sfswq
       Sfswq  = 0.0
       Ice    = Pkswq
     else
       Sfswq = Sfswq + vaporf
       Ice   = Ice + vaporf
     end if
   end if
!  
! ----------------------------------------------------------------------
! Done with iteration etc, now Update the liquid water content of the
! surface layer  
! ----------------------------------------------------------------------
!
   MaxLiqWat = LIQ_WC * Sfswq
!
   if  (sfwater .gt. MaxLiqWat) then
     snowm = sfwater - MaxLiqWat
     sfwater  = MaxLiqWat
   else
     snowm = 0.0
   end if
!  
! ----------------------------------------------------------------------
! Refreeze liquid water in the pack. variable 'Refreeze' is the
! heat released to the snow pack if all liquid water were refrozen.  
! if Refreeze < PackCC then all water IS refrozen  PackCC always <=0.0
! WORK IN PROGRESS: This energy is NOT added to MeltEnergy, since this
!  does not involve energy transported to the pixel. Instead heat from
!  the snow pack is used to refreeze water
! ----------------------------------------------------------------------
!  
! ----------------------------------------------------------------------
! add surface layer outflow to pack liquid water
! ----------------------------------------------------------------------
!
   pkwater = pkwater + snowm
!
   Refreeze = pkwater * Lf * rhoh2o
!
! ----------------------------------------------------------------------
! calculate energy released to freeze
! ----------------------------------------------------------------------
!
   if (PkCC .lt. -Refreeze) then      ! cold content not fully depleted
     Pkswq = Pkswq + pkwater         ! refreeze all water and update
     Ice   = Ice + pkwater
     pkwater = 0.0
     if (Pkswq .gt. 0.0) then
!            PkCC = Pkswq * CH_ICE * Tpk + Refreeze ! original vic equation
       PkCC = PkCC + Refreeze
       Tpk  = PkCC / (CH_ICE * Pkswq)
     else 
       Tpk = 0.0
     end if
   else  
!
! ----------------------------------------------------------------------
! cold content has been either exactly satisfied or exceeded. If
! PackCC = refreeze then pack is ripe and all pack water is
! refrozen, else if energy released in refreezing exceeds PackCC 
! then exactly the right amount of water is refrozen to satify PackCC.
! The refrozen water is added to PackSwq and Ice
! ----------------------------------------------------------------------
!
     Tpk      = 0.0
     DeltPkswq = -PkCC/(Lf * rhoh2o)  ! the amount of pkwater -> Pkswq
!
     pkwater = pkwater - DeltPkSwq
     Pkswq = Pkswq + DeltPkSwq
     Ice = Ice + DeltPkSwq 
   end if
!  
! ----------------------------------------------------------------------
! Update the liquid water content of the pack
! ----------------------------------------------------------------------
!
   MaxLiqWat = LIQ_WC * Pkswq
   if (pkwater .gt. MaxLiqWat) then
     snowm = pkwater - MaxLiqWat
     pkwater = MaxLiqWat
   else
     snowm = 0.0
   end if
!
! ----------------------------------------------------------------------
! Update snow properties
! ----------------------------------------------------------------------
!
   Ice  = Pkswq + Sfswq
!
   if (Ice .gt. MAX_SS) then
     SfCC   = CH_ICE * Tsf * Sfswq
     PkCC   = CH_ICE * Tpk * Pkswq
     if (Sfswq .gt. MAX_SS) then
       tmp = SfCC * (Sfswq-MAX_SS)/Sfswq
       PkCC  = PkCC + tmp
       SfCC  = SfCC - tmp
       Pkswq = Pkswq + Sfswq - MAX_SS
       Sfswq = MAX_SS
     else
       tmp = PkCC * (MAX_SS - Sfswq) / Pkswq
       PkCC  = PkCC - tmp
       SfCC  = SfCC + tmp
       Pkswq = Pkswq + Sfswq - MAX_SS
       Sfswq = MAX_SS
     end if
!
     Tpk = PkCC / (CH_ICE * Pkswq)
     Tsf = SfCC / (CH_ICE * Sfswq)
   else 
     Pkswq = 0.0
     Tpk = 0.0
   end if
!
   swq = Ice + pkwater + sfwater
!
! ----------------------------------------------------------------------
! convert units ( c -> K)
! ----------------------------------------------------------------------
!
   Tair  = Tair + t0c             ! c -> k
   Tgrnd = Tgrnd + t0c            ! c -> k
   Tsf   = Tsf + t0c              ! c -> k
   Tpk   = Tpk + t0c              ! c -> k
!
#endif
   return
   end subroutine vic_snow_melt
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vic_snow_canopy(ccsnow, tmpintstr, releasmss, drip)
!-------------------------------------------------------------------------------
!
! subprogram: vic_snow_canopy
!
! abstract: calculates mass release of snow from canopy
!
! program history log:
!   2003-06~09  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from 4.0.3 uw vic 
!                        (vic_snow_canopy.c and part of dhsvm)
!
!-------------------------------------------------------------------------------
#include <vartyp.h>
! ---------------------- local variables -------------------------------
   real, parameter  ::  min_is = 0.005    ! min_interception_storage (m)
   real             ::  threshold, maxrelease, tempdrip, tmprelmass
!-----------------------------------------------------------------------
!  
! ------------------- input/output variables ---------------------------
   real             ::  ccsnow       ! canopy intercepted snow (m)
   real             ::  tmpintstr    ! temporary intercepted storage (m)
   real             ::  releasmss    ! released snow mass from leaf (m)
   real             ::  drip         ! leaf drip (m)
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! if the amount of snow in the canopy is greater than some minimum
! value, min_is, then calculate mass release and drip
! ----------------------------------------------------------------------
!
55 if (ccsnow .gt. min_is) then
     threshold  = 0.10 * ccsnow
     maxrelease = 0.17 * ccsnow
!    
! ----------------------------------------------------------------------
! if the amount of vic_snow_melt after interception, vic_snow_melt, 
!  is >= the theshhold then there is mass release.  
!  if vic_snow_melt is < the treshhold
!  then there is no mass release but that water remains in
!  tempintstorage which will be augmented during the next
!  compute period
! ----------------------------------------------------------------------
!    
     if (tmpintstr .ge. threshold) then
       drip  = drip  + threshold
       ccsnow = ccsnow - threshold
       tmpintstr = tmpintstr - threshold
!
       if (ccsnow .lt. min_is) then
         tmprelmass = 0.0
       else
         tmprelmass = min((ccsnow-min_is),maxrelease)
       end if
!
       releasmss = releasmss + tmprelmass
       ccsnow = ccsnow - tmprelmass
       go to 55
     else
       tempdrip = min(tmpintstr, ccsnow)
       drip  = drip + tempdrip
       ccsnow = ccsnow - tempdrip
     end if
   else
!
! ----------------------------------------------------------------------
! interceptedsnow < min_is) if the amount of snow in
!   the canopy is less than some minimum value, min_is,
!   then only melt can occur and there is no mass release.
! ----------------------------------------------------------------------
!
     tempdrip = min(tmpintstr, ccsnow)
     drip = drip + tempdrip
     ccsnow = ccsnow - tempdrip
     tmpintstr = 0.0
   end if
!
   return
   end subroutine vic_snow_canopy
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function snow_density(new_snow, tair, swq, depth, dtime, tsf)
!-------------------------------------------------------------------------------
!
! abstract: this function computes the snow density based on the day of 
!           the year. density information comes from a plot of seasonal 
!           variation of typical snow densities found in bras 
!           (figure 6.10, p 258). the equation was developed by regressing 
!           against the curve for southern manitoba, so this routine should 
!           be modified if used outside the plains of south central canada, 
!           and the north central us. 
!
! program history log:
!   2003-06   ji chen   modified from 4.0.3 uw vic 
!   2003-07   ji chen
!   2003-08   ji chen
!   2003-09   ji chen
!
!-------------------------------------------------------------------------------
#include <vartyp.h>
! ------------------- input variables ----------------------------------
   real  ::  dtime      ! model time step
   real  ::  new_snow   ! amount of new snow (m/time step)
   real  ::  tair       ! air temperature (k)
   real  ::  swq        ! snow water equivalent (m)
   real  ::  depth      ! snow depth (m)
   real  ::  tsf        ! snow surface temperature (k)
! ----------------------------------------------------------------------
!
! ------------------- output variables ---------------------------------
   real  ::  snow_density ! (kg/m3)
! ----------------------------------------------------------------------
!
! ---------------------- local variables -------------------------------
!
!       new_sdens   ::  new snow density(kg/m3)
!       max_change  ::  maximum fraction of snowpack depth change
!
   real, parameter  ::  new_sdens=50., max_change=0.9
!  
!        eta0       ::  viscosity of snow at t = 0c and density = 0
!                   ::  used in calculation of true viscosity (ns/m2)
!        c5         ::  constant used in snow viscosity calculation,
!                   ::  taken from snthrm.89 (/c)
!        c6         ::  constant used in snow viscosity calculation
!        g          ::  gravitational accelleration (m/(s^2))
!
   real, parameter  ::  eta0 = 3.6e6, c5 = 0.08, c6 = 0.021, g=9.81
!
   real, parameter  ::  t0c=273.15, rhoh2o=1.e3
!
   real             ::  ta, depth_new
   real             ::  density, pdepth, density_new
   real             ::  overburden, viscosity, ddepth
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! compaction of snow pack by new snow fall *** bras pg. 257 !
! ----------------------------------------------------------------------
!
   if(new_snow .gt. 0) then
!
! ----------------------------------------------------------------------
! estimate density of new snow based on air temperature 
! ----------------------------------------------------------------------
!
     ta = (tair-t0c) * 9. / 5. + 32.
     if(ta .gt. 0) then
       density_new = new_sdens + 1000.                                       &
                * (ta / 100.) * (ta / 100.)
     else
       density_new = new_sdens
     end if
!
     if(depth .gt. 0.) then
!
! ----------------------------------------------------------------------
! compact current snowpack by weight of new snowfall
! ----------------------------------------------------------------------
!
       pdepth = ((new_snow / 0.0254)*(depth / 0.0254)) /                      &
            (swq/0.0254)* (((depth/0.0254)/10.)**0.35) * 0.0254
!  
! ----------------------------------------------------------------------
! pdepth cannot be greater than depth 
! ----------------------------------------------------------------------
!
       if (pdepth .ge. depth) then
         pdepth = max_change * depth
       end if
!
       depth_new = 1000.0 * new_snow / density_new
       depth = depth - pdepth + depth_new
       swq   = swq + new_snow
       if(depth.gt.0) then
         density = 1000. * swq / depth
       else
         density = density_new
       endif
!
     else
!
! ----------------------------------------------------------------------
! no snowpack present, so snow density equals that of new snow
! ----------------------------------------------------------------------
!
       density = density_new
       swq = swq + new_snow
     end if
!
   else
     if(depth.gt.0) then
       density = 1000. * swq / depth
     else
       density = density_new
     endif
   end if
!
! ----------------------------------------------------------------------
! densification of the snow pack due to aging
! based on snthrm89 r. jordan 1991 - used in barts dhsvm code
! ----------------------------------------------------------------------
!
   if(swq.gt.0) then
     depth       = 1000. * swq / density
     overburden  = 0.5 * g * rhoh2o * swq
     viscosity   = eta0 * exp(-c5 * (tsf-t0c) + c6 * density)
     ddepth      = -overburden / viscosity * depth * dtime
     if(abs(ddepth).gt.0.05*depth*dtime/86400.0) then
       ddepth = - 0.05*depth*dtime/86400.0
     endif
     depth       = depth + ddepth
!      
     if(depth .le. 0) depth = (1-max_change)*depth
!      
     density     = 1000. * swq / depth
   endif
!
   snow_density = density
!
   return
   end function snow_density
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vic_soil_moisture_solver(msl,    nsl,    dphnd,    dph,          &
                                     smxnd,  expnd,    bubnd,    sld,          &
                                       bkd,    qrt,     fsmc,     tn,          &
                                      smnd,  icend,    kapnd,   csnd)
!-------------------------------------------------------------------------------
!
! subroutine: vic_soil_moisture_solver
!
! abstract: this subroutine determines the moisture and ice contents
!        of each soil thermal node based on the current node temperature
!        and layer moisture content. thermal conductivity and volumetric
!        heat capacity are then estimated for each node based on the
!        division of moisture contents. soil thermal conductivity
!        calculated using johansens method.
!
! reference: farouki, o.t., "thermal properties of soils" 1986
! chapter 7: methods for calculating the thermal conductivity of soils
!     h.b.h. - refers to the handbook of hydrology.
!
! prgmmr:            ji chen   
! org:               ecpc/crd/sio/ucsd 
! date:              june, july, august & september 2003
! prgm history:      modified from 4.0.3 uw vic c version
!             (distribute_node_moisture_properties in soil_conduction.c)
! program history: vic_soil_moisture_solver.f is originally a part of 
!                   uw_vic (distribute_node_moisture_properties in 
!                           soil_conduction.c)
!
! program history log:
!   2003-09-01  ji chen                modified from 4.0.3 uw vic (runoff.c)
!   2008-08-01  kyeong-hee seol        debugged, scm option
!
!-------------------------------------------------------------------------------
#include <vartyp.h>
! ------------------- input variables ----------------------------------
   integer  ::  msl           ! number of soil moisture layers
   integer  ::  nsl           ! number of soil thermal nodes
   real     ::  dphnd(nsl)    ! thermal node thickness (m)
   real     ::  dph(msl)      ! soil layer thickness (m)
   real     ::  smxnd(nsl)    ! thermal node maximum moisture content (mm/mm)
   real     ::  expnd(nsl)    ! thermal node exponential (n/a)
   real     ::  bubnd(nsl)    ! thermal node bubbling pressure (cm)
   real     ::  sld(msl)      ! soil density (kg m-3)
   real     ::  bkd(msl)      ! soil layer bulk density (kg m-3)
   real     ::  qrt(msl)      ! soil qrt content (fraction)
   real     ::  fsmc(msl)     ! soil layer moisture (mm)
   real     ::  tn(nsl)       ! thermal node temperature (k)
! ----------------------------------------------------------------------
!
! ------------------- output variables ---------------------------------
   real     ::  smnd(nsl)     ! thermal node moisture content (mm/mm)
   real     ::  icend(nsl)    ! thermal node ice content (mm/mm)
   real     ::  kapnd(nsl)    ! thermal node thermal conductivity (w m-1 k-1)
   real     ::  csnd(nsl)     ! thermal node heat capacity (j m-3 k-1)
! ----------------------------------------------------------------------
!
! ---------------------- local variables -------------------------------
   real     ::  lsum, zsum
   real     ::  sfc_max_unfwat, soil_condvty, heat_capcty ! function names
   integer  ::  lidx, nidx
   logical  ::  past_bottom
   real, parameter  ::  t0c = 273.15
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! initialize the computation
! ----------------------------------------------------------------------
!
   past_bottom = .false.
   lidx = 1
   lsum = 0.
   zsum = 0.
!
   do nidx = 1,nsl
     zsum = zsum + dphnd(nidx)
     if(zsum.gt.lsum.and..not.past_bottom) then
       lsum = lsum + dph(lidx)
       lidx = lidx + 1
       if( lidx .ge. msl ) then
         past_bottom = .true.
         lidx = msl
       end if
!
       do while(zsum.gt.lsum)
         lsum = lsum + dph(lidx)
         lidx = lidx + 1
         if( lidx .ge. msl ) then
           past_bottom = .true.
           lidx = msl
         end if
       end do
     end if
!
!  node on layer boundary                  
!
     if(zsum.eq.lsum.and.nidx.ne.1                                             &
       .and.lidx.ne.1)then
       smnd(nidx) = (fsmc(lidx-1)/dph(lidx-1)                                  &
                   + fsmc(lidx)/dph(lidx))/1000.0/2.
     else                            ! node completely in layer
       smnd(nidx) = fsmc(lidx)/dph(lidx)/1000.0
     end if
!
     if(tn(nidx).lt.t0c) then        ! compute moisture/ice contents
       icend(nidx) = smnd(nidx) -                                              &
       sfc_max_unfwat(tn(nidx),smxnd(nidx),                                    &
                  bubnd(nidx),expnd(nidx))
!
       if(icend(nidx).lt.0) icend(nidx)=0
!
! ----------------------------------------------------------------------
! compute thermal conductivity
! ----------------------------------------------------------------------
!
       kapnd(nidx) = soil_condvty(smnd(nidx),                                  &
                 smnd(nidx)-icend(nidx),                                       &
                 sld(lidx),bkd(lidx),qrt(lidx))
     else                          !compute moisture and ice contents
       icend(nidx)   = 0
!
! ---------------------------------------------------------------------
! compute thermal conductivity
! ----------------------------------------------------------------------
!
       kapnd(nidx)= soil_condvty(smnd(nidx),                                   &
                 smnd(nidx), sld(lidx), bkd(lidx), qrt(lidx))
!
     end if
!
! ----------------------------------------------------------------------
! compute volumetric heat capacity 
! ----------------------------------------------------------------------
!
     csnd(nidx) = heat_capcty(bkd(lidx)/sld(lidx),                             &
                   smnd(nidx)-icend(nidx),icend(nidx))
   end do
!
   return
   end subroutine vic_soil_moisture_solver
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function heat_capcty(soil_fract, water_fract, ice_fract)
!-------------------------------------------------------------------------------
!
! source file:       heat_capcty.f
!
! purpose:           this function calculates the soil volumetric heat 
!                    capacity based on the fractional volume of its 
!                    component parts.
!
! program history log:
!     03-11-01  ji chen
!     modified from 4.0.3 uw vic c version
!     (volumetric_heat_capacity in soil_conduction.c)
!
! notes:
!       constant values are volumetric heat capacities in j/m^3/k
!       soil value is for clay or quartz - assumed for all other types
!
!-------------------------------------------------------------------------------
#include <vartyp.h>
!
! input variables 
!
   real  ::  soil_fract    ! fraction of soilvolume composed of actual soil
   real  ::  water_fract   ! fraction of soilvolume composed of liquid water
   real  ::  ice_fract     ! fraction of soilvolume composed of ice
!
! output variables 
!
   real  ::  heat_capcty   ! heat capacity (j m-3 k-1)
!
! local variables 
!
   real  ::  cs
   real  ::  organic_fract
!-----------------------------------------------------------------------
!
   organic_fract = 0.0
!
   cs = 2.0e6 * (soil_fract - organic_fract)
   cs = cs + 4.2e6 * water_fract
   cs = cs + 1.9e6 * ice_fract
   cs = cs + 2.7e6 * organic_fract
   cs = cs + 1.3e3 * (1.-(soil_fract+water_fract+                              &
        ice_fract+organic_fract))
!
   heat_capcty = cs
!
   return
   end function heat_capcty 
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function soil_condvty(moist, wu, rhosoil, rhobulk, quartz)
!-------------------------------------------------------------------------------
!
! abstract:
!   - soil thermal conductivity using johansens method
!
! program history log: 
!   2003-06,07,08,09  ji chen     modified from 4.0.3 uw vic 
!                    (soil_conductivity in soil_conduction.c)
!
! references: 
!     farouki, o.t., "thermal properties of soils" 1986
!     chapter 7: methods for calculating the thermal conductivity of
!     soils.     h.b.h. - refers to the handbook of hydrology.
!
!  porosity = n = porosity
!  ratio = sr = fractionaldegree of saturation
!  all k values are conductivity in w/mk
!  wu is the fractional volume of unfrozen water
!
!-------------------------------------------------------------------------------
#include <vartyp.h>
! ------------------- input variables ----------------------------------
   real  ::  moist         ! total moisture content (mm/mm)
   real  ::  wu            ! liquid water content (mm/mm)
   real  ::  rhosoil       ! soil density (kg m-3)
   real  ::  rhobulk       ! soil bulk density (kg m-3)
   real  ::  quartz        ! soil quartz content (fraction)
! ----------------------------------------------------------------------
!
! ------------------- output variables ---------------------------------
   real  ::  soil_condvty  ! soil thermal conductivity (w m-1 k-1)
!-----------------------------------------------------------------------
!
! ---------------------- local variables -------------------------------
   real  ::  ke
   real  ::  ksat
   real  ::  ks            ! thermal conductivity of solid (w/mk)
                           ! function of quartz content
   real  ::  kdry
   real  ::  sr            ! fractional degree of saturation
   real  ::  k
   real  ::  porosity
!
!               ki  ::  thermal conductivity of ice (w/mk)
!               kw  ::  thermal conductivity of water (w/mk)
!
   real, parameter  ::  ki=2.2, kw = 0.57
! ----------------------------------------------------------------------
!
   kdry = (0.135*rhobulk+64.7)/(rhosoil - 0.947*rhobulk)
!
   if(moist.gt.0.) then
     porosity = 1.0 - rhobulk / rhosoil
     sr = moist/porosity
     ks = (7.7**quartz) * (2.2**(1.0-quartz))
!
     if(wu.eq.moist) then             !soil unfrozen
       ksat = (ks**(1.0-porosity)) * (kw**porosity)
       ke = 0.7 * log10(sr) + 1.0
     else                             !soil frozen
       ksat = (ks**(1.0-porosity))*(ki**(porosity-wu))*(kw**wu)
       ke = sr
     end if
!
     k = (ksat-kdry)*ke+kdry
     if(k.lt.kdry) k=kdry
   else
     k = kdry
   end if
!
   soil_condvty = k
!
   return
   end function soil_condvty
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vic_surface_flux(msl,    nsl,  dtime,  month,                    &
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
                      ftnd,  falbd,   lath,   senh,                            &
                     grndh,   frtf, fsnowm, fsnowe,                            &
#ifdef SMP_RA2SFC
                    gdelth,  netsw, sncvfr)
#else
                    gdelth, sncvfr)
#endif
!-------------------------------------------------------------------------------
!
! subprogram: surface_fluxes
!
! abstract: this routine computes surface fluxes, and solves 
!   the snow accumulation and ablation algorithm. solutions
!   are for the current snow band and vegetation type
!
! program history log:
!   2003-06~09  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from 4.0.3 uw vic (surface_fluxes.c)
!
!-------------------------------------------------------------------------------
   use vic_veglib
#include <vartyp.h>
! ------------------- input variables ----------------------------------
!
! -- model basic parameters
!
   integer  ::  msl              ! number of soil layer
   integer  ::  nsl              ! number of soil thermal nodes
   integer  ::  month            ! current month
   real     ::  dtime            ! time step (second)
!
! -- atmosphere variables
!
   real     ::  rfall            ! rainfall (m/time step)
   real     ::  sfall            ! snowfall (m/time step)
   real     ::  pgcm             ! pressure (pa) (vic: use elev(m) to get p)
   real     ::  tgcm             ! air temperature (k)
   real     ::  flwds            ! long wave (w/m^2)
   real     ::  sols             ! short wave (w/m^2)
#ifdef SMP_RA2SFC
   real     ::  netsw            ! net short wave (w/m^2)
#endif
   real     ::  vpair            ! actual vapor pressure of air (pa)
   real     ::  vpd              ! vapor pressure deficit (pa)
   real     ::  rhoair           ! air density (kg/m^3)
   real     ::  u(3)             ! wind speed (m/s)
   real     ::  ra(3)            ! vic_exch resistance (s/m)
!
! -- soil parameters
!
   real     ::  binf             ! vic infiltration parameter
   real     ::  dph(msl)         ! soil layer thickness (m)
   real     ::  qrt(msl)         ! quartz content of soil (fraction)
   real     ::  bkd(msl)         ! bulk density of soil layer (kg/m^3)
   real     ::  sld(msl)         ! soil particle density (kg/m^3)
   real     ::  wcr(msl)         ! ~70% of field capacity (mm)
   real     ::  wpw(msl)         ! wilting point soil moisture (mm)
   real     ::  smr(msl)         ! residual moisture (mm)
   real     ::  dphnd(nsl)       ! soil node thickness (m)
   real     ::  smxnd(nsl)       ! maximum soil moisture at soil node (m3/m3)
   real     ::  expnd(nsl)       ! parameter for ksat with soil moisture (n/a)
   real     ::  bubnd(nsl)       ! bubbling pressure of soil (cm)
   real     ::  alpnd(nsl)       ! coef for computing soiltemp profile
   real     ::  betnd(nsl)       ! coef for computing soiltemp profile
   real     ::  gamnd(nsl)       ! coef for computing soiltemp profile
!
! -- vegetation parameters
!
   integer  ::  vtype            ! vegetation class
   real     ::  frt(msl)         ! root content (fraction)
   real     ::  fflai            ! leaf area index
!
! -- land surface variable
!
   integer  ::  lstsn         ! surface snow age (number of model time step)
! ----------------------------------------------------------------------
!
! ------------------- modified variables -------------------------------
!
   real     ::  fswq             ! snow water equivalent (m h2o)
   real     ::  frsn             ! snow density (kg/m^3)
   real     ::  ftsf             ! snow surface temperature (k)
   real     ::  ftpk             ! snow pack temperature (k)
   real     ::  fsfw             ! snow surface water equivalent (m h2o)
   real     ::  fpkw             ! snow pack water equivalent (m h2o)
   real     ::  fcwt             ! canopy intercepted water (m h2o)
   real     ::  fcsn             ! canopy intercepted snow (m h2o)
   real     ::  fsmc(msl)        ! soil moisture (liquid + ice) (mm)
   real     ::  fsic(msl)        ! soil ice (mm)
   real     ::  ftnd(nsl)        ! soil temperature profile (k)
! ----------------------------------------------------------------------
!
! ------------------- output variables ---------------------------------
!
   real     ::  falbd            ! land surface albedo (fraction)
   real     ::  lath             ! latent heat flux (w/m^2) (+: to surface)
   real     ::  senh             ! sensible heat flux (w/m^2) (+: to surface)
   real     ::  grndh            ! ground heat flux (w/m^2) (+: to surface)
   real     ::  frtf             ! rain throughfall (m/time step)
   real     ::  fsnowm           ! snow melt (m/time step)
   real     ::  fsnowe           ! snow evap (w/m2)
   real     ::  gdelth           ! ground heat storage (w/m^2) (+: to surface)
   real     ::  sncvfr           ! snow cover fraction (fraction)
! ----------------------------------------------------------------------
!
! ---------------------- local variables -------------------------------
!
   real     ::  falbds           ! snow related surface albedo (fraction)
   real     ::  falbdv           ! without snow surface albedo (fraction)
   real     ::  tgrnd            ! ground surface temperature (k)
   real     ::  laths            ! latent heat flux over snow surface (w/m^2)
   real     ::  senhs            ! sensible heat over snow surface (w/m^2)
   real     ::  netshts          ! net shortwave radiation for snow ground (w/m2)
   real     ::  netshtv          ! net shortwave radiation without snow (w/m2)
   integer  ::  ost              ! vegetation overstory (1: with, 0: no)
   real     ::  snowage          ! snow age (day)
   real     ::  snow_albedo      ! function name
   integer, parameter  ::  secpday=86400
! ----------------------------------------------------------------------
!
#ifdef DBGVIC
   write(6,*) 'vic debug -- enter vic_surface_flux',fsmc,fsic
!
#endif
! ----------------------------------------------------------------------
! initialize variables
! ----------------------------------------------------------------------
!
   ost = veg_ost(vtype)
   lath  = 0.0
   senh  = 0.0
   grndh = 0.0
   frtf  = 0.0
   fsnowm = 0.
   fsnowe = 0.
!
   tgrnd  = ftnd(1)       ! ground surface temperature
   laths  = 0.
   senhs  = 0.
   sncvfr = 0.
   gdelth = 0.
!
   falbds = 0.85          ! new snow albedo
   falbdv = veg_alb(month, vtype)
!
   if(fswq.gt.0.or.sfall.gt.0.0)then
!
! ----------------------------------------------------------------------
! compute snow pack albedo 
! ----------------------------------------------------------------------
!
     if(lstsn.lt.0) then            ! due to snow existing on ground
       write(6,*) 'warning in vic_surface_flux for setting lstsn =',lstsn
       write(6,*) ' in vic_surface_flux fswq, sfall',fswq,sfall
       lstsn = 0
     endif
!      
     snowage = float(lstsn)*dtime/float(secpday)
     falbds  = snow_albedo(snowage, ftsf)
   end if
!
   netshts = (1. - falbds)*sols
   netshtv = (1. - falbdv)*sols
!
#ifdef SMP_RA2SFC
! use the net radiation which is passed from rad_diurnal_cycle
!
   print*,' vic_surface_flux org netsw',netshtv
   print*,'         new netsw',netsw

   netshtv = netsw

   print*,'vic_surface_flux albedo',falbdv
   print*,'vic_surface_flux DNSW UPSW NSW',sols,falbdv*sols,netshtv
!
#endif
! ----------------------------------------------------------------------
! solve for overstory canopy hydrology, snow pack accumulation & ablation
! ----------------------------------------------------------------------
!
   if(fswq.gt.0.or.sfall.gt.0.or.(fcsn.gt.0.and.ost.eq.1)) then
!
#ifdef DBGVIC
     write(6,*) 'vic debug in vic_surface_flux before vic_snow_solver'
     write(6,*) ' fswq=',fswq,' fcsn=',fcsn,' sfall=',sfall
     write(6,*) ' bef vic_snow_solver fsmc,fsic',fsmc,fsic
!
#endif
     call vic_snow_solver(msl,     dtime,    month,     rfall,                 &
                        sfall,      pgcm,     tgcm,     flwds,                 &
                      netshts,   netshtv,   rhoair,     vpair,                 &
                          vpd,        ra,        u,       wcr,                 &
                          wpw,     vtype,      ost,       frt,                 &
                        fflai,     tgrnd,     fsmc,      fsic,                 &
                         fswq,      frsn,     ftsf,      ftpk,                 &
                         fsfw,      fpkw,     fcwt,      fcsn,                 &
                         frtf,    fsnowm,   sncvfr,     laths,                 &
                        senhs)

     fsnowe = - laths   ! positive upward for fsnowe (w/m2)
!
#ifdef DBGVIC
     write(6,*) ' after vic_snow_solver fsmc,fsic',fsmc,fsic
#endif
   else
!
#ifdef DBGVIC
     write(6,*) 'vic in vic_surface_flux no snow cover '
#endif
     frtf = rfall
!
   end if
!
   falbd = falbds*sncvfr + (1.-sncvfr)*falbdv
!
! ----------------------------------------------------------------------
! solve energy balance components for ground
! ----------------------------------------------------------------------
!
#ifdef DBGVIC
   write(6,*)  'vic debug in vic_surface_flux before vic_surface_temp ra=',ra
   write(6,*) 'in vic_surface_flux fflai',fflai,' vtype',vtype
!
#endif
   call vic_surface_temp(dtime,  month,                                        &
                  frtf,   pgcm,   tgcm,  flwds,                                &
               netshtv,  vpair,    vpd, rhoair,                                &
                  u(1),  ra(1),   binf,    dph,                                &
                   qrt,    bkd,    sld,    wcr,                                &
                   wpw,    smr,  dphnd,  smxnd,                                &
                 expnd,  bubnd,  alpnd,  betnd,                                &
                 gamnd,  vtype,    frt,  fflai,                                &
                sncvfr,   fswq,   frsn,   ftsf,                                &
                  fcwt,   fcsn,   fsmc,   fsic,                                &
                  ftnd,   lath,   senh,  grndh,                                &
                gdelth)
!
#ifdef DBGVIC
   write(6,*)  'in vic_surface_flux lath=',lath,' laths=',laths
   write(6,*)  'in vic_surface_flux senh=',senh,' senhs=',senhs
!
#endif
   lath = (1.0-sncvfr)*lath + laths
   senh = (1.0-sncvfr)*senh + senhs
!
   if(abs(lath).gt.800) then
     write(6,*)'** warning for lath in vic_surface_flux ',lath, laths, sncvfr
   endif
!
   if(abs(senh).gt.800) then
     write(6,*)'** warning for senh in vic_surface_flux ',senh, senhs, sncvfr
!         stop
   endif
!
#ifdef DBGVIC
   write(6,*)  'vic debug -- end of vic_surface_flux'
!
#endif
   return
   end subroutine vic_surface_flux
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vic_surface_temp(dtime,  month,                                  &
                      frTF,   pgcm,   Tair,   longw,                           &
                    netsht,  vpair,    vpd,  rhoair,                           &
                      wind,     ra,   binf,     dph,                           &
                       qrt,    bkd,    sld,     wcr,                           &
                       wpw,    smr,  dphnd,   smxnd,                           &
                     expnd,  bubnd,  alpnd,   betnd,                           &
                     gamnd,  vtype,    frt,   fflai,                           &
                    sncvfr,   fswq,   frsn,    ftsf,                           &
                      fcwt,   fcsn,   fsmc,    fsic,                           &
                      ftnd,   latH,   senH,   grndH,                           &
                    gdelth)
!-------------------------------------------------------------------------------
!
! subprogram: vic_surface_temp
!
! abstract: calculates the surface temperature, in the
!   case of no snow cover.  Evaporation is computed using the
!   previous ground heat flux, and then used to comput latent 
!   heat in the energy balance routine.  Surface temperature
!   is found using the Frt Brent method (Numerical Recipies).
!
! program history log:
!   2003-06~09  ji chen (ecpc/crd/sio/ucsd) 
!                        modified from 4.0.3 uw vic (calc_surf_energy_bal.c)
!
!-------------------------------------------------------------------------------
   use varsfc, only : nsoil_,lsoil_ 
   use vic_veglib
#include <define.h>
!
#ifdef VICLSM1
#include <vartyp.h>
#include <vic_surfeb.h>
!-------------------------------------------------------------------------------
! ------------------- input variables ----------------------------------
! -- model basic parameter
!
   integer  ::  month    ! current month
   real     ::  dtime       ! time step (second)
!
! -- atmosphere data
!
   real     ::  pgcm        ! pressure (pa)
   real     ::  Tair        ! air temperature (k)
   real     ::  rhoair      ! air density (kg/m3)
   real     ::  vpair       ! Actual vapor pressure of air (Pa)
   real     ::  vpd         ! vapor pressure deficit (Pa)
   real     ::  ra          ! vic_exch resistance (s/m)
   real     ::  wind        ! wind speed (m/s)
   real     ::  netsht      ! net shortwave radiation (w/m2)
   real     ::  longw       ! longwave radiation (w/m2)
!
! -- soil parameters
!
   real     ::  binf        ! vic infiltration parameter (n/a)
   real     ::  dph(msl)    ! soil moisture layer thickness (m)
   real     ::  wcr(msl)    ! ~70% of field capacity (mm)
   real     ::  wpw(msl)    ! wilting point soil moisture (mm)
   real     ::  smr(msl)    ! residual moisture (mm)
   real     ::  sld(msl)    ! soil particle density (kg/m^3)
   real     ::  bkd(msl)    ! bulk density of soil layer (kg/m^3)
   real     ::  qrt(msl)    ! quartz content of soil (fraction)
   real     ::  dphnd(nsl)  ! soil thermal node thicknesses (m)
   real     ::  alpnd(nsl)  ! coef for computing soiltemp profile
   real     ::  betnd(nsl)  ! coef for computing soiltemp profile
   real     ::  gamnd(nsl)  ! coef for computing soiltemp profile
   real     ::  smxnd(nsl)  ! maximum soil moisture (m3/m3)
   real     ::  bubnd(nsl)  ! bubbling pressure of soil (cm)
   real     ::  expnd(nsl)  ! parameter for ksat with soil moisture (n/a)
   real     ::  kapnd(nsl)  ! soil layer thermal conductivity (w/m/k)
   real     ::  csnd(nsl)   ! soil layer heat capacity (j/m^3/k)
   real     ::  smnd(nsl)   ! soil moisture at node (m3/m3)
   real     ::  icend(nsl)  ! soil ice at node (m3/m3)
!
! -- vegetation parameters
!
   integer  ::  vtype       ! vegetation type
   real     ::  fflai       ! leaf area index (fraction)
   real     ::  frt(msl)    ! root content (fraction)
!
! -- land surface variables
!
   real     ::  sncvfr      ! snow cover fraction (fraction)
   real     ::  fcsn        ! snow on vegetation (m)
   real     ::  ftsf        ! snow surface temperature (k)
   real     ::  fswq        ! snow water equivalent (m)
   real     ::  frsn        ! snow density (kg/m3)
   real     ::  snowd       ! snow depth (m)
!
! ----------------------------------------------------------------------
!
! ------------------- modified variables -------------------------------
!
   real     ::  fcwt        ! dew and rain trapped on vegetation (m)
   real     ::  frtf        ! rain throughfall (m/time step)
   real     ::  fsmc(msl)   ! soil moisture (liquid + ice) (mm)
   real     ::  fsic(msl)   ! soil ice (mm)
   real     ::  ftnd(nsl)   ! soil tmperature profiles (k)
! ----------------------------------------------------------------------
!
! ------------------- output variables ---------------------------------
!
   real     ::  lath        ! latent heat flux (w/m2)
   real     ::  senh        ! sensible heat flux (w/m2)
   real     ::  grndh       ! ground heat flux (w/m2)
   real     ::  gdelth      ! ground heat storage (w/m2)
!-----------------------------------------------------------------------
!
! ---------------------- local variables -------------------------------
!
   real     ::  d                ! vegetation displacement (m)
   real     ::  z0               ! surface roughness (m)
   real     ::  t0               ! land surface temperature (k)
   real     ::  d1               ! first soil layer thickness (m)
   real     ::  t_upper          ! up-boundary of temperature (k)
   real     ::  t_lower          ! low-boundary of temperature (k)
   real     ::  tsurf            ! surface temperature (k)
   real     ::  stabilityc       ! function name
   real     ::  root_brent       ! function name
   integer  ::  n, m             ! loop index
!
!     rhoh2o          ::  water density (kg/m^3)
!
   real, parameter    ::  surf_dt = 20.0,rhoh2o=1.e3, huge_resist = 1.e2
!
   character(len=80)  ::  ctrfct  ! determine function
! ----------------------------------------------------------------------
!
#ifdef DBGVIC
   print *,'vic debug -- enter vic_surface_temp',vtype,month
!
#endif
! ----------------------------------------------------------------------
! correct vic_exch resistance for stability conditions
! ----------------------------------------------------------------------
!
!      if(vtype.gt.0) then
!         d  = veg_d(month, vtype)
!         z0 = veg_rough(month, vtype)
!      else
!         d  = 0
!         z0 = z0_soil
!      end if
!
   d = 0.0
!
!      if(fswq.gt.0.0) then
!         t0 = ftsf
!         z0 = 0.03
!      else
!
! because only ra is used for computing sensible and latent over
! no snow ground
!
   t0 = ftnd(1)
   z0 = 0.01
!
!      end if
!
#ifdef DBGVIC
   print *,'vic debug in vic_surface_temp before computing stability ',ra
#endif
!
   if (wind .gt. 0.0) then
     ra = ra/stabilityc(2.0, d, t0, tair, wind, z0)
   else
     ra = huge_resist
   end if  
!
#ifdef DBGVIC
   print *,'vic debug in vic_surface_temp after computing stability ',ra
!
#endif
! ----------------------------------------------------------------------
! find surface temperature using frt brent method
! ----------------------------------------------------------------------
!
#ifdef DBGVIC
   print *,'vic debug in vic_surface_temp before vic_soil_moisture_solver'
   print *,'vic_soil_moisture_solver variables',                               &
                       msl,   nsl, dphnd,   dph,                               &
                     smxnd, expnd, bubnd,   sld,                               &
                       bkd,   qrt,  fsmc,  ftnd
#endif
   call vic_soil_moisture_solver(msl,   nsl, dphnd,   dph,                     &
                               smxnd, expnd, bubnd,   sld,                     &
                                 bkd,   qrt,  fsmc,  ftnd,                     &
                                smnd, icend, kapnd,  csnd)

#ifdef DBGVIC
   print *, 'vic debug in vic_surface_temp after vic_soil_moisture_solver'
   print *,'vic_soil_moisture_solver variables',                               &
                 smnd, icend, kapnd,  csnd
#endif
!
! ----------------------------------------------------------------------
! calculate snow dph (h.b.h. 7.2.1)
! ----------------------------------------------------------------------
!
#ifdef DBGVIC
   print *,'check frsn ',frsn
#endif
!
   if(sncvfr.gt.1.e-6) then
     snowd = rhoh2o * fswq / frsn / sncvfr ! in the units of m
   else
     snowd = 0.0
   endif
!
   d1 = dph(1)
!
! ----------------------------------------------------------------------
! added for temporary backwards compatability
! ----------------------------------------------------------------------
!
#ifdef DBGVIC
   print *,'vic debug in vic_surface_temp before root_brent'
   print *,'dtime',dtime,pgcm,tair,rhoair,vpd,ra,wind,longw,netsht,
  1        binf,dph(1),wcr(1),wpw(1),smr(1),sld(1),bkd(1),
  2        qrt(1),frt(1),fsmc(1),fsic(1),dphnd(1),alpnd(1),
  3        betnd(1),gamnd(1),smxnd(1),bubnd(1),expnd(1),
  4        kapnd(1),csnd(1),smnd(1),icend(1),ftnd(1)
   print *,'vtpe ',vtype,fflai,frtf,fcwt,fcsn,ftsf,fswq,
  1        snowd
#endif
   t_upper = t0 + surf_dt
   t_lower = t0 - surf_dt
!
! ----------------------------------------------------------------------
! initialize vic_surfeb common block
! ----------------------------------------------------------------------
!
   edtime = dtime
   esncvfr= sncvfr
   epgcm  = pgcm
   etair  = tair
   erhoair= rhoair
   evpair = vpair
   evpd   = vpd
   era    = ra
   ewind  = wind
   elongw = longw
   enetsht= netsht
   ebinf  = binf
!   
   do m = 1,msl
     edph(m) = dph(m)
     ewcr(m) = wcr(m)
     ewpw(m) = wpw(m)
     esmr(m) = smr(m)
     esld(m) = sld(m)
     ebkd(m) = bkd(m)
     eqrt(m) = qrt(m)
     efrt(m) = frt(m)
     efsmc(m) = fsmc(m)
     efsic(m) = fsic(m)
   end do
!
   do n = 1,nsl
     edphnd(n) = dphnd(n)
     ealpnd(n) = alpnd(n)
     ebetnd(n) = betnd(n)
     egamnd(n) = gamnd(n)
     esmxnd(n) = smxnd(n)
     ebubnd(n) = bubnd(n)
     eexpnd(n) = expnd(n)
     ekapnd(n) = kapnd(n)
     ecsnd(n)  = csnd(n)
     esmnd(n)  = smnd(n)
     eicend(n) = icend(n)
     eftnd(n) = ftnd(n)
   end do
!
   evtype = vtype
   efflai = fflai
!   
   efrtf = frtf
   efcwt = fcwt
   efcsn = fcsn
   eftsf = ftsf
   efswq = fswq
   efrsn = frsn
   esnowd = snowd
!
   ctrfct = 'vic_funct_energy'
   tsurf = root_brent(t_upper, t_lower, ctrfct)
!
   fcwt = mfcwt
   frtf = mfrtf
!   
   do m = 1,msl
     fsmc(m) = mfsmc(m)
   end do
!
   do n = 1,nsl
     ftnd(n) = mftnd(n)
   end do
!
   lath = olath
   senh = osenh
   grndh = ogrndh
   gdelth = odelth
!
#ifdef DBGVIC
   print *,'vic debug -- end of vic_surface_temp  tsurf (k) ',tsurf
#endif
#endif
   return
   end subroutine vic_surface_temp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function snow_albedo(snowage, tsf)
!-------------------------------------------------------------------------------
!
! abstract: this function computes the snow pack surface albedo based 
!           on snow age and season, using the tables generated in 
!           snow_table_albedo.
!
! program history log:
!   2003-06  ji chen   modified from 4.0.3 uw vic 
!   2003-07  ji chen
!
!-------------------------------------------------------------------------------
#include <vartyp.h>
! ------------------- input variables ----------------------------------
   real             ::  snowage     ! surface snow age (day)
   real             ::  tsf         ! temperature of snow pack surface layer (k)
! ----------------------------------------------------------------------
!
! ------------------- output variables ---------------------------------
   real             ::  snow_albedo ! snow albedo
! ----------------------------------------------------------------------
!
! ---------------------- local variables -------------------------------
   real, parameter  ::  t0c = 273.15
   real, parameter  ::  new_snow_alb = 0.85
   real, parameter  ::  snow_alb_accum_a=0.94, snow_alb_accum_b=0.58
   real, parameter  ::  snow_alb_thaw_a=0.82, snow_alb_thaw_b=0.46
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! new snow 
! ----------------------------------------------------------------------
!
   if(snowage .eq. 0.0) then
     snow_albedo = new_snow_alb
!
! ----------------------------------------------------------------------
! aged snow: accumulation season 
! ----------------------------------------------------------------------
!
   else
     if(tsf .lt. t0c) then
       snow_albedo = new_snow_alb*snow_alb_accum_a**                           &
                     (snowage**snow_alb_accum_b)
!
! ----------------------------------------------------------------------
! melt season
! ----------------------------------------------------------------------
!
     else
       snow_albedo = new_snow_alb*snow_alb_thaw_a**                            &
                     (snowage**snow_alb_thaw_b)
     end if
   end if
!
   return
   end function snow_albedo
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vic_stop_run
!-------------------------------------------------------------------------------
!
! subprogram: vic_stop_run
!
! pbstract: terminate job
!
! program history log:
!   2003-11
!
!-------------------------------------------------------------------------------
#ifdef MP
#ifdef RMP
   call rmpabort
#else
   call mpabort
#endif
#else
   call abort
#endif
   end subroutine vic_stop_run
!-------------------------------------------------------------------------------
