#include <define.h>
#ifndef LWRMDC
   subroutine rad_lw_nasa_driver
   end subroutine rad_lw_nasa_driver
#else
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_driver(ipts,hlw,slwnc,tlwuc,                         &
#ifdef CLR
                 slwn0,tlwu0,                                                  &
#endif
                 pl2,temp,rh2o,qo3,taulw,cld0)
!-------------------------------------------------------------------------------
!
! subroutine: rad_lw_nasa_driver
!
! [06/07/00] fanglin yang
!     1) add "save" statement for data blocks.
!     2) change all "real*8" in longwave subroutine to "real".
!     3) fix nt=1 in rad_lw_nasa.F, and rename h11,h21,h81,c1 and o1
!        as hh11,hh21,hh81,cc1 and oo1 to avoid memeory confusion
!        for shared-memeroy computers
!
! [05/01/00] fanglin yang
!     computes ir fluxes due to water vapor, co2, o3, trace gases
!     (n2o, ch4, cfc11, cfc12, cfc22, co2-minor),clouds, and aerosols.
! --input variables
!     pl2    : sigma-layer pressure (mb), k=1-->toa
!     temp   : integer-layer temperature (k),k=1--> top layer,
!              with temp(lp1) = surface temperature
!     rh2o   : water vapor mixing ratio (g/g), k=1-->top layer
!     qo3    : o3 mixing ratio (g/g), k=1-->top layer
!     taulw  : cloud optical depth, k=1-->top layer
!     cld0   : cloud fraction, k=1-->top layer
! --output variables
!     hlw    : longwave cooling rate(k/day), negative
!     slwnc  : surface net lw flux (upward +) (1.e-3 w/m^2)
!     tlwuc  : top upward(+) lw flux (1.e-3 w/m^2)
!     slwn0  : clear-sky surface net lw flux(upward +) (1.e-3 w/m^2)
!     tlwu0  : clear-sky top upward(+) lw flux (1.e-3 w/m^2)
!
! structure: 
!     [rad_lw_nasa_driver] - [rad_lw_nasa]
!
! program history log:
!   2000-01-01  fanglin yang           implementation from uiuc chou scheme
!   2000-07-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel
   use rdparm
   use comconsts, only   : ictin,icbin
!-------------------------------------------------------------------------------
!
! -- input
!
   real                 ::  pl2(imax,lp1), temp(imax,lp1), rh2o(imax,l)
   real                 ::  qo3(imax,l),   taulw(imax,l), cld0(imax,l)
!
! -- output, radiative fluxex are in cm.g.s. system
!
   real                 ::  hlw(imax,l), slwnc(imax), tlwuc(imax)
   real                 ::  slwn0(imax), tlwu0(imax)
!
! --temporary, all radiative fluxes are in w/m^2, downward(+)
!
   real                 ::  pir(imax,lp1), ta(imax,l), wa(imax,l), oa(imax,l)
   real                 ::  ccld(imax,l),  taucl(imax,l), taual(imax,l,10)
   real                 ::  tg(imax)
   real                 ::  flxd(imax,lp1), flxu(imax,lp1), flx(imax,lp1)
   real                 ::  flxdc(imax,lp1),flxuc(imax,lp1),flc(imax,lp1)
   real                 ::  dfdts(imax,lp1),st4(imax)
   real                 ::  n2o
   integer              ::  n2oppm
!
! * if model top is above 30mb, please use high=.true.
!
   logical              ::  high,trace
!c data high/.true./
!c data trace/.false./
   high=.true.
   trace=.false.  !without ch4, cfcs, and n2o
!
! -- find cloud tops
!   level index separating high and middle clouds (ict)
!   level index separating middle and low clouds   (icb)
!  ict       =15 for 28 layer   !about 440mb for mls standard atmos.
!  icb       =19 for 28 layer   !about 702mb for mls standard atmos.
!
   ict=ictin
   icb=icbin
!
   co2ppm    =348.0
   ch4ppm    =1.75
   n2oppm    =0.28
   f22ppm    =0.2e-03      !ppmv
   f11ppm    =0.2145E-3    !cfc11 concentration in ppm
   f12ppm    =0.3711E-3    !cfc12 concentration in ppm
!
! -----exchange unit
!
   co2       =co2ppm*1.e-6
   ch4       =ch4ppm*1.e-6
   n2o       =n2oppm*1.e-6
   cfc11     =f11ppm*1.e-6
   cfc12     =f12ppm*1.e-6
   cfc22     =f22ppm*1.e-6
!
! -- exchange temporary variables
!
   do i = 1,ipts
     do k = 1,l
       pir(i,k)     =pl2(i,k)
       ta(i,k)      =temp(i,k)
       wa(i,k)      =rh2o(i,k)
       oa(i,k)      =qo3(i,k)
       ccld(i,k)    =cld0(i,k)
       taucl(i,k)   =taulw(i,k)
       do kk=1,10
         taual(i,k,kk)=0.    !aerosol optical depth, set to zero
       enddo
     enddo
     tg(i)       =temp(i,lp1)
     pir(i,lp1)  =pl2(i,lp1)
   enddo
!
   call rad_lw_nasa (ipts,imax,l,pir,ta,wa,oa,tg,co2,                          &
        n2o,ch4,cfc11,cfc12,cfc22,                                             &
        taucl,ccld,ict,icb,taual,                                              &
        high,trace,flxd,flxu,flx,flxdc,flxuc,flc,dfdts,st4)
!
!     ------------------------------------------------------------
!     ------------------------------------------------------------
!     constant for calculating heating rate (k/day)
!     hrate=daysecond*gravity/cpc=24*3600*9.80616/1004.0
!     where cpc is the specific heat of dry air (j/deg/kg)
!
   hrate = 843.8767171
!
   do i = 1,ipts
     do k = 1,l
       hlw(i,k)=0.01*hrate*(flx(i,k)-flx(i,k+1))/(pir(i,k+1)-pir(i,k))
     enddo
     slwnc(i) =-1000.*flx(i,lp1)
     slwn0(i) =-1000.*flc(i,lp1)
     tlwuc(i) =-1000.*flxu(i,1)
     tlwu0(i) =-1000.*flxuc(i,1)
   enddo
!
   return
   end subroutine rad_lw_nasa_driver
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa (ipts,m,np,prse,ta,wa,oa,ts,co2,                     &
              n2o,ch4,cfc11,cfc12,cfc22,                                       &
              taucl,fcld,ict,icb,taual,                                        &
              high,trace,flxd,flxu,flx,flxdc,flxuc,flc,dfdts,st4)
!-------------------------------------------------------------------------------
   use funct_rad_lw, only : rad_lw_nasa_b10exps, rad_lw_nasa_b10kdis,          &
                            rad_lw_nasa_cfcexps, rad_lw_nasa_cfckdis,          &
                            rad_lw_nasa_ch4exps, rad_lw_nasa_ch4kdis,          &
                            rad_lw_nasa_co2exps, rad_lw_nasa_co2kdis,          &
                            rad_lw_nasa_column,  rad_lw_nasa_comexps,          &
                            rad_lw_nasa_comkdis, rad_lw_nasa_conexps,          &
                            rad_lw_nasa_h2oexps, rad_lw_nasa_h2okdis,          &
                            rad_lw_nasa_n2oexps, rad_lw_nasa_n2okdis,          &
                            rad_lw_nasa_tablup
!-------------------------------------------------------------------------------
!
! this routine computes ir fluxes due to water vapor, co2, o3,           
!   trace gases (n2o, ch4, cfc11, cfc12, cfc22, co2-minor),           
!   clouds, and aerosols.
!  
! this is a vectorized code.  it computes fluxes simultaneously for
!   (m) soundings. 
!
! some detailed descriptions of the radiation routine are given in
!   chou and suarez (1994).
!
! there are options for computing cooling rate profiles.
!
!   if high = .true., transmission functions in the co2, o3, and the
!   three water vapor bands with strong absorption are computed using
!   table look-up.  cooling rates are computed accurately from the
!   surface up to 0.01 mb.
!   if high = .false., transmission functions are computed using the
!   k-distribution method with linear pressure scaling.  cooling rates
!   are not calculated accurately for pressures less than 20 mb.
!   the computation is faster with high=.false. than with high=.true.
!   if trace = .true., absorption due to n2o, ch4, cfcs, and the 
!   two minor co2 bands in the window region are included.
!   if trace = .false., absorption in minor bands are not included.
!
! the ir spectrum is divided into nine bands:
!   
!   band     wavenumber (/cm)   absorber         method
!
!    1           0 - 340           h2o            k/t
!    2         340 - 540           h2o            k/t
!    3         540 - 800       h2o,cont,co2,n2o   k,s,k/t,k
!    4         800 - 980       h2o,cont           k,s
!                              co2,f11,f12,f22 k,k,k,k,k
!    5         980 - 1100      h2o,cont,o3        k,s,t
!                              co2,f11            k,k
!    6        1100 - 1215      h2o,cont           k,s
!                              n2o,ch4,f12,f22 k,k,k,k,k
!    7        1215 - 1380      h2o,cont           k,s
!                              n2o,ch4            k,k
!    8        1380 - 1900          h2o            k/t
!    9        1900 - 3000          h2o            k 
!
! in addition, a narrow band in the 15 micrometer region is added to
!    compute flux reduction due to n2o
!
!    10        540 - 620       h2o,cont,co2,n2o   k,s,k,k
!
! note : "h2o" for h2o line absorption
!        "cont" for h2o continuum absorption
!        "k" for k-distribution method
!        "s" for one-parameter temperature scaling
!        "t" for table look-up
!
! band 3 (540-800/cm) is further divided into
!   3 sub-bands :
!
!   subband   wavenumber (/cm)
!
!    1          540 - 620
!    2          620 - 720
!    3          720 - 800
!
!---- input parameters                               units    size
!
!   number of soundings in zonal direction (m)        n/d      1
!   number of atmospheric layers (np)                 n/d      1
!   level pressure (prse)                             mb      m*(np+1)
!   layer temperature (ta)                            k       m*np
!   layer specific humidity (wa)                      g/g     m*np
!   layer ozone mixing ratio by mass (oa)             g/g     m*np
!   surface temperature (ts)                          k        m  
!   co2 mixing ratio by volumn (co2)                  pppv     1
!   n2o mixing ratio by volumn (n2o)                  pppv     1
!   ch4 mixing ratio by volumn (ch4)                  pppv     1
!   cfc11 mixing ratio by volumn (cfc11)              pppv     1
!   cfc12 mixing ratio by volumn (cfc12)              pppv     1
!   cfc22 mixing ratio by volumn (cfc22)              pppv     1
!   cloud optical thickness (taucl)                   n/d     m*np
!   cloud amount (fcld)                             fraction  m*np
!   level index separating high and middle            n/d      1
!           clouds (ict)
!   level index separating middle and low             n/d      1
!           clouds (icb)
!   aerosol optical thickness (taual)                 n/d     m*np*10
!   high (see explanation above)                               1
!   trace (see explanation above)                              1
!
! pre-computed tables used in table look-up for transmittance calculations:
!
!   c1 , c2, c3: for co2 (band 3)
!   o1 , o2, o3: for  o3 (band 5)
!   h11,h12,h13: for h2o (band 1)
!   h21,h22,h23: for h2o (band 2)
!   h81,h82,h83: for h2o (band 8)
!
!---- output parameters
!
!   net downward flux, all-sky   (flx)             w/m**2     m*(np+1)
!   net downward flux, clear-sky (flc)             w/m**2     m*(np+1)
!   downward flux, all-sky   (flxd)                w/m**2     m*(np+1)
!   upward flux, all-sky   (flxu)                  w/m**2     m*(np+1)
!   downward flux, clear-sky   (flxdc)             w/m**2     m*(np+1)
!   upward flux, clear-sky   (flxuc)               w/m**2     m*(np+1)
!   sensitivity of net downward flux  
!       to surface temperature (dfdts)             w/m**2/k   m*(np+1)
!   emission by the surface (st4)                  w/m**2     m 
! 
! notes: 
!
!   (1)  water vapor continuum absorption is included in 540-1380 /cm.
!   (2)  scattering by clouds and aerosols is not included.
!   (3)  clouds are assumed "gray" bodies.
!   (4)  the diffuse cloud and aerosol transmissions are computed 
!        to be exp(-1.66*tau).
!   (5)  if there are no clouds, flx=flc.
!   (6)  plevel(1) is the pressure at the top of the model atmosphere, and
!        plevel(np+1) is the surface pressure.
!   (7)  downward flux is positive, and upward flux is negative.
!   (8)  dfdts is always negative because upward flux is defined as negative.
!   (9)  for questions and coding errors, please contact with ming-dah chou,           &
!        code 913, nasa/goddard space flight center, greenbelt, md 20771.
!        phone: 301-286-4012, fax: 301-286-1759,           &
!        e-mail: chou@climate.gsfc.nasa.gov
!
!-----parameters defining the size of the pre-computed tables for transmittance
!     calculations using table look-up.
!
!     "nx" is the number of intervals in pressure
!     "no" is the number of intervals in o3 amount
!     "nc" is the number of intervals in co2 amount
!     "nh" is the number of intervals in h2o amount
!     "nt" is the number of copies to be made from the pre-computed
!          transmittance tables to reduce "memory-bank conflict"
!          in parallel machines and, hence, enhancing the speed of
!          computations using table look-up. 
!          if such advantage does not exist, "nt" can be set to 1.
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer,parameter    ::  nx=26,no=21,nc=24,nh=31,nt=1
   integer              ::  ipts
   integer              ::  i,k,ip,iw,it,ib,ik,iq,isb,k1,k2
!
!---- input parameters ------
!
   integer              ::  m,np,ict,icb,jlat
   real                 ::  taucl(m,np),fcld(m,np),prse(m,np+1)
   real                 ::  ta(m,np),wa(m,np),oa(m,np),ts(m)
   real                 ::  taual(m,np,10)
   real                 ::  co2,n2o,ch4,cfc11,cfc12,cfc22,tofday

   logical              ::  high,trace
!
!---- output parameters ------
!
   real                 ::  flx(m,np+1),flc(m,np+1),dfdts(m,np+1)
   real                 ::  st4(m),flcd(m,np+1),flcu(m,np+1)
   real                 ::  flxu(m,np+1),flxd(m,np+1)
   real                 ::  flxuc(m,np+1),flxdc(m,np+1)
!
!---- static data -----
!
   real                 ::  cb(5,10)
   real                 ::  xkw(9),aw(9),bw(9),pm(9),fkw(6,9),gkw(6,3),xke(9)
   integer              ::  ne(9),mw(9)
!
!---- temporary arrays -----
!
   real                 ::  pa(m,np),dt(m,np)
   real                 ::  sh2o(m,np+1),swpre(m,np+1),swtem(m,np+1)
   real                 ::  sco3(m,np+1),scopre(m,np+1),scotem(m,np+1)
   real                 ::  dh2o(m,np),dcont(m,np),dco2(m,np),do3(m,np)
   real                 ::  dn2o(m,np),dch4(m,np)
   real                 ::  df11(m,np),df12(m,np),df22(m,np)
   real                 ::  th2o(m,6),tcon(m,3),tco2(m,6,2)
   real                 ::  tn2o(m,4),tch4(m,4),tcom(m,2)
   real                 ::  tf11(m),tf12(m),tf22(m)
   real                 ::  h2oexp(m,np,6), conexp(m,np,3), co2exp(m,np,6,2)
   real                 ::  n2oexp(m,np,4), ch4exp(m,np,4), comexp(m,np,2)
   real                 ::  f11exp(m,np), f12exp(m,np), f22exp(m,np)
   real                 ::  clr(m,0:np+1),fclr(m)
   real                 ::  blayer(m,0:np+1),dbs(m)
   real                 ::  clrlw(m),clrmd(m),clrhi(m)
   real                 ::  trant(m),tranal(m)
   real                 ::  rflx(m,np+1),rflc(m,np+1)
!
   logical              ::  oznbnd,co2bnd,h2otbl,conbnd,n2obnd
   logical              ::  ch4bnd,combnd,f11bnd,f12bnd,f22bnd,b10bnd
!
   real                 ::  c1 (nx,nc,nt),c2 (nx,nc,nt),c3 (nx,nc,nt)
   real                 ::  o1 (nx,no,nt),o2 (nx,no,nt),o3 (nx,no,nt)
   real                 ::  h11(nx,nh,nt),h12(nx,nh,nt),h13(nx,nh,nt)
   real                 ::  h21(nx,nh,nt),h22(nx,nh,nt),h23(nx,nh,nt)
   real                 ::  h81(nx,nh,nt),h82(nx,nh,nt),h83(nx,nh,nt)
   real                 ::  cc1 (nx,nc,nt),oo1 (nx,no,nt)
   real                 ::  hh11(nx,nh,nt),hh21(nx,nh,nt),hh81(nx,nh,nt)
!
   real                 ::  dp,xx,trantc, w1,p1,dwe,dpe
   real                 ::  a1,b1,fk1,a2,b2,fk2
   real                 ::  tofdayrd
   integer              ::  jlatrd,ibrd,k1rd,k2rd
!
!-----the following coefficients (equivalent to table 2 of
!     chou and suarez, 1995) are for computing spectrally 
!     integtrated planck fluxes using eq. (22)
!
    data cb/                                                                   &
         -2.6844e-1,-8.8994e-2, 1.5676e-3,-2.9349e-6, 2.2233e-9,               &
          3.7315e+1,-7.4758e-1, 4.6151e-3,-6.3260e-6, 3.5647e-9,               &
          3.7187e+1,-3.9085e-1,-6.1072e-4, 1.4534e-5,-1.6863e-8,               &
         -4.1928e+1, 1.0027e+0,-8.5789e-3, 2.9199e-5,-2.5654e-8,               &
         -4.9163e+1, 9.8457e-1,-7.0968e-3, 2.0478e-5,-1.5514e-8,               &
         -4.7107e+1, 8.9130e-1,-5.9735e-3, 1.5596e-5,-9.5911e-9,               &
         -5.4041e+1, 9.5332e-1,-5.7784e-3, 1.2555e-5,-2.9377e-9,               &
         -6.9233e+0,-1.5878e-1, 3.9160e-3,-2.4496e-5, 4.9301e-8,               &
          1.1483e+2,-2.2376e+0, 1.6394e-2,-5.3672e-5, 6.6456e-8,               &
          1.9668e+1,-3.1161e-1, 1.2886e-3, 3.6835e-7,-1.6212e-9/
   save cb
!
!-----xkw  are the absorption coefficients for the first
!     k-distribution function due to water vapor line absorption
!     (tables 4 and 7).  units are cm**2/g    
!
   data xkw / 29.55  , 4.167e-1, 1.328e-2, 5.250e-4,                           &
              5.25e-4, 8.938e-3, 4.719e-2, 1.320e-0, 5.250e-4/
   save xkw
!
!-----xke are the absorption coefficients for the first
!     k-distribution function due to water vapor continuum absorption
!     (table 6).  units are cm**2/g
!
   data xke /  0.00,   0.00,   27.40,   15.8,                                  &
               9.40,   7.75,     0.0,    0.0,    0.0/
   save xke
!
!-----mw are the ratios between neighboring absorption coefficients
!     for water vapor line absorption (tables 4 and 7).
!
   data mw /6,6,8,6,6,8,9,6,16/
   save mw
!
!-----aw and bw (table 3) are the coefficients for temperature scaling
!     in eq. (25).
!
   data aw/ 0.0021, 0.0140, 0.0167, 0.0302,                                    &
            0.0307, 0.0195, 0.0152, 0.0008, 0.0096/
   data bw/ -1.01e-5, 5.57e-5, 8.54e-5, 2.96e-4,                               &
             2.86e-4, 1.108e-4, 7.608e-5, -3.52e-6, 1.64e-5/
!
   data pm/ 1.0, 1.0, 1.0, 1.0, 1.0, 0.77, 0.5, 1.0, 1.0/
   save aw, bw, pm
!
!-----fkw is the planck-weighted k-distribution function due to h2o
!     line absorption given in table 4 of chou and suarez (1995).
!     the k-distribution function for the third band, fkw(*,3), is not used
!
   data fkw / 0.2747,0.2717,0.2752,0.1177,0.0352,0.0255,                       &
              0.1521,0.3974,0.1778,0.1826,0.0374,0.0527,                       &
              6*1.00,                                                          &
              0.4654,0.2991,0.1343,0.0646,0.0226,0.0140,                       &
              0.5543,0.2723,0.1131,0.0443,0.0160,0.0000,                       &
              0.5955,0.2693,0.0953,0.0335,0.0064,0.0000,                       &
              0.1958,0.3469,0.3147,0.1013,0.0365,0.0048,                       &
              0.0740,0.1636,0.4174,0.1783,0.1101,0.0566,                       &
              0.1437,0.2197,0.3185,0.2351,0.0647,0.0183/
   save fkw
!
!-----gkw is the planck-weighted k-distribution function due to h2o
!     line absorption in the 3 subbands (800-720,620-720,540-620 /cm)
!     of band 3 given in table 7.  note that the order of the sub-bands
!     is reversed.
!
   data gkw/  0.1782,0.0593,0.0215,0.0068,0.0022,0.0000,                       &
              0.0923,0.1675,0.0923,0.0187,0.0178,0.0000,                       &
              0.0000,0.1083,0.1581,0.0455,0.0274,0.0041/
   save gkw
!
!-----ne is the number of terms used in each band to compute water vapor
!     continuum transmittance (table 6).
!
   data ne /0,0,3,1,1,1,0,0,0/
   save ne
!
!-----include tables used in the table look-up for co2 (band 3), 
!     o3 (band 5), and h2o (bands 1, 2, and 7) transmission functions.
!
!c    include "lookup_tran"
!------data from lookup_tran
!
   data ((h11(ip,iw,1),iw=1,31),ip= 1, 1)/                                     &
     0.6160e-04,  0.9815e-04,  0.1474e-03,  0.2092e-03,  0.2823e-03,           &
     0.3662e-03,  0.4615e-03,  0.5710e-03,  0.6998e-03,  0.8554e-03,           &
     0.1049e-02,  0.1295e-02,  0.1612e-02,  0.2021e-02,  0.2548e-02,           &
     0.3230e-02,  0.4123e-02,  0.5306e-02,  0.6887e-02,  0.9021e-02,           &
     0.1193e-01,  0.1590e-01,  0.2135e-01,  0.2885e-01,  0.3914e-01,           &
     0.5317e-01,  0.7223e-01,  0.9800e-01,  0.1326e+00,  0.1783e+00,           &
     0.2373e+00/
   data ((h12(ip,iw,1),iw=1,31),ip= 1, 1)/                                     &
    -0.2021e-06, -0.3628e-06, -0.5891e-06, -0.8735e-06, -0.1204e-05,           &
    -0.1579e-05, -0.2002e-05, -0.2494e-05, -0.3093e-05, -0.3852e-05,           &
    -0.4835e-05, -0.6082e-05, -0.7591e-05, -0.9332e-05, -0.1128e-04,           &
    -0.1347e-04, -0.1596e-04, -0.1890e-04, -0.2241e-04, -0.2672e-04,           &
    -0.3208e-04, -0.3884e-04, -0.4747e-04, -0.5854e-04, -0.7272e-04,           &
    -0.9092e-04, -0.1146e-03, -0.1458e-03, -0.1877e-03, -0.2435e-03,           &
    -0.3159e-03/
   data ((h13(ip,iw,1),iw=1,31),ip= 1, 1)/                                     &
     0.5907e-09,  0.8541e-09,  0.1095e-08,  0.1272e-08,  0.1297e-08,           &
     0.1105e-08,  0.6788e-09, -0.5585e-10, -0.1147e-08, -0.2746e-08,           &
    -0.5001e-08, -0.7715e-08, -0.1037e-07, -0.1227e-07, -0.1287e-07,           &
    -0.1175e-07, -0.8517e-08, -0.2920e-08,  0.4786e-08,  0.1407e-07,           &
     0.2476e-07,  0.3781e-07,  0.5633e-07,  0.8578e-07,  0.1322e-06,           &
     0.2013e-06,  0.3006e-06,  0.4409e-06,  0.6343e-06,  0.8896e-06,           &
     0.1216e-05/
   data ((h11(ip,iw,1),iw=1,31),ip= 2, 2)/                                     &
     0.6166e-04,  0.9828e-04,  0.1477e-03,  0.2097e-03,  0.2833e-03,           &
     0.3680e-03,  0.4648e-03,  0.5768e-03,  0.7101e-03,  0.8736e-03,           &
     0.1080e-02,  0.1348e-02,  0.1700e-02,  0.2162e-02,  0.2767e-02,           &
     0.3563e-02,  0.4621e-02,  0.6039e-02,  0.7953e-02,  0.1056e-01,           &
     0.1412e-01,  0.1901e-01,  0.2574e-01,  0.3498e-01,  0.4763e-01,           &
     0.6484e-01,  0.8815e-01,  0.1196e+00,  0.1614e+00,  0.2157e+00,           &
     0.2844e+00/
   data ((h12(ip,iw,1),iw=1,31),ip= 2, 2)/                                     &
    -0.2017e-06, -0.3620e-06, -0.5878e-06, -0.8713e-06, -0.1201e-05,           &
    -0.1572e-05, -0.1991e-05, -0.2476e-05, -0.3063e-05, -0.3808e-05,           &
    -0.4776e-05, -0.6011e-05, -0.7516e-05, -0.9272e-05, -0.1127e-04,           &
    -0.1355e-04, -0.1620e-04, -0.1936e-04, -0.2321e-04, -0.2797e-04,           &
    -0.3399e-04, -0.4171e-04, -0.5172e-04, -0.6471e-04, -0.8150e-04,           &
    -0.1034e-03, -0.1321e-03, -0.1705e-03, -0.2217e-03, -0.2889e-03,           &
    -0.3726e-03/
   data ((h13(ip,iw,1),iw=1,31),ip= 2, 2)/                                     &
     0.5894e-09,  0.8519e-09,  0.1092e-08,  0.1267e-08,  0.1289e-08,           &
     0.1093e-08,  0.6601e-09, -0.7831e-10, -0.1167e-08, -0.2732e-08,           &
    -0.4864e-08, -0.7334e-08, -0.9581e-08, -0.1097e-07, -0.1094e-07,           &
    -0.8999e-08, -0.4669e-08,  0.2391e-08,  0.1215e-07,  0.2424e-07,           &
     0.3877e-07,  0.5711e-07,  0.8295e-07,  0.1218e-06,  0.1793e-06,           &
     0.2621e-06,  0.3812e-06,  0.5508e-06,  0.7824e-06,  0.1085e-05,           &
     0.1462e-05/
   data ((h11(ip,iw,1),iw=1,31),ip= 3, 3)/                                     &
     0.6175e-04,  0.9849e-04,  0.1481e-03,  0.2106e-03,  0.2849e-03,           &
     0.3708e-03,  0.4698e-03,  0.5857e-03,  0.7257e-03,  0.9006e-03,           &
     0.1126e-02,  0.1425e-02,  0.1823e-02,  0.2353e-02,  0.3059e-02,           &
     0.4002e-02,  0.5270e-02,  0.6984e-02,  0.9316e-02,  0.1251e-01,           &
     0.1689e-01,  0.2292e-01,  0.3123e-01,  0.4262e-01,  0.5814e-01,           &
     0.7921e-01,  0.1077e+00,  0.1458e+00,  0.1957e+00,  0.2595e+00,           &
     0.3380e+00/
   data ((h12(ip,iw,1),iw=1,31),ip= 3, 3)/                                     &
    -0.2011e-06, -0.3609e-06, -0.5859e-06, -0.8680e-06, -0.1195e-05,           &
    -0.1563e-05, -0.1975e-05, -0.2450e-05, -0.3024e-05, -0.3755e-05,           &
    -0.4711e-05, -0.5941e-05, -0.7455e-05, -0.9248e-05, -0.1132e-04,           &
    -0.1373e-04, -0.1659e-04, -0.2004e-04, -0.2431e-04, -0.2966e-04,           &
    -0.3653e-04, -0.4549e-04, -0.5724e-04, -0.7259e-04, -0.9265e-04,           &
    -0.1191e-03, -0.1543e-03, -0.2013e-03, -0.2633e-03, -0.3421e-03,           &
    -0.4350e-03/
   data ((h13(ip,iw,1),iw=1,31),ip= 3, 3)/                                     &
     0.5872e-09,  0.8484e-09,  0.1087e-08,  0.1259e-08,  0.1279e-08,           &
     0.1077e-08,  0.6413e-09, -0.9334e-10, -0.1161e-08, -0.2644e-08,           &
    -0.4588e-08, -0.6709e-08, -0.8474e-08, -0.9263e-08, -0.8489e-08,           &
    -0.5553e-08,  0.1203e-09,  0.9035e-08,  0.2135e-07,  0.3689e-07,           &
     0.5610e-07,  0.8097e-07,  0.1155e-06,  0.1649e-06,  0.2350e-06,           &
     0.3353e-06,  0.4806e-06,  0.6858e-06,  0.9617e-06,  0.1315e-05,           &
     0.1741e-05/
   data ((h11(ip,iw,1),iw=1,31),ip= 4, 4)/                                     &
     0.6189e-04,  0.9882e-04,  0.1488e-03,  0.2119e-03,  0.2873e-03,           &
     0.3752e-03,  0.4776e-03,  0.5993e-03,  0.7490e-03,  0.9403e-03,           &
     0.1192e-02,  0.1531e-02,  0.1990e-02,  0.2610e-02,  0.3446e-02,           &
     0.4576e-02,  0.6109e-02,  0.8196e-02,  0.1105e-01,  0.1498e-01,           &
     0.2039e-01,  0.2785e-01,  0.3809e-01,  0.5209e-01,  0.7112e-01,           &
     0.9688e-01,  0.1315e+00,  0.1773e+00,  0.2363e+00,  0.3100e+00,           &
     0.3976e+00/
   data ((h12(ip,iw,1),iw=1,31),ip= 4, 4)/                                     &
    -0.2001e-06, -0.3592e-06, -0.5829e-06, -0.8631e-06, -0.1187e-05,           &
    -0.1549e-05, -0.1953e-05, -0.2415e-05, -0.2975e-05, -0.3694e-05,           &
    -0.4645e-05, -0.5882e-05, -0.7425e-05, -0.9279e-05, -0.1147e-04,           &
    -0.1406e-04, -0.1717e-04, -0.2100e-04, -0.2580e-04, -0.3191e-04,           &
    -0.3989e-04, -0.5042e-04, -0.6432e-04, -0.8261e-04, -0.1068e-03,           &
    -0.1389e-03, -0.1820e-03, -0.2391e-03, -0.3127e-03, -0.4021e-03,           &
    -0.5002e-03/
   data ((h13(ip,iw,1),iw=1,31),ip= 4, 4)/                                     &
     0.5838e-09,  0.8426e-09,  0.1081e-08,  0.1249e-08,  0.1267e-08,           &
     0.1062e-08,  0.6313e-09, -0.8241e-10, -0.1094e-08, -0.2436e-08,           &
    -0.4100e-08, -0.5786e-08, -0.6992e-08, -0.7083e-08, -0.5405e-08,           &
    -0.1259e-08,  0.6099e-08,  0.1732e-07,  0.3276e-07,  0.5256e-07,           &
     0.7756e-07,  0.1103e-06,  0.1547e-06,  0.2159e-06,  0.3016e-06,           &
     0.4251e-06,  0.6033e-06,  0.8499e-06,  0.1175e-05,  0.1579e-05,           &
     0.2044e-05/
   data ((h11(ip,iw,1),iw=1,31),ip= 5, 5)/                                     &
     0.6211e-04,  0.9932e-04,  0.1499e-03,  0.2140e-03,  0.2911e-03,           &
     0.3820e-03,  0.4895e-03,  0.6196e-03,  0.7834e-03,  0.9973e-03,           &
     0.1285e-02,  0.1678e-02,  0.2216e-02,  0.2951e-02,  0.3953e-02,           &
     0.5319e-02,  0.7186e-02,  0.9743e-02,  0.1326e-01,  0.1811e-01,           &
     0.2479e-01,  0.3400e-01,  0.4662e-01,  0.6379e-01,  0.8708e-01,           &
     0.1185e+00,  0.1603e+00,  0.2147e+00,  0.2835e+00,  0.3667e+00,           &
     0.4620e+00/
   data ((h12(ip,iw,1),iw=1,31),ip= 5, 5)/                                     &
    -0.1987e-06, -0.3565e-06, -0.5784e-06, -0.8557e-06, -0.1175e-05,           &
    -0.1530e-05, -0.1923e-05, -0.2372e-05, -0.2919e-05, -0.3631e-05,           &
    -0.4587e-05, -0.5848e-05, -0.7442e-05, -0.9391e-05, -0.1173e-04,           &
    -0.1455e-04, -0.1801e-04, -0.2232e-04, -0.2779e-04, -0.3489e-04,           &
    -0.4428e-04, -0.5678e-04, -0.7333e-04, -0.9530e-04, -0.1246e-03,           &
    -0.1639e-03, -0.2164e-03, -0.2848e-03, -0.3697e-03, -0.4665e-03,           &
    -0.5646e-03/
   data ((h13(ip,iw,1),iw=1,31),ip= 5, 5)/                                     &
     0.5785e-09,  0.8338e-09,  0.1071e-08,  0.1239e-08,  0.1256e-08,           &
     0.1057e-08,  0.6480e-09, -0.1793e-10, -0.9278e-09, -0.2051e-08,           &
    -0.3337e-08, -0.4514e-08, -0.5067e-08, -0.4328e-08, -0.1545e-08,           &
     0.4100e-08,  0.1354e-07,  0.2762e-07,  0.4690e-07,  0.7190e-07,           &
     0.1040e-06,  0.1459e-06,  0.2014e-06,  0.2764e-06,  0.3824e-06,           &
     0.5359e-06,  0.7532e-06,  0.1047e-05,  0.1424e-05,  0.1873e-05,           &
     0.2356e-05/
   data ((h11(ip,iw,1),iw=1,31),ip= 6, 6)/                                     &
     0.6246e-04,  0.1001e-03,  0.1515e-03,  0.2171e-03,  0.2970e-03,           &
     0.3924e-03,  0.5072e-03,  0.6495e-03,  0.8329e-03,  0.1078e-02,           &
     0.1413e-02,  0.1876e-02,  0.2516e-02,  0.3399e-02,  0.4612e-02,           &
     0.6276e-02,  0.8562e-02,  0.1171e-01,  0.1605e-01,  0.2205e-01,           &
     0.3032e-01,  0.4167e-01,  0.5717e-01,  0.7821e-01,  0.1067e+00,           &
     0.1447e+00,  0.1948e+00,  0.2586e+00,  0.3372e+00,  0.4290e+00,           &
     0.5295e+00/
   data ((h12(ip,iw,1),iw=1,31),ip= 6, 6)/                                     &
    -0.1964e-06, -0.3526e-06, -0.5717e-06, -0.8451e-06, -0.1158e-05,           &
    -0.1504e-05, -0.1886e-05, -0.2322e-05, -0.2861e-05, -0.3576e-05,           &
    -0.4552e-05, -0.5856e-05, -0.7529e-05, -0.9609e-05, -0.1216e-04,           &
    -0.1528e-04, -0.1916e-04, -0.2408e-04, -0.3043e-04, -0.3880e-04,           &
    -0.4997e-04, -0.6488e-04, -0.8474e-04, -0.1113e-03, -0.1471e-03,           &
    -0.1950e-03, -0.2583e-03, -0.3384e-03, -0.4326e-03, -0.5319e-03,           &
    -0.6244e-03/
   data ((h13(ip,iw,1),iw=1,31),ip= 6, 6)/                                     &
     0.5713e-09,  0.8263e-09,  0.1060e-08,  0.1226e-08,  0.1252e-08,           &
     0.1076e-08,  0.7149e-09,  0.1379e-09, -0.6043e-09, -0.1417e-08,           &
    -0.2241e-08, -0.2830e-08, -0.2627e-08, -0.8950e-09,  0.3231e-08,           &
     0.1075e-07,  0.2278e-07,  0.4037e-07,  0.6439e-07,  0.9576e-07,           &
     0.1363e-06,  0.1886e-06,  0.2567e-06,  0.3494e-06,  0.4821e-06,           &
     0.6719e-06,  0.9343e-06,  0.1280e-05,  0.1705e-05,  0.2184e-05,           &
     0.2651e-05/
   data ((h11(ip,iw,1),iw=1,31),ip= 7, 7)/                                     &
     0.6299e-04,  0.1013e-03,  0.1541e-03,  0.2220e-03,  0.3058e-03,           &
     0.4078e-03,  0.5334e-03,  0.6928e-03,  0.9032e-03,  0.1190e-02,           &
     0.1587e-02,  0.2140e-02,  0.2912e-02,  0.3982e-02,  0.5460e-02,           &
     0.7500e-02,  0.1031e-01,  0.1420e-01,  0.1958e-01,  0.2700e-01,           &
     0.3721e-01,  0.5118e-01,  0.7019e-01,  0.9593e-01,  0.1305e+00,           &
     0.1763e+00,  0.2354e+00,  0.3091e+00,  0.3969e+00,  0.4952e+00,           &
     0.5978e+00/
   data ((h12(ip,iw,1),iw=1,31),ip= 7, 7)/                                     &
    -0.1932e-06, -0.3467e-06, -0.5623e-06, -0.8306e-06, -0.1136e-05,           &
    -0.1472e-05, -0.1842e-05, -0.2269e-05, -0.2807e-05, -0.3539e-05,           &
    -0.4553e-05, -0.5925e-05, -0.7710e-05, -0.9968e-05, -0.1278e-04,           &
    -0.1629e-04, -0.2073e-04, -0.2644e-04, -0.3392e-04, -0.4390e-04,           &
    -0.5727e-04, -0.7516e-04, -0.9916e-04, -0.1315e-03, -0.1752e-03,           &
    -0.2333e-03, -0.3082e-03, -0.3988e-03, -0.4982e-03, -0.5947e-03,           &
    -0.6764e-03/
   data ((h13(ip,iw,1),iw=1,31),ip= 7, 7)/                                     &
     0.5612e-09,  0.8116e-09,  0.1048e-08,  0.1222e-08,  0.1270e-08,           &
     0.1141e-08,  0.8732e-09,  0.4336e-09, -0.6548e-10, -0.4774e-09,           &
    -0.7556e-09, -0.6577e-09,  0.4377e-09,  0.3359e-08,  0.9159e-08,           &
     0.1901e-07,  0.3422e-07,  0.5616e-07,  0.8598e-07,  0.1251e-06,           &
     0.1752e-06,  0.2392e-06,  0.3228e-06,  0.4389e-06,  0.6049e-06,           &
     0.8370e-06,  0.1150e-05,  0.1547e-05,  0.2012e-05,  0.2493e-05,           &
     0.2913e-05/
   data ((h11(ip,iw,1),iw=1,31),ip= 8, 8)/                                     &
     0.6378e-04,  0.1032e-03,  0.1579e-03,  0.2293e-03,  0.3190e-03,           &
     0.4305e-03,  0.5712e-03,  0.7541e-03,  0.1001e-02,  0.1342e-02,           &
     0.1819e-02,  0.2489e-02,  0.3427e-02,  0.4734e-02,  0.6548e-02,           &
     0.9059e-02,  0.1254e-01,  0.1735e-01,  0.2401e-01,  0.3318e-01,           &
     0.4577e-01,  0.6292e-01,  0.8620e-01,  0.1176e+00,  0.1594e+00,           &
     0.2139e+00,  0.2827e+00,  0.3660e+00,  0.4614e+00,  0.5634e+00,           &
     0.6649e+00/
   data ((h12(ip,iw,1),iw=1,31),ip= 8, 8)/                                     &
    -0.1885e-06, -0.3385e-06, -0.5493e-06, -0.8114e-06, -0.1109e-05,           &
    -0.1436e-05, -0.1796e-05, -0.2219e-05, -0.2770e-05, -0.3535e-05,           &
    -0.4609e-05, -0.6077e-05, -0.8016e-05, -0.1051e-04, -0.1367e-04,           &
    -0.1768e-04, -0.2283e-04, -0.2955e-04, -0.3849e-04, -0.5046e-04,           &
    -0.6653e-04, -0.8813e-04, -0.1173e-03, -0.1569e-03, -0.2100e-03,           &
    -0.2794e-03, -0.3656e-03, -0.4637e-03, -0.5629e-03, -0.6512e-03,           &
    -0.7167e-03/
   data ((h13(ip,iw,1),iw=1,31),ip= 8, 8)/                                     &
     0.5477e-09,  0.8000e-09,  0.1039e-08,  0.1234e-08,  0.1331e-08,           &
     0.1295e-08,  0.1160e-08,  0.9178e-09,  0.7535e-09,  0.8301e-09,           &
     0.1184e-08,  0.2082e-08,  0.4253e-08,  0.8646e-08,  0.1650e-07,           &
     0.2920e-07,  0.4834e-07,  0.7564e-07,  0.1125e-06,  0.1606e-06,           &
     0.2216e-06,  0.2992e-06,  0.4031e-06,  0.5493e-06,  0.7549e-06,           &
     0.1035e-05,  0.1400e-05,  0.1843e-05,  0.2327e-05,  0.2774e-05,           &
     0.3143e-05/
   data ((h11(ip,iw,1),iw=1,31),ip= 9, 9)/                                     &
     0.6495e-04,  0.1059e-03,  0.1635e-03,  0.2400e-03,  0.3381e-03,           &
     0.4631e-03,  0.6246e-03,  0.8394e-03,  0.1134e-02,  0.1546e-02,           &
     0.2126e-02,  0.2944e-02,  0.4093e-02,  0.5699e-02,  0.7934e-02,           &
     0.1104e-01,  0.1535e-01,  0.2131e-01,  0.2956e-01,  0.4089e-01,           &
     0.5636e-01,  0.7739e-01,  0.1058e+00,  0.1439e+00,  0.1939e+00,           &
     0.2578e+00,  0.3364e+00,  0.4283e+00,  0.5290e+00,  0.6314e+00,           &
     0.7292e+00/
   data ((h12(ip,iw,1),iw=1,31),ip= 9, 9)/                                     &
    -0.1822e-06, -0.3274e-06, -0.5325e-06, -0.7881e-06, -0.1079e-05,           &
    -0.1398e-05, -0.1754e-05, -0.2184e-05, -0.2763e-05, -0.3581e-05,           &
    -0.4739e-05, -0.6341e-05, -0.8484e-05, -0.1128e-04, -0.1490e-04,           &
    -0.1955e-04, -0.2561e-04, -0.3364e-04, -0.4438e-04, -0.5881e-04,           &
    -0.7822e-04, -0.1045e-03, -0.1401e-03, -0.1884e-03, -0.2523e-03,           &
    -0.3335e-03, -0.4289e-03, -0.5296e-03, -0.6231e-03, -0.6980e-03,           &
    -0.7406e-03/
   data ((h13(ip,iw,1),iw=1,31),ip= 9, 9)/                                     &
     0.5334e-09,  0.7859e-09,  0.1043e-08,  0.1279e-08,  0.1460e-08,           &
     0.1560e-08,  0.1618e-08,  0.1657e-08,  0.1912e-08,  0.2569e-08,           &
     0.3654e-08,  0.5509e-08,  0.8964e-08,  0.1518e-07,  0.2560e-07,           &
     0.4178e-07,  0.6574e-07,  0.9958e-07,  0.1449e-06,  0.2031e-06,           &
     0.2766e-06,  0.3718e-06,  0.5022e-06,  0.6849e-06,  0.9360e-06,           &
     0.1268e-05,  0.1683e-05,  0.2157e-05,  0.2625e-05,  0.3020e-05,           &
     0.3364e-05/
   data ((h11(ip,iw,1),iw=1,31),ip=10,10)/                                     &
     0.6661e-04,  0.1098e-03,  0.1716e-03,  0.2554e-03,  0.3653e-03,           &
     0.5088e-03,  0.6986e-03,  0.9557e-03,  0.1313e-02,  0.1816e-02,           &
     0.2527e-02,  0.3532e-02,  0.4947e-02,  0.6929e-02,  0.9694e-02,           &
     0.1354e-01,  0.1888e-01,  0.2628e-01,  0.3647e-01,  0.5043e-01,           &
     0.6941e-01,  0.9514e-01,  0.1297e+00,  0.1755e+00,  0.2347e+00,           &
     0.3084e+00,  0.3962e+00,  0.4947e+00,  0.5974e+00,  0.6973e+00,           &
     0.7898e+00/
   data ((h12(ip,iw,1),iw=1,31),ip=10,10)/                                     &
    -0.1742e-06, -0.3134e-06, -0.5121e-06, -0.7619e-06, -0.1048e-05,           &
    -0.1364e-05, -0.1725e-05, -0.2177e-05, -0.2801e-05, -0.3694e-05,           &
    -0.4969e-05, -0.6748e-05, -0.9161e-05, -0.1236e-04, -0.1655e-04,           &
    -0.2203e-04, -0.2927e-04, -0.3894e-04, -0.5192e-04, -0.6936e-04,           &
    -0.9294e-04, -0.1250e-03, -0.1686e-03, -0.2271e-03, -0.3027e-03,           &
    -0.3944e-03, -0.4951e-03, -0.5928e-03, -0.6755e-03, -0.7309e-03,           &
    -0.7417e-03/
   data ((h13(ip,iw,1),iw=1,31),ip=10,10)/                                     &
     0.5179e-09,  0.7789e-09,  0.1071e-08,  0.1382e-08,  0.1690e-08,           &
     0.1979e-08,  0.2297e-08,  0.2704e-08,  0.3466e-08,  0.4794e-08,           &
     0.6746e-08,  0.9739e-08,  0.1481e-07,  0.2331e-07,  0.3679e-07,           &
     0.5726e-07,  0.8716e-07,  0.1289e-06,  0.1837e-06,  0.2534e-06,           &
     0.3424e-06,  0.4609e-06,  0.6245e-06,  0.8495e-06,  0.1151e-05,           &
     0.1536e-05,  0.1991e-05,  0.2468e-05,  0.2891e-05,  0.3245e-05,           &
     0.3580e-05/
   data ((h11(ip,iw,1),iw=1,31),ip=11,11)/                                     &
     0.6888e-04,  0.1151e-03,  0.1828e-03,  0.2766e-03,  0.4028e-03,           &
     0.5713e-03,  0.7987e-03,  0.1111e-02,  0.1548e-02,  0.2167e-02,           &
     0.3044e-02,  0.4285e-02,  0.6035e-02,  0.8490e-02,  0.1192e-01,           &
     0.1670e-01,  0.2333e-01,  0.3249e-01,  0.4506e-01,  0.6220e-01,           &
     0.8546e-01,  0.1168e+00,  0.1587e+00,  0.2131e+00,  0.2820e+00,           &
     0.3653e+00,  0.4609e+00,  0.5630e+00,  0.6645e+00,  0.7599e+00,           &
     0.8458e+00/
   data ((h12(ip,iw,1),iw=1,31),ip=11,11)/                                     &
    -0.1647e-06, -0.2974e-06, -0.4900e-06, -0.7358e-06, -0.1022e-05,           &
    -0.1344e-05, -0.1721e-05, -0.2212e-05, -0.2901e-05, -0.3896e-05,           &
    -0.5327e-05, -0.7342e-05, -0.1011e-04, -0.1382e-04, -0.1875e-04,           &
    -0.2530e-04, -0.3403e-04, -0.4573e-04, -0.6145e-04, -0.8264e-04,           &
    -0.1114e-03, -0.1507e-03, -0.2039e-03, -0.2737e-03, -0.3607e-03,           &
    -0.4599e-03, -0.5604e-03, -0.6497e-03, -0.7161e-03, -0.7443e-03,           &
    -0.7133e-03/
   data ((h13(ip,iw,1),iw=1,31),ip=11,11)/                                     &
     0.5073e-09,  0.7906e-09,  0.1134e-08,  0.1560e-08,  0.2046e-08,           &
     0.2589e-08,  0.3254e-08,  0.4107e-08,  0.5481e-08,  0.7602e-08,           &
     0.1059e-07,  0.1501e-07,  0.2210e-07,  0.3334e-07,  0.5055e-07,           &
     0.7629e-07,  0.1134e-06,  0.1642e-06,  0.2298e-06,  0.3133e-06,           &
     0.4225e-06,  0.5709e-06,  0.7739e-06,  0.1047e-05,  0.1401e-05,           &
     0.1833e-05,  0.2308e-05,  0.2753e-05,  0.3125e-05,  0.3467e-05,           &
     0.3748e-05/
   data ((h11(ip,iw,1),iw=1,31),ip=12,12)/                                     &
     0.7179e-04,  0.1221e-03,  0.1975e-03,  0.3049e-03,  0.4529e-03,           &
     0.6547e-03,  0.9312e-03,  0.1315e-02,  0.1855e-02,  0.2620e-02,           &
     0.3705e-02,  0.5243e-02,  0.7414e-02,  0.1046e-01,  0.1472e-01,           &
     0.2065e-01,  0.2888e-01,  0.4019e-01,  0.5566e-01,  0.7668e-01,           &
     0.1051e+00,  0.1432e+00,  0.1932e+00,  0.2571e+00,  0.3358e+00,           &
     0.4278e+00,  0.5285e+00,  0.6310e+00,  0.7289e+00,  0.8184e+00,           &
     0.8954e+00/
   data ((h12(ip,iw,1),iw=1,31),ip=12,12)/                                     &
    -0.1548e-06, -0.2808e-06, -0.4683e-06, -0.7142e-06, -0.1008e-05,           &
    -0.1347e-05, -0.1758e-05, -0.2306e-05, -0.3083e-05, -0.4214e-05,           &
    -0.5851e-05, -0.8175e-05, -0.1140e-04, -0.1577e-04, -0.2166e-04,           &
    -0.2955e-04, -0.4014e-04, -0.5434e-04, -0.7343e-04, -0.9931e-04,           &
    -0.1346e-03, -0.1826e-03, -0.2467e-03, -0.3283e-03, -0.4246e-03,           &
    -0.5264e-03, -0.6211e-03, -0.6970e-03, -0.7402e-03, -0.7316e-03,           &
    -0.6486e-03/
   data ((h13(ip,iw,1),iw=1,31),ip=12,12)/                                     &
     0.5078e-09,  0.8244e-09,  0.1255e-08,  0.1826e-08,  0.2550e-08,           &
     0.3438e-08,  0.4532e-08,  0.5949e-08,  0.8041e-08,  0.1110e-07,           &
     0.1534e-07,  0.2157e-07,  0.3116e-07,  0.4570e-07,  0.6747e-07,           &
     0.9961e-07,  0.1451e-06,  0.2061e-06,  0.2843e-06,  0.3855e-06,           &
     0.5213e-06,  0.7060e-06,  0.9544e-06,  0.1280e-05,  0.1684e-05,           &
     0.2148e-05,  0.2609e-05,  0.3002e-05,  0.3349e-05,  0.3670e-05,           &
     0.3780e-05/
   data ((h11(ip,iw,1),iw=1,31),ip=13,13)/                                     &
     0.7529e-04,  0.1306e-03,  0.2160e-03,  0.3410e-03,  0.5175e-03,           &
     0.7626e-03,  0.1103e-02,  0.1577e-02,  0.2246e-02,  0.3196e-02,           &
     0.4544e-02,  0.6455e-02,  0.9152e-02,  0.1294e-01,  0.1823e-01,           &
     0.2560e-01,  0.3577e-01,  0.4971e-01,  0.6870e-01,  0.9443e-01,           &
     0.1290e+00,  0.1748e+00,  0.2340e+00,  0.3078e+00,  0.3956e+00,           &
     0.4942e+00,  0.5969e+00,  0.6970e+00,  0.7896e+00,  0.8714e+00,           &
     0.9364e+00/
   data ((h12(ip,iw,1),iw=1,31),ip=13,13)/                                     &
    -0.1461e-06, -0.2663e-06, -0.4512e-06, -0.7027e-06, -0.1014e-05,           &
    -0.1387e-05, -0.1851e-05, -0.2478e-05, -0.3373e-05, -0.4682e-05,           &
    -0.6588e-05, -0.9311e-05, -0.1311e-04, -0.1834e-04, -0.2544e-04,           &
    -0.3502e-04, -0.4789e-04, -0.6515e-04, -0.8846e-04, -0.1202e-03,           &
    -0.1635e-03, -0.2217e-03, -0.2975e-03, -0.3897e-03, -0.4913e-03,           &
    -0.5902e-03, -0.6740e-03, -0.7302e-03, -0.7415e-03, -0.6858e-03,           &
    -0.5447e-03/
   data ((h13(ip,iw,1),iw=1,31),ip=13,13)/                                     &
     0.5236e-09,  0.8873e-09,  0.1426e-08,  0.2193e-08,  0.3230e-08,           &
     0.4555e-08,  0.6200e-08,  0.8298e-08,  0.1126e-07,  0.1544e-07,           &
     0.2130e-07,  0.2978e-07,  0.4239e-07,  0.6096e-07,  0.8829e-07,           &
     0.1280e-06,  0.1830e-06,  0.2555e-06,  0.3493e-06,  0.4740e-06,           &
     0.6431e-06,  0.8701e-06,  0.1169e-05,  0.1547e-05,  0.1992e-05,           &
     0.2460e-05,  0.2877e-05,  0.3230e-05,  0.3569e-05,  0.3782e-05,           &
     0.3591e-05/
   data ((h11(ip,iw,1),iw=1,31),ip=14,14)/                                     &
     0.7911e-04,  0.1403e-03,  0.2378e-03,  0.3847e-03,  0.5973e-03,           &
     0.8978e-03,  0.1319e-02,  0.1909e-02,  0.2741e-02,  0.3923e-02,           &
     0.5600e-02,  0.7977e-02,  0.1133e-01,  0.1604e-01,  0.2262e-01,           &
     0.3174e-01,  0.4429e-01,  0.6143e-01,  0.8469e-01,  0.1161e+00,           &
     0.1579e+00,  0.2124e+00,  0.2813e+00,  0.3647e+00,  0.4603e+00,           &
     0.5625e+00,  0.6641e+00,  0.7596e+00,  0.8456e+00,  0.9170e+00,           &
     0.9670e+00/
   data ((h12(ip,iw,1),iw=1,31),ip=14,14)/                                     &
    -0.1402e-06, -0.2569e-06, -0.4428e-06, -0.7076e-06, -0.1051e-05,           &
    -0.1478e-05, -0.2019e-05, -0.2752e-05, -0.3802e-05, -0.5343e-05,           &
    -0.7594e-05, -0.1082e-04, -0.1536e-04, -0.2166e-04, -0.3028e-04,           &
    -0.4195e-04, -0.5761e-04, -0.7867e-04, -0.1072e-03, -0.1462e-03,           &
    -0.1990e-03, -0.2687e-03, -0.3559e-03, -0.4558e-03, -0.5572e-03,           &
    -0.6476e-03, -0.7150e-03, -0.7439e-03, -0.7133e-03, -0.6015e-03,           &
    -0.4089e-03/
   data ((h13(ip,iw,1),iw=1,31),ip=14,14)/                                     &
     0.5531e-09,  0.9757e-09,  0.1644e-08,  0.2650e-08,  0.4074e-08,           &
     0.5957e-08,  0.8314e-08,  0.1128e-07,  0.1528e-07,  0.2087e-07,           &
     0.2874e-07,  0.4002e-07,  0.5631e-07,  0.7981e-07,  0.1139e-06,           &
     0.1621e-06,  0.2275e-06,  0.3136e-06,  0.4280e-06,  0.5829e-06,           &
     0.7917e-06,  0.1067e-05,  0.1419e-05,  0.1844e-05,  0.2310e-05,           &
     0.2747e-05,  0.3113e-05,  0.3455e-05,  0.3739e-05,  0.3715e-05,           &
     0.3125e-05/
   data ((h11(ip,iw,1),iw=1,31),ip=15,15)/                                     &
     0.8292e-04,  0.1503e-03,  0.2614e-03,  0.4345e-03,  0.6915e-03,           &
     0.1061e-02,  0.1584e-02,  0.2320e-02,  0.3359e-02,  0.4831e-02,           &
     0.6920e-02,  0.9877e-02,  0.1405e-01,  0.1990e-01,  0.2805e-01,           &
     0.3933e-01,  0.5477e-01,  0.7579e-01,  0.1042e+00,  0.1423e+00,           &
     0.1924e+00,  0.2564e+00,  0.3351e+00,  0.4271e+00,  0.5280e+00,           &
     0.6306e+00,  0.7286e+00,  0.8182e+00,  0.8952e+00,  0.9530e+00,           &
     0.9864e+00/
   data ((h12(ip,iw,1),iw=1,31),ip=15,15)/                                     &
    -0.1378e-06, -0.2542e-06, -0.4461e-06, -0.7333e-06, -0.1125e-05,           &
    -0.1630e-05, -0.2281e-05, -0.3159e-05, -0.4410e-05, -0.6246e-05,           &
    -0.8933e-05, -0.1280e-04, -0.1826e-04, -0.2589e-04, -0.3639e-04,           &
    -0.5059e-04, -0.6970e-04, -0.9552e-04, -0.1307e-03, -0.1784e-03,           &
    -0.2422e-03, -0.3237e-03, -0.4203e-03, -0.5227e-03, -0.6184e-03,           &
    -0.6953e-03, -0.7395e-03, -0.7315e-03, -0.6487e-03, -0.4799e-03,           &
    -0.2625e-03/
   data ((h13(ip,iw,1),iw=1,31),ip=15,15)/                                     &
     0.5891e-09,  0.1074e-08,  0.1885e-08,  0.3167e-08,  0.5051e-08,           &
     0.7631e-08,  0.1092e-07,  0.1500e-07,  0.2032e-07,  0.2769e-07,           &
     0.3810e-07,  0.5279e-07,  0.7361e-07,  0.1032e-06,  0.1450e-06,           &
     0.2026e-06,  0.2798e-06,  0.3832e-06,  0.5242e-06,  0.7159e-06,           &
     0.9706e-06,  0.1299e-05,  0.1701e-05,  0.2159e-05,  0.2612e-05,           &
     0.2998e-05,  0.3341e-05,  0.3661e-05,  0.3775e-05,  0.3393e-05,           &
     0.2384e-05/
   data ((h11(ip,iw,1),iw=1,31),ip=16,16)/                                     &
     0.8637e-04,  0.1598e-03,  0.2853e-03,  0.4876e-03,  0.7970e-03,           &
     0.1251e-02,  0.1901e-02,  0.2820e-02,  0.4118e-02,  0.5955e-02,           &
     0.8557e-02,  0.1224e-01,  0.1742e-01,  0.2467e-01,  0.3476e-01,           &
     0.4864e-01,  0.6759e-01,  0.9332e-01,  0.1280e+00,  0.1738e+00,           &
     0.2330e+00,  0.3069e+00,  0.3949e+00,  0.4935e+00,  0.5964e+00,           &
     0.6965e+00,  0.7893e+00,  0.8713e+00,  0.9363e+00,  0.9780e+00,           &
     0.9961e+00/
   data ((h12(ip,iw,1),iw=1,31),ip=16,16)/                                     &
    -0.1383e-06, -0.2577e-06, -0.4608e-06, -0.7793e-06, -0.1237e-05,           &
    -0.1850e-05, -0.2652e-05, -0.3728e-05, -0.5244e-05, -0.7451e-05,           &
    -0.1067e-04, -0.1532e-04, -0.2193e-04, -0.3119e-04, -0.4395e-04,           &
    -0.6126e-04, -0.8466e-04, -0.1164e-03, -0.1596e-03, -0.2177e-03,           &
    -0.2933e-03, -0.3855e-03, -0.4874e-03, -0.5870e-03, -0.6718e-03,           &
    -0.7290e-03, -0.7411e-03, -0.6859e-03, -0.5450e-03, -0.3353e-03,           &
    -0.1363e-03/
   data ((h13(ip,iw,1),iw=1,31),ip=16,16)/                                     &
     0.6217e-09,  0.1165e-08,  0.2116e-08,  0.3685e-08,  0.6101e-08,           &
     0.9523e-08,  0.1400e-07,  0.1959e-07,  0.2668e-07,  0.3629e-07,           &
     0.4982e-07,  0.6876e-07,  0.9523e-07,  0.1321e-06,  0.1825e-06,           &
     0.2505e-06,  0.3420e-06,  0.4677e-06,  0.6416e-06,  0.8760e-06,           &
     0.1183e-05,  0.1565e-05,  0.2010e-05,  0.2472e-05,  0.2882e-05,           &
     0.3229e-05,  0.3564e-05,  0.3777e-05,  0.3589e-05,  0.2786e-05,           &
     0.1487e-05/
   data ((h11(ip,iw,1),iw=1,31),ip=17,17)/                                     &
     0.8925e-04,  0.1682e-03,  0.3074e-03,  0.5403e-03,  0.9085e-03,           &
     0.1463e-02,  0.2268e-02,  0.3416e-02,  0.5040e-02,  0.7333e-02,           &
     0.1057e-01,  0.1515e-01,  0.2157e-01,  0.3055e-01,  0.4297e-01,           &
     0.6001e-01,  0.8323e-01,  0.1146e+00,  0.1565e+00,  0.2111e+00,           &
     0.2801e+00,  0.3636e+00,  0.4594e+00,  0.5618e+00,  0.6636e+00,           &
     0.7592e+00,  0.8454e+00,  0.9169e+00,  0.9669e+00,  0.9923e+00,           &
     0.9995e+00/
   data ((h12(ip,iw,1),iw=1,31),ip=17,17)/                                     &
    -0.1405e-06, -0.2649e-06, -0.4829e-06, -0.8398e-06, -0.1379e-05,           &
    -0.2132e-05, -0.3138e-05, -0.4487e-05, -0.6353e-05, -0.9026e-05,           &
    -0.1290e-04, -0.1851e-04, -0.2650e-04, -0.3772e-04, -0.5319e-04,           &
    -0.7431e-04, -0.1031e-03, -0.1422e-03, -0.1951e-03, -0.2648e-03,           &
    -0.3519e-03, -0.4518e-03, -0.5537e-03, -0.6449e-03, -0.7133e-03,           &
    -0.7432e-03, -0.7133e-03, -0.6018e-03, -0.4092e-03, -0.1951e-03,           &
    -0.5345e-04/
   data ((h13(ip,iw,1),iw=1,31),ip=17,17)/                                     &
     0.6457e-09,  0.1235e-08,  0.2303e-08,  0.4149e-08,  0.7120e-08,           &
     0.1152e-07,  0.1749e-07,  0.2508e-07,  0.3462e-07,  0.4718e-07,           &
     0.6452e-07,  0.8874e-07,  0.1222e-06,  0.1675e-06,  0.2276e-06,           &
     0.3076e-06,  0.4174e-06,  0.5714e-06,  0.7837e-06,  0.1067e-05,           &
     0.1428e-05,  0.1859e-05,  0.2327e-05,  0.2760e-05,  0.3122e-05,           &
     0.3458e-05,  0.3739e-05,  0.3715e-05,  0.3126e-05,  0.1942e-05,           &
     0.6977e-06/
   data ((h11(ip,iw,1),iw=1,31),ip=18,18)/                                     &
     0.9151e-04,  0.1750e-03,  0.3265e-03,  0.5891e-03,  0.1020e-02,           &
     0.1688e-02,  0.2679e-02,  0.4109e-02,  0.6138e-02,  0.9002e-02,           &
     0.1304e-01,  0.1871e-01,  0.2667e-01,  0.3773e-01,  0.5299e-01,           &
     0.7386e-01,  0.1022e+00,  0.1403e+00,  0.1905e+00,  0.2546e+00,           &
     0.3335e+00,  0.4258e+00,  0.5269e+00,  0.6297e+00,  0.7280e+00,           &
     0.8178e+00,  0.8950e+00,  0.9529e+00,  0.9864e+00,  0.9983e+00,           &
     0.1000e+01/
   data ((h12(ip,iw,1),iw=1,31),ip=18,18)/                                     &
    -0.1431e-06, -0.2731e-06, -0.5072e-06, -0.9057e-06, -0.1537e-05,           &
    -0.2460e-05, -0.3733e-05, -0.5449e-05, -0.7786e-05, -0.1106e-04,           &
    -0.1574e-04, -0.2249e-04, -0.3212e-04, -0.4564e-04, -0.6438e-04,           &
    -0.9019e-04, -0.1256e-03, -0.1737e-03, -0.2378e-03, -0.3196e-03,           &
    -0.4163e-03, -0.5191e-03, -0.6154e-03, -0.6931e-03, -0.7384e-03,           &
    -0.7313e-03, -0.6492e-03, -0.4805e-03, -0.2629e-03, -0.8897e-04,           &
    -0.1432e-04/
   data ((h13(ip,iw,1),iw=1,31),ip=18,18)/                                     &
     0.6607e-09,  0.1282e-08,  0.2441e-08,  0.4522e-08,  0.8027e-08,           &
     0.1348e-07,  0.2122e-07,  0.3139e-07,  0.4435e-07,  0.6095e-07,           &
     0.8319e-07,  0.1139e-06,  0.1557e-06,  0.2107e-06,  0.2819e-06,           &
     0.3773e-06,  0.5107e-06,  0.6982e-06,  0.9542e-06,  0.1290e-05,           &
     0.1703e-05,  0.2170e-05,  0.2628e-05,  0.3013e-05,  0.3352e-05,           &
     0.3669e-05,  0.3780e-05,  0.3397e-05,  0.2386e-05,  0.1062e-05,           &
     0.2216e-06/
   data ((h11(ip,iw,1),iw=1,31),ip=19,19)/                                     &
     0.9320e-04,  0.1803e-03,  0.3422e-03,  0.6317e-03,  0.1124e-02,           &
     0.1915e-02,  0.3121e-02,  0.4890e-02,  0.7421e-02,  0.1100e-01,           &
     0.1602e-01,  0.2306e-01,  0.3288e-01,  0.4647e-01,  0.6515e-01,           &
     0.9066e-01,  0.1252e+00,  0.1710e+00,  0.2304e+00,  0.3045e+00,           &
     0.3928e+00,  0.4918e+00,  0.5951e+00,  0.6956e+00,  0.7887e+00,           &
     0.8709e+00,  0.9361e+00,  0.9780e+00,  0.9961e+00,  0.9999e+00,           &
     0.1000e+01/
   data ((h12(ip,iw,1),iw=1,31),ip=19,19)/                                     &
    -0.1454e-06, -0.2805e-06, -0.5296e-06, -0.9685e-06, -0.1695e-05,           &
    -0.2812e-05, -0.4412e-05, -0.6606e-05, -0.9573e-05, -0.1363e-04,           &
    -0.1932e-04, -0.2743e-04, -0.3897e-04, -0.5520e-04, -0.7787e-04,           &
    -0.1094e-03, -0.1529e-03, -0.2117e-03, -0.2880e-03, -0.3809e-03,           &
    -0.4834e-03, -0.5836e-03, -0.6692e-03, -0.7275e-03, -0.7408e-03,           &
    -0.6865e-03, -0.5459e-03, -0.3360e-03, -0.1365e-03, -0.2935e-04,           &
    -0.2173e-05/
   data ((h13(ip,iw,1),iw=1,31),ip=19,19)/                                     &
     0.6693e-09,  0.1312e-08,  0.2538e-08,  0.4802e-08,  0.8778e-08,           &
     0.1528e-07,  0.2501e-07,  0.3836e-07,  0.5578e-07,  0.7806e-07,           &
     0.1069e-06,  0.1456e-06,  0.1970e-06,  0.2631e-06,  0.3485e-06,           &
     0.4642e-06,  0.6268e-06,  0.8526e-06,  0.1157e-05,  0.1545e-05,           &
     0.2002e-05,  0.2478e-05,  0.2897e-05,  0.3245e-05,  0.3578e-05,           &
     0.3789e-05,  0.3598e-05,  0.2792e-05,  0.1489e-05,  0.4160e-06,           &
     0.3843e-07/
   data ((h11(ip,iw,1),iw=1,31),ip=20,20)/                                     &
     0.9443e-04,  0.1843e-03,  0.3545e-03,  0.6669e-03,  0.1217e-02,           &
     0.2133e-02,  0.3577e-02,  0.5742e-02,  0.8880e-02,  0.1333e-01,           &
     0.1958e-01,  0.2830e-01,  0.4039e-01,  0.5705e-01,  0.7988e-01,           &
     0.1110e+00,  0.1526e+00,  0.2072e+00,  0.2764e+00,  0.3604e+00,           &
     0.4567e+00,  0.5597e+00,  0.6620e+00,  0.7581e+00,  0.8447e+00,           &
     0.9165e+00,  0.9668e+00,  0.9923e+00,  0.9995e+00,  0.1000e+01,           &
     0.1000e+01/
   data ((h12(ip,iw,1),iw=1,31),ip=20,20)/                                     &
    -0.1472e-06, -0.2866e-06, -0.5485e-06, -0.1024e-05, -0.1842e-05,           &
    -0.3160e-05, -0.5136e-05, -0.7922e-05, -0.1171e-04, -0.1682e-04,           &
    -0.2381e-04, -0.3355e-04, -0.4729e-04, -0.6673e-04, -0.9417e-04,           &
    -0.1327e-03, -0.1858e-03, -0.2564e-03, -0.3449e-03, -0.4463e-03,           &
    -0.5495e-03, -0.6420e-03, -0.7116e-03, -0.7427e-03, -0.7139e-03,           &
    -0.6031e-03, -0.4104e-03, -0.1957e-03, -0.5358e-04, -0.6176e-05,           &
    -0.1347e-06/
   data ((h13(ip,iw,1),iw=1,31),ip=20,20)/                                     &
     0.6750e-09,  0.1332e-08,  0.2602e-08,  0.5003e-08,  0.9367e-08,           &
     0.1684e-07,  0.2863e-07,  0.4566e-07,  0.6865e-07,  0.9861e-07,           &
     0.1368e-06,  0.1856e-06,  0.2479e-06,  0.3274e-06,  0.4315e-06,           &
     0.5739e-06,  0.7710e-06,  0.1040e-05,  0.1394e-05,  0.1829e-05,           &
     0.2309e-05,  0.2759e-05,  0.3131e-05,  0.3472e-05,  0.3755e-05,           &
     0.3730e-05,  0.3138e-05,  0.1948e-05,  0.6994e-06,  0.1022e-06,           &
     0.2459e-08/
   data ((h11(ip,iw,1),iw=1,31),ip=21,21)/                                     &
     0.9530e-04,  0.1872e-03,  0.3638e-03,  0.6949e-03,  0.1295e-02,           &
     0.2330e-02,  0.4022e-02,  0.6633e-02,  0.1049e-01,  0.1601e-01,           &
     0.2378e-01,  0.3456e-01,  0.4941e-01,  0.6980e-01,  0.9765e-01,           &
     0.1353e+00,  0.1852e+00,  0.2494e+00,  0.3287e+00,  0.4216e+00,           &
     0.5234e+00,  0.6272e+00,  0.7261e+00,  0.8166e+00,  0.8943e+00,           &
     0.9526e+00,  0.9863e+00,  0.9983e+00,  0.1000e+01,  0.1000e+01,           &
     0.1000e+01/
   data ((h12(ip,iw,1),iw=1,31),ip=21,21)/                                     &
    -0.1487e-06, -0.2912e-06, -0.5636e-06, -0.1069e-05, -0.1969e-05,           &
    -0.3483e-05, -0.5858e-05, -0.9334e-05, -0.1416e-04, -0.2067e-04,           &
    -0.2936e-04, -0.4113e-04, -0.5750e-04, -0.8072e-04, -0.1139e-03,           &
    -0.1606e-03, -0.2246e-03, -0.3076e-03, -0.4067e-03, -0.5121e-03,           &
    -0.6110e-03, -0.6909e-03, -0.7378e-03, -0.7321e-03, -0.6509e-03,           &
    -0.4825e-03, -0.2641e-03, -0.8936e-04, -0.1436e-04, -0.5966e-06,           &
     0.0000e+00/
   data ((h13(ip,iw,1),iw=1,31),ip=21,21)/                                     &
     0.6777e-09,  0.1344e-08,  0.2643e-08,  0.5138e-08,  0.9798e-08,           &
     0.1809e-07,  0.3185e-07,  0.5285e-07,  0.8249e-07,  0.1222e-06,           &
     0.1730e-06,  0.2351e-06,  0.3111e-06,  0.4078e-06,  0.5366e-06,           &
     0.7117e-06,  0.9495e-06,  0.1266e-05,  0.1667e-05,  0.2132e-05,           &
     0.2600e-05,  0.3001e-05,  0.3354e-05,  0.3679e-05,  0.3796e-05,           &
     0.3414e-05,  0.2399e-05,  0.1067e-05,  0.2222e-06,  0.1075e-07,           &
     0.0000e+00/
   data ((h11(ip,iw,1),iw=1,31),ip=22,22)/                                     &
     0.9589e-04,  0.1893e-03,  0.3706e-03,  0.7162e-03,  0.1357e-02,           &
     0.2500e-02,  0.4434e-02,  0.7523e-02,  0.1220e-01,  0.1900e-01,           &
     0.2860e-01,  0.4190e-01,  0.6015e-01,  0.8508e-01,  0.1189e+00,           &
     0.1643e+00,  0.2233e+00,  0.2976e+00,  0.3865e+00,  0.4865e+00,           &
     0.5909e+00,  0.6925e+00,  0.7866e+00,  0.8696e+00,  0.9355e+00,           &
     0.9778e+00,  0.9961e+00,  0.9999e+00,  0.1000e+01,  0.1000e+01,           &
     0.1000e+01/
   data ((h12(ip,iw,1),iw=1,31),ip=22,22)/                                     &
    -0.1496e-06, -0.2947e-06, -0.5749e-06, -0.1105e-05, -0.2074e-05,           &
    -0.3763e-05, -0.6531e-05, -0.1076e-04, -0.1682e-04, -0.2509e-04,           &
    -0.3605e-04, -0.5049e-04, -0.7012e-04, -0.9787e-04, -0.1378e-03,           &
    -0.1939e-03, -0.2695e-03, -0.3641e-03, -0.4703e-03, -0.5750e-03,           &
    -0.6648e-03, -0.7264e-03, -0.7419e-03, -0.6889e-03, -0.5488e-03,           &
    -0.3382e-03, -0.1375e-03, -0.2951e-04, -0.2174e-05,  0.0000e+00,           &
     0.0000e+00/
   data ((h13(ip,iw,1),iw=1,31),ip=22,22)/                                     &
     0.6798e-09,  0.1350e-08,  0.2667e-08,  0.5226e-08,  0.1010e-07,           &
     0.1903e-07,  0.3455e-07,  0.5951e-07,  0.9658e-07,  0.1479e-06,           &
     0.2146e-06,  0.2951e-06,  0.3903e-06,  0.5101e-06,  0.6693e-06,           &
     0.8830e-06,  0.1168e-05,  0.1532e-05,  0.1968e-05,  0.2435e-05,           &
     0.2859e-05,  0.3222e-05,  0.3572e-05,  0.3797e-05,  0.3615e-05,           &
     0.2811e-05,  0.1500e-05,  0.4185e-06,  0.3850e-07,  0.0000e+00,           &
     0.0000e+00/
   data ((h11(ip,iw,1),iw=1,31),ip=23,23)/                                     &
     0.9629e-04,  0.1907e-03,  0.3755e-03,  0.7321e-03,  0.1406e-02,           &
     0.2639e-02,  0.4796e-02,  0.8368e-02,  0.1394e-01,  0.2221e-01,           &
     0.3400e-01,  0.5037e-01,  0.7279e-01,  0.1033e+00,  0.1442e+00,           &
     0.1983e+00,  0.2673e+00,  0.3516e+00,  0.4489e+00,  0.5533e+00,           &
     0.6572e+00,  0.7547e+00,  0.8425e+00,  0.9153e+00,  0.9663e+00,           &
     0.9922e+00,  0.9995e+00,  0.1000e+01,  0.1000e+01,  0.1000e+01,           &
     0.1000e+01/
   data ((h12(ip,iw,1),iw=1,31),ip=23,23)/                                     &
    -0.1503e-06, -0.2971e-06, -0.5832e-06, -0.1131e-05, -0.2154e-05,           &
    -0.3992e-05, -0.7122e-05, -0.1211e-04, -0.1954e-04, -0.2995e-04,           &
    -0.4380e-04, -0.6183e-04, -0.8577e-04, -0.1191e-03, -0.1668e-03,           &
    -0.2333e-03, -0.3203e-03, -0.4237e-03, -0.5324e-03, -0.6318e-03,           &
    -0.7075e-03, -0.7429e-03, -0.7168e-03, -0.6071e-03, -0.4139e-03,           &
    -0.1976e-03, -0.5410e-04, -0.6215e-05, -0.1343e-06,  0.0000e+00,           &
     0.0000e+00/
   data ((h13(ip,iw,1),iw=1,31),ip=23,23)/                                     &
     0.6809e-09,  0.1356e-08,  0.2683e-08,  0.5287e-08,  0.1030e-07,           &
     0.1971e-07,  0.3665e-07,  0.6528e-07,  0.1100e-06,  0.1744e-06,           &
     0.2599e-06,  0.3650e-06,  0.4887e-06,  0.6398e-06,  0.8358e-06,           &
     0.1095e-05,  0.1429e-05,  0.1836e-05,  0.2286e-05,  0.2716e-05,           &
     0.3088e-05,  0.3444e-05,  0.3748e-05,  0.3740e-05,  0.3157e-05,           &
     0.1966e-05,  0.7064e-06,  0.1030e-06,  0.2456e-08,  0.0000e+00,           &
     0.0000e+00/
   data ((h11(ip,iw,1),iw=1,31),ip=24,24)/                                     &
     0.9655e-04,  0.1917e-03,  0.3789e-03,  0.7436e-03,  0.1442e-02,           &
     0.2748e-02,  0.5100e-02,  0.9128e-02,  0.1563e-01,  0.2553e-01,           &
     0.3987e-01,  0.5994e-01,  0.8746e-01,  0.1246e+00,  0.1739e+00,           &
     0.2376e+00,  0.3170e+00,  0.4107e+00,  0.5141e+00,  0.6198e+00,           &
     0.7208e+00,  0.8130e+00,  0.8923e+00,  0.9517e+00,  0.9860e+00,           &
     0.9983e+00,  0.1000e+01,  0.1000e+01,  0.1000e+01,  0.1000e+01,           &
     0.1000e+01/
   data ((h12(ip,iw,1),iw=1,31),ip=24,24)/                                     &
    -0.1508e-06, -0.2989e-06, -0.5892e-06, -0.1151e-05, -0.2216e-05,           &
    -0.4175e-05, -0.7619e-05, -0.1333e-04, -0.2217e-04, -0.3497e-04,           &
    -0.5238e-04, -0.7513e-04, -0.1049e-03, -0.1455e-03, -0.2021e-03,           &
    -0.2790e-03, -0.3757e-03, -0.4839e-03, -0.5902e-03, -0.6794e-03,           &
    -0.7344e-03, -0.7341e-03, -0.6557e-03, -0.4874e-03, -0.2674e-03,           &
    -0.9059e-04, -0.1455e-04, -0.5986e-06,  0.0000e+00,  0.0000e+00,           &
     0.0000e+00/
   data ((h13(ip,iw,1),iw=1,31),ip=24,24)/                                     &
     0.6812e-09,  0.1356e-08,  0.2693e-08,  0.5328e-08,  0.1045e-07,           &
     0.2021e-07,  0.3826e-07,  0.6994e-07,  0.1218e-06,  0.1997e-06,           &
     0.3069e-06,  0.4428e-06,  0.6064e-06,  0.8015e-06,  0.1043e-05,           &
     0.1351e-05,  0.1733e-05,  0.2168e-05,  0.2598e-05,  0.2968e-05,           &
     0.3316e-05,  0.3662e-05,  0.3801e-05,  0.3433e-05,  0.2422e-05,           &
     0.1081e-05,  0.2256e-06,  0.1082e-07,  0.0000e+00,  0.0000e+00,           &
     0.0000e+00/
   data ((h11(ip,iw,1),iw=1,31),ip=25,25)/                                     &
     0.9672e-04,  0.1923e-03,  0.3813e-03,  0.7518e-03,  0.1469e-02,           &
     0.2832e-02,  0.5343e-02,  0.9779e-02,  0.1719e-01,  0.2882e-01,           &
     0.4606e-01,  0.7051e-01,  0.1042e+00,  0.1493e+00,  0.2081e+00,           &
     0.2824e+00,  0.3720e+00,  0.4736e+00,  0.5803e+00,  0.6844e+00,           &
     0.7810e+00,  0.8663e+00,  0.9339e+00,  0.9772e+00,  0.9960e+00,           &
     0.9999e+00,  0.1000e+01,  0.1000e+01,  0.1000e+01,  0.1000e+01,           &
     0.1000e+01/
   data ((h12(ip,iw,1),iw=1,31),ip=25,25)/                                     &
    -0.1511e-06, -0.3001e-06, -0.5934e-06, -0.1166e-05, -0.2263e-05,           &
    -0.4319e-05, -0.8028e-05, -0.1438e-04, -0.2460e-04, -0.3991e-04,           &
    -0.6138e-04, -0.9005e-04, -0.1278e-03, -0.1778e-03, -0.2447e-03,           &
    -0.3313e-03, -0.4342e-03, -0.5424e-03, -0.6416e-03, -0.7146e-03,           &
    -0.7399e-03, -0.6932e-03, -0.5551e-03, -0.3432e-03, -0.1398e-03,           &
    -0.3010e-04, -0.2229e-05,  0.0000e+00,  0.0000e+00,  0.0000e+00,           &
     0.0000e+00/
   data ((h13(ip,iw,1),iw=1,31),ip=25,25)/                                     &
     0.6815e-09,  0.1358e-08,  0.2698e-08,  0.5355e-08,  0.1054e-07,           &
     0.2056e-07,  0.3942e-07,  0.7349e-07,  0.1315e-06,  0.2226e-06,           &
     0.3537e-06,  0.5266e-06,  0.7407e-06,  0.9958e-06,  0.1296e-05,           &
     0.1657e-05,  0.2077e-05,  0.2512e-05,  0.2893e-05,  0.3216e-05,           &
     0.3562e-05,  0.3811e-05,  0.3644e-05,  0.2841e-05,  0.1524e-05,           &
     0.4276e-06,  0.3960e-07,  0.0000e+00,  0.0000e+00,  0.0000e+00,           &
     0.0000e+00/
   data ((h11(ip,iw,1),iw=1,31),ip=26,26)/                                     &
     0.9682e-04,  0.1928e-03,  0.3829e-03,  0.7576e-03,  0.1489e-02,           &
     0.2894e-02,  0.5533e-02,  0.1031e-01,  0.1856e-01,  0.3195e-01,           &
     0.5237e-01,  0.8187e-01,  0.1227e+00,  0.1771e+00,  0.2468e+00,           &
     0.3321e+00,  0.4312e+00,  0.5384e+00,  0.6455e+00,  0.7463e+00,           &
     0.8372e+00,  0.9126e+00,  0.9652e+00,  0.9919e+00,  0.9994e+00,           &
     0.1000e+01,  0.1000e+01,  0.1000e+01,  0.1000e+01,  0.1000e+01,           &
     0.1000e+01/
   data ((h12(ip,iw,1),iw=1,31),ip=26,26)/                                     &
    -0.1513e-06, -0.3009e-06, -0.5966e-06, -0.1176e-05, -0.2299e-05,           &
    -0.4430e-05, -0.8352e-05, -0.1526e-04, -0.2674e-04, -0.4454e-04,           &
    -0.7042e-04, -0.1062e-03, -0.1540e-03, -0.2163e-03, -0.2951e-03,           &
    -0.3899e-03, -0.4948e-03, -0.5983e-03, -0.6846e-03, -0.7332e-03,           &
    -0.7182e-03, -0.6142e-03, -0.4209e-03, -0.2014e-03, -0.5530e-04,           &
    -0.6418e-05, -0.1439e-06,  0.0000e+00,  0.0000e+00,  0.0000e+00,           &
     0.0000e+00/
   data ((h13(ip,iw,1),iw=1,31),ip=26,26)/                                     &
     0.6817e-09,  0.1359e-08,  0.2702e-08,  0.5374e-08,  0.1061e-07,           &
     0.2079e-07,  0.4022e-07,  0.7610e-07,  0.1392e-06,  0.2428e-06,           &
     0.3992e-06,  0.6149e-06,  0.8893e-06,  0.1220e-05,  0.1599e-05,           &
     0.2015e-05,  0.2453e-05,  0.2853e-05,  0.3173e-05,  0.3488e-05,           &
     0.3792e-05,  0.3800e-05,  0.3210e-05,  0.2002e-05,  0.7234e-06,           &
     0.1068e-06,  0.2646e-08,  0.0000e+00,  0.0000e+00,  0.0000e+00,           &
     0.0000e+00/
   data ((h21(ip,iw,1),iw=1,31),ip= 1, 1)/                                     &
     0.3920e-05,  0.7617e-05,  0.1455e-04,  0.2706e-04,  0.4855e-04,           &
     0.8315e-04,  0.1349e-03,  0.2063e-03,  0.2984e-03,  0.4109e-03,           &
     0.5422e-03,  0.6896e-03,  0.8537e-03,  0.1041e-02,  0.1262e-02,           &
     0.1534e-02,  0.1870e-02,  0.2286e-02,  0.2803e-02,  0.3444e-02,           &
     0.4242e-02,  0.5244e-02,  0.6511e-02,  0.8138e-02,  0.1027e-01,           &
     0.1312e-01,  0.1697e-01,  0.2222e-01,  0.2941e-01,  0.3923e-01,           &
     0.5258e-01/
   data ((h22(ip,iw,1),iw=1,31),ip= 1, 1)/                                     &
    -0.5622e-07, -0.1071e-06, -0.1983e-06, -0.3533e-06, -0.5991e-06,           &
    -0.9592e-06, -0.1444e-05, -0.2049e-05, -0.2764e-05, -0.3577e-05,           &
    -0.4469e-05, -0.5467e-05, -0.6654e-05, -0.8137e-05, -0.1002e-04,           &
    -0.1237e-04, -0.1528e-04, -0.1884e-04, -0.2310e-04, -0.2809e-04,           &
    -0.3396e-04, -0.4098e-04, -0.4960e-04, -0.6058e-04, -0.7506e-04,           &
    -0.9451e-04, -0.1207e-03, -0.1558e-03, -0.2026e-03, -0.2648e-03,           &
    -0.3468e-03/
   data ((h23(ip,iw,1),iw=1,31),ip= 1, 1)/                                     &
    -0.2195e-09, -0.4031e-09, -0.7043e-09, -0.1153e-08, -0.1737e-08,           &
    -0.2395e-08, -0.3020e-08, -0.3549e-08, -0.4034e-08, -0.4421e-08,           &
    -0.4736e-08, -0.5681e-08, -0.8289e-08, -0.1287e-07, -0.1873e-07,           &
    -0.2523e-07, -0.3223e-07, -0.3902e-07, -0.4409e-07, -0.4699e-07,           &
    -0.4782e-07, -0.4705e-07, -0.4657e-07, -0.4885e-07, -0.5550e-07,           &
    -0.6619e-07, -0.7656e-07, -0.8027e-07, -0.7261e-07, -0.4983e-07,           &
    -0.1101e-07/
   data ((h21(ip,iw,1),iw=1,31),ip= 2, 2)/                                     &
     0.3920e-05,  0.7617e-05,  0.1455e-04,  0.2706e-04,  0.4856e-04,           &
     0.8318e-04,  0.1349e-03,  0.2065e-03,  0.2986e-03,  0.4114e-03,           &
     0.5431e-03,  0.6912e-03,  0.8566e-03,  0.1046e-02,  0.1272e-02,           &
     0.1552e-02,  0.1902e-02,  0.2342e-02,  0.2899e-02,  0.3605e-02,           &
     0.4501e-02,  0.5648e-02,  0.7129e-02,  0.9063e-02,  0.1163e-01,           &
     0.1509e-01,  0.1980e-01,  0.2626e-01,  0.3510e-01,  0.4717e-01,           &
     0.6351e-01/
   data ((h22(ip,iw,1),iw=1,31),ip= 2, 2)/                                     &
    -0.5622e-07, -0.1071e-06, -0.1983e-06, -0.3534e-06, -0.5992e-06,           &
    -0.9594e-06, -0.1445e-05, -0.2050e-05, -0.2766e-05, -0.3580e-05,           &
    -0.4476e-05, -0.5479e-05, -0.6677e-05, -0.8179e-05, -0.1009e-04,           &
    -0.1251e-04, -0.1553e-04, -0.1928e-04, -0.2384e-04, -0.2930e-04,           &
    -0.3588e-04, -0.4393e-04, -0.5403e-04, -0.6714e-04, -0.8458e-04,           &
    -0.1082e-03, -0.1400e-03, -0.1829e-03, -0.2401e-03, -0.3157e-03,           &
    -0.4147e-03/
   data ((h23(ip,iw,1),iw=1,31),ip= 2, 2)/                                     &
    -0.2195e-09, -0.4032e-09, -0.7046e-09, -0.1153e-08, -0.1738e-08,           &
    -0.2395e-08, -0.3021e-08, -0.3550e-08, -0.4035e-08, -0.4423e-08,           &
    -0.4740e-08, -0.5692e-08, -0.8314e-08, -0.1292e-07, -0.1882e-07,           &
    -0.2536e-07, -0.3242e-07, -0.3927e-07, -0.4449e-07, -0.4767e-07,           &
    -0.4889e-07, -0.4857e-07, -0.4860e-07, -0.5132e-07, -0.5847e-07,           &
    -0.6968e-07, -0.8037e-07, -0.8400e-07, -0.7521e-07, -0.4830e-07,           &
    -0.7562e-09/
   data ((h21(ip,iw,1),iw=1,31),ip= 3, 3)/                                     &
     0.3920e-05,  0.7617e-05,  0.1455e-04,  0.2707e-04,  0.4858e-04,           &
     0.8322e-04,  0.1350e-03,  0.2066e-03,  0.2990e-03,  0.4122e-03,           &
     0.5444e-03,  0.6937e-03,  0.8611e-03,  0.1054e-02,  0.1287e-02,           &
     0.1579e-02,  0.1949e-02,  0.2424e-02,  0.3037e-02,  0.3829e-02,           &
     0.4854e-02,  0.6190e-02,  0.7942e-02,  0.1026e-01,  0.1338e-01,           &
     0.1761e-01,  0.2340e-01,  0.3134e-01,  0.4223e-01,  0.5703e-01,           &
     0.7694e-01/
   data ((h22(ip,iw,1),iw=1,31),ip= 3, 3)/                                     &
    -0.5622e-07, -0.1071e-06, -0.1983e-06, -0.3535e-06, -0.5994e-06,           &
    -0.9599e-06, -0.1446e-05, -0.2052e-05, -0.2769e-05, -0.3586e-05,           &
    -0.4487e-05, -0.5499e-05, -0.6712e-05, -0.8244e-05, -0.1021e-04,           &
    -0.1272e-04, -0.1591e-04, -0.1992e-04, -0.2489e-04, -0.3097e-04,           &
    -0.3845e-04, -0.4782e-04, -0.5982e-04, -0.7558e-04, -0.9674e-04,           &
    -0.1254e-03, -0.1644e-03, -0.2167e-03, -0.2863e-03, -0.3777e-03,           &
    -0.4959e-03/
   data ((h23(ip,iw,1),iw=1,31),ip= 3, 3)/                                     &
    -0.2196e-09, -0.4033e-09, -0.7048e-09, -0.1154e-08, -0.1739e-08,           &
    -0.2396e-08, -0.3022e-08, -0.3551e-08, -0.4036e-08, -0.4425e-08,           &
    -0.4746e-08, -0.5710e-08, -0.8354e-08, -0.1300e-07, -0.1894e-07,           &
    -0.2554e-07, -0.3265e-07, -0.3958e-07, -0.4502e-07, -0.4859e-07,           &
    -0.5030e-07, -0.5053e-07, -0.5104e-07, -0.5427e-07, -0.6204e-07,           &
    -0.7388e-07, -0.8477e-07, -0.8760e-07, -0.7545e-07, -0.4099e-07,           &
     0.2046e-07/
   data ((h21(ip,iw,1),iw=1,31),ip= 4, 4)/                                     &
     0.3919e-05,  0.7618e-05,  0.1455e-04,  0.2708e-04,  0.4860e-04,           &
     0.8327e-04,  0.1352e-03,  0.2070e-03,  0.2997e-03,  0.4134e-03,           &
     0.5466e-03,  0.6976e-03,  0.8682e-03,  0.1067e-02,  0.1310e-02,           &
     0.1619e-02,  0.2020e-02,  0.2543e-02,  0.3232e-02,  0.4136e-02,           &
     0.5328e-02,  0.6906e-02,  0.9004e-02,  0.1182e-01,  0.1562e-01,           &
     0.2081e-01,  0.2794e-01,  0.3773e-01,  0.5111e-01,  0.6919e-01,           &
     0.9329e-01/
   data ((h22(ip,iw,1),iw=1,31),ip= 4, 4)/                                     &
    -0.5623e-07, -0.1071e-06, -0.1984e-06, -0.3536e-06, -0.5997e-06,           &
    -0.9606e-06, -0.1447e-05, -0.2055e-05, -0.2775e-05, -0.3596e-05,           &
    -0.4504e-05, -0.5529e-05, -0.6768e-05, -0.8345e-05, -0.1039e-04,           &
    -0.1304e-04, -0.1645e-04, -0.2082e-04, -0.2633e-04, -0.3322e-04,           &
    -0.4187e-04, -0.5292e-04, -0.6730e-04, -0.8640e-04, -0.1122e-03,           &
    -0.1472e-03, -0.1948e-03, -0.2585e-03, -0.3428e-03, -0.4523e-03,           &
    -0.5915e-03/
   data ((h23(ip,iw,1),iw=1,31),ip= 4, 4)/                                     &
    -0.2196e-09, -0.4034e-09, -0.7050e-09, -0.1154e-08, -0.1740e-08,           &
    -0.2398e-08, -0.3024e-08, -0.3552e-08, -0.4037e-08, -0.4428e-08,           &
    -0.4756e-08, -0.5741e-08, -0.8418e-08, -0.1310e-07, -0.1910e-07,           &
    -0.2575e-07, -0.3293e-07, -0.3998e-07, -0.4572e-07, -0.4980e-07,           &
    -0.5211e-07, -0.5287e-07, -0.5390e-07, -0.5782e-07, -0.6650e-07,           &
    -0.7892e-07, -0.8940e-07, -0.8980e-07, -0.7119e-07, -0.2452e-07,           &
     0.5823e-07/
   data ((h21(ip,iw,1),iw=1,31),ip= 5, 5)/                                     &
     0.3919e-05,  0.7618e-05,  0.1455e-04,  0.2709e-04,  0.4863e-04,           &
     0.8337e-04,  0.1354e-03,  0.2075e-03,  0.3007e-03,  0.4152e-03,           &
     0.5500e-03,  0.7038e-03,  0.8793e-03,  0.1086e-02,  0.1345e-02,           &
     0.1680e-02,  0.2123e-02,  0.2712e-02,  0.3500e-02,  0.4552e-02,           &
     0.5958e-02,  0.7844e-02,  0.1038e-01,  0.1381e-01,  0.1847e-01,           &
     0.2487e-01,  0.3366e-01,  0.4572e-01,  0.6209e-01,  0.8406e-01,           &
     0.1130e+00/
   data ((h22(ip,iw,1),iw=1,31),ip= 5, 5)/                                     &
    -0.5623e-07, -0.1071e-06, -0.1985e-06, -0.3538e-06, -0.6002e-06,           &
    -0.9618e-06, -0.1450e-05, -0.2059e-05, -0.2783e-05, -0.3611e-05,           &
    -0.4531e-05, -0.5577e-05, -0.6855e-05, -0.8499e-05, -0.1066e-04,           &
    -0.1351e-04, -0.1723e-04, -0.2207e-04, -0.2829e-04, -0.3621e-04,           &
    -0.4636e-04, -0.5954e-04, -0.7690e-04, -0.1002e-03, -0.1317e-03,           &
    -0.1746e-03, -0.2326e-03, -0.3099e-03, -0.4111e-03, -0.5407e-03,           &
    -0.7020e-03/
   data ((h23(ip,iw,1),iw=1,31),ip= 5, 5)/                                     &
    -0.2197e-09, -0.4037e-09, -0.7054e-09, -0.1155e-08, -0.1741e-08,           &
    -0.2401e-08, -0.3027e-08, -0.3553e-08, -0.4039e-08, -0.4431e-08,           &
    -0.4775e-08, -0.5784e-08, -0.8506e-08, -0.1326e-07, -0.1931e-07,           &
    -0.2600e-07, -0.3324e-07, -0.4048e-07, -0.4666e-07, -0.5137e-07,           &
    -0.5428e-07, -0.5558e-07, -0.5730e-07, -0.6228e-07, -0.7197e-07,           &
    -0.8455e-07, -0.9347e-07, -0.8867e-07, -0.5945e-07,  0.5512e-08,           &
     0.1209e-06/
   data ((h21(ip,iw,1),iw=1,31),ip= 6, 6)/                                     &
     0.3919e-05,  0.7619e-05,  0.1456e-04,  0.2710e-04,  0.4868e-04,           &
     0.8351e-04,  0.1357e-03,  0.2082e-03,  0.3022e-03,  0.4181e-03,           &
     0.5554e-03,  0.7134e-03,  0.8963e-03,  0.1117e-02,  0.1397e-02,           &
     0.1770e-02,  0.2270e-02,  0.2947e-02,  0.3865e-02,  0.5107e-02,           &
     0.6787e-02,  0.9062e-02,  0.1215e-01,  0.1634e-01,  0.2209e-01,           &
     0.2999e-01,  0.4083e-01,  0.5563e-01,  0.7559e-01,  0.1021e+00,           &
     0.1364e+00/
   data ((h22(ip,iw,1),iw=1,31),ip= 6, 6)/                                     &
    -0.5624e-07, -0.1072e-06, -0.1986e-06, -0.3541e-06, -0.6010e-06,           &
    -0.9636e-06, -0.1453e-05, -0.2067e-05, -0.2796e-05, -0.3634e-05,           &
    -0.4572e-05, -0.5652e-05, -0.6987e-05, -0.8733e-05, -0.1107e-04,           &
    -0.1418e-04, -0.1832e-04, -0.2378e-04, -0.3092e-04, -0.4017e-04,           &
    -0.5221e-04, -0.6806e-04, -0.8916e-04, -0.1176e-03, -0.1562e-03,           &
    -0.2087e-03, -0.2793e-03, -0.3724e-03, -0.4928e-03, -0.6440e-03,           &
    -0.8270e-03/
   data ((h23(ip,iw,1),iw=1,31),ip= 6, 6)/                                     &
    -0.2198e-09, -0.4040e-09, -0.7061e-09, -0.1156e-08, -0.1744e-08,           &
    -0.2405e-08, -0.3032e-08, -0.3556e-08, -0.4040e-08, -0.4444e-08,           &
    -0.4800e-08, -0.5848e-08, -0.8640e-08, -0.1346e-07, -0.1957e-07,           &
    -0.2627e-07, -0.3357e-07, -0.4114e-07, -0.4793e-07, -0.5330e-07,           &
    -0.5676e-07, -0.5873e-07, -0.6152e-07, -0.6783e-07, -0.7834e-07,           &
    -0.9023e-07, -0.9530e-07, -0.8162e-07, -0.3634e-07,  0.5638e-07,           &
     0.2189e-06/
   data ((h21(ip,iw,1),iw=1,31),ip= 7, 7)/                                     &
     0.3919e-05,  0.7620e-05,  0.1457e-04,  0.2713e-04,  0.4877e-04,           &
     0.8374e-04,  0.1363e-03,  0.2095e-03,  0.3047e-03,  0.4227e-03,           &
     0.5637e-03,  0.7282e-03,  0.9224e-03,  0.1162e-02,  0.1475e-02,           &
     0.1898e-02,  0.2475e-02,  0.3267e-02,  0.4354e-02,  0.5839e-02,           &
     0.7866e-02,  0.1063e-01,  0.1440e-01,  0.1957e-01,  0.2666e-01,           &
     0.3641e-01,  0.4975e-01,  0.6784e-01,  0.9201e-01,  0.1236e+00,           &
     0.1639e+00/
   data ((h22(ip,iw,1),iw=1,31),ip= 7, 7)/                                     &
    -0.5626e-07, -0.1072e-06, -0.1987e-06, -0.3545e-06, -0.6022e-06,           &
    -0.9665e-06, -0.1460e-05, -0.2078e-05, -0.2817e-05, -0.3671e-05,           &
    -0.4637e-05, -0.5767e-05, -0.7188e-05, -0.9080e-05, -0.1165e-04,           &
    -0.1513e-04, -0.1981e-04, -0.2609e-04, -0.3441e-04, -0.4534e-04,           &
    -0.5978e-04, -0.7897e-04, -0.1047e-03, -0.1396e-03, -0.1870e-03,           &
    -0.2510e-03, -0.3363e-03, -0.4475e-03, -0.5888e-03, -0.7621e-03,           &
    -0.9647e-03/
   data ((h23(ip,iw,1),iw=1,31),ip= 7, 7)/                                     &
    -0.2200e-09, -0.4045e-09, -0.7071e-09, -0.1159e-08, -0.1748e-08,           &
    -0.2411e-08, -0.3040e-08, -0.3561e-08, -0.4046e-08, -0.4455e-08,           &
    -0.4839e-08, -0.5941e-08, -0.8815e-08, -0.1371e-07, -0.1983e-07,           &
    -0.2652e-07, -0.3400e-07, -0.4207e-07, -0.4955e-07, -0.5554e-07,           &
    -0.5966e-07, -0.6261e-07, -0.6688e-07, -0.7454e-07, -0.8521e-07,           &
    -0.9470e-07, -0.9275e-07, -0.6525e-07,  0.3686e-08,  0.1371e-06,           &
     0.3623e-06/
   data ((h21(ip,iw,1),iw=1,31),ip= 8, 8)/                                     &
     0.3919e-05,  0.7622e-05,  0.1458e-04,  0.2717e-04,  0.4890e-04,           &
     0.8410e-04,  0.1372e-03,  0.2114e-03,  0.3085e-03,  0.4298e-03,           &
     0.5765e-03,  0.7508e-03,  0.9618e-03,  0.1229e-02,  0.1586e-02,           &
     0.2077e-02,  0.2757e-02,  0.3697e-02,  0.5000e-02,  0.6794e-02,           &
     0.9260e-02,  0.1264e-01,  0.1728e-01,  0.2365e-01,  0.3242e-01,           &
     0.4444e-01,  0.6080e-01,  0.8278e-01,  0.1118e+00,  0.1491e+00,           &
     0.1956e+00/
   data ((h22(ip,iw,1),iw=1,31),ip= 8, 8)/                                     &
    -0.5628e-07, -0.1073e-06, -0.1990e-06, -0.3553e-06, -0.6042e-06,           &
    -0.9710e-06, -0.1469e-05, -0.2096e-05, -0.2849e-05, -0.3728e-05,           &
    -0.4738e-05, -0.5942e-05, -0.7490e-05, -0.9586e-05, -0.1247e-04,           &
    -0.1644e-04, -0.2184e-04, -0.2916e-04, -0.3898e-04, -0.5205e-04,           &
    -0.6948e-04, -0.9285e-04, -0.1244e-03, -0.1672e-03, -0.2251e-03,           &
    -0.3028e-03, -0.4051e-03, -0.5365e-03, -0.6998e-03, -0.8940e-03,           &
    -0.1112e-02/
   data ((h23(ip,iw,1),iw=1,31),ip= 8, 8)/                                     &
    -0.2204e-09, -0.4052e-09, -0.7088e-09, -0.1162e-08, -0.1755e-08,           &
    -0.2422e-08, -0.3053e-08, -0.3572e-08, -0.4052e-08, -0.4474e-08,           &
    -0.4898e-08, -0.6082e-08, -0.9046e-08, -0.1400e-07, -0.2009e-07,           &
    -0.2683e-07, -0.3463e-07, -0.4334e-07, -0.5153e-07, -0.5811e-07,           &
    -0.6305e-07, -0.6749e-07, -0.7346e-07, -0.8208e-07, -0.9173e-07,           &
    -0.9603e-07, -0.8264e-07, -0.3505e-07,  0.6878e-07,  0.2586e-06,           &
     0.5530e-06/
   data ((h21(ip,iw,1),iw=1,31),ip= 9, 9)/                                     &
     0.3919e-05,  0.7625e-05,  0.1459e-04,  0.2723e-04,  0.4910e-04,           &
     0.8465e-04,  0.1385e-03,  0.2143e-03,  0.3143e-03,  0.4407e-03,           &
     0.5960e-03,  0.7850e-03,  0.1020e-02,  0.1326e-02,  0.1744e-02,           &
     0.2324e-02,  0.3135e-02,  0.4267e-02,  0.5845e-02,  0.8031e-02,           &
     0.1105e-01,  0.1520e-01,  0.2093e-01,  0.2880e-01,  0.3962e-01,           &
     0.5440e-01,  0.7435e-01,  0.1009e+00,  0.1353e+00,  0.1787e+00,           &
     0.2317e+00/
   data ((h22(ip,iw,1),iw=1,31),ip= 9, 9)/                                     &
    -0.5630e-07, -0.1074e-06, -0.1994e-06, -0.3564e-06, -0.6072e-06,           &
    -0.9779e-06, -0.1484e-05, -0.2124e-05, -0.2900e-05, -0.3817e-05,           &
    -0.4891e-05, -0.6205e-05, -0.7931e-05, -0.1031e-04, -0.1362e-04,           &
    -0.1821e-04, -0.2454e-04, -0.3320e-04, -0.4493e-04, -0.6068e-04,           &
    -0.8186e-04, -0.1104e-03, -0.1491e-03, -0.2016e-03, -0.2722e-03,           &
    -0.3658e-03, -0.4873e-03, -0.6404e-03, -0.8253e-03, -0.1037e-02,           &
    -0.1267e-02/
   data ((h23(ip,iw,1),iw=1,31),ip= 9, 9)/                                     &
    -0.2207e-09, -0.4061e-09, -0.7117e-09, -0.1169e-08, -0.1767e-08,           &
    -0.2439e-08, -0.3074e-08, -0.3588e-08, -0.4062e-08, -0.4510e-08,           &
    -0.4983e-08, -0.6261e-08, -0.9324e-08, -0.1430e-07, -0.2036e-07,           &
    -0.2725e-07, -0.3561e-07, -0.4505e-07, -0.5384e-07, -0.6111e-07,           &
    -0.6731e-07, -0.7355e-07, -0.8112e-07, -0.8978e-07, -0.9616e-07,           &
    -0.9157e-07, -0.6114e-07,  0.1622e-07,  0.1694e-06,  0.4277e-06,           &
     0.7751e-06/
   data ((h21(ip,iw,1),iw=1,31),ip=10,10)/                                     &
     0.3919e-05,  0.7629e-05,  0.1462e-04,  0.2733e-04,  0.4940e-04,           &
     0.8549e-04,  0.1406e-03,  0.2188e-03,  0.3232e-03,  0.4572e-03,           &
     0.6253e-03,  0.8354e-03,  0.1105e-02,  0.1464e-02,  0.1960e-02,           &
     0.2657e-02,  0.3637e-02,  0.5014e-02,  0.6942e-02,  0.9622e-02,           &
     0.1333e-01,  0.1846e-01,  0.2553e-01,  0.3526e-01,  0.4859e-01,           &
     0.6667e-01,  0.9084e-01,  0.1225e+00,  0.1629e+00,  0.2127e+00,           &
     0.2721e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=10,10)/                                     &
    -0.5636e-07, -0.1076e-06, -0.2000e-06, -0.3582e-06, -0.6119e-06,           &
    -0.9888e-06, -0.1507e-05, -0.2168e-05, -0.2978e-05, -0.3952e-05,           &
    -0.5122e-05, -0.6592e-05, -0.8565e-05, -0.1132e-04, -0.1518e-04,           &
    -0.2060e-04, -0.2811e-04, -0.3848e-04, -0.5261e-04, -0.7173e-04,           &
    -0.9758e-04, -0.1326e-03, -0.1801e-03, -0.2442e-03, -0.3296e-03,           &
    -0.4415e-03, -0.5840e-03, -0.7591e-03, -0.9636e-03, -0.1189e-02,           &
    -0.1427e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=10,10)/                                     &
    -0.2214e-09, -0.4080e-09, -0.7156e-09, -0.1178e-08, -0.1784e-08,           &
    -0.2466e-08, -0.3105e-08, -0.3617e-08, -0.4087e-08, -0.4563e-08,           &
    -0.5110e-08, -0.6492e-08, -0.9643e-08, -0.1461e-07, -0.2069e-07,           &
    -0.2796e-07, -0.3702e-07, -0.4717e-07, -0.5662e-07, -0.6484e-07,           &
    -0.7271e-07, -0.8079e-07, -0.8928e-07, -0.9634e-07, -0.9625e-07,           &
    -0.7776e-07, -0.2242e-07,  0.9745e-07,  0.3152e-06,  0.6388e-06,           &
     0.9992e-06/
   data ((h21(ip,iw,1),iw=1,31),ip=11,11)/                                     &
     0.3919e-05,  0.7639e-05,  0.1466e-04,  0.2748e-04,  0.4986e-04,           &
     0.8675e-04,  0.1437e-03,  0.2255e-03,  0.3365e-03,  0.4816e-03,           &
     0.6682e-03,  0.9082e-03,  0.1224e-02,  0.1653e-02,  0.2253e-02,           &
     0.3099e-02,  0.4296e-02,  0.5984e-02,  0.8352e-02,  0.1165e-01,           &
     0.1623e-01,  0.2257e-01,  0.3132e-01,  0.4333e-01,  0.5968e-01,           &
     0.8166e-01,  0.1107e+00,  0.1481e+00,  0.1947e+00,  0.2509e+00,           &
     0.3170e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=11,11)/                                     &
    -0.5645e-07, -0.1079e-06, -0.2009e-06, -0.3610e-06, -0.6190e-06,           &
    -0.1005e-05, -0.1541e-05, -0.2235e-05, -0.3096e-05, -0.4155e-05,           &
    -0.5463e-05, -0.7150e-05, -0.9453e-05, -0.1269e-04, -0.1728e-04,           &
    -0.2375e-04, -0.3279e-04, -0.4531e-04, -0.6246e-04, -0.8579e-04,           &
    -0.1175e-03, -0.1604e-03, -0.2186e-03, -0.2964e-03, -0.3990e-03,           &
    -0.5311e-03, -0.6957e-03, -0.8916e-03, -0.1112e-02, -0.1346e-02,           &
    -0.1590e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=11,11)/                                     &
    -0.2225e-09, -0.4104e-09, -0.7217e-09, -0.1192e-08, -0.1811e-08,           &
    -0.2509e-08, -0.3155e-08, -0.3668e-08, -0.4138e-08, -0.4650e-08,           &
    -0.5296e-08, -0.6785e-08, -0.9991e-08, -0.1494e-07, -0.2122e-07,           &
    -0.2911e-07, -0.3895e-07, -0.4979e-07, -0.6002e-07, -0.6964e-07,           &
    -0.7935e-07, -0.8887e-07, -0.9699e-07, -0.9967e-07, -0.8883e-07,           &
    -0.4988e-07,  0.4156e-07,  0.2197e-06,  0.5081e-06,  0.8667e-06,           &
     0.1212e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=12,12)/                                     &
     0.3920e-05,  0.7649e-05,  0.1472e-04,  0.2768e-04,  0.5051e-04,           &
     0.8856e-04,  0.1481e-03,  0.2352e-03,  0.3557e-03,  0.5169e-03,           &
     0.7295e-03,  0.1010e-02,  0.1388e-02,  0.1909e-02,  0.2640e-02,           &
     0.3678e-02,  0.5151e-02,  0.7231e-02,  0.1015e-01,  0.1423e-01,           &
     0.1990e-01,  0.2775e-01,  0.3856e-01,  0.5333e-01,  0.7328e-01,           &
     0.9980e-01,  0.1343e+00,  0.1778e+00,  0.2308e+00,  0.2936e+00,           &
     0.3667e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=12,12)/                                     &
    -0.5658e-07, -0.1083e-06, -0.2023e-06, -0.3650e-06, -0.6295e-06,           &
    -0.1030e-05, -0.1593e-05, -0.2334e-05, -0.3273e-05, -0.4455e-05,           &
    -0.5955e-05, -0.7935e-05, -0.1067e-04, -0.1455e-04, -0.2007e-04,           &
    -0.2788e-04, -0.3885e-04, -0.5409e-04, -0.7503e-04, -0.1036e-03,           &
    -0.1425e-03, -0.1951e-03, -0.2660e-03, -0.3598e-03, -0.4817e-03,           &
    -0.6355e-03, -0.8218e-03, -0.1035e-02, -0.1267e-02, -0.1508e-02,           &
    -0.1755e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=12,12)/                                     &
    -0.2241e-09, -0.4142e-09, -0.7312e-09, -0.1214e-08, -0.1854e-08,           &
    -0.2578e-08, -0.3250e-08, -0.3765e-08, -0.4238e-08, -0.4809e-08,           &
    -0.5553e-08, -0.7132e-08, -0.1035e-07, -0.1538e-07, -0.2211e-07,           &
    -0.3079e-07, -0.4142e-07, -0.5303e-07, -0.6437e-07, -0.7566e-07,           &
    -0.8703e-07, -0.9700e-07, -0.1025e-06, -0.9718e-07, -0.6973e-07,           &
    -0.5265e-09,  0.1413e-06,  0.3895e-06,  0.7321e-06,  0.1085e-05,           &
     0.1449e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=13,13)/                                     &
     0.3920e-05,  0.7668e-05,  0.1480e-04,  0.2798e-04,  0.5143e-04,           &
     0.9106e-04,  0.1542e-03,  0.2488e-03,  0.3827e-03,  0.5661e-03,           &
     0.8143e-03,  0.1150e-02,  0.1608e-02,  0.2246e-02,  0.3147e-02,           &
     0.4428e-02,  0.6248e-02,  0.8823e-02,  0.1244e-01,  0.1750e-01,           &
     0.2452e-01,  0.3424e-01,  0.4756e-01,  0.6564e-01,  0.8981e-01,           &
     0.1215e+00,  0.1619e+00,  0.2118e+00,  0.2713e+00,  0.3409e+00,           &
     0.4215e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=13,13)/                                     &
    -0.5677e-07, -0.1090e-06, -0.2043e-06, -0.3709e-06, -0.6448e-06,           &
    -0.1066e-05, -0.1669e-05, -0.2480e-05, -0.3532e-05, -0.4887e-05,           &
    -0.6650e-05, -0.9017e-05, -0.1232e-04, -0.1701e-04, -0.2372e-04,           &
    -0.3325e-04, -0.4664e-04, -0.6528e-04, -0.9095e-04, -0.1260e-03,           &
    -0.1737e-03, -0.2381e-03, -0.3238e-03, -0.4359e-03, -0.5789e-03,           &
    -0.7549e-03, -0.9607e-03, -0.1188e-02, -0.1427e-02, -0.1674e-02,           &
    -0.1914e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=13,13)/                                     &
    -0.2262e-09, -0.4191e-09, -0.7441e-09, -0.1244e-08, -0.1916e-08,           &
    -0.2687e-08, -0.3407e-08, -0.3947e-08, -0.4432e-08, -0.5059e-08,           &
    -0.5896e-08, -0.7538e-08, -0.1079e-07, -0.1606e-07, -0.2346e-07,           &
    -0.3305e-07, -0.4459e-07, -0.5720e-07, -0.6999e-07, -0.8300e-07,           &
    -0.9536e-07, -0.1039e-06, -0.1036e-06, -0.8526e-07, -0.3316e-07,           &
     0.7909e-07,  0.2865e-06,  0.6013e-06,  0.9580e-06,  0.1303e-05,           &
     0.1792e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=14,14)/                                     &
     0.3923e-05,  0.7691e-05,  0.1491e-04,  0.2836e-04,  0.5261e-04,           &
     0.9432e-04,  0.1623e-03,  0.2667e-03,  0.4187e-03,  0.6323e-03,           &
     0.9278e-03,  0.1334e-02,  0.1897e-02,  0.2685e-02,  0.3800e-02,           &
     0.5388e-02,  0.7646e-02,  0.1084e-01,  0.1533e-01,  0.2160e-01,           &
     0.3031e-01,  0.4231e-01,  0.5866e-01,  0.8064e-01,  0.1097e+00,           &
     0.1471e+00,  0.1938e+00,  0.2501e+00,  0.3163e+00,  0.3932e+00,           &
     0.4810e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=14,14)/                                     &
    -0.5703e-07, -0.1098e-06, -0.2071e-06, -0.3788e-06, -0.6657e-06,           &
    -0.1116e-05, -0.1776e-05, -0.2687e-05, -0.3898e-05, -0.5493e-05,           &
    -0.7607e-05, -0.1048e-04, -0.1450e-04, -0.2024e-04, -0.2845e-04,           &
    -0.4014e-04, -0.5658e-04, -0.7947e-04, -0.1110e-03, -0.1541e-03,           &
    -0.2125e-03, -0.2907e-03, -0.3936e-03, -0.5261e-03, -0.6912e-03,           &
    -0.8880e-03, -0.1109e-02, -0.1346e-02, -0.1591e-02, -0.1837e-02,           &
    -0.2054e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=14,14)/                                     &
    -0.2288e-09, -0.4265e-09, -0.7627e-09, -0.1289e-08, -0.2011e-08,           &
    -0.2861e-08, -0.3673e-08, -0.4288e-08, -0.4812e-08, -0.5475e-08,           &
    -0.6365e-08, -0.8052e-08, -0.1142e-07, -0.1711e-07, -0.2533e-07,           &
    -0.3597e-07, -0.4862e-07, -0.6259e-07, -0.7708e-07, -0.9150e-07,           &
    -0.1035e-06, -0.1079e-06, -0.9742e-07, -0.5928e-07,  0.2892e-07,           &
     0.1998e-06,  0.4789e-06,  0.8298e-06,  0.1172e-05,  0.1583e-05,           &
     0.2329e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=15,15)/                                     &
     0.3928e-05,  0.7724e-05,  0.1504e-04,  0.2883e-04,  0.5403e-04,           &
     0.9826e-04,  0.1722e-03,  0.2891e-03,  0.4644e-03,  0.7172e-03,           &
     0.1074e-02,  0.1572e-02,  0.2267e-02,  0.3245e-02,  0.4632e-02,           &
     0.6606e-02,  0.9411e-02,  0.1338e-01,  0.1896e-01,  0.2674e-01,           &
     0.3751e-01,  0.5227e-01,  0.7222e-01,  0.9875e-01,  0.1332e+00,           &
     0.1768e+00,  0.2299e+00,  0.2927e+00,  0.3660e+00,  0.4504e+00,           &
     0.5444e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=15,15)/                                     &
    -0.5736e-07, -0.1109e-06, -0.2106e-06, -0.3890e-06, -0.6928e-06,           &
    -0.1181e-05, -0.1917e-05, -0.2965e-05, -0.4396e-05, -0.6315e-05,           &
    -0.8891e-05, -0.1242e-04, -0.1736e-04, -0.2442e-04, -0.3454e-04,           &
    -0.4892e-04, -0.6916e-04, -0.9735e-04, -0.1362e-03, -0.1891e-03,           &
    -0.2602e-03, -0.3545e-03, -0.4768e-03, -0.6310e-03, -0.8179e-03,           &
    -0.1032e-02, -0.1265e-02, -0.1508e-02, -0.1757e-02, -0.1989e-02,           &
    -0.2159e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=15,15)/                                     &
    -0.2321e-09, -0.4350e-09, -0.7861e-09, -0.1347e-08, -0.2144e-08,           &
    -0.3120e-08, -0.4107e-08, -0.4892e-08, -0.5511e-08, -0.6164e-08,           &
    -0.7054e-08, -0.8811e-08, -0.1240e-07, -0.1862e-07, -0.2773e-07,           &
    -0.3957e-07, -0.5371e-07, -0.6939e-07, -0.8564e-07, -0.1006e-06,           &
    -0.1100e-06, -0.1066e-06, -0.8018e-07, -0.1228e-07,  0.1263e-06,           &
     0.3678e-06,  0.7022e-06,  0.1049e-05,  0.1411e-05,  0.2015e-05,           &
     0.3099e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=16,16)/                                     &
     0.2064e-05,  0.4104e-05,  0.8141e-05,  0.1607e-04,  0.3145e-04,           &
     0.6080e-04,  0.1153e-03,  0.2129e-03,  0.3802e-03,  0.6538e-03,           &
     0.1078e-02,  0.1699e-02,  0.2576e-02,  0.3805e-02,  0.5540e-02,           &
     0.8004e-02,  0.1149e-01,  0.1640e-01,  0.2329e-01,  0.3290e-01,           &
     0.4622e-01,  0.6441e-01,  0.8888e-01,  0.1210e+00,  0.1622e+00,           &
     0.2132e+00,  0.2745e+00,  0.3470e+00,  0.4317e+00,  0.5286e+00,           &
     0.6335e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=16,16)/                                     &
    -0.3122e-07, -0.6175e-07, -0.1214e-06, -0.2361e-06, -0.4518e-06,           &
    -0.8438e-06, -0.1524e-05, -0.2643e-05, -0.4380e-05, -0.6922e-05,           &
    -0.1042e-04, -0.1504e-04, -0.2125e-04, -0.2987e-04, -0.4200e-04,           &
    -0.5923e-04, -0.8383e-04, -0.1186e-03, -0.1670e-03, -0.2328e-03,           &
    -0.3204e-03, -0.4347e-03, -0.5802e-03, -0.7595e-03, -0.9703e-03,           &
    -0.1205e-02, -0.1457e-02, -0.1720e-02, -0.1980e-02, -0.2191e-02,           &
    -0.2290e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=16,16)/                                     &
    -0.1376e-09, -0.2699e-09, -0.5220e-09, -0.9897e-09, -0.1819e-08,           &
    -0.3186e-08, -0.5224e-08, -0.7896e-08, -0.1090e-07, -0.1349e-07,           &
    -0.1443e-07, -0.1374e-07, -0.1386e-07, -0.1673e-07, -0.2237e-07,           &
    -0.3248e-07, -0.5050e-07, -0.7743e-07, -0.1097e-06, -0.1369e-06,           &
    -0.1463e-06, -0.1268e-06, -0.6424e-07,  0.5941e-07,  0.2742e-06,           &
     0.5924e-06,  0.9445e-06,  0.1286e-05,  0.1819e-05,  0.2867e-05,           &
     0.4527e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=17,17)/                                     &
     0.2452e-05,  0.4881e-05,  0.9693e-05,  0.1918e-04,  0.3768e-04,           &
     0.7322e-04,  0.1399e-03,  0.2607e-03,  0.4694e-03,  0.8102e-03,           &
     0.1334e-02,  0.2098e-02,  0.3180e-02,  0.4710e-02,  0.6879e-02,           &
     0.9952e-02,  0.1429e-01,  0.2038e-01,  0.2894e-01,  0.4085e-01,           &
     0.5722e-01,  0.7939e-01,  0.1088e+00,  0.1468e+00,  0.1945e+00,           &
     0.2524e+00,  0.3210e+00,  0.4017e+00,  0.4948e+00,  0.5978e+00,           &
     0.7040e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=17,17)/                                     &
    -0.3547e-07, -0.7029e-07, -0.1386e-06, -0.2709e-06, -0.5218e-06,           &
    -0.9840e-06, -0.1799e-05, -0.3156e-05, -0.5272e-05, -0.8357e-05,           &
    -0.1260e-04, -0.1827e-04, -0.2598e-04, -0.3667e-04, -0.5169e-04,           &
    -0.7312e-04, -0.1037e-03, -0.1467e-03, -0.2060e-03, -0.2857e-03,           &
    -0.3907e-03, -0.5257e-03, -0.6940e-03, -0.8954e-03, -0.1124e-02,           &
    -0.1371e-02, -0.1632e-02, -0.1897e-02, -0.2131e-02, -0.2275e-02,           &
    -0.2265e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=17,17)/                                     &
    -0.1482e-09, -0.2910e-09, -0.5667e-09, -0.1081e-08, -0.2005e-08,           &
    -0.3554e-08, -0.5902e-08, -0.8925e-08, -0.1209e-07, -0.1448e-07,           &
    -0.1536e-07, -0.1565e-07, -0.1763e-07, -0.2088e-07, -0.2564e-07,           &
    -0.3635e-07, -0.5791e-07, -0.8907e-07, -0.1213e-06, -0.1418e-06,           &
    -0.1397e-06, -0.1000e-06, -0.4427e-08,  0.1713e-06,  0.4536e-06,           &
     0.8086e-06,  0.1153e-05,  0.1588e-05,  0.2437e-05,  0.3905e-05,           &
     0.5874e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=18,18)/                                     &
     0.2878e-05,  0.5733e-05,  0.1140e-04,  0.2258e-04,  0.4446e-04,           &
     0.8668e-04,  0.1664e-03,  0.3120e-03,  0.5656e-03,  0.9812e-03,           &
     0.1621e-02,  0.2559e-02,  0.3902e-02,  0.5812e-02,  0.8522e-02,           &
     0.1235e-01,  0.1773e-01,  0.2530e-01,  0.3590e-01,  0.5058e-01,           &
     0.7059e-01,  0.9736e-01,  0.1323e+00,  0.1767e+00,  0.2312e+00,           &
     0.2962e+00,  0.3728e+00,  0.4619e+00,  0.5623e+00,  0.6685e+00,           &
     0.7716e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=18,18)/                                     &
    -0.4064e-07, -0.8066e-07, -0.1593e-06, -0.3124e-06, -0.6049e-06,           &
    -0.1148e-05, -0.2118e-05, -0.3751e-05, -0.6314e-05, -0.1006e-04,           &
    -0.1526e-04, -0.2232e-04, -0.3196e-04, -0.4525e-04, -0.6394e-04,           &
    -0.9058e-04, -0.1284e-03, -0.1812e-03, -0.2533e-03, -0.3493e-03,           &
    -0.4740e-03, -0.6315e-03, -0.8228e-03, -0.1044e-02, -0.1286e-02,           &
    -0.1544e-02, -0.1810e-02, -0.2061e-02, -0.2243e-02, -0.2291e-02,           &
    -0.2152e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=18,18)/                                     &
    -0.1630e-09, -0.3213e-09, -0.6266e-09, -0.1201e-08, -0.2248e-08,           &
    -0.4030e-08, -0.6770e-08, -0.1033e-07, -0.1392e-07, -0.1640e-07,           &
    -0.1768e-07, -0.1932e-07, -0.2229e-07, -0.2508e-07, -0.2940e-07,           &
    -0.4200e-07, -0.6717e-07, -0.1002e-06, -0.1286e-06, -0.1402e-06,           &
    -0.1216e-06, -0.5487e-07,  0.8418e-07,  0.3246e-06,  0.6610e-06,           &
     0.1013e-05,  0.1394e-05,  0.2073e-05,  0.3337e-05,  0.5175e-05,           &
     0.7255e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=19,19)/                                     &
     0.3282e-05,  0.6539e-05,  0.1300e-04,  0.2578e-04,  0.5085e-04,           &
     0.9936e-04,  0.1914e-03,  0.3608e-03,  0.6586e-03,  0.1153e-02,           &
     0.1924e-02,  0.3074e-02,  0.4742e-02,  0.7126e-02,  0.1051e-01,           &
     0.1526e-01,  0.2195e-01,  0.3133e-01,  0.4441e-01,  0.6238e-01,           &
     0.8664e-01,  0.1186e+00,  0.1597e+00,  0.2108e+00,  0.2723e+00,           &
     0.3450e+00,  0.4300e+00,  0.5272e+00,  0.6324e+00,  0.7378e+00,           &
     0.8330e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=19,19)/                                     &
    -0.4629e-07, -0.9195e-07, -0.1819e-06, -0.3572e-06, -0.6936e-06,           &
    -0.1323e-05, -0.2456e-05, -0.4385e-05, -0.7453e-05, -0.1200e-04,           &
    -0.1843e-04, -0.2731e-04, -0.3943e-04, -0.5606e-04, -0.7936e-04,           &
    -0.1123e-03, -0.1588e-03, -0.2231e-03, -0.3101e-03, -0.4247e-03,           &
    -0.5713e-03, -0.7522e-03, -0.9651e-03, -0.1202e-02, -0.1456e-02,           &
    -0.1721e-02, -0.1983e-02, -0.2196e-02, -0.2296e-02, -0.2224e-02,           &
    -0.1952e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=19,19)/                                     &
    -0.1827e-09, -0.3607e-09, -0.7057e-09, -0.1359e-08, -0.2552e-08,           &
    -0.4615e-08, -0.7854e-08, -0.1218e-07, -0.1670e-07, -0.2008e-07,           &
    -0.2241e-07, -0.2516e-07, -0.2796e-07, -0.3015e-07, -0.3506e-07,           &
    -0.4958e-07, -0.7627e-07, -0.1070e-06, -0.1289e-06, -0.1286e-06,           &
    -0.8843e-07,  0.1492e-07,  0.2118e-06,  0.5155e-06,  0.8665e-06,           &
     0.1220e-05,  0.1765e-05,  0.2825e-05,  0.4498e-05,  0.6563e-05,           &
     0.8422e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=20,20)/                                     &
     0.3608e-05,  0.7189e-05,  0.1430e-04,  0.2837e-04,  0.5601e-04,           &
     0.1097e-03,  0.2120e-03,  0.4017e-03,  0.7396e-03,  0.1311e-02,           &
     0.2226e-02,  0.3624e-02,  0.5685e-02,  0.8652e-02,  0.1286e-01,           &
     0.1878e-01,  0.2709e-01,  0.3868e-01,  0.5471e-01,  0.7654e-01,           &
     0.1057e+00,  0.1435e+00,  0.1912e+00,  0.2492e+00,  0.3181e+00,           &
     0.3990e+00,  0.4926e+00,  0.5960e+00,  0.7027e+00,  0.8027e+00,           &
     0.8852e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=20,20)/                                     &
    -0.5164e-07, -0.1026e-06, -0.2031e-06, -0.3994e-06, -0.7771e-06,           &
    -0.1488e-05, -0.2776e-05, -0.5001e-05, -0.8610e-05, -0.1411e-04,           &
    -0.2209e-04, -0.3328e-04, -0.4860e-04, -0.6954e-04, -0.9861e-04,           &
    -0.1393e-03, -0.1961e-03, -0.2738e-03, -0.3778e-03, -0.5132e-03,           &
    -0.6831e-03, -0.8868e-03, -0.1118e-02, -0.1368e-02, -0.1632e-02,           &
    -0.1899e-02, -0.2136e-02, -0.2282e-02, -0.2273e-02, -0.2067e-02,           &
    -0.1679e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=20,20)/                                     &
    -0.2058e-09, -0.4066e-09, -0.7967e-09, -0.1539e-08, -0.2904e-08,           &
    -0.5293e-08, -0.9116e-08, -0.1447e-07, -0.2058e-07, -0.2608e-07,           &
    -0.3053e-07, -0.3418e-07, -0.3619e-07, -0.3766e-07, -0.4313e-07,           &
    -0.5817e-07, -0.8299e-07, -0.1072e-06, -0.1195e-06, -0.1031e-06,           &
    -0.3275e-07,  0.1215e-06,  0.3835e-06,  0.7220e-06,  0.1062e-05,           &
     0.1504e-05,  0.2367e-05,  0.3854e-05,  0.5842e-05,  0.7875e-05,           &
     0.9082e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=21,21)/                                     &
     0.3827e-05,  0.7626e-05,  0.1518e-04,  0.3013e-04,  0.5955e-04,           &
     0.1169e-03,  0.2268e-03,  0.4325e-03,  0.8046e-03,  0.1449e-02,           &
     0.2512e-02,  0.4185e-02,  0.6707e-02,  0.1038e-01,  0.1561e-01,           &
     0.2298e-01,  0.3328e-01,  0.4754e-01,  0.6708e-01,  0.9340e-01,           &
     0.1280e+00,  0.1722e+00,  0.2267e+00,  0.2918e+00,  0.3688e+00,           &
     0.4584e+00,  0.5594e+00,  0.6662e+00,  0.7700e+00,  0.8597e+00,           &
     0.9266e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=21,21)/                                     &
    -0.5584e-07, -0.1110e-06, -0.2198e-06, -0.4329e-06, -0.8444e-06,           &
    -0.1623e-05, -0.3049e-05, -0.5551e-05, -0.9714e-05, -0.1627e-04,           &
    -0.2609e-04, -0.4015e-04, -0.5955e-04, -0.8603e-04, -0.1223e-03,           &
    -0.1724e-03, -0.2413e-03, -0.3346e-03, -0.4578e-03, -0.6155e-03,           &
    -0.8087e-03, -0.1033e-02, -0.1279e-02, -0.1540e-02, -0.1811e-02,           &
    -0.2065e-02, -0.2251e-02, -0.2301e-02, -0.2163e-02, -0.1828e-02,           &
    -0.1365e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=21,21)/                                     &
    -0.2274e-09, -0.4498e-09, -0.8814e-09, -0.1708e-08, -0.3247e-08,           &
    -0.5972e-08, -0.1045e-07, -0.1707e-07, -0.2545e-07, -0.3440e-07,           &
    -0.4259e-07, -0.4822e-07, -0.5004e-07, -0.5061e-07, -0.5485e-07,           &
    -0.6687e-07, -0.8483e-07, -0.9896e-07, -0.9646e-07, -0.5557e-07,           &
     0.5765e-07,  0.2752e-06,  0.5870e-06,  0.9188e-06,  0.1291e-05,           &
     0.1971e-05,  0.3251e-05,  0.5115e-05,  0.7221e-05,  0.8825e-05,           &
     0.9032e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=22,22)/                                     &
     0.3938e-05,  0.7849e-05,  0.1563e-04,  0.3105e-04,  0.6149e-04,           &
     0.1210e-03,  0.2361e-03,  0.4537e-03,  0.8543e-03,  0.1565e-02,           &
     0.2774e-02,  0.4739e-02,  0.7781e-02,  0.1229e-01,  0.1877e-01,           &
     0.2791e-01,  0.4064e-01,  0.5812e-01,  0.8180e-01,  0.1132e+00,           &
     0.1539e+00,  0.2047e+00,  0.2662e+00,  0.3391e+00,  0.4247e+00,           &
     0.5226e+00,  0.6287e+00,  0.7351e+00,  0.8311e+00,  0.9066e+00,           &
     0.9569e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=22,22)/                                     &
    -0.5833e-07, -0.1160e-06, -0.2300e-06, -0.4540e-06, -0.8885e-06,           &
    -0.1718e-05, -0.3256e-05, -0.6010e-05, -0.1072e-04, -0.1838e-04,           &
    -0.3026e-04, -0.4772e-04, -0.7223e-04, -0.1057e-03, -0.1512e-03,           &
    -0.2128e-03, -0.2961e-03, -0.4070e-03, -0.5514e-03, -0.7320e-03,           &
    -0.9467e-03, -0.1187e-02, -0.1446e-02, -0.1717e-02, -0.1985e-02,           &
    -0.2203e-02, -0.2307e-02, -0.2237e-02, -0.1965e-02, -0.1532e-02,           &
    -0.1044e-02/
   data ((h23(ip,iw,1),iw=1,31),ip=22,22)/                                     &
    -0.2426e-09, -0.4805e-09, -0.9447e-09, -0.1841e-08, -0.3519e-08,           &
    -0.6565e-08, -0.1172e-07, -0.1979e-07, -0.3095e-07, -0.4443e-07,           &
    -0.5821e-07, -0.6868e-07, -0.7282e-07, -0.7208e-07, -0.7176e-07,           &
    -0.7562e-07, -0.8110e-07, -0.7934e-07, -0.5365e-07,  0.2483e-07,           &
     0.1959e-06,  0.4731e-06,  0.7954e-06,  0.1123e-05,  0.1652e-05,           &
     0.2711e-05,  0.4402e-05,  0.6498e-05,  0.8392e-05,  0.9154e-05,           &
     0.8261e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=23,23)/                                     &
     0.3973e-05,  0.7921e-05,  0.1578e-04,  0.3139e-04,  0.6229e-04,           &
     0.1231e-03,  0.2414e-03,  0.4680e-03,  0.8923e-03,  0.1663e-02,           &
     0.3011e-02,  0.5268e-02,  0.8869e-02,  0.1434e-01,  0.2230e-01,           &
     0.3359e-01,  0.4924e-01,  0.7057e-01,  0.9908e-01,  0.1363e+00,           &
     0.1834e+00,  0.2411e+00,  0.3100e+00,  0.3913e+00,  0.4856e+00,           &
     0.5901e+00,  0.6981e+00,  0.7994e+00,  0.8832e+00,  0.9424e+00,           &
     0.9773e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=23,23)/                                     &
    -0.5929e-07, -0.1180e-06, -0.2344e-06, -0.4638e-06, -0.9118e-06,           &
    -0.1775e-05, -0.3401e-05, -0.6375e-05, -0.1160e-04, -0.2039e-04,           &
    -0.3444e-04, -0.5575e-04, -0.8641e-04, -0.1287e-03, -0.1856e-03,           &
    -0.2615e-03, -0.3618e-03, -0.4928e-03, -0.6594e-03, -0.8621e-03,           &
    -0.1095e-02, -0.1349e-02, -0.1618e-02, -0.1894e-02, -0.2140e-02,           &
    -0.2293e-02, -0.2290e-02, -0.2085e-02, -0.1696e-02, -0.1212e-02,           &
    -0.7506e-03/
   data ((h23(ip,iw,1),iw=1,31),ip=23,23)/                                     &
    -0.2496e-09, -0.4954e-09, -0.9780e-09, -0.1915e-08, -0.3697e-08,           &
    -0.6991e-08, -0.1279e-07, -0.2231e-07, -0.3653e-07, -0.5541e-07,           &
    -0.7688e-07, -0.9614e-07, -0.1069e-06, -0.1065e-06, -0.9866e-07,           &
    -0.8740e-07, -0.7192e-07, -0.4304e-07,  0.1982e-07,  0.1525e-06,           &
     0.3873e-06,  0.6947e-06,  0.1000e-05,  0.1409e-05,  0.2253e-05,           &
     0.3739e-05,  0.5744e-05,  0.7812e-05,  0.9067e-05,  0.8746e-05,           &
     0.6940e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=24,24)/                                     &
     0.3978e-05,  0.7933e-05,  0.1581e-04,  0.3149e-04,  0.6260e-04,           &
     0.1240e-03,  0.2445e-03,  0.4777e-03,  0.9213e-03,  0.1743e-02,           &
     0.3217e-02,  0.5756e-02,  0.9929e-02,  0.1644e-01,  0.2612e-01,           &
     0.3995e-01,  0.5910e-01,  0.8497e-01,  0.1190e+00,  0.1626e+00,           &
     0.2165e+00,  0.2813e+00,  0.3583e+00,  0.4485e+00,  0.5505e+00,           &
     0.6590e+00,  0.7646e+00,  0.8562e+00,  0.9245e+00,  0.9676e+00,           &
     0.9897e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=24,24)/                                     &
    -0.5950e-07, -0.1185e-06, -0.2358e-06, -0.4678e-06, -0.9235e-06,           &
    -0.1810e-05, -0.3503e-05, -0.6664e-05, -0.1236e-04, -0.2223e-04,           &
    -0.3849e-04, -0.6396e-04, -0.1017e-03, -0.1545e-03, -0.2257e-03,           &
    -0.3192e-03, -0.4399e-03, -0.5933e-03, -0.7824e-03, -0.1005e-02,           &
    -0.1251e-02, -0.1516e-02, -0.1794e-02, -0.2060e-02, -0.2257e-02,           &
    -0.2318e-02, -0.2186e-02, -0.1853e-02, -0.1386e-02, -0.9021e-03,           &
    -0.5050e-03/
   data ((h23(ip,iw,1),iw=1,31),ip=24,24)/                                     &
    -0.2515e-09, -0.5001e-09, -0.9904e-09, -0.1951e-08, -0.3800e-08,           &
    -0.7288e-08, -0.1362e-07, -0.2452e-07, -0.4184e-07, -0.6663e-07,           &
    -0.9770e-07, -0.1299e-06, -0.1533e-06, -0.1584e-06, -0.1425e-06,           &
    -0.1093e-06, -0.5972e-07,  0.1426e-07,  0.1347e-06,  0.3364e-06,           &
     0.6209e-06,  0.9169e-06,  0.1243e-05,  0.1883e-05,  0.3145e-05,           &
     0.5011e-05,  0.7136e-05,  0.8785e-05,  0.9048e-05,  0.7686e-05,           &
     0.5368e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=25,25)/                                     &
     0.3978e-05,  0.7934e-05,  0.1582e-04,  0.3153e-04,  0.6275e-04,           &
     0.1246e-03,  0.2466e-03,  0.4847e-03,  0.9432e-03,  0.1808e-02,           &
     0.3392e-02,  0.6189e-02,  0.1092e-01,  0.1852e-01,  0.3009e-01,           &
     0.4683e-01,  0.7008e-01,  0.1012e+00,  0.1416e+00,  0.1920e+00,           &
     0.2531e+00,  0.3258e+00,  0.4115e+00,  0.5102e+00,  0.6180e+00,           &
     0.7267e+00,  0.8252e+00,  0.9030e+00,  0.9550e+00,  0.9838e+00,           &
     0.9965e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=25,25)/                                     &
    -0.5953e-07, -0.1187e-06, -0.2363e-06, -0.4697e-06, -0.9304e-06,           &
    -0.1833e-05, -0.3578e-05, -0.6889e-05, -0.1299e-04, -0.2384e-04,           &
    -0.4227e-04, -0.7201e-04, -0.1174e-03, -0.1825e-03, -0.2710e-03,           &
    -0.3861e-03, -0.5313e-03, -0.7095e-03, -0.9203e-03, -0.1158e-02,           &
    -0.1417e-02, -0.1691e-02, -0.1968e-02, -0.2200e-02, -0.2320e-02,           &
    -0.2263e-02, -0.1998e-02, -0.1563e-02, -0.1068e-02, -0.6300e-03,           &
    -0.3174e-03/
   data ((h23(ip,iw,1),iw=1,31),ip=25,25)/                                     &
    -0.2520e-09, -0.5016e-09, -0.9963e-09, -0.1971e-08, -0.3867e-08,           &
    -0.7500e-08, -0.1427e-07, -0.2634e-07, -0.4656e-07, -0.7753e-07,           &
    -0.1196e-06, -0.1683e-06, -0.2113e-06, -0.2309e-06, -0.2119e-06,           &
    -0.1518e-06, -0.5200e-07,  0.9274e-07,  0.2973e-06,  0.5714e-06,           &
     0.8687e-06,  0.1152e-05,  0.1621e-05,  0.2634e-05,  0.4310e-05,           &
     0.6418e-05,  0.8347e-05,  0.9162e-05,  0.8319e-05,  0.6209e-05,           &
     0.3844e-05/
   data ((h21(ip,iw,1),iw=1,31),ip=26,26)/                                     &
     0.3976e-05,  0.7932e-05,  0.1582e-04,  0.3154e-04,  0.6284e-04,           &
     0.1250e-03,  0.2479e-03,  0.4896e-03,  0.9592e-03,  0.1857e-02,           &
     0.3532e-02,  0.6558e-02,  0.1181e-01,  0.2048e-01,  0.3402e-01,           &
     0.5399e-01,  0.8188e-01,  0.1190e+00,  0.1664e+00,  0.2243e+00,           &
     0.2933e+00,  0.3747e+00,  0.4695e+00,  0.5755e+00,  0.6858e+00,           &
     0.7903e+00,  0.8773e+00,  0.9391e+00,  0.9758e+00,  0.9934e+00,           &
     0.9996e+00/
   data ((h22(ip,iw,1),iw=1,31),ip=26,26)/                                     &
    -0.5954e-07, -0.1187e-06, -0.2366e-06, -0.4709e-06, -0.9349e-06,           &
    -0.1849e-05, -0.3632e-05, -0.7058e-05, -0.1350e-04, -0.2521e-04,           &
    -0.4564e-04, -0.7958e-04, -0.1329e-03, -0.2114e-03, -0.3200e-03,           &
    -0.4611e-03, -0.6353e-03, -0.8409e-03, -0.1072e-02, -0.1324e-02,           &
    -0.1592e-02, -0.1871e-02, -0.2127e-02, -0.2297e-02, -0.2312e-02,           &
    -0.2123e-02, -0.1738e-02, -0.1247e-02, -0.7744e-03, -0.4117e-03,           &
    -0.1850e-03/
   data ((h23(ip,iw,1),iw=1,31),ip=26,26)/                                     &
    -0.2522e-09, -0.5025e-09, -0.9997e-09, -0.1983e-08, -0.3912e-08,           &
    -0.7650e-08, -0.1474e-07, -0.2777e-07, -0.5055e-07, -0.8745e-07,           &
    -0.1414e-06, -0.2095e-06, -0.2790e-06, -0.3241e-06, -0.3135e-06,           &
    -0.2269e-06, -0.5896e-07,  0.1875e-06,  0.4996e-06,  0.8299e-06,           &
     0.1115e-05,  0.1467e-05,  0.2236e-05,  0.3672e-05,  0.5668e-05,           &
     0.7772e-05,  0.9094e-05,  0.8827e-05,  0.7041e-05,  0.4638e-05,           &
     0.2539e-05/
   data ((h81(ip,iw,1),iw=1,31),ip= 1, 1)/                                     &
     0.1342e-04,  0.2640e-04,  0.5137e-04,  0.9829e-04,  0.1832e-03,           &
     0.3284e-03,  0.5585e-03,  0.8907e-03,  0.1327e-02,  0.1856e-02,           &
     0.2461e-02,  0.3138e-02,  0.3898e-02,  0.4766e-02,  0.5783e-02,           &
     0.7007e-02,  0.8521e-02,  0.1041e-01,  0.1278e-01,  0.1570e-01,           &
     0.1928e-01,  0.2361e-01,  0.2884e-01,  0.3520e-01,  0.4305e-01,           &
     0.5286e-01,  0.6531e-01,  0.8127e-01,  0.1019e+00,  0.1287e+00,           &
     0.1632e+00/
   data ((h82(ip,iw,1),iw=1,31),ip= 1, 1)/                                     &
    -0.5685e-08, -0.1331e-07, -0.3249e-07, -0.8137e-07, -0.2048e-06,           &
    -0.4973e-06, -0.1118e-05, -0.2246e-05, -0.3982e-05, -0.6290e-05,           &
    -0.9040e-05, -0.1215e-04, -0.1567e-04, -0.1970e-04, -0.2449e-04,           &
    -0.3046e-04, -0.3798e-04, -0.4725e-04, -0.5831e-04, -0.7123e-04,           &
    -0.8605e-04, -0.1028e-03, -0.1212e-03, -0.1413e-03, -0.1635e-03,           &
    -0.1884e-03, -0.2160e-03, -0.2461e-03, -0.2778e-03, -0.3098e-03,           &
    -0.3411e-03/
   data ((h83(ip,iw,1),iw=1,31),ip= 1, 1)/                                     &
     0.2169e-10,  0.5237e-10,  0.1296e-09,  0.3204e-09,  0.7665e-09,           &
     0.1691e-08,  0.3222e-08,  0.5110e-08,  0.6779e-08,  0.7681e-08,           &
     0.7378e-08,  0.5836e-08,  0.3191e-08, -0.1491e-08, -0.1022e-07,           &
    -0.2359e-07, -0.3957e-07, -0.5553e-07, -0.6927e-07, -0.7849e-07,           &
    -0.8139e-07, -0.7853e-07, -0.7368e-07, -0.7220e-07, -0.7780e-07,           &
    -0.9091e-07, -0.1038e-06, -0.9929e-07, -0.5422e-07,  0.5379e-07,           &
     0.2350e-06/
   data ((h81(ip,iw,1),iw=1,31),ip= 2, 2)/                                     &
     0.1342e-04,  0.2640e-04,  0.5137e-04,  0.9829e-04,  0.1832e-03,           &
     0.3284e-03,  0.5586e-03,  0.8908e-03,  0.1327e-02,  0.1856e-02,           &
     0.2462e-02,  0.3140e-02,  0.3902e-02,  0.4774e-02,  0.5797e-02,           &
     0.7033e-02,  0.8570e-02,  0.1050e-01,  0.1294e-01,  0.1599e-01,           &
     0.1979e-01,  0.2448e-01,  0.3024e-01,  0.3738e-01,  0.4633e-01,           &
     0.5766e-01,  0.7218e-01,  0.9096e-01,  0.1154e+00,  0.1471e+00,           &
     0.1880e+00/
   data ((h82(ip,iw,1),iw=1,31),ip= 2, 2)/                                     &
    -0.5684e-08, -0.1331e-07, -0.3248e-07, -0.8133e-07, -0.2047e-06,           &
    -0.4971e-06, -0.1117e-05, -0.2245e-05, -0.3981e-05, -0.6287e-05,           &
    -0.9035e-05, -0.1215e-04, -0.1565e-04, -0.1967e-04, -0.2444e-04,           &
    -0.3036e-04, -0.3780e-04, -0.4694e-04, -0.5779e-04, -0.7042e-04,           &
    -0.8491e-04, -0.1013e-03, -0.1196e-03, -0.1399e-03, -0.1625e-03,           &
    -0.1879e-03, -0.2163e-03, -0.2474e-03, -0.2803e-03, -0.3140e-03,           &
    -0.3478e-03/
   data ((h83(ip,iw,1),iw=1,31),ip= 2, 2)/                                     &
     0.2168e-10,  0.5242e-10,  0.1295e-09,  0.3201e-09,  0.7662e-09,           &
     0.1690e-08,  0.3220e-08,  0.5106e-08,  0.6776e-08,  0.7673e-08,           &
     0.7362e-08,  0.5808e-08,  0.3138e-08, -0.1595e-08, -0.1041e-07,           &
    -0.2390e-07, -0.4010e-07, -0.5636e-07, -0.7045e-07, -0.7972e-07,           &
    -0.8178e-07, -0.7677e-07, -0.6876e-07, -0.6381e-07, -0.6583e-07,           &
    -0.7486e-07, -0.8229e-07, -0.7017e-07, -0.1497e-07,  0.1051e-06,           &
     0.2990e-06/
   data ((h81(ip,iw,1),iw=1,31),ip= 3, 3)/                                     &
     0.1342e-04,  0.2640e-04,  0.5137e-04,  0.9830e-04,  0.1832e-03,           &
     0.3285e-03,  0.5587e-03,  0.8911e-03,  0.1328e-02,  0.1857e-02,           &
     0.2465e-02,  0.3144e-02,  0.3909e-02,  0.4787e-02,  0.5820e-02,           &
     0.7075e-02,  0.8646e-02,  0.1064e-01,  0.1319e-01,  0.1643e-01,           &
     0.2053e-01,  0.2569e-01,  0.3215e-01,  0.4028e-01,  0.5059e-01,           &
     0.6380e-01,  0.8088e-01,  0.1031e+00,  0.1321e+00,  0.1698e+00,           &
     0.2179e+00/
   data ((h82(ip,iw,1),iw=1,31),ip= 3, 3)/                                     &
    -0.5682e-08, -0.1330e-07, -0.3247e-07, -0.8129e-07, -0.2046e-06,           &
    -0.4968e-06, -0.1117e-05, -0.2244e-05, -0.3978e-05, -0.6283e-05,           &
    -0.9027e-05, -0.1213e-04, -0.1563e-04, -0.1963e-04, -0.2436e-04,           &
    -0.3021e-04, -0.3754e-04, -0.4649e-04, -0.5709e-04, -0.6940e-04,           &
    -0.8359e-04, -0.9986e-04, -0.1182e-03, -0.1388e-03, -0.1620e-03,           &
    -0.1882e-03, -0.2175e-03, -0.2498e-03, -0.2843e-03, -0.3203e-03,           &
    -0.3573e-03/
   data ((h83(ip,iw,1),iw=1,31),ip= 3, 3)/                                     &
     0.2167e-10,  0.5238e-10,  0.1294e-09,  0.3198e-09,  0.7656e-09,           &
     0.1688e-08,  0.3217e-08,  0.5104e-08,  0.6767e-08,  0.7661e-08,           &
     0.7337e-08,  0.5764e-08,  0.3051e-08, -0.1752e-08, -0.1068e-07,           &
    -0.2436e-07, -0.4081e-07, -0.5740e-07, -0.7165e-07, -0.8046e-07,           &
    -0.8082e-07, -0.7289e-07, -0.6141e-07, -0.5294e-07, -0.5134e-07,           &
    -0.5552e-07, -0.5609e-07, -0.3464e-07,  0.3275e-07,  0.1669e-06,           &
     0.3745e-06/
   data ((h81(ip,iw,1),iw=1,31),ip= 4, 4)/                                     &
     0.1342e-04,  0.2640e-04,  0.5138e-04,  0.9830e-04,  0.1832e-03,           &
     0.3286e-03,  0.5589e-03,  0.8915e-03,  0.1329e-02,  0.1859e-02,           &
     0.2468e-02,  0.3150e-02,  0.3920e-02,  0.4806e-02,  0.5855e-02,           &
     0.7140e-02,  0.8763e-02,  0.1085e-01,  0.1356e-01,  0.1707e-01,           &
     0.2158e-01,  0.2736e-01,  0.3471e-01,  0.4407e-01,  0.5608e-01,           &
     0.7161e-01,  0.9184e-01,  0.1183e+00,  0.1528e+00,  0.1973e+00,           &
     0.2537e+00/
   data ((h82(ip,iw,1),iw=1,31),ip= 4, 4)/                                     &
    -0.5680e-08, -0.1329e-07, -0.3243e-07, -0.8121e-07, -0.2044e-06,           &
    -0.4963e-06, -0.1115e-05, -0.2242e-05, -0.3974e-05, -0.6276e-05,           &
    -0.9015e-05, -0.1211e-04, -0.1559e-04, -0.1956e-04, -0.2423e-04,           &
    -0.2999e-04, -0.3716e-04, -0.4588e-04, -0.5618e-04, -0.6818e-04,           &
    -0.8218e-04, -0.9847e-04, -0.1171e-03, -0.1382e-03, -0.1621e-03,           &
    -0.1892e-03, -0.2197e-03, -0.2535e-03, -0.2902e-03, -0.3293e-03,           &
    -0.3700e-03/
   data ((h83(ip,iw,1),iw=1,31),ip= 4, 4)/                                     &
     0.2166e-10,  0.5229e-10,  0.1294e-09,  0.3193e-09,  0.7644e-09,           &
     0.1686e-08,  0.3213e-08,  0.5092e-08,  0.6753e-08,  0.7640e-08,           &
     0.7302e-08,  0.5696e-08,  0.2917e-08, -0.1984e-08, -0.1108e-07,           &
    -0.2497e-07, -0.4171e-07, -0.5849e-07, -0.7254e-07, -0.8017e-07,           &
    -0.7802e-07, -0.6662e-07, -0.5153e-07, -0.3961e-07, -0.3387e-07,           &
    -0.3219e-07, -0.2426e-07,  0.8700e-08,  0.9027e-07,  0.2400e-06,           &
     0.4623e-06/
   data ((h81(ip,iw,1),iw=1,31),ip= 5, 5)/                                     &
     0.1342e-04,  0.2640e-04,  0.5138e-04,  0.9832e-04,  0.1833e-03,           &
     0.3287e-03,  0.5591e-03,  0.8922e-03,  0.1330e-02,  0.1862e-02,           &
     0.2473e-02,  0.3159e-02,  0.3937e-02,  0.4837e-02,  0.5911e-02,           &
     0.7240e-02,  0.8943e-02,  0.1117e-01,  0.1411e-01,  0.1798e-01,           &
     0.2304e-01,  0.2961e-01,  0.3807e-01,  0.4896e-01,  0.6308e-01,           &
     0.8148e-01,  0.1056e+00,  0.1371e+00,  0.1780e+00,  0.2303e+00,           &
     0.2958e+00/
   data ((h82(ip,iw,1),iw=1,31),ip= 5, 5)/                                     &
    -0.5675e-08, -0.1328e-07, -0.3239e-07, -0.8110e-07, -0.2040e-06,           &
    -0.4954e-06, -0.1114e-05, -0.2238e-05, -0.3968e-05, -0.6265e-05,           &
    -0.8996e-05, -0.1208e-04, -0.1553e-04, -0.1945e-04, -0.2404e-04,           &
    -0.2966e-04, -0.3663e-04, -0.4508e-04, -0.5508e-04, -0.6686e-04,           &
    -0.8082e-04, -0.9732e-04, -0.1165e-03, -0.1382e-03, -0.1630e-03,           &
    -0.1913e-03, -0.2234e-03, -0.2593e-03, -0.2989e-03, -0.3417e-03,           &
    -0.3857e-03/
   data ((h83(ip,iw,1),iw=1,31),ip= 5, 5)/                                     &
     0.2163e-10,  0.5209e-10,  0.1291e-09,  0.3186e-09,  0.7626e-09,           &
     0.1682e-08,  0.3203e-08,  0.5078e-08,  0.6730e-08,  0.7606e-08,           &
     0.7246e-08,  0.5592e-08,  0.2735e-08, -0.2325e-08, -0.1162e-07,           &
    -0.2576e-07, -0.4268e-07, -0.5938e-07, -0.7262e-07, -0.7827e-07,           &
    -0.7297e-07, -0.5786e-07, -0.3930e-07, -0.2373e-07, -0.1295e-07,           &
    -0.3728e-08,  0.1465e-07,  0.6114e-07,  0.1590e-06,  0.3257e-06,           &
     0.5622e-06/
   data ((h81(ip,iw,1),iw=1,31),ip= 6, 6)/                                     &
     0.1342e-04,  0.2640e-04,  0.5139e-04,  0.9834e-04,  0.1833e-03,           &
     0.3288e-03,  0.5596e-03,  0.8932e-03,  0.1332e-02,  0.1866e-02,           &
     0.2481e-02,  0.3174e-02,  0.3963e-02,  0.4885e-02,  0.5997e-02,           &
     0.7394e-02,  0.9215e-02,  0.1164e-01,  0.1489e-01,  0.1924e-01,           &
     0.2501e-01,  0.3258e-01,  0.4243e-01,  0.5523e-01,  0.7195e-01,           &
     0.9387e-01,  0.1226e+00,  0.1601e+00,  0.2084e+00,  0.2695e+00,           &
     0.3446e+00/
   data ((h82(ip,iw,1),iw=1,31),ip= 6, 6)/                                     &
    -0.5671e-08, -0.1326e-07, -0.3234e-07, -0.8091e-07, -0.2035e-06,           &
    -0.4941e-06, -0.1111e-05, -0.2232e-05, -0.3958e-05, -0.6247e-05,           &
    -0.8966e-05, -0.1202e-04, -0.1544e-04, -0.1929e-04, -0.2377e-04,           &
    -0.2921e-04, -0.3593e-04, -0.4409e-04, -0.5385e-04, -0.6555e-04,           &
    -0.7965e-04, -0.9656e-04, -0.1163e-03, -0.1390e-03, -0.1649e-03,           &
    -0.1947e-03, -0.2288e-03, -0.2675e-03, -0.3109e-03, -0.3575e-03,           &
    -0.4039e-03/
   data ((h83(ip,iw,1),iw=1,31),ip= 6, 6)/                                     &
     0.2155e-10,  0.5188e-10,  0.1288e-09,  0.3175e-09,  0.7599e-09,           &
     0.1675e-08,  0.3190e-08,  0.5059e-08,  0.6699e-08,  0.7551e-08,           &
     0.7154e-08,  0.5435e-08,  0.2452e-08, -0.2802e-08, -0.1235e-07,           &
    -0.2668e-07, -0.4353e-07, -0.5962e-07, -0.7134e-07, -0.7435e-07,           &
    -0.6551e-07, -0.4676e-07, -0.2475e-07, -0.4876e-08,  0.1235e-07,           &
     0.3092e-07,  0.6192e-07,  0.1243e-06,  0.2400e-06,  0.4247e-06,           &
     0.6755e-06/
   data ((h81(ip,iw,1),iw=1,31),ip= 7, 7)/                                     &
     0.1342e-04,  0.2640e-04,  0.5139e-04,  0.9837e-04,  0.1834e-03,           &
     0.3291e-03,  0.5603e-03,  0.8949e-03,  0.1336e-02,  0.1873e-02,           &
     0.2494e-02,  0.3198e-02,  0.4005e-02,  0.4960e-02,  0.6130e-02,           &
     0.7628e-02,  0.9618e-02,  0.1232e-01,  0.1599e-01,  0.2097e-01,           &
     0.2764e-01,  0.3646e-01,  0.4803e-01,  0.6319e-01,  0.8312e-01,           &
     0.1093e+00,  0.1436e+00,  0.1881e+00,  0.2448e+00,  0.3153e+00,           &
     0.4001e+00/
   data ((h82(ip,iw,1),iw=1,31),ip= 7, 7)/                                     &
    -0.5665e-08, -0.1322e-07, -0.3224e-07, -0.8063e-07, -0.2027e-06,           &
    -0.4921e-06, -0.1106e-05, -0.2223e-05, -0.3942e-05, -0.6220e-05,           &
    -0.8920e-05, -0.1194e-04, -0.1530e-04, -0.1905e-04, -0.2337e-04,           &
    -0.2860e-04, -0.3505e-04, -0.4296e-04, -0.5259e-04, -0.6439e-04,           &
    -0.7884e-04, -0.9635e-04, -0.1170e-03, -0.1407e-03, -0.1681e-03,           &
    -0.1998e-03, -0.2366e-03, -0.2790e-03, -0.3265e-03, -0.3763e-03,           &
    -0.4235e-03/
   data ((h83(ip,iw,1),iw=1,31),ip= 7, 7)/                                     &
     0.2157e-10,  0.5178e-10,  0.1283e-09,  0.3162e-09,  0.7558e-09,           &
     0.1665e-08,  0.3169e-08,  0.5027e-08,  0.6645e-08,  0.7472e-08,           &
     0.7017e-08,  0.5212e-08,  0.2059e-08, -0.3443e-08, -0.1321e-07,           &
    -0.2754e-07, -0.4389e-07, -0.5869e-07, -0.6825e-07, -0.6819e-07,           &
    -0.5570e-07, -0.3343e-07, -0.7592e-08,  0.1778e-07,  0.4322e-07,           &
     0.7330e-07,  0.1190e-06,  0.1993e-06,  0.3348e-06,  0.5376e-06,           &
     0.8030e-06/
   data ((h81(ip,iw,1),iw=1,31),ip= 8, 8)/                                     &
     0.1342e-04,  0.2641e-04,  0.5141e-04,  0.9842e-04,  0.1836e-03,           &
     0.3295e-03,  0.5614e-03,  0.8975e-03,  0.1341e-02,  0.1884e-02,           &
     0.2514e-02,  0.3234e-02,  0.4071e-02,  0.5075e-02,  0.6332e-02,           &
     0.7976e-02,  0.1021e-01,  0.1328e-01,  0.1751e-01,  0.2329e-01,           &
     0.3109e-01,  0.4146e-01,  0.5517e-01,  0.7325e-01,  0.9710e-01,           &
     0.1284e+00,  0.1693e+00,  0.2217e+00,  0.2876e+00,  0.3679e+00,           &
     0.4616e+00/
   data ((h82(ip,iw,1),iw=1,31),ip= 8, 8)/                                     &
    -0.5652e-08, -0.1318e-07, -0.3210e-07, -0.8018e-07, -0.2014e-06,           &
    -0.4888e-06, -0.1099e-05, -0.2210e-05, -0.3918e-05, -0.6179e-05,           &
    -0.8849e-05, -0.1182e-04, -0.1509e-04, -0.1871e-04, -0.2284e-04,           &
    -0.2782e-04, -0.3403e-04, -0.4177e-04, -0.5145e-04, -0.6354e-04,           &
    -0.7853e-04, -0.9681e-04, -0.1185e-03, -0.1437e-03, -0.1729e-03,           &
    -0.2072e-03, -0.2475e-03, -0.2942e-03, -0.3457e-03, -0.3973e-03,           &
    -0.4434e-03/
   data ((h83(ip,iw,1),iw=1,31),ip= 8, 8)/                                     &
     0.2153e-10,  0.5151e-10,  0.1273e-09,  0.3136e-09,  0.7488e-09,           &
     0.1649e-08,  0.3142e-08,  0.4980e-08,  0.6559e-08,  0.7346e-08,           &
     0.6813e-08,  0.4884e-08,  0.1533e-08, -0.4209e-08, -0.1409e-07,           &
    -0.2801e-07, -0.4320e-07, -0.5614e-07, -0.6312e-07, -0.5976e-07,           &
    -0.4369e-07, -0.1775e-07,  0.1280e-07,  0.4534e-07,  0.8106e-07,           &
     0.1246e-06,  0.1874e-06,  0.2873e-06,  0.4433e-06,  0.6651e-06,           &
     0.9477e-06/
   data ((h81(ip,iw,1),iw=1,31),ip= 9, 9)/                                     &
     0.1342e-04,  0.2641e-04,  0.5143e-04,  0.9849e-04,  0.1838e-03,           &
     0.3302e-03,  0.5631e-03,  0.9016e-03,  0.1350e-02,  0.1901e-02,           &
     0.2546e-02,  0.3292e-02,  0.4171e-02,  0.5250e-02,  0.6634e-02,           &
     0.8486e-02,  0.1104e-01,  0.1460e-01,  0.1955e-01,  0.2635e-01,           &
     0.3555e-01,  0.4786e-01,  0.6421e-01,  0.8588e-01,  0.1145e+00,           &
     0.1519e+00,  0.2003e+00,  0.2615e+00,  0.3372e+00,  0.4269e+00,           &
     0.5281e+00/
   data ((h82(ip,iw,1),iw=1,31),ip= 9, 9)/                                     &
    -0.5629e-08, -0.1310e-07, -0.3186e-07, -0.7948e-07, -0.1995e-06,           &
    -0.4837e-06, -0.1088e-05, -0.2188e-05, -0.3880e-05, -0.6115e-05,           &
    -0.8743e-05, -0.1165e-04, -0.1480e-04, -0.1824e-04, -0.2216e-04,           &
    -0.2691e-04, -0.3293e-04, -0.4067e-04, -0.5057e-04, -0.6314e-04,           &
    -0.7885e-04, -0.9813e-04, -0.1212e-03, -0.1482e-03, -0.1799e-03,           &
    -0.2175e-03, -0.2622e-03, -0.3135e-03, -0.3678e-03, -0.4193e-03,           &
    -0.4627e-03/
   data ((h83(ip,iw,1),iw=1,31),ip= 9, 9)/                                     &
     0.2121e-10,  0.5076e-10,  0.1257e-09,  0.3091e-09,  0.7379e-09,           &
     0.1623e-08,  0.3097e-08,  0.4904e-08,  0.6453e-08,  0.7168e-08,           &
     0.6534e-08,  0.4458e-08,  0.8932e-09, -0.5026e-08, -0.1469e-07,           &
    -0.2765e-07, -0.4103e-07, -0.5169e-07, -0.5585e-07, -0.4913e-07,           &
    -0.2954e-07,  0.6372e-09,  0.3738e-07,  0.7896e-07,  0.1272e-06,           &
     0.1867e-06,  0.2682e-06,  0.3895e-06,  0.5672e-06,  0.8091e-06,           &
     0.1114e-05/
   data ((h81(ip,iw,1),iw=1,31),ip=10,10)/                                     &
     0.1343e-04,  0.2642e-04,  0.5147e-04,  0.9861e-04,  0.1842e-03,           &
     0.3312e-03,  0.5658e-03,  0.9080e-03,  0.1364e-02,  0.1928e-02,           &
     0.2595e-02,  0.3379e-02,  0.4325e-02,  0.5513e-02,  0.7077e-02,           &
     0.9215e-02,  0.1220e-01,  0.1640e-01,  0.2226e-01,  0.3032e-01,           &
     0.4127e-01,  0.5598e-01,  0.7559e-01,  0.1016e+00,  0.1359e+00,           &
     0.1804e+00,  0.2372e+00,  0.3081e+00,  0.3935e+00,  0.4916e+00,           &
     0.5975e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=10,10)/                                     &
    -0.5597e-08, -0.1300e-07, -0.3148e-07, -0.7838e-07, -0.1964e-06,           &
    -0.4759e-06, -0.1071e-05, -0.2155e-05, -0.3822e-05, -0.6019e-05,           &
    -0.8586e-05, -0.1139e-04, -0.1439e-04, -0.1764e-04, -0.2134e-04,           &
    -0.2591e-04, -0.3188e-04, -0.3978e-04, -0.5011e-04, -0.6334e-04,           &
    -0.7998e-04, -0.1006e-03, -0.1253e-03, -0.1547e-03, -0.1895e-03,           &
    -0.2315e-03, -0.2811e-03, -0.3363e-03, -0.3917e-03, -0.4413e-03,           &
    -0.4809e-03/
   data ((h83(ip,iw,1),iw=1,31),ip=10,10)/                                     &
     0.2109e-10,  0.5017e-10,  0.1235e-09,  0.3021e-09,  0.7217e-09,           &
     0.1585e-08,  0.3028e-08,  0.4796e-08,  0.6285e-08,  0.6910e-08,           &
     0.6178e-08,  0.3945e-08,  0.2436e-09, -0.5632e-08, -0.1464e-07,           &
    -0.2596e-07, -0.3707e-07, -0.4527e-07, -0.4651e-07, -0.3644e-07,           &
    -0.1296e-07,  0.2250e-07,  0.6722e-07,  0.1202e-06,  0.1831e-06,           &
     0.2605e-06,  0.3627e-06,  0.5062e-06,  0.7064e-06,  0.9725e-06,           &
     0.1304e-05/
   data ((h81(ip,iw,1),iw=1,31),ip=11,11)/                                     &
     0.1343e-04,  0.2644e-04,  0.5152e-04,  0.9878e-04,  0.1847e-03,           &
     0.3328e-03,  0.5701e-03,  0.9179e-03,  0.1385e-02,  0.1969e-02,           &
     0.2671e-02,  0.3513e-02,  0.4554e-02,  0.5900e-02,  0.7713e-02,           &
     0.1023e-01,  0.1379e-01,  0.1879e-01,  0.2578e-01,  0.3542e-01,           &
     0.4854e-01,  0.6622e-01,  0.8983e-01,  0.1211e+00,  0.1619e+00,           &
     0.2145e+00,  0.2807e+00,  0.3614e+00,  0.4559e+00,  0.5603e+00,           &
     0.6675e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=11,11)/                                     &
    -0.5538e-08, -0.1280e-07, -0.3089e-07, -0.7667e-07, -0.1917e-06,           &
    -0.4642e-06, -0.1045e-05, -0.2106e-05, -0.3736e-05, -0.5878e-05,           &
    -0.8363e-05, -0.1104e-04, -0.1387e-04, -0.1692e-04, -0.2044e-04,           &
    -0.2493e-04, -0.3101e-04, -0.3926e-04, -0.5020e-04, -0.6429e-04,           &
    -0.8213e-04, -0.1044e-03, -0.1314e-03, -0.1637e-03, -0.2027e-03,           &
    -0.2498e-03, -0.3042e-03, -0.3617e-03, -0.4163e-03, -0.4625e-03,           &
    -0.4969e-03/
   data ((h83(ip,iw,1),iw=1,31),ip=11,11)/                                     &
     0.2067e-10,  0.4903e-10,  0.1200e-09,  0.2917e-09,  0.6965e-09,           &
     0.1532e-08,  0.2925e-08,  0.4632e-08,  0.6054e-08,  0.6590e-08,           &
     0.5746e-08,  0.3436e-08, -0.2251e-09, -0.5703e-08, -0.1344e-07,           &
    -0.2256e-07, -0.3120e-07, -0.3690e-07, -0.3520e-07, -0.2164e-07,           &
     0.6510e-08,  0.4895e-07,  0.1037e-06,  0.1702e-06,  0.2502e-06,           &
     0.3472e-06,  0.4710e-06,  0.6379e-06,  0.8633e-06,  0.1159e-05,           &
     0.1514e-05/
   data ((h81(ip,iw,1),iw=1,31),ip=12,12)/                                     &
     0.1344e-04,  0.2647e-04,  0.5161e-04,  0.9906e-04,  0.1856e-03,           &
     0.3353e-03,  0.5765e-03,  0.9331e-03,  0.1417e-02,  0.2032e-02,           &
     0.2785e-02,  0.3712e-02,  0.4892e-02,  0.6456e-02,  0.8604e-02,           &
     0.1162e-01,  0.1590e-01,  0.2192e-01,  0.3033e-01,  0.4193e-01,           &
     0.5774e-01,  0.7906e-01,  0.1075e+00,  0.1449e+00,  0.1934e+00,           &
     0.2549e+00,  0.3309e+00,  0.4213e+00,  0.5232e+00,  0.6308e+00,           &
     0.7348e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=12,12)/                                     &
    -0.5476e-08, -0.1257e-07, -0.3008e-07, -0.7418e-07, -0.1848e-06,           &
    -0.4468e-06, -0.1006e-05, -0.2032e-05, -0.3611e-05, -0.5679e-05,           &
    -0.8058e-05, -0.1059e-04, -0.1324e-04, -0.1612e-04, -0.1956e-04,           &
    -0.2411e-04, -0.3046e-04, -0.3925e-04, -0.5098e-04, -0.6619e-04,           &
    -0.8562e-04, -0.1100e-03, -0.1399e-03, -0.1761e-03, -0.2202e-03,           &
    -0.2726e-03, -0.3306e-03, -0.3885e-03, -0.4404e-03, -0.4820e-03,           &
    -0.5082e-03/
   data ((h83(ip,iw,1),iw=1,31),ip=12,12)/                                     &
     0.2041e-10,  0.4771e-10,  0.1149e-09,  0.2782e-09,  0.6614e-09,           &
     0.1451e-08,  0.2778e-08,  0.4401e-08,  0.5736e-08,  0.6189e-08,           &
     0.5315e-08,  0.3087e-08, -0.2518e-09, -0.4806e-08, -0.1071e-07,           &
    -0.1731e-07, -0.2346e-07, -0.2659e-07, -0.2184e-07, -0.4261e-08,           &
     0.2975e-07,  0.8112e-07,  0.1484e-06,  0.2308e-06,  0.3296e-06,           &
     0.4475e-06,  0.5942e-06,  0.7859e-06,  0.1041e-05,  0.1369e-05,           &
     0.1726e-05/
   data ((h81(ip,iw,1),iw=1,31),ip=13,13)/                                     &
     0.1345e-04,  0.2650e-04,  0.5173e-04,  0.9947e-04,  0.1868e-03,           &
     0.3389e-03,  0.5861e-03,  0.9559e-03,  0.1466e-02,  0.2125e-02,           &
     0.2955e-02,  0.4004e-02,  0.5376e-02,  0.7235e-02,  0.9824e-02,           &
     0.1348e-01,  0.1866e-01,  0.2596e-01,  0.3613e-01,  0.5016e-01,           &
     0.6929e-01,  0.9504e-01,  0.1292e+00,  0.1738e+00,  0.2309e+00,           &
     0.3021e+00,  0.3880e+00,  0.4868e+00,  0.5936e+00,  0.7003e+00,           &
     0.7964e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=13,13)/                                     &
    -0.5362e-08, -0.1223e-07, -0.2895e-07, -0.7071e-07, -0.1748e-06,           &
    -0.4219e-06, -0.9516e-06, -0.1928e-05, -0.3436e-05, -0.5409e-05,           &
    -0.7666e-05, -0.1005e-04, -0.1254e-04, -0.1533e-04, -0.1880e-04,           &
    -0.2358e-04, -0.3038e-04, -0.3988e-04, -0.5264e-04, -0.6934e-04,           &
    -0.9083e-04, -0.1179e-03, -0.1515e-03, -0.1927e-03, -0.2424e-03,           &
    -0.2994e-03, -0.3591e-03, -0.4155e-03, -0.4634e-03, -0.4982e-03,           &
    -0.5096e-03/
   data ((h83(ip,iw,1),iw=1,31),ip=13,13)/                                     &
     0.1976e-10,  0.4551e-10,  0.1086e-09,  0.2601e-09,  0.6126e-09,           &
     0.1345e-08,  0.2583e-08,  0.4112e-08,  0.5365e-08,  0.5796e-08,           &
     0.5031e-08,  0.3182e-08,  0.5970e-09, -0.2547e-08, -0.6172e-08,           &
    -0.1017e-07, -0.1388e-07, -0.1430e-07, -0.6118e-08,  0.1624e-07,           &
     0.5791e-07,  0.1205e-06,  0.2025e-06,  0.3032e-06,  0.4225e-06,           &
     0.5619e-06,  0.7322e-06,  0.9528e-06,  0.1243e-05,  0.1592e-05,           &
     0.1904e-05/
   data ((h81(ip,iw,1),iw=1,31),ip=14,14)/                                     &
     0.1346e-04,  0.2655e-04,  0.5192e-04,  0.1001e-03,  0.1887e-03,           &
     0.3442e-03,  0.6001e-03,  0.9892e-03,  0.1536e-02,  0.2262e-02,           &
     0.3200e-02,  0.4420e-02,  0.6053e-02,  0.8300e-02,  0.1146e-01,           &
     0.1592e-01,  0.2224e-01,  0.3112e-01,  0.4348e-01,  0.6051e-01,           &
     0.8369e-01,  0.1147e+00,  0.1556e+00,  0.2084e+00,  0.2749e+00,           &
     0.3561e+00,  0.4511e+00,  0.5562e+00,  0.6644e+00,  0.7655e+00,           &
     0.8499e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=14,14)/                                     &
    -0.5210e-08, -0.1172e-07, -0.2731e-07, -0.6598e-07, -0.1615e-06,           &
    -0.3880e-06, -0.8769e-06, -0.1787e-05, -0.3204e-05, -0.5066e-05,           &
    -0.7197e-05, -0.9451e-05, -0.1185e-04, -0.1465e-04, -0.1831e-04,           &
    -0.2346e-04, -0.3088e-04, -0.4132e-04, -0.5545e-04, -0.7410e-04,           &
    -0.9820e-04, -0.1288e-03, -0.1670e-03, -0.2140e-03, -0.2692e-03,           &
    -0.3293e-03, -0.3886e-03, -0.4417e-03, -0.4840e-03, -0.5073e-03,           &
    -0.4944e-03/
   data ((h83(ip,iw,1),iw=1,31),ip=14,14)/                                     &
     0.1880e-10,  0.4271e-10,  0.9966e-10,  0.2352e-09,  0.5497e-09,           &
     0.1205e-08,  0.2334e-08,  0.3765e-08,  0.4993e-08,  0.5532e-08,           &
     0.5148e-08,  0.4055e-08,  0.2650e-08,  0.1326e-08,  0.2019e-09,           &
    -0.1124e-08, -0.2234e-08,  0.2827e-09,  0.1247e-07,  0.4102e-07,           &
     0.9228e-07,  0.1682e-06,  0.2676e-06,  0.3885e-06,  0.5286e-06,           &
     0.6904e-06,  0.8871e-06,  0.1142e-05,  0.1466e-05,  0.1800e-05,           &
     0.2004e-05/
   data ((h81(ip,iw,1),iw=1,31),ip=15,15)/                                     &
     0.1348e-04,  0.2663e-04,  0.5216e-04,  0.1009e-03,  0.1912e-03,           &
     0.3515e-03,  0.6196e-03,  0.1036e-02,  0.1637e-02,  0.2456e-02,           &
     0.3544e-02,  0.4996e-02,  0.6973e-02,  0.9722e-02,  0.1360e-01,           &
     0.1908e-01,  0.2681e-01,  0.3766e-01,  0.5272e-01,  0.7343e-01,           &
     0.1015e+00,  0.1388e+00,  0.1874e+00,  0.2492e+00,  0.3257e+00,           &
     0.4165e+00,  0.5191e+00,  0.6275e+00,  0.7324e+00,  0.8235e+00,           &
     0.8939e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=15,15)/                                     &
    -0.5045e-08, -0.1113e-07, -0.2540e-07, -0.6008e-07, -0.1449e-06,           &
    -0.3457e-06, -0.7826e-06, -0.1609e-05, -0.2920e-05, -0.4665e-05,           &
    -0.6691e-05, -0.8868e-05, -0.1127e-04, -0.1422e-04, -0.1820e-04,           &
    -0.2389e-04, -0.3213e-04, -0.4380e-04, -0.5975e-04, -0.8092e-04,           &
    -0.1083e-03, -0.1433e-03, -0.1873e-03, -0.2402e-03, -0.2997e-03,           &
    -0.3607e-03, -0.4178e-03, -0.4662e-03, -0.4994e-03, -0.5028e-03,           &
    -0.4563e-03/
   data ((h83(ip,iw,1),iw=1,31),ip=15,15)/                                     &
     0.1804e-10,  0.3983e-10,  0.9045e-10,  0.2080e-09,  0.4786e-09,           &
     0.1046e-08,  0.2052e-08,  0.3413e-08,  0.4704e-08,  0.5565e-08,           &
     0.5887e-08,  0.5981e-08,  0.6202e-08,  0.6998e-08,  0.8493e-08,           &
     0.1002e-07,  0.1184e-07,  0.1780e-07,  0.3483e-07,  0.7122e-07,           &
     0.1341e-06,  0.2259e-06,  0.3446e-06,  0.4866e-06,  0.6486e-06,           &
     0.8343e-06,  0.1063e-05,  0.1356e-05,  0.1690e-05,  0.1951e-05,           &
     0.2005e-05/
   data ((h81(ip,iw,1),iw=1,31),ip=16,16)/                                     &
     0.1351e-04,  0.2672e-04,  0.5247e-04,  0.1019e-03,  0.1944e-03,           &
     0.3609e-03,  0.6451e-03,  0.1098e-02,  0.1772e-02,  0.2718e-02,           &
     0.4009e-02,  0.5767e-02,  0.8190e-02,  0.1158e-01,  0.1636e-01,           &
     0.2311e-01,  0.3260e-01,  0.4586e-01,  0.6425e-01,  0.8940e-01,           &
     0.1232e+00,  0.1678e+00,  0.2251e+00,  0.2967e+00,  0.3831e+00,           &
     0.4825e+00,  0.5901e+00,  0.6976e+00,  0.7946e+00,  0.8725e+00,           &
     0.9285e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=16,16)/                                     &
    -0.4850e-08, -0.1045e-07, -0.2334e-07, -0.5367e-07, -0.1265e-06,           &
    -0.2980e-06, -0.6750e-06, -0.1406e-05, -0.2601e-05, -0.4239e-05,           &
    -0.6201e-05, -0.8389e-05, -0.1091e-04, -0.1413e-04, -0.1859e-04,           &
    -0.2500e-04, -0.3432e-04, -0.4761e-04, -0.6595e-04, -0.9030e-04,           &
    -0.1219e-03, -0.1624e-03, -0.2126e-03, -0.2708e-03, -0.3327e-03,           &
    -0.3926e-03, -0.4458e-03, -0.4871e-03, -0.5045e-03, -0.4777e-03,           &
    -0.3954e-03/
   data ((h83(ip,iw,1),iw=1,31),ip=16,16)/                                     &
     0.1717e-10,  0.3723e-10,  0.8093e-10,  0.1817e-09,  0.4100e-09,           &
     0.8932e-09,  0.1791e-08,  0.3126e-08,  0.4634e-08,  0.6095e-08,           &
     0.7497e-08,  0.9170e-08,  0.1136e-07,  0.1453e-07,  0.1892e-07,           &
     0.2369e-07,  0.2909e-07,  0.3922e-07,  0.6232e-07,  0.1083e-06,           &
     0.1847e-06,  0.2943e-06,  0.4336e-06,  0.5970e-06,  0.7815e-06,           &
     0.9959e-06,  0.1263e-05,  0.1583e-05,  0.1880e-05,  0.2009e-05,           &
     0.1914e-05/
   data ((h81(ip,iw,1),iw=1,31),ip=17,17)/                                     &
     0.1353e-04,  0.2681e-04,  0.5283e-04,  0.1031e-03,  0.1981e-03,           &
     0.3721e-03,  0.6761e-03,  0.1176e-02,  0.1943e-02,  0.3055e-02,           &
     0.4612e-02,  0.6766e-02,  0.9760e-02,  0.1396e-01,  0.1986e-01,           &
     0.2818e-01,  0.3984e-01,  0.5609e-01,  0.7850e-01,  0.1090e+00,           &
     0.1495e+00,  0.2023e+00,  0.2691e+00,  0.3508e+00,  0.4465e+00,           &
     0.5524e+00,  0.6613e+00,  0.7633e+00,  0.8485e+00,  0.9119e+00,           &
     0.9543e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=17,17)/                                     &
    -0.4673e-08, -0.9862e-08, -0.2135e-07, -0.4753e-07, -0.1087e-06,           &
    -0.2512e-06, -0.5671e-06, -0.1199e-05, -0.2281e-05, -0.3842e-05,           &
    -0.5804e-05, -0.8110e-05, -0.1088e-04, -0.1452e-04, -0.1961e-04,           &
    -0.2696e-04, -0.3768e-04, -0.5311e-04, -0.7444e-04, -0.1028e-03,           &
    -0.1397e-03, -0.1865e-03, -0.2427e-03, -0.3047e-03, -0.3667e-03,           &
    -0.4237e-03, -0.4712e-03, -0.5003e-03, -0.4921e-03, -0.4286e-03,           &
    -0.3188e-03/
   data ((h83(ip,iw,1),iw=1,31),ip=17,17)/                                     &
     0.1653e-10,  0.3436e-10,  0.7431e-10,  0.1605e-09,  0.3548e-09,           &
     0.7723e-09,  0.1595e-08,  0.2966e-08,  0.4849e-08,  0.7169e-08,           &
     0.1003e-07,  0.1366e-07,  0.1825e-07,  0.2419e-07,  0.3186e-07,           &
     0.4068e-07,  0.5064e-07,  0.6618e-07,  0.9684e-07,  0.1536e-06,           &
     0.2450e-06,  0.3730e-06,  0.5328e-06,  0.7184e-06,  0.9291e-06,           &
     0.1180e-05,  0.1484e-05,  0.1798e-05,  0.1992e-05,  0.1968e-05,           &
     0.1736e-05/
   data ((h81(ip,iw,1),iw=1,31),ip=18,18)/                                     &
     0.1356e-04,  0.2691e-04,  0.5318e-04,  0.1043e-03,  0.2020e-03,           &
     0.3841e-03,  0.7100e-03,  0.1263e-02,  0.2144e-02,  0.3464e-02,           &
     0.5359e-02,  0.8019e-02,  0.1174e-01,  0.1694e-01,  0.2426e-01,           &
     0.3451e-01,  0.4883e-01,  0.6870e-01,  0.9593e-01,  0.1326e+00,           &
     0.1809e+00,  0.2428e+00,  0.3196e+00,  0.4112e+00,  0.5147e+00,           &
     0.6239e+00,  0.7297e+00,  0.8217e+00,  0.8928e+00,  0.9421e+00,           &
     0.9726e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=18,18)/                                     &
    -0.4532e-08, -0.9395e-08, -0.1978e-07, -0.4272e-07, -0.9442e-07,           &
    -0.2124e-06, -0.4747e-06, -0.1017e-05, -0.2003e-05, -0.3524e-05,           &
    -0.5567e-05, -0.8108e-05, -0.1127e-04, -0.1547e-04, -0.2138e-04,           &
    -0.2996e-04, -0.4251e-04, -0.6059e-04, -0.8563e-04, -0.1190e-03,           &
    -0.1623e-03, -0.2156e-03, -0.2767e-03, -0.3403e-03, -0.4006e-03,           &
    -0.4530e-03, -0.4912e-03, -0.4995e-03, -0.4563e-03, -0.3592e-03,           &
    -0.2383e-03/
   data ((h83(ip,iw,1),iw=1,31),ip=18,18)/                                     &
     0.1593e-10,  0.3276e-10,  0.6896e-10,  0.1476e-09,  0.3190e-09,           &
     0.6944e-09,  0.1474e-08,  0.2935e-08,  0.5300e-08,  0.8697e-08,           &
     0.1336e-07,  0.1946e-07,  0.2707e-07,  0.3637e-07,  0.4800e-07,           &
     0.6187e-07,  0.7806e-07,  0.1008e-06,  0.1404e-06,  0.2089e-06,           &
     0.3153e-06,  0.4613e-06,  0.6416e-06,  0.8506e-06,  0.1095e-05,           &
     0.1387e-05,  0.1708e-05,  0.1956e-05,  0.2003e-05,  0.1836e-05,           &
     0.1483e-05/
   data ((h81(ip,iw,1),iw=1,31),ip=19,19)/                                     &
     0.1358e-04,  0.2700e-04,  0.5350e-04,  0.1054e-03,  0.2057e-03,           &
     0.3955e-03,  0.7434e-03,  0.1353e-02,  0.2361e-02,  0.3928e-02,           &
     0.6237e-02,  0.9531e-02,  0.1416e-01,  0.2063e-01,  0.2968e-01,           &
     0.4233e-01,  0.5990e-01,  0.8412e-01,  0.1170e+00,  0.1608e+00,           &
     0.2177e+00,  0.2894e+00,  0.3764e+00,  0.4768e+00,  0.5855e+00,           &
     0.6941e+00,  0.7921e+00,  0.8710e+00,  0.9276e+00,  0.9641e+00,           &
     0.9846e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=19,19)/                                     &
    -0.4448e-08, -0.9085e-08, -0.1877e-07, -0.3946e-07, -0.8472e-07,           &
    -0.1852e-06, -0.4074e-06, -0.8791e-06, -0.1789e-05, -0.3314e-05,           &
    -0.5521e-05, -0.8425e-05, -0.1215e-04, -0.1711e-04, -0.2407e-04,           &
    -0.3421e-04, -0.4905e-04, -0.7032e-04, -0.9985e-04, -0.1394e-03,           &
    -0.1897e-03, -0.2491e-03, -0.3132e-03, -0.3763e-03, -0.4332e-03,           &
    -0.4786e-03, -0.5005e-03, -0.4775e-03, -0.3970e-03, -0.2794e-03,           &
    -0.1652e-03/
   data ((h83(ip,iw,1),iw=1,31),ip=19,19)/                                     &
     0.1566e-10,  0.3219e-10,  0.6635e-10,  0.1400e-09,  0.2999e-09,           &
     0.6513e-09,  0.1406e-08,  0.2953e-08,  0.5789e-08,  0.1037e-07,           &
     0.1709e-07,  0.2623e-07,  0.3777e-07,  0.5159e-07,  0.6823e-07,           &
     0.8864e-07,  0.1134e-06,  0.1461e-06,  0.1960e-06,  0.2761e-06,           &
     0.3962e-06,  0.5583e-06,  0.7580e-06,  0.9957e-06,  0.1282e-05,           &
     0.1607e-05,  0.1898e-05,  0.2020e-05,  0.1919e-05,  0.1623e-05,           &
     0.1171e-05/
   data ((h81(ip,iw,1),iw=1,31),ip=20,20)/                                     &
     0.1360e-04,  0.2707e-04,  0.5375e-04,  0.1063e-03,  0.2087e-03,           &
     0.4052e-03,  0.7731e-03,  0.1437e-02,  0.2575e-02,  0.4414e-02,           &
     0.7212e-02,  0.1128e-01,  0.1704e-01,  0.2509e-01,  0.3631e-01,           &
     0.5188e-01,  0.7338e-01,  0.1028e+00,  0.1422e+00,  0.1940e+00,           &
     0.2603e+00,  0.3422e+00,  0.4387e+00,  0.5459e+00,  0.6563e+00,           &
     0.7597e+00,  0.8461e+00,  0.9105e+00,  0.9536e+00,  0.9791e+00,           &
     0.9920e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=20,20)/                                     &
    -0.4403e-08, -0.8896e-08, -0.1818e-07, -0.3751e-07, -0.7880e-07,           &
    -0.1683e-06, -0.3640e-06, -0.7852e-06, -0.1640e-05, -0.3191e-05,           &
    -0.5634e-05, -0.9046e-05, -0.1355e-04, -0.1953e-04, -0.2786e-04,           &
    -0.3995e-04, -0.5752e-04, -0.8256e-04, -0.1174e-03, -0.1638e-03,           &
    -0.2211e-03, -0.2854e-03, -0.3507e-03, -0.4116e-03, -0.4633e-03,           &
    -0.4966e-03, -0.4921e-03, -0.4309e-03, -0.3215e-03, -0.2016e-03,           &
    -0.1061e-03/
   data ((h83(ip,iw,1),iw=1,31),ip=20,20)/                                     &
     0.1551e-10,  0.3147e-10,  0.6419e-10,  0.1356e-09,  0.2860e-09,           &
     0.6178e-09,  0.1353e-08,  0.2934e-08,  0.6095e-08,  0.1174e-07,           &
     0.2067e-07,  0.3346e-07,  0.5014e-07,  0.7024e-07,  0.9377e-07,           &
     0.1226e-06,  0.1592e-06,  0.2056e-06,  0.2678e-06,  0.3584e-06,           &
     0.4892e-06,  0.6651e-06,  0.8859e-06,  0.1160e-05,  0.1488e-05,           &
     0.1814e-05,  0.2010e-05,  0.1984e-05,  0.1748e-05,  0.1338e-05,           &
     0.8445e-06/
   data ((h81(ip,iw,1),iw=1,31),ip=21,21)/                                     &
     0.1361e-04,  0.2711e-04,  0.5393e-04,  0.1070e-03,  0.2110e-03,           &
     0.4129e-03,  0.7973e-03,  0.1508e-02,  0.2769e-02,  0.4886e-02,           &
     0.8229e-02,  0.1322e-01,  0.2037e-01,  0.3038e-01,  0.4427e-01,           &
     0.6342e-01,  0.8963e-01,  0.1250e+00,  0.1720e+00,  0.2327e+00,           &
     0.3089e+00,  0.4007e+00,  0.5052e+00,  0.6163e+00,  0.7241e+00,           &
     0.8179e+00,  0.8905e+00,  0.9408e+00,  0.9720e+00,  0.9887e+00,           &
     0.9961e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=21,21)/                                     &
    -0.4379e-08, -0.8801e-08, -0.1782e-07, -0.3642e-07, -0.7536e-07,           &
    -0.1581e-06, -0.3366e-06, -0.7227e-06, -0.1532e-05, -0.3106e-05,           &
    -0.5810e-05, -0.9862e-05, -0.1540e-04, -0.2279e-04, -0.3292e-04,           &
    -0.4738e-04, -0.6817e-04, -0.9765e-04, -0.1384e-03, -0.1918e-03,           &
    -0.2551e-03, -0.3226e-03, -0.3876e-03, -0.4452e-03, -0.4883e-03,           &
    -0.5005e-03, -0.4598e-03, -0.3633e-03, -0.2416e-03, -0.1349e-03,           &
    -0.6278e-04/
   data ((h83(ip,iw,1),iw=1,31),ip=21,21)/                                     &
     0.1542e-10,  0.3111e-10,  0.6345e-10,  0.1310e-09,  0.2742e-09,           &
     0.5902e-09,  0.1289e-08,  0.2826e-08,  0.6103e-08,  0.1250e-07,           &
     0.2355e-07,  0.4041e-07,  0.6347e-07,  0.9217e-07,  0.1256e-06,           &
     0.1658e-06,  0.2175e-06,  0.2824e-06,  0.3607e-06,  0.4614e-06,           &
     0.6004e-06,  0.7880e-06,  0.1034e-05,  0.1349e-05,  0.1698e-05,           &
     0.1965e-05,  0.2021e-05,  0.1857e-05,  0.1500e-05,  0.1015e-05,           &
     0.5467e-06/
   data ((h81(ip,iw,1),iw=1,31),ip=22,22)/                                     &
     0.1362e-04,  0.2714e-04,  0.5405e-04,  0.1074e-03,  0.2127e-03,           &
     0.4187e-03,  0.8159e-03,  0.1565e-02,  0.2933e-02,  0.5314e-02,           &
     0.9222e-02,  0.1525e-01,  0.2406e-01,  0.3649e-01,  0.5367e-01,           &
     0.7717e-01,  0.1090e+00,  0.1514e+00,  0.2067e+00,  0.2770e+00,           &
     0.3633e+00,  0.4640e+00,  0.5742e+00,  0.6853e+00,  0.7859e+00,           &
     0.8670e+00,  0.9253e+00,  0.9629e+00,  0.9842e+00,  0.9942e+00,           &
     0.9983e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=22,22)/                                     &
    -0.4366e-08, -0.8749e-08, -0.1761e-07, -0.3578e-07, -0.7322e-07,           &
    -0.1517e-06, -0.3189e-06, -0.6785e-06, -0.1446e-05, -0.3014e-05,           &
    -0.5933e-05, -0.1069e-04, -0.1755e-04, -0.2683e-04, -0.3936e-04,           &
    -0.5675e-04, -0.8137e-04, -0.1160e-03, -0.1630e-03, -0.2223e-03,           &
    -0.2899e-03, -0.3589e-03, -0.4230e-03, -0.4755e-03, -0.5031e-03,           &
    -0.4834e-03, -0.4036e-03, -0.2849e-03, -0.1687e-03, -0.8356e-04,           &
    -0.3388e-04/
   data ((h83(ip,iw,1),iw=1,31),ip=22,22)/                                     &
     0.1536e-10,  0.3086e-10,  0.6248e-10,  0.1288e-09,  0.2664e-09,           &
     0.5637e-09,  0.1222e-08,  0.2680e-08,  0.5899e-08,  0.1262e-07,           &
     0.2527e-07,  0.4621e-07,  0.7678e-07,  0.1165e-06,  0.1640e-06,           &
     0.2199e-06,  0.2904e-06,  0.3783e-06,  0.4787e-06,  0.5925e-06,           &
     0.7377e-06,  0.9389e-06,  0.1216e-05,  0.1560e-05,  0.1879e-05,           &
     0.2025e-05,  0.1940e-05,  0.1650e-05,  0.1194e-05,  0.6981e-06,           &
     0.3103e-06/
   data ((h81(ip,iw,1),iw=1,31),ip=23,23)/                                     &
     0.1362e-04,  0.2716e-04,  0.5413e-04,  0.1077e-03,  0.2139e-03,           &
     0.4228e-03,  0.8295e-03,  0.1609e-02,  0.3063e-02,  0.5676e-02,           &
     0.1013e-01,  0.1726e-01,  0.2799e-01,  0.4332e-01,  0.6452e-01,           &
     0.9328e-01,  0.1317e+00,  0.1820e+00,  0.2467e+00,  0.3270e+00,           &
     0.4228e+00,  0.5308e+00,  0.6434e+00,  0.7499e+00,  0.8395e+00,           &
     0.9065e+00,  0.9515e+00,  0.9782e+00,  0.9916e+00,  0.9973e+00,           &
     0.9993e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=23,23)/                                     &
    -0.4359e-08, -0.8720e-08, -0.1749e-07, -0.3527e-07, -0.7175e-07,           &
    -0.1473e-06, -0.3062e-06, -0.6451e-06, -0.1372e-05, -0.2902e-05,           &
    -0.5936e-05, -0.1133e-04, -0.1971e-04, -0.3143e-04, -0.4715e-04,           &
    -0.6833e-04, -0.9759e-04, -0.1379e-03, -0.1907e-03, -0.2542e-03,           &
    -0.3239e-03, -0.3935e-03, -0.4559e-03, -0.4991e-03, -0.5009e-03,           &
    -0.4414e-03, -0.3306e-03, -0.2077e-03, -0.1093e-03, -0.4754e-04,           &
    -0.1642e-04/
   data ((h83(ip,iw,1),iw=1,31),ip=23,23)/                                     &
     0.1531e-10,  0.3070e-10,  0.6184e-10,  0.1257e-09,  0.2578e-09,           &
     0.5451e-09,  0.1159e-08,  0.2526e-08,  0.5585e-08,  0.1225e-07,           &
     0.2576e-07,  0.5017e-07,  0.8855e-07,  0.1417e-06,  0.2078e-06,           &
     0.2858e-06,  0.3802e-06,  0.4946e-06,  0.6226e-06,  0.7572e-06,           &
     0.9137e-06,  0.1133e-05,  0.1438e-05,  0.1772e-05,  0.1994e-05,           &
     0.1994e-05,  0.1779e-05,  0.1375e-05,  0.8711e-06,  0.4273e-06,           &
     0.1539e-06/
   data ((h81(ip,iw,1),iw=1,31),ip=24,24)/                                     &
     0.1362e-04,  0.2717e-04,  0.5417e-04,  0.1079e-03,  0.2146e-03,           &
     0.4256e-03,  0.8393e-03,  0.1641e-02,  0.3163e-02,  0.5968e-02,           &
     0.1092e-01,  0.1915e-01,  0.3195e-01,  0.5066e-01,  0.7669e-01,           &
     0.1117e+00,  0.1578e+00,  0.2173e+00,  0.2919e+00,  0.3824e+00,           &
     0.4867e+00,  0.5992e+00,  0.7100e+00,  0.8076e+00,  0.8838e+00,           &
     0.9370e+00,  0.9702e+00,  0.9880e+00,  0.9959e+00,  0.9988e+00,           &
     0.9998e+00/
   data ((h82(ip,iw,1),iw=1,31),ip=24,24)/                                     &
    -0.4354e-08, -0.8703e-08, -0.1742e-07, -0.3499e-07, -0.7074e-07,           &
    -0.1441e-06, -0.2971e-06, -0.6195e-06, -0.1309e-05, -0.2780e-05,           &
    -0.5823e-05, -0.1165e-04, -0.2152e-04, -0.3616e-04, -0.5604e-04,           &
    -0.8230e-04, -0.1173e-03, -0.1635e-03, -0.2211e-03, -0.2868e-03,           &
    -0.3567e-03, -0.4260e-03, -0.4844e-03, -0.5097e-03, -0.4750e-03,           &
    -0.3779e-03, -0.2522e-03, -0.1409e-03, -0.6540e-04, -0.2449e-04,           &
    -0.6948e-05/
   data ((h83(ip,iw,1),iw=1,31),ip=24,24)/                                     &
     0.1529e-10,  0.3060e-10,  0.6142e-10,  0.1241e-09,  0.2535e-09,           &
     0.5259e-09,  0.1107e-08,  0.2383e-08,  0.5243e-08,  0.1161e-07,           &
     0.2523e-07,  0.5188e-07,  0.9757e-07,  0.1657e-06,  0.2553e-06,           &
     0.3629e-06,  0.4878e-06,  0.6323e-06,  0.7923e-06,  0.9575e-06,           &
     0.1139e-05,  0.1381e-05,  0.1687e-05,  0.1952e-05,  0.2029e-05,           &
     0.1890e-05,  0.1552e-05,  0.1062e-05,  0.5728e-06,  0.2280e-06,           &
     0.6762e-07/
   data ((h81(ip,iw,1),iw=1,31),ip=25,25)/                                     &
     0.1363e-04,  0.2718e-04,  0.5420e-04,  0.1080e-03,  0.2151e-03,           &
     0.4275e-03,  0.8461e-03,  0.1664e-02,  0.3237e-02,  0.6195e-02,           &
     0.1156e-01,  0.2080e-01,  0.3573e-01,  0.5818e-01,  0.8981e-01,           &
     0.1323e+00,  0.1874e+00,  0.2570e+00,  0.3423e+00,  0.4425e+00,           &
     0.5534e+00,  0.6669e+00,  0.7714e+00,  0.8568e+00,  0.9191e+00,           &
     0.9597e+00,  0.9828e+00,  0.9938e+00,  0.9981e+00,  0.9996e+00,           &
     0.1000e+01/
   data ((h82(ip,iw,1),iw=1,31),ip=25,25)/                                     &
    -0.4352e-08, -0.8693e-08, -0.1738e-07, -0.3483e-07, -0.7006e-07,           &
    -0.1423e-06, -0.2905e-06, -0.6008e-06, -0.1258e-05, -0.2663e-05,           &
    -0.5638e-05, -0.1165e-04, -0.2270e-04, -0.4044e-04, -0.6554e-04,           &
    -0.9855e-04, -0.1407e-03, -0.1928e-03, -0.2534e-03, -0.3197e-03,           &
    -0.3890e-03, -0.4563e-03, -0.5040e-03, -0.4998e-03, -0.4249e-03,           &
    -0.3025e-03, -0.1794e-03, -0.8860e-04, -0.3575e-04, -0.1122e-04,           &
    -0.2506e-05/
   data ((h83(ip,iw,1),iw=1,31),ip=25,25)/                                     &
     0.1527e-10,  0.3053e-10,  0.6115e-10,  0.1230e-09,  0.2492e-09,           &
     0.5149e-09,  0.1068e-08,  0.2268e-08,  0.4932e-08,  0.1089e-07,           &
     0.2408e-07,  0.5156e-07,  0.1028e-06,  0.1859e-06,  0.3028e-06,           &
     0.4476e-06,  0.6124e-06,  0.7932e-06,  0.9879e-06,  0.1194e-05,           &
     0.1417e-05,  0.1673e-05,  0.1929e-05,  0.2064e-05,  0.1997e-05,           &
     0.1725e-05,  0.1267e-05,  0.7464e-06,  0.3312e-06,  0.1066e-06,           &
     0.2718e-07/
   data ((h81(ip,iw,1),iw=1,31),ip=26,26)/                                     &
     0.1363e-04,  0.2719e-04,  0.5422e-04,  0.1081e-03,  0.2154e-03,           &
     0.4287e-03,  0.8506e-03,  0.1680e-02,  0.3291e-02,  0.6364e-02,           &
     0.1206e-01,  0.2218e-01,  0.3912e-01,  0.6544e-01,  0.1033e+00,           &
     0.1544e+00,  0.2198e+00,  0.3008e+00,  0.3970e+00,  0.5060e+00,           &
     0.6209e+00,  0.7311e+00,  0.8254e+00,  0.8972e+00,  0.9462e+00,           &
     0.9757e+00,  0.9907e+00,  0.9970e+00,  0.9992e+00,  0.9999e+00,           &
     0.1000e+01/
   data ((h82(ip,iw,1),iw=1,31),ip=26,26)/                                     &
    -0.4351e-08, -0.8688e-08, -0.1736e-07, -0.3473e-07, -0.6966e-07,           &
    -0.1405e-06, -0.2857e-06, -0.5867e-06, -0.1218e-05, -0.2563e-05,           &
    -0.5435e-05, -0.1144e-04, -0.2321e-04, -0.4379e-04, -0.7487e-04,           &
    -0.1163e-03, -0.1670e-03, -0.2250e-03, -0.2876e-03, -0.3535e-03,           &
    -0.4215e-03, -0.4826e-03, -0.5082e-03, -0.4649e-03, -0.3564e-03,           &
    -0.2264e-03, -0.1188e-03, -0.5128e-04, -0.1758e-04, -0.4431e-05,           &
    -0.7275e-06/
   data ((h83(ip,iw,1),iw=1,31),ip=26,26)/                                     &
     0.1525e-10,  0.3048e-10,  0.6097e-10,  0.1223e-09,  0.2466e-09,           &
     0.5021e-09,  0.1032e-08,  0.2195e-08,  0.4688e-08,  0.1027e-07,           &
     0.2279e-07,  0.4999e-07,  0.1046e-06,  0.2009e-06,  0.3460e-06,           &
     0.5335e-06,  0.7478e-06,  0.9767e-06,  0.1216e-05,  0.1469e-05,           &
     0.1735e-05,  0.1977e-05,  0.2121e-05,  0.2103e-05,  0.1902e-05,           &
     0.1495e-05,  0.9541e-06,  0.4681e-06,  0.1672e-06,  0.4496e-07,           &
     0.9859e-08/
   data ((c1(ip,iw,1),iw=1,24),ip= 1, 1)/                                      &
     0.1444e-03,  0.2378e-03,  0.3644e-03,  0.5245e-03,  0.7311e-03,           &
     0.1015e-02,  0.1416e-02,  0.1975e-02,  0.2715e-02,  0.3651e-02,           &
     0.4813e-02,  0.6264e-02,  0.8084e-02,  0.1037e-01,  0.1327e-01,           &
     0.1693e-01,  0.2155e-01,  0.2734e-01,  0.3453e-01,  0.4337e-01,           &
     0.5414e-01,  0.6726e-01,  0.8335e-01,  0.1032e+00/
   data ((c2(ip,iw,1),iw=1,24),ip= 1, 1)/                                      &
    -0.2878e-06, -0.6295e-06, -0.1270e-05, -0.2323e-05, -0.3910e-05,           &
    -0.6235e-05, -0.9527e-05, -0.1391e-04, -0.1947e-04, -0.2652e-04,           &
    -0.3568e-04, -0.4777e-04, -0.6363e-04, -0.8414e-04, -0.1102e-03,           &
    -0.1425e-03, -0.1822e-03, -0.2301e-03, -0.2869e-03, -0.3534e-03,           &
    -0.4305e-03, -0.5188e-03, -0.6173e-03, -0.7229e-03/
   data ((c3(ip,iw,1),iw=1,24),ip= 1, 1)/                                      &
    -0.2369e-09, -0.5453e-09, -0.1434e-08, -0.3513e-08, -0.7277e-08,           &
    -0.1283e-07, -0.1983e-07, -0.2838e-07, -0.4075e-07, -0.6115e-07,           &
    -0.9201e-07, -0.1325e-06, -0.1844e-06, -0.2534e-06, -0.3430e-06,           &
    -0.4550e-06, -0.5896e-06, -0.7434e-06, -0.9121e-06, -0.1096e-05,           &
    -0.1298e-05, -0.1508e-05, -0.1689e-05, -0.1773e-05/
   data ((c1(ip,iw,1),iw=1,24),ip= 2, 2)/                                      &
     0.1444e-03,  0.2379e-03,  0.3646e-03,  0.5250e-03,  0.7320e-03,           &
     0.1016e-02,  0.1419e-02,  0.1981e-02,  0.2726e-02,  0.3669e-02,           &
     0.4847e-02,  0.6325e-02,  0.8193e-02,  0.1056e-01,  0.1358e-01,           &
     0.1742e-01,  0.2228e-01,  0.2843e-01,  0.3611e-01,  0.4561e-01,           &
     0.5730e-01,  0.7169e-01,  0.8952e-01,  0.1117e+00/
   data ((c2(ip,iw,1),iw=1,24),ip= 2, 2)/                                      &
    -0.2875e-06, -0.6289e-06, -0.1269e-05, -0.2321e-05, -0.3907e-05,           &
    -0.6229e-05, -0.9516e-05, -0.1389e-04, -0.1944e-04, -0.2646e-04,           &
    -0.3559e-04, -0.4762e-04, -0.6340e-04, -0.8384e-04, -0.1098e-03,           &
    -0.1423e-03, -0.1824e-03, -0.2309e-03, -0.2887e-03, -0.3566e-03,           &
    -0.4355e-03, -0.5260e-03, -0.6274e-03, -0.7367e-03/
   data ((c3(ip,iw,1),iw=1,24),ip= 2, 2)/                                      &
    -0.2382e-09, -0.5476e-09, -0.1437e-08, -0.3520e-08, -0.7291e-08,           &
    -0.1286e-07, -0.1989e-07, -0.2849e-07, -0.4093e-07, -0.6144e-07,           &
    -0.9241e-07, -0.1330e-06, -0.1852e-06, -0.2542e-06, -0.3435e-06,           &
    -0.4549e-06, -0.5891e-06, -0.7427e-06, -0.9114e-06, -0.1096e-05,           &
    -0.1298e-05, -0.1509e-05, -0.1689e-05, -0.1772e-05/
   data ((c1(ip,iw,1),iw=1,24),ip= 3, 3)/                                      &
     0.1445e-03,  0.2381e-03,  0.3650e-03,  0.5258e-03,  0.7335e-03,           &
     0.1019e-02,  0.1424e-02,  0.1990e-02,  0.2742e-02,  0.3698e-02,           &
     0.4899e-02,  0.6418e-02,  0.8353e-02,  0.1083e-01,  0.1401e-01,           &
     0.1807e-01,  0.2326e-01,  0.2985e-01,  0.3814e-01,  0.4849e-01,           &
     0.6134e-01,  0.7733e-01,  0.9732e-01,  0.1223e+00/
   data ((c2(ip,iw,1),iw=1,24),ip= 3, 3)/                                      &
    -0.2871e-06, -0.6279e-06, -0.1267e-05, -0.2318e-05, -0.3901e-05,           &
    -0.6219e-05, -0.9500e-05, -0.1386e-04, -0.1939e-04, -0.2639e-04,           &
    -0.3547e-04, -0.4743e-04, -0.6313e-04, -0.8354e-04, -0.1096e-03,           &
    -0.1425e-03, -0.1831e-03, -0.2325e-03, -0.2916e-03, -0.3613e-03,           &
    -0.4425e-03, -0.5361e-03, -0.6412e-03, -0.7550e-03/
   data ((c3(ip,iw,1),iw=1,24),ip= 3, 3)/                                      &
    -0.2404e-09, -0.5516e-09, -0.1443e-08, -0.3531e-08, -0.7310e-08,           &
    -0.1290e-07, -0.1997e-07, -0.2865e-07, -0.4121e-07, -0.6185e-07,           &
    -0.9295e-07, -0.1337e-06, -0.1860e-06, -0.2549e-06, -0.3438e-06,           &
    -0.4547e-06, -0.5888e-06, -0.7424e-06, -0.9114e-06, -0.1097e-05,           &
    -0.1300e-05, -0.1511e-05, -0.1691e-05, -0.1771e-05/
   data ((c1(ip,iw,1),iw=1,24),ip= 4, 4)/                                      &
     0.1446e-03,  0.2383e-03,  0.3656e-03,  0.5270e-03,  0.7358e-03,           &
     0.1023e-02,  0.1432e-02,  0.2003e-02,  0.2766e-02,  0.3743e-02,           &
     0.4979e-02,  0.6557e-02,  0.8587e-02,  0.1121e-01,  0.1459e-01,           &
     0.1895e-01,  0.2455e-01,  0.3170e-01,  0.4076e-01,  0.5217e-01,           &
     0.6649e-01,  0.8448e-01,  0.1071e+00,  0.1353e+00/
   data ((c2(ip,iw,1),iw=1,24),ip= 4, 4)/                                      &
    -0.2865e-06, -0.6264e-06, -0.1264e-05, -0.2313e-05, -0.3892e-05,           &
    -0.6204e-05, -0.9474e-05, -0.1382e-04, -0.1932e-04, -0.2628e-04,           &
    -0.3530e-04, -0.4720e-04, -0.6286e-04, -0.8333e-04, -0.1097e-03,           &
    -0.1430e-03, -0.1845e-03, -0.2352e-03, -0.2960e-03, -0.3681e-03,           &
    -0.4523e-03, -0.5496e-03, -0.6595e-03, -0.7790e-03/
   data ((c3(ip,iw,1),iw=1,24),ip= 4, 4)/                                      &
    -0.2433e-09, -0.5580e-09, -0.1454e-08, -0.3550e-08, -0.7345e-08,           &
    -0.1296e-07, -0.2010e-07, -0.2889e-07, -0.4159e-07, -0.6239e-07,           &
    -0.9360e-07, -0.1345e-06, -0.1867e-06, -0.2554e-06, -0.3439e-06,           &
    -0.4548e-06, -0.5890e-06, -0.7429e-06, -0.9127e-06, -0.1099e-05,           &
    -0.1304e-05, -0.1516e-05, -0.1693e-05, -0.1770e-05/
   data ((c1(ip,iw,1),iw=1,24),ip= 5, 5)/                                     &
     0.1448e-03,  0.2387e-03,  0.3666e-03,  0.5290e-03,  0.7394e-03,           &
     0.1030e-02,  0.1444e-02,  0.2025e-02,  0.2805e-02,  0.3811e-02,           &
     0.5099e-02,  0.6761e-02,  0.8919e-02,  0.1173e-01,  0.1538e-01,           &
     0.2011e-01,  0.2622e-01,  0.3408e-01,  0.4412e-01,  0.5686e-01,           &
     0.7302e-01,  0.9346e-01,  0.1192e+00,  0.1512e+00/
   data ((c2(ip,iw,1),iw=1,24),ip= 5, 5)/                                      &
    -0.2854e-06, -0.6241e-06, -0.1260e-05, -0.2305e-05, -0.3878e-05,           &
    -0.6181e-05, -0.9435e-05, -0.1375e-04, -0.1922e-04, -0.2613e-04,           &
    -0.3510e-04, -0.4696e-04, -0.6265e-04, -0.8334e-04, -0.1102e-03,           &
    -0.1443e-03, -0.1870e-03, -0.2393e-03, -0.3025e-03, -0.3775e-03,           &
    -0.4656e-03, -0.5678e-03, -0.6835e-03, -0.8098e-03/
   data ((c3(ip,iw,1),iw=1,24),ip= 5, 5)/                                      &
    -0.2485e-09, -0.5672e-09, -0.1468e-08, -0.3577e-08, -0.7395e-08,           &
    -0.1306e-07, -0.2029e-07, -0.2923e-07, -0.4213e-07, -0.6307e-07,           &
    -0.9429e-07, -0.1351e-06, -0.1872e-06, -0.2557e-06, -0.3442e-06,           &
    -0.4553e-06, -0.5899e-06, -0.7447e-06, -0.9163e-06, -0.1105e-05,           &
    -0.1311e-05, -0.1522e-05, -0.1697e-05, -0.1766e-05/
   data ((c1(ip,iw,1),iw=1,24),ip= 6, 6)/                                      &
     0.1450e-03,  0.2394e-03,  0.3680e-03,  0.5320e-03,  0.7452e-03,           &
     0.1040e-02,  0.1462e-02,  0.2058e-02,  0.2863e-02,  0.3914e-02,           &
     0.5275e-02,  0.7051e-02,  0.9378e-02,  0.1243e-01,  0.1641e-01,           &
     0.2162e-01,  0.2838e-01,  0.3713e-01,  0.4839e-01,  0.6282e-01,           &
     0.8125e-01,  0.1046e+00,  0.1340e+00,  0.1702e+00/
   data ((c2(ip,iw,1),iw=1,24),ip= 6, 6)/                                      &
    -0.2837e-06, -0.6203e-06, -0.1252e-05, -0.2292e-05, -0.3857e-05,           &
    -0.6146e-05, -0.9377e-05, -0.1366e-04, -0.1909e-04, -0.2596e-04,           &
    -0.3489e-04, -0.4676e-04, -0.6263e-04, -0.8372e-04, -0.1113e-03,           &
    -0.1465e-03, -0.1908e-03, -0.2454e-03, -0.3115e-03, -0.3903e-03,           &
    -0.4833e-03, -0.5915e-03, -0.7144e-03, -0.8483e-03/
   data ((c3(ip,iw,1),iw=1,24),ip= 6, 6)/                                      &
    -0.2557e-09, -0.5835e-09, -0.1495e-08, -0.3619e-08, -0.7472e-08,           &
    -0.1320e-07, -0.2055e-07, -0.2970e-07, -0.4281e-07, -0.6379e-07,           &
    -0.9492e-07, -0.1356e-06, -0.1876e-06, -0.2561e-06, -0.3451e-06,           &
    -0.4568e-06, -0.5924e-06, -0.7489e-06, -0.9229e-06, -0.1114e-05,           &
    -0.1321e-05, -0.1531e-05, -0.1700e-05, -0.1756e-05/
   data ((c1(ip,iw,1),iw=1,24),ip= 7, 7)/                                      &
     0.1454e-03,  0.2404e-03,  0.3704e-03,  0.5368e-03,  0.7542e-03,           &
     0.1056e-02,  0.1491e-02,  0.2108e-02,  0.2951e-02,  0.4066e-02,           &
     0.5529e-02,  0.7457e-02,  0.1000e-01,  0.1336e-01,  0.1777e-01,           &
     0.2357e-01,  0.3115e-01,  0.4102e-01,  0.5382e-01,  0.7035e-01,           &
     0.9154e-01,  0.1184e+00,  0.1519e+00,  0.1926e+00/
   data ((c2(ip,iw,1),iw=1,24),ip= 7, 7)/                                      &
    -0.2811e-06, -0.6147e-06, -0.1241e-05, -0.2272e-05, -0.3825e-05,           &
    -0.6093e-05, -0.9294e-05, -0.1354e-04, -0.1893e-04, -0.2577e-04,           &
    -0.3472e-04, -0.4672e-04, -0.6292e-04, -0.8464e-04, -0.1132e-03,           &
    -0.1500e-03, -0.1964e-03, -0.2539e-03, -0.3237e-03, -0.4074e-03,           &
    -0.5065e-03, -0.6222e-03, -0.7532e-03, -0.8955e-03/
   data ((c3(ip,iw,1),iw=1,24),ip= 7, 7)/                                      &
    -0.2673e-09, -0.6022e-09, -0.1531e-08, -0.3682e-08, -0.7588e-08,           &
    -0.1342e-07, -0.2092e-07, -0.3030e-07, -0.4357e-07, -0.6448e-07,           &
    -0.9533e-07, -0.1359e-06, -0.1879e-06, -0.2570e-06, -0.3470e-06,           &
    -0.4598e-06, -0.5971e-06, -0.7563e-06, -0.9334e-06, -0.1127e-05,           &
    -0.1336e-05, -0.1543e-05, -0.1701e-05, -0.1738e-05/
   data ((c1(ip,iw,1),iw=1,24),ip= 8, 8)/                                      &
     0.1460e-03,  0.2420e-03,  0.3740e-03,  0.5443e-03,  0.7681e-03,           &
     0.1081e-02,  0.1535e-02,  0.2184e-02,  0.3082e-02,  0.4286e-02,           &
     0.5885e-02,  0.8011e-02,  0.1083e-01,  0.1458e-01,  0.1954e-01,           &
     0.2608e-01,  0.3468e-01,  0.4595e-01,  0.6068e-01,  0.7979e-01,           &
     0.1043e+00,  0.1351e+00,  0.1731e+00,  0.2184e+00/
   data ((c2(ip,iw,1),iw=1,24),ip= 8, 8)/                                      &
    -0.2770e-06, -0.6058e-06, -0.1224e-05, -0.2243e-05, -0.3777e-05,           &
    -0.6018e-05, -0.9181e-05, -0.1339e-04, -0.1876e-04, -0.2561e-04,           &
    -0.3467e-04, -0.4695e-04, -0.6370e-04, -0.8635e-04, -0.1163e-03,           &
    -0.1551e-03, -0.2043e-03, -0.2654e-03, -0.3400e-03, -0.4299e-03,           &
    -0.5366e-03, -0.6609e-03, -0.8010e-03, -0.9521e-03/
   data ((c3(ip,iw,1),iw=1,24),ip= 8, 8)/                                      &
    -0.2845e-09, -0.6367e-09, -0.1591e-08, -0.3779e-08, -0.7752e-08,           &
    -0.1370e-07, -0.2140e-07, -0.3101e-07, -0.4437e-07, -0.6498e-07,           &
    -0.9551e-07, -0.1360e-06, -0.1887e-06, -0.2590e-06, -0.3503e-06,           &
    -0.4650e-06, -0.6052e-06, -0.7681e-06, -0.9488e-06, -0.1145e-05,           &
    -0.1354e-05, -0.1555e-05, -0.1698e-05, -0.1701e-05/
   data ((c1(ip,iw,1),iw=1,24),ip= 9, 9)/                                      &
     0.1470e-03,  0.2444e-03,  0.3796e-03,  0.5557e-03,  0.7894e-03,           &
     0.1119e-02,  0.1600e-02,  0.2296e-02,  0.3271e-02,  0.4596e-02,           &
     0.6373e-02,  0.8752e-02,  0.1193e-01,  0.1617e-01,  0.2181e-01,           &
     0.2928e-01,  0.3916e-01,  0.5220e-01,  0.6930e-01,  0.9150e-01,           &
     0.1198e+00,  0.1551e+00,  0.1978e+00,  0.2475e+00/
   data ((c2(ip,iw,1),iw=1,24),ip= 9, 9)/                                      &
    -0.2709e-06, -0.5925e-06, -0.1198e-05, -0.2199e-05, -0.3708e-05,           &
    -0.5913e-05, -0.9036e-05, -0.1322e-04, -0.1861e-04, -0.2556e-04,           &
    -0.3486e-04, -0.4760e-04, -0.6516e-04, -0.8909e-04, -0.1210e-03,           &
    -0.1624e-03, -0.2151e-03, -0.2808e-03, -0.3615e-03, -0.4590e-03,           &
    -0.5747e-03, -0.7088e-03, -0.8585e-03, -0.1018e-02/
   data ((c3(ip,iw,1),iw=1,24),ip= 9, 9)/                                      &
    -0.3119e-09, -0.6873e-09, -0.1676e-08, -0.3915e-08, -0.7973e-08,           &
    -0.1407e-07, -0.2200e-07, -0.3179e-07, -0.4505e-07, -0.6522e-07,           &
    -0.9551e-07, -0.1364e-06, -0.1902e-06, -0.2622e-06, -0.3557e-06,           &
    -0.4735e-06, -0.6178e-06, -0.7853e-06, -0.9700e-06, -0.1169e-05,           &
    -0.1377e-05, -0.1567e-05, -0.1683e-05, -0.1631e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=10,10)/                                      &
     0.1484e-03,  0.2482e-03,  0.3882e-03,  0.5732e-03,  0.8217e-03,           &
     0.1176e-02,  0.1696e-02,  0.2457e-02,  0.3536e-02,  0.5021e-02,           &
     0.7028e-02,  0.9728e-02,  0.1336e-01,  0.1821e-01,  0.2470e-01,           &
     0.3334e-01,  0.4482e-01,  0.6004e-01,  0.8003e-01,  0.1059e+00,           &
     0.1385e+00,  0.1785e+00,  0.2259e+00,  0.2794e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=10,10)/                                      &
    -0.2620e-06, -0.5728e-06, -0.1161e-05, -0.2136e-05, -0.3612e-05,           &
    -0.5777e-05, -0.8867e-05, -0.1306e-04, -0.1854e-04, -0.2572e-04,           &
    -0.3543e-04, -0.4887e-04, -0.6756e-04, -0.9321e-04, -0.1275e-03,           &
    -0.1723e-03, -0.2295e-03, -0.3011e-03, -0.3894e-03, -0.4960e-03,           &
    -0.6220e-03, -0.7665e-03, -0.9259e-03, -0.1095e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=10,10)/                                      &
    -0.3497e-09, -0.7575e-09, -0.1791e-08, -0.4091e-08, -0.8254e-08,           &
    -0.1452e-07, -0.2265e-07, -0.3255e-07, -0.4557e-07, -0.6530e-07,           &
    -0.9557e-07, -0.1374e-06, -0.1930e-06, -0.2675e-06, -0.3644e-06,           &
    -0.4866e-06, -0.6362e-06, -0.8090e-06, -0.9981e-06, -0.1199e-05,           &
    -0.1402e-05, -0.1574e-05, -0.1643e-05, -0.1506e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=11,11)/                                      &
     0.1506e-03,  0.2538e-03,  0.4010e-03,  0.5992e-03,  0.8693e-03,           &
     0.1258e-02,  0.1834e-02,  0.2683e-02,  0.3901e-02,  0.5591e-02,           &
     0.7890e-02,  0.1100e-01,  0.1519e-01,  0.2082e-01,  0.2837e-01,           &
     0.3845e-01,  0.5192e-01,  0.6981e-01,  0.9323e-01,  0.1232e+00,           &
     0.1606e+00,  0.2055e+00,  0.2571e+00,  0.3133e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=11,11)/                                      &
    -0.2494e-06, -0.5450e-06, -0.1108e-05, -0.2049e-05, -0.3485e-05,           &
    -0.5612e-05, -0.8695e-05, -0.1296e-04, -0.1865e-04, -0.2620e-04,           &
    -0.3654e-04, -0.5099e-04, -0.7121e-04, -0.9909e-04, -0.1365e-03,           &
    -0.1855e-03, -0.2483e-03, -0.3274e-03, -0.4248e-03, -0.5421e-03,           &
    -0.6794e-03, -0.8346e-03, -0.1003e-02, -0.1181e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=11,11)/                                      &
    -0.3973e-09, -0.8499e-09, -0.1944e-08, -0.4319e-08, -0.8596e-08,           &
    -0.1501e-07, -0.2332e-07, -0.3326e-07, -0.4600e-07, -0.6541e-07,           &
    -0.9592e-07, -0.1392e-06, -0.1975e-06, -0.2757e-06, -0.3774e-06,           &
    -0.5057e-06, -0.6614e-06, -0.8401e-06, -0.1034e-05, -0.1235e-05,           &
    -0.1427e-05, -0.1565e-05, -0.1560e-05, -0.1310e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=12,12)/                                      &
     0.1537e-03,  0.2618e-03,  0.4196e-03,  0.6368e-03,  0.9377e-03,           &
     0.1374e-02,  0.2024e-02,  0.2990e-02,  0.4387e-02,  0.6339e-02,           &
     0.9008e-02,  0.1263e-01,  0.1752e-01,  0.2410e-01,  0.3296e-01,           &
     0.4485e-01,  0.6077e-01,  0.8185e-01,  0.1092e+00,  0.1439e+00,           &
     0.1861e+00,  0.2356e+00,  0.2906e+00,  0.3483e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=12,12)/                                      &
    -0.2326e-06, -0.5076e-06, -0.1037e-05, -0.1935e-05, -0.3329e-05,           &
    -0.5433e-05, -0.8553e-05, -0.1298e-04, -0.1902e-04, -0.2716e-04,           &
    -0.3840e-04, -0.5421e-04, -0.7643e-04, -0.1071e-03, -0.1485e-03,           &
    -0.2028e-03, -0.2727e-03, -0.3608e-03, -0.4689e-03, -0.5981e-03,           &
    -0.7471e-03, -0.9128e-03, -0.1091e-02, -0.1276e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=12,12)/                                      &
    -0.4576e-09, -0.9621e-09, -0.2122e-08, -0.4569e-08, -0.8941e-08,           &
    -0.1550e-07, -0.2398e-07, -0.3404e-07, -0.4663e-07, -0.6577e-07,           &
    -0.9677e-07, -0.1421e-06, -0.2042e-06, -0.2878e-06, -0.3962e-06,           &
    -0.5318e-06, -0.6946e-06, -0.8798e-06, -0.1078e-05, -0.1275e-05,           &
    -0.1445e-05, -0.1523e-05, -0.1413e-05, -0.1030e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=13,13)/                                      &
     0.1579e-03,  0.2729e-03,  0.4454e-03,  0.6893e-03,  0.1033e-02,           &
     0.1533e-02,  0.2280e-02,  0.3396e-02,  0.5018e-02,  0.7299e-02,           &
     0.1043e-01,  0.1468e-01,  0.2044e-01,  0.2820e-01,  0.3867e-01,           &
     0.5278e-01,  0.7166e-01,  0.9648e-01,  0.1283e+00,  0.1679e+00,           &
     0.2150e+00,  0.2685e+00,  0.3256e+00,  0.3836e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=13,13)/                                      &
    -0.2120e-06, -0.4613e-06, -0.9487e-06, -0.1798e-05, -0.3153e-05,           &
    -0.5260e-05, -0.8480e-05, -0.1319e-04, -0.1976e-04, -0.2876e-04,           &
    -0.4125e-04, -0.5885e-04, -0.8357e-04, -0.1178e-03, -0.1640e-03,           &
    -0.2250e-03, -0.3036e-03, -0.4023e-03, -0.5226e-03, -0.6643e-03,           &
    -0.8250e-03, -0.1001e-02, -0.1187e-02, -0.1378e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=13,13)/                                      &
    -0.5223e-09, -0.1082e-08, -0.2314e-08, -0.4826e-08, -0.9255e-08,           &
    -0.1596e-07, -0.2472e-07, -0.3507e-07, -0.4774e-07, -0.6677e-07,           &
    -0.9841e-07, -0.1466e-06, -0.2141e-06, -0.3049e-06, -0.4219e-06,           &
    -0.5662e-06, -0.7370e-06, -0.9290e-06, -0.1129e-05, -0.1315e-05,           &
    -0.1441e-05, -0.1427e-05, -0.1184e-05, -0.6731e-06/
   data ((c1(ip,iw,1),iw=1,24),ip=14,14)/                                      &
     0.1634e-03,  0.2873e-03,  0.4793e-03,  0.7589e-03,  0.1159e-02,           &
     0.1742e-02,  0.2614e-02,  0.3916e-02,  0.5817e-02,  0.8503e-02,           &
     0.1220e-01,  0.1723e-01,  0.2403e-01,  0.3324e-01,  0.4572e-01,           &
     0.6253e-01,  0.8488e-01,  0.1140e+00,  0.1507e+00,  0.1952e+00,           &
     0.2468e+00,  0.3031e+00,  0.3612e+00,  0.4186e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=14,14)/                                      &
    -0.1896e-06, -0.4100e-06, -0.8500e-06, -0.1645e-05, -0.2972e-05,           &
    -0.5116e-05, -0.8508e-05, -0.1363e-04, -0.2097e-04, -0.3115e-04,           &
    -0.4531e-04, -0.6517e-04, -0.9298e-04, -0.1315e-03, -0.1837e-03,           &
    -0.2528e-03, -0.3417e-03, -0.4525e-03, -0.5860e-03, -0.7405e-03,           &
    -0.9123e-03, -0.1097e-02, -0.1290e-02, -0.1482e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=14,14)/                                      &
    -0.5738e-09, -0.1185e-08, -0.2465e-08, -0.5015e-08, -0.9502e-08,           &
    -0.1640e-07, -0.2567e-07, -0.3675e-07, -0.4990e-07, -0.6893e-07,           &
    -0.1014e-06, -0.1533e-06, -0.2278e-06, -0.3278e-06, -0.4553e-06,           &
    -0.6098e-06, -0.7899e-06, -0.9876e-06, -0.1184e-05, -0.1342e-05,           &
    -0.1393e-05, -0.1254e-05, -0.8703e-06, -0.2593e-06/
   data ((c1(ip,iw,1),iw=1,24),ip=15,15)/                                      &
     0.1697e-03,  0.3043e-03,  0.5205e-03,  0.8458e-03,  0.1319e-02,           &
     0.2009e-02,  0.3034e-02,  0.4563e-02,  0.6800e-02,  0.9975e-02,           &
     0.1436e-01,  0.2033e-01,  0.2842e-01,  0.3941e-01,  0.5431e-01,           &
     0.7432e-01,  0.1007e+00,  0.1344e+00,  0.1762e+00,  0.2255e+00,           &
     0.2806e+00,  0.3385e+00,  0.3964e+00,  0.4532e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=15,15)/                                      &
    -0.1689e-06, -0.3610e-06, -0.7531e-06, -0.1495e-05, -0.2803e-05,           &
    -0.5013e-05, -0.8646e-05, -0.1433e-04, -0.2270e-04, -0.3446e-04,           &
    -0.5077e-04, -0.7342e-04, -0.1049e-03, -0.1486e-03, -0.2080e-03,           &
    -0.2866e-03, -0.3873e-03, -0.5115e-03, -0.6586e-03, -0.8256e-03,           &
    -0.1008e-02, -0.1201e-02, -0.1398e-02, -0.1583e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=15,15)/                                      &
    -0.6063e-09, -0.1245e-08, -0.2559e-08, -0.5110e-08, -0.9652e-08,           &
    -0.1685e-07, -0.2693e-07, -0.3933e-07, -0.5373e-07, -0.7330e-07,           &
    -0.1066e-06, -0.1625e-06, -0.2453e-06, -0.3564e-06, -0.4963e-06,           &
    -0.6636e-06, -0.8532e-06, -0.1054e-05, -0.1236e-05, -0.1339e-05,           &
    -0.1280e-05, -0.9947e-06, -0.4862e-06,  0.1889e-06/
   data ((c1(ip,iw,1),iw=1,24),ip=16,16)/                                      &
     0.1763e-03,  0.3227e-03,  0.5664e-03,  0.9465e-03,  0.1511e-02,           &
     0.2335e-02,  0.3549e-02,  0.5348e-02,  0.7977e-02,  0.1172e-01,           &
     0.1692e-01,  0.2402e-01,  0.3369e-01,  0.4685e-01,  0.6465e-01,           &
     0.8835e-01,  0.1192e+00,  0.1579e+00,  0.2046e+00,  0.2579e+00,           &
     0.3153e+00,  0.3736e+00,  0.4308e+00,  0.4875e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=16,16)/                                      &
    -0.1533e-06, -0.3224e-06, -0.6725e-06, -0.1365e-05, -0.2657e-05,           &
    -0.4953e-05, -0.8878e-05, -0.1525e-04, -0.2490e-04, -0.3866e-04,           &
    -0.5766e-04, -0.8372e-04, -0.1197e-03, -0.1695e-03, -0.2372e-03,           &
    -0.3266e-03, -0.4401e-03, -0.5783e-03, -0.7390e-03, -0.9180e-03,           &
    -0.1110e-02, -0.1309e-02, -0.1503e-02, -0.1676e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=16,16)/                                      &
    -0.6109e-09, -0.1253e-08, -0.2555e-08, -0.5100e-08, -0.9703e-08,           &
    -0.1729e-07, -0.2846e-07, -0.4291e-07, -0.5988e-07, -0.8135e-07,           &
    -0.1158e-06, -0.1750e-06, -0.2657e-06, -0.3890e-06, -0.5438e-06,           &
    -0.7262e-06, -0.9256e-06, -0.1122e-05, -0.1269e-05, -0.1284e-05,           &
    -0.1086e-05, -0.6553e-06, -0.5588e-07,  0.6660e-06/
   data ((c1(ip,iw,1),iw=1,24),ip=17,17)/                                      &
     0.1824e-03,  0.3403e-03,  0.6126e-03,  0.1054e-02,  0.1727e-02,           &
     0.2716e-02,  0.4160e-02,  0.6275e-02,  0.9353e-02,  0.1375e-01,           &
     0.1990e-01,  0.2837e-01,  0.3995e-01,  0.5571e-01,  0.7690e-01,           &
     0.1048e+00,  0.1404e+00,  0.1841e+00,  0.2353e+00,  0.2916e+00,           &
     0.3500e+00,  0.4078e+00,  0.4646e+00,  0.5215e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=17,17)/                                      &
    -0.1436e-06, -0.2978e-06, -0.6179e-06, -0.1270e-05, -0.2546e-05,           &
    -0.4927e-05, -0.9164e-05, -0.1628e-04, -0.2740e-04, -0.4353e-04,           &
    -0.6579e-04, -0.9606e-04, -0.1375e-03, -0.1943e-03, -0.2712e-03,           &
    -0.3722e-03, -0.4993e-03, -0.6515e-03, -0.8255e-03, -0.1016e-02,           &
    -0.1217e-02, -0.1417e-02, -0.1602e-02, -0.1757e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=17,17)/                                      &
    -0.6015e-09, -0.1225e-08, -0.2489e-08, -0.4978e-08, -0.9608e-08,           &
    -0.1756e-07, -0.2993e-07, -0.4715e-07, -0.6853e-07, -0.9454e-07,           &
    -0.1314e-06, -0.1920e-06, -0.2877e-06, -0.4222e-06, -0.5940e-06,           &
    -0.7936e-06, -0.1002e-05, -0.1181e-05, -0.1262e-05, -0.1156e-05,           &
    -0.8086e-06, -0.2604e-06,  0.4056e-06,  0.1164e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=18,18)/                                      &
     0.1875e-03,  0.3555e-03,  0.6547e-03,  0.1159e-02,  0.1954e-02,           &
     0.3140e-02,  0.4861e-02,  0.7345e-02,  0.1093e-01,  0.1606e-01,           &
     0.2330e-01,  0.3339e-01,  0.4729e-01,  0.6613e-01,  0.9118e-01,           &
     0.1236e+00,  0.1642e+00,  0.2127e+00,  0.2674e+00,  0.3254e+00,           &
     0.3836e+00,  0.4408e+00,  0.4978e+00,  0.5551e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=18,18)/                                      &
    -0.1390e-06, -0.2850e-06, -0.5878e-06, -0.1213e-05, -0.2477e-05,           &
    -0.4934e-05, -0.9466e-05, -0.1731e-04, -0.2990e-04, -0.4860e-04,           &
    -0.7471e-04, -0.1102e-03, -0.1582e-03, -0.2230e-03, -0.3098e-03,           &
    -0.4228e-03, -0.5632e-03, -0.7291e-03, -0.9158e-03, -0.1117e-02,           &
    -0.1323e-02, -0.1519e-02, -0.1689e-02, -0.1824e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=18,18)/                                      &
    -0.5875e-09, -0.1187e-08, -0.2392e-08, -0.4777e-08, -0.9321e-08,           &
    -0.1745e-07, -0.3091e-07, -0.5121e-07, -0.7879e-07, -0.1130e-06,           &
    -0.1555e-06, -0.2161e-06, -0.3110e-06, -0.4527e-06, -0.6400e-06,           &
    -0.8571e-06, -0.1069e-05, -0.1209e-05, -0.1193e-05, -0.9472e-06,           &
    -0.4712e-06,  0.1597e-06,  0.8892e-06,  0.1671e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=19,19)/                                      &
     0.1914e-03,  0.3677e-03,  0.6900e-03,  0.1252e-02,  0.2172e-02,           &
     0.3581e-02,  0.5632e-02,  0.8553e-02,  0.1271e-01,  0.1866e-01,           &
     0.2715e-01,  0.3914e-01,  0.5576e-01,  0.7817e-01,  0.1076e+00,           &
     0.1449e+00,  0.1903e+00,  0.2429e+00,  0.2999e+00,  0.3582e+00,           &
     0.4157e+00,  0.4727e+00,  0.5303e+00,  0.5878e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=19,19)/                                      &
    -0.1371e-06, -0.2799e-06, -0.5756e-06, -0.1190e-05, -0.2455e-05,           &
    -0.4982e-05, -0.9770e-05, -0.1824e-04, -0.3214e-04, -0.5332e-04,           &
    -0.8371e-04, -0.1255e-03, -0.1814e-03, -0.2555e-03, -0.3528e-03,           &
    -0.4774e-03, -0.6306e-03, -0.8095e-03, -0.1008e-02, -0.1218e-02,           &
    -0.1424e-02, -0.1610e-02, -0.1763e-02, -0.1875e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=19,19)/                                      &
    -0.5739e-09, -0.1144e-08, -0.2283e-08, -0.4524e-08, -0.8848e-08,           &
    -0.1686e-07, -0.3096e-07, -0.5411e-07, -0.8876e-07, -0.1343e-06,           &
    -0.1880e-06, -0.2511e-06, -0.3405e-06, -0.4800e-06, -0.6756e-06,           &
    -0.9041e-06, -0.1105e-05, -0.1184e-05, -0.1049e-05, -0.6701e-06,           &
    -0.1061e-06,  0.5899e-06,  0.1376e-05,  0.2160e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=20,20)/                                      &
     0.1942e-03,  0.3767e-03,  0.7173e-03,  0.1328e-02,  0.2364e-02,           &
     0.4005e-02,  0.6433e-02,  0.9873e-02,  0.1471e-01,  0.2160e-01,           &
     0.3150e-01,  0.4565e-01,  0.6537e-01,  0.9183e-01,  0.1260e+00,           &
     0.1683e+00,  0.2183e+00,  0.2738e+00,  0.3318e+00,  0.3894e+00,           &
     0.4463e+00,  0.5038e+00,  0.5619e+00,  0.6191e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=20,20)/                                      &
    -0.1368e-06, -0.2791e-06, -0.5744e-06, -0.1192e-05, -0.2475e-05,           &
    -0.5078e-05, -0.1008e-04, -0.1903e-04, -0.3392e-04, -0.5726e-04,           &
    -0.9202e-04, -0.1409e-03, -0.2063e-03, -0.2912e-03, -0.3997e-03,           &
    -0.5357e-03, -0.7004e-03, -0.8911e-03, -0.1100e-02, -0.1314e-02,           &
    -0.1514e-02, -0.1687e-02, -0.1822e-02, -0.1908e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=20,20)/                                      &
    -0.5617e-09, -0.1109e-08, -0.2184e-08, -0.4264e-08, -0.8251e-08,           &
    -0.1582e-07, -0.3000e-07, -0.5536e-07, -0.9656e-07, -0.1549e-06,           &
    -0.2254e-06, -0.2999e-06, -0.3855e-06, -0.5114e-06, -0.6982e-06,           &
    -0.9171e-06, -0.1082e-05, -0.1078e-05, -0.8280e-06, -0.3561e-06,           &
     0.2636e-06,  0.1015e-05,  0.1836e-05,  0.2588e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=21,21)/                                      &
     0.1962e-03,  0.3831e-03,  0.7370e-03,  0.1386e-02,  0.2520e-02,           &
     0.4379e-02,  0.7206e-02,  0.1125e-01,  0.1692e-01,  0.2493e-01,           &
     0.3644e-01,  0.5299e-01,  0.7611e-01,  0.1070e+00,  0.1463e+00,           &
     0.1936e+00,  0.2474e+00,  0.3046e+00,  0.3621e+00,  0.4188e+00,           &
     0.4758e+00,  0.5339e+00,  0.5921e+00,  0.6482e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=21,21)/                                      &
    -0.1371e-06, -0.2803e-06, -0.5786e-06, -0.1206e-05, -0.2522e-05,           &
    -0.5211e-05, -0.1040e-04, -0.1969e-04, -0.3522e-04, -0.6021e-04,           &
    -0.9890e-04, -0.1552e-03, -0.2312e-03, -0.3285e-03, -0.4496e-03,           &
    -0.5973e-03, -0.7728e-03, -0.9733e-03, -0.1188e-02, -0.1399e-02,           &
    -0.1590e-02, -0.1748e-02, -0.1862e-02, -0.1922e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=21,21)/                                      &
    -0.5521e-09, -0.1082e-08, -0.2100e-08, -0.4025e-08, -0.7661e-08,           &
    -0.1462e-07, -0.2837e-07, -0.5493e-07, -0.1013e-06, -0.1716e-06,           &
    -0.2627e-06, -0.3609e-06, -0.4548e-06, -0.5617e-06, -0.7128e-06,           &
    -0.8814e-06, -0.9684e-06, -0.8687e-06, -0.5417e-06, -0.3260e-07,           &
     0.6298e-06,  0.1417e-05,  0.2233e-05,  0.2884e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=22,22)/                                      &
     0.1974e-03,  0.3872e-03,  0.7505e-03,  0.1427e-02,  0.2638e-02,           &
     0.4682e-02,  0.7896e-02,  0.1261e-01,  0.1927e-01,  0.2865e-01,           &
     0.4207e-01,  0.6125e-01,  0.8802e-01,  0.1236e+00,  0.1681e+00,           &
     0.2201e+00,  0.2766e+00,  0.3342e+00,  0.3906e+00,  0.4468e+00,           &
     0.5043e+00,  0.5628e+00,  0.6204e+00,  0.6743e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=22,22)/                                      &
    -0.1377e-06, -0.2818e-06, -0.5836e-06, -0.1222e-05, -0.2570e-05,           &
    -0.5339e-05, -0.1069e-04, -0.2019e-04, -0.3606e-04, -0.6211e-04,           &
    -0.1038e-03, -0.1666e-03, -0.2537e-03, -0.3653e-03, -0.5012e-03,           &
    -0.6623e-03, -0.8484e-03, -0.1055e-02, -0.1268e-02, -0.1470e-02,           &
    -0.1649e-02, -0.1790e-02, -0.1882e-02, -0.1920e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=22,22)/                                      &
    -0.5464e-09, -0.1064e-08, -0.2050e-08, -0.3874e-08, -0.7265e-08,           &
    -0.1373e-07, -0.2688e-07, -0.5376e-07, -0.1035e-06, -0.1833e-06,           &
    -0.2955e-06, -0.4267e-06, -0.5474e-06, -0.6443e-06, -0.7332e-06,           &
    -0.7933e-06, -0.7465e-06, -0.5424e-06, -0.1941e-06,  0.3078e-06,           &
     0.9845e-06,  0.1771e-05,  0.2515e-05,  0.2962e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=23,23)/                                      &
     0.1981e-03,  0.3899e-03,  0.7594e-03,  0.1455e-02,  0.2722e-02,           &
     0.4913e-02,  0.8469e-02,  0.1387e-01,  0.2168e-01,  0.3275e-01,           &
     0.4842e-01,  0.7058e-01,  0.1012e+00,  0.1414e+00,  0.1908e+00,           &
     0.2467e+00,  0.3050e+00,  0.3620e+00,  0.4176e+00,  0.4739e+00,           &
     0.5318e+00,  0.5900e+00,  0.6460e+00,  0.6969e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=23,23)/                                      &
    -0.1374e-06, -0.2821e-06, -0.5849e-06, -0.1227e-05, -0.2587e-05,           &
    -0.5392e-05, -0.1082e-04, -0.2043e-04, -0.3640e-04, -0.6291e-04,           &
    -0.1064e-03, -0.1742e-03, -0.2713e-03, -0.3983e-03, -0.5523e-03,           &
    -0.7302e-03, -0.9280e-03, -0.1135e-02, -0.1338e-02, -0.1526e-02,           &
    -0.1690e-02, -0.1811e-02, -0.1880e-02, -0.1904e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=23,23)/                                      &
    -0.5505e-09, -0.1066e-08, -0.2044e-08, -0.3856e-08, -0.7199e-08,           &
    -0.1353e-07, -0.2638e-07, -0.5316e-07, -0.1044e-06, -0.1908e-06,           &
    -0.3207e-06, -0.4858e-06, -0.6471e-06, -0.7512e-06, -0.7676e-06,           &
    -0.6661e-06, -0.4252e-06, -0.1067e-06,  0.2334e-06,  0.6913e-06,           &
     0.1330e-05,  0.2055e-05,  0.2617e-05,  0.2760e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=24,24)/                                      &
     0.1986e-03,  0.3917e-03,  0.7654e-03,  0.1475e-02,  0.2782e-02,           &
     0.5086e-02,  0.8931e-02,  0.1497e-01,  0.2401e-01,  0.3707e-01,           &
     0.5549e-01,  0.8107e-01,  0.1157e+00,  0.1604e+00,  0.2138e+00,           &
     0.2723e+00,  0.3311e+00,  0.3877e+00,  0.4433e+00,  0.5001e+00,           &
     0.5580e+00,  0.6152e+00,  0.6683e+00,  0.7158e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=24,24)/                                      &
    -0.1372e-06, -0.2802e-06, -0.5797e-06, -0.1212e-05, -0.2549e-05,           &
    -0.5312e-05, -0.1068e-04, -0.2023e-04, -0.3609e-04, -0.6255e-04,           &
    -0.1068e-03, -0.1777e-03, -0.2830e-03, -0.4252e-03, -0.6004e-03,           &
    -0.8005e-03, -0.1011e-02, -0.1213e-02, -0.1399e-02, -0.1569e-02,           &
    -0.1712e-02, -0.1809e-02, -0.1858e-02, -0.1883e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=24,24)/                                      &
    -0.5589e-09, -0.1071e-08, -0.2083e-08, -0.3965e-08, -0.7465e-08,           &
    -0.1408e-07, -0.2720e-07, -0.5414e-07, -0.1061e-06, -0.1958e-06,           &
    -0.3365e-06, -0.5263e-06, -0.7236e-06, -0.8420e-06, -0.7920e-06,           &
    -0.5166e-06, -0.5606e-07,  0.4032e-06,  0.7697e-06,  0.1160e-05,           &
     0.1692e-05,  0.2251e-05,  0.2483e-05,  0.2295e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=25,25)/                                      &
     0.1990e-03,  0.3929e-03,  0.7699e-03,  0.1490e-02,  0.2829e-02,           &
     0.5226e-02,  0.9314e-02,  0.1593e-01,  0.2618e-01,  0.4144e-01,           &
     0.6312e-01,  0.9271e-01,  0.1317e+00,  0.1803e+00,  0.2364e+00,           &
     0.2956e+00,  0.3541e+00,  0.4111e+00,  0.4680e+00,  0.5255e+00,           &
     0.5829e+00,  0.6376e+00,  0.6869e+00,  0.7311e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=25,25)/                                      &
    -0.1366e-06, -0.2772e-06, -0.5687e-06, -0.1179e-05, -0.2458e-05,           &
    -0.5093e-05, -0.1024e-04, -0.1950e-04, -0.3503e-04, -0.6107e-04,           &
    -0.1053e-03, -0.1782e-03, -0.2902e-03, -0.4470e-03, -0.6455e-03,           &
    -0.8715e-03, -0.1095e-02, -0.1290e-02, -0.1458e-02, -0.1604e-02,           &
    -0.1716e-02, -0.1783e-02, -0.1821e-02, -0.1862e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=25,25)/                                      &
    -0.5574e-09, -0.1099e-08, -0.2136e-08, -0.4136e-08, -0.7893e-08,           &
    -0.1505e-07, -0.2894e-07, -0.5650e-07, -0.1084e-06, -0.1986e-06,           &
    -0.3416e-06, -0.5381e-06, -0.7439e-06, -0.8545e-06, -0.7482e-06,           &
    -0.3541e-06,  0.2762e-06,  0.9306e-06,  0.1416e-05,  0.1754e-05,           &
     0.2105e-05,  0.2343e-05,  0.2113e-05,  0.1689e-05/
   data ((c1(ip,iw,1),iw=1,24),ip=26,26)/                                      &
     0.1993e-03,  0.3940e-03,  0.7739e-03,  0.1503e-02,  0.2871e-02,           &
     0.5350e-02,  0.9655e-02,  0.1679e-01,  0.2817e-01,  0.4563e-01,           &
     0.7094e-01,  0.1052e+00,  0.1488e+00,  0.2006e+00,  0.2575e+00,           &
     0.3155e+00,  0.3732e+00,  0.4317e+00,  0.4909e+00,  0.5496e+00,           &
     0.6058e+00,  0.6569e+00,  0.7019e+00,  0.7436e+00/
   data ((c2(ip,iw,1),iw=1,24),ip=26,26)/                                      &
    -0.1351e-06, -0.2730e-06, -0.5547e-06, -0.1136e-05, -0.2337e-05,           &
    -0.4789e-05, -0.9589e-05, -0.1836e-04, -0.3337e-04, -0.5890e-04,           &
    -0.1031e-03, -0.1781e-03, -0.2974e-03, -0.4696e-03, -0.6912e-03,           &
    -0.9400e-03, -0.1174e-02, -0.1368e-02, -0.1521e-02, -0.1635e-02,           &
    -0.1706e-02, -0.1742e-02, -0.1779e-02, -0.1846e-02/
   data ((c3(ip,iw,1),iw=1,24),ip=26,26)/                                      &
    -0.5597e-09, -0.1126e-08, -0.2190e-08, -0.4298e-08, -0.8333e-08,           &
    -0.1610e-07, -0.3092e-07, -0.5912e-07, -0.1108e-06, -0.1982e-06,           &
    -0.3334e-06, -0.5129e-06, -0.6865e-06, -0.7484e-06, -0.5999e-06,           &
    -0.1961e-06,  0.5028e-06,  0.1395e-05,  0.2125e-05,  0.2489e-05,           &
     0.2584e-05,  0.2312e-05,  0.1587e-05,  0.1122e-05/
   data ((o1(ip,iw,1),iw=1,21),ip= 1, 1)/                                      &
     0.7117e-05,  0.1419e-04,  0.2828e-04,  0.5629e-04,  0.1117e-03,           &
     0.2207e-03,  0.4326e-03,  0.8367e-03,  0.1586e-02,  0.2911e-02,           &
     0.5106e-02,  0.8433e-02,  0.1302e-01,  0.1896e-01,  0.2658e-01,           &
     0.3657e-01,  0.4988e-01,  0.6765e-01,  0.9092e-01,  0.1204e+00,           &
     0.1559e+00/
   data ((o2(ip,iw,1),iw=1,21),ip= 1, 1)/                                      &
     0.6531e-10,  0.5926e-10, -0.1646e-09, -0.1454e-08, -0.7376e-08,           &
    -0.2968e-07, -0.1071e-06, -0.3584e-06, -0.1125e-05, -0.3289e-05,           &
    -0.8760e-05, -0.2070e-04, -0.4259e-04, -0.7691e-04, -0.1264e-03,           &
    -0.1957e-03, -0.2895e-03, -0.4107e-03, -0.5588e-03, -0.7300e-03,           &
    -0.9199e-03/
   data ((o3(ip,iw,1),iw=1,21),ip= 1, 1)/                                      &
    -0.2438e-10, -0.4826e-10, -0.9474e-10, -0.1828e-09, -0.3406e-09,           &
    -0.6223e-09, -0.1008e-08, -0.1412e-08, -0.1244e-08,  0.8485e-09,           &
     0.6343e-08,  0.1201e-07,  0.2838e-08, -0.4024e-07, -0.1257e-06,           &
    -0.2566e-06, -0.4298e-06, -0.6184e-06, -0.7657e-06, -0.8153e-06,           &
    -0.7552e-06/
   data ((o1(ip,iw,1),iw=1,21),ip= 2, 2)/                                      &
     0.7117e-05,  0.1419e-04,  0.2828e-04,  0.5629e-04,  0.1117e-03,           &
     0.2207e-03,  0.4326e-03,  0.8367e-03,  0.1586e-02,  0.2912e-02,           &
     0.5107e-02,  0.8435e-02,  0.1303e-01,  0.1897e-01,  0.2660e-01,           &
     0.3660e-01,  0.4995e-01,  0.6777e-01,  0.9114e-01,  0.1207e+00,           &
     0.1566e+00/
   data ((o2(ip,iw,1),iw=1,21),ip= 2, 2)/                                      &
     0.6193e-10,  0.5262e-10, -0.1774e-09, -0.1478e-08, -0.7416e-08,           &
    -0.2985e-07, -0.1071e-06, -0.3584e-06, -0.1124e-05, -0.3287e-05,           &
    -0.8753e-05, -0.2069e-04, -0.4256e-04, -0.7686e-04, -0.1264e-03,           &
    -0.1956e-03, -0.2893e-03, -0.4103e-03, -0.5580e-03, -0.7285e-03,           &
    -0.9171e-03/
   data ((o3(ip,iw,1),iw=1,21),ip= 2, 2)/                                      &
    -0.2436e-10, -0.4822e-10, -0.9466e-10, -0.1827e-09, -0.3404e-09,           &
    -0.6220e-09, -0.1008e-08, -0.1414e-08, -0.1247e-08,  0.8360e-09,           &
     0.6312e-08,  0.1194e-07,  0.2753e-08, -0.4040e-07, -0.1260e-06,           &
    -0.2571e-06, -0.4307e-06, -0.6202e-06, -0.7687e-06, -0.8204e-06,           &
    -0.7636e-06/
   data ((o1(ip,iw,1),iw=1,21),ip= 3, 3)/                                      &
     0.7117e-05,  0.1419e-04,  0.2828e-04,  0.5628e-04,  0.1117e-03,           &
     0.2207e-03,  0.4326e-03,  0.8367e-03,  0.1586e-02,  0.2912e-02,           &
     0.5109e-02,  0.8439e-02,  0.1303e-01,  0.1899e-01,  0.2664e-01,           &
     0.3666e-01,  0.5005e-01,  0.6795e-01,  0.9147e-01,  0.1213e+00,           &
     0.1576e+00/
   data ((o2(ip,iw,1),iw=1,21),ip= 3, 3)/                                      &
     0.5658e-10,  0.4212e-10, -0.1977e-09, -0.1516e-08, -0.7481e-08,           &
    -0.2995e-07, -0.1072e-06, -0.3583e-06, -0.1123e-05, -0.3283e-05,           &
    -0.8744e-05, -0.2067e-04, -0.4252e-04, -0.7679e-04, -0.1262e-03,           &
    -0.1953e-03, -0.2889e-03, -0.4096e-03, -0.5567e-03, -0.7263e-03,           &
    -0.9130e-03/
   data ((o3(ip,iw,1),iw=1,21),ip= 3, 3)/                                      &
    -0.2433e-10, -0.4815e-10, -0.9453e-10, -0.1825e-09, -0.3400e-09,           &
    -0.6215e-09, -0.1007e-08, -0.1415e-08, -0.1253e-08,  0.8143e-09,           &
     0.6269e-08,  0.1186e-07,  0.2604e-08, -0.4067e-07, -0.1264e-06,           &
    -0.2579e-06, -0.4321e-06, -0.6229e-06, -0.7732e-06, -0.8277e-06,           &
    -0.7752e-06/
   data ((o1(ip,iw,1),iw=1,21),ip= 4, 4)/                                      &
     0.7116e-05,  0.1419e-04,  0.2828e-04,  0.5628e-04,  0.1117e-03,           &
     0.2207e-03,  0.4326e-03,  0.8368e-03,  0.1586e-02,  0.2913e-02,           &
     0.5111e-02,  0.8444e-02,  0.1305e-01,  0.1902e-01,  0.2669e-01,           &
     0.3676e-01,  0.5022e-01,  0.6825e-01,  0.9199e-01,  0.1222e+00,           &
     0.1591e+00/
   data ((o2(ip,iw,1),iw=1,21),ip= 4, 4)/                                      &
     0.4814e-10,  0.2552e-10, -0.2298e-09, -0.1576e-08, -0.7579e-08,           &
    -0.3009e-07, -0.1074e-06, -0.3581e-06, -0.1122e-05, -0.3278e-05,           &
    -0.8729e-05, -0.2063e-04, -0.4245e-04, -0.7667e-04, -0.1260e-03,           &
    -0.1950e-03, -0.2883e-03, -0.4086e-03, -0.5549e-03, -0.7229e-03,           &
    -0.9071e-03/
   data ((o3(ip,iw,1),iw=1,21),ip= 4, 4)/                                      &
    -0.2428e-10, -0.4805e-10, -0.9433e-10, -0.1821e-09, -0.3394e-09,           &
    -0.6206e-09, -0.1008e-08, -0.1416e-08, -0.1261e-08,  0.7860e-09,           &
     0.6188e-08,  0.1171e-07,  0.2389e-08, -0.4109e-07, -0.1271e-06,           &
    -0.2591e-06, -0.4344e-06, -0.6267e-06, -0.7797e-06, -0.8378e-06,           &
    -0.7901e-06/
   data ((o1(ip,iw,1),iw=1,21),ip= 5, 5)/                                      &
     0.7116e-05,  0.1419e-04,  0.2827e-04,  0.5627e-04,  0.1117e-03,           &
     0.2207e-03,  0.4326e-03,  0.8368e-03,  0.1586e-02,  0.2914e-02,           &
     0.5114e-02,  0.8454e-02,  0.1307e-01,  0.1906e-01,  0.2677e-01,           &
     0.3690e-01,  0.5048e-01,  0.6872e-01,  0.9281e-01,  0.1236e+00,           &
     0.1615e+00/
   data ((o2(ip,iw,1),iw=1,21),ip= 5, 5)/                                      &
     0.3482e-10, -0.6492e-12, -0.2805e-09, -0.1671e-08, -0.7740e-08,           &
    -0.3032e-07, -0.1076e-06, -0.3582e-06, -0.1120e-05, -0.3270e-05,           &
    -0.8704e-05, -0.2058e-04, -0.4235e-04, -0.7649e-04, -0.1257e-03,           &
    -0.1945e-03, -0.2874e-03, -0.4070e-03, -0.5521e-03, -0.7181e-03,           &
    -0.8990e-03/
   data ((o3(ip,iw,1),iw=1,21),ip= 5, 5)/                                      &
    -0.2419e-10, -0.4788e-10, -0.9401e-10, -0.1815e-09, -0.3385e-09,           &
    -0.6192e-09, -0.1006e-08, -0.1417e-08, -0.1273e-08,  0.7404e-09,           &
     0.6068e-08,  0.1148e-07,  0.2021e-08, -0.4165e-07, -0.1281e-06,           &
    -0.2609e-06, -0.4375e-06, -0.6323e-06, -0.7887e-06, -0.8508e-06,           &
    -0.8067e-06/
   data ((o1(ip,iw,1),iw=1,21),ip= 6, 6)/                                      &
     0.7114e-05,  0.1419e-04,  0.2827e-04,  0.5627e-04,  0.1117e-03,           &
     0.2207e-03,  0.4325e-03,  0.8369e-03,  0.1587e-02,  0.2916e-02,           &
     0.5120e-02,  0.8468e-02,  0.1310e-01,  0.1913e-01,  0.2690e-01,           &
     0.3714e-01,  0.5090e-01,  0.6944e-01,  0.9407e-01,  0.1258e+00,           &
     0.1651e+00/
   data ((o2(ip,iw,1),iw=1,21),ip= 6, 6)/                                      &
     0.1388e-10, -0.4180e-10, -0.3601e-09, -0.1820e-08, -0.7993e-08,           &
    -0.3068e-07, -0.1081e-06, -0.3580e-06, -0.1117e-05, -0.3257e-05,           &
    -0.8667e-05, -0.2049e-04, -0.4218e-04, -0.7620e-04, -0.1253e-03,           &
    -0.1937e-03, -0.2860e-03, -0.4047e-03, -0.5481e-03, -0.7115e-03,           &
    -0.8885e-03/
   data ((o3(ip,iw,1),iw=1,21),ip= 6, 6)/                                      &
    -0.2406e-10, -0.4762e-10, -0.9351e-10, -0.1806e-09, -0.3370e-09,           &
    -0.6170e-09, -0.1004e-08, -0.1417e-08, -0.1297e-08,  0.6738e-09,           &
     0.5895e-08,  0.1113e-07,  0.1466e-08, -0.4265e-07, -0.1298e-06,           &
    -0.2636e-06, -0.4423e-06, -0.6402e-06, -0.8005e-06, -0.8658e-06,           &
    -0.8222e-06/
   data ((o1(ip,iw,1),iw=1,21),ip= 7, 7)/                                      &
     0.7113e-05,  0.1418e-04,  0.2826e-04,  0.5625e-04,  0.1117e-03,           &
     0.2206e-03,  0.4325e-03,  0.8371e-03,  0.1588e-02,  0.2918e-02,           &
     0.5128e-02,  0.8491e-02,  0.1315e-01,  0.1923e-01,  0.2710e-01,           &
     0.3750e-01,  0.5154e-01,  0.7056e-01,  0.9600e-01,  0.1290e+00,           &
     0.1703e+00/
   data ((o2(ip,iw,1),iw=1,21),ip= 7, 7)/                                      &
    -0.1889e-10, -0.1062e-09, -0.4847e-09, -0.2053e-08, -0.8389e-08,           &
    -0.3140e-07, -0.1089e-06, -0.3577e-06, -0.1112e-05, -0.3236e-05,           &
    -0.8607e-05, -0.2035e-04, -0.4192e-04, -0.7576e-04, -0.1245e-03,           &
    -0.1925e-03, -0.2840e-03, -0.4013e-03, -0.5427e-03, -0.7029e-03,           &
    -0.8756e-03/
   data ((o3(ip,iw,1),iw=1,21),ip= 7, 7)/                                      &
    -0.2385e-10, -0.4722e-10, -0.9273e-10, -0.1791e-09, -0.3348e-09,           &
    -0.6121e-09, -0.9974e-09, -0.1422e-08, -0.1326e-08,  0.5603e-09,           &
     0.5604e-08,  0.1061e-07,  0.6106e-09, -0.4398e-07, -0.1321e-06,           &
    -0.2676e-06, -0.4490e-06, -0.6507e-06, -0.8145e-06, -0.8801e-06,           &
    -0.8311e-06/
   data ((o1(ip,iw,1),iw=1,21),ip= 8, 8)/                                      &
     0.7110e-05,  0.1418e-04,  0.2825e-04,  0.5623e-04,  0.1116e-03,           &
     0.2206e-03,  0.4325e-03,  0.8373e-03,  0.1589e-02,  0.2923e-02,           &
     0.5141e-02,  0.8526e-02,  0.1324e-01,  0.1940e-01,  0.2741e-01,           &
     0.3807e-01,  0.5253e-01,  0.7227e-01,  0.9889e-01,  0.1338e+00,           &
     0.1777e+00/
   data ((o2(ip,iw,1),iw=1,21),ip= 8, 8)/                                      &
    -0.6983e-10, -0.2063e-09, -0.6785e-09, -0.2416e-08, -0.9000e-08,           &
    -0.3243e-07, -0.1100e-06, -0.3574e-06, -0.1104e-05, -0.3205e-05,           &
    -0.8516e-05, -0.2014e-04, -0.4151e-04, -0.7508e-04, -0.1234e-03,           &
    -0.1907e-03, -0.2811e-03, -0.3966e-03, -0.5355e-03, -0.6924e-03,           &
    -0.8613e-03/
   data ((o3(ip,iw,1),iw=1,21),ip= 8, 8)/                                      &
    -0.2353e-10, -0.4659e-10, -0.9153e-10, -0.1769e-09, -0.3313e-09,           &
    -0.6054e-09, -0.9899e-09, -0.1430e-08, -0.1375e-08,  0.3874e-09,           &
     0.5171e-08,  0.9807e-08, -0.7345e-09, -0.4604e-07, -0.1356e-06,           &
    -0.2731e-06, -0.4577e-06, -0.6632e-06, -0.8284e-06, -0.8894e-06,           &
    -0.8267e-06/
   data ((o1(ip,iw,1),iw=1,21),ip= 9, 9)/                                      &
     0.7105e-05,  0.1417e-04,  0.2823e-04,  0.5620e-04,  0.1116e-03,           &
     0.2205e-03,  0.4325e-03,  0.8376e-03,  0.1591e-02,  0.2929e-02,           &
     0.5162e-02,  0.8581e-02,  0.1336e-01,  0.1966e-01,  0.2790e-01,           &
     0.3894e-01,  0.5404e-01,  0.7484e-01,  0.1031e+00,  0.1405e+00,           &
     0.1880e+00/
   data ((o2(ip,iw,1),iw=1,21),ip= 9, 9)/                                      &
    -0.1481e-09, -0.3601e-09, -0.9762e-09, -0.2973e-08, -0.1014e-07,           &
    -0.3421e-07, -0.1121e-06, -0.3569e-06, -0.1092e-05, -0.3156e-05,           &
    -0.8375e-05, -0.1981e-04, -0.4090e-04, -0.7405e-04, -0.1218e-03,           &
    -0.1881e-03, -0.2770e-03, -0.3906e-03, -0.5269e-03, -0.6810e-03,           &
    -0.8471e-03/
   data ((o3(ip,iw,1),iw=1,21),ip= 9, 9)/                                      &
    -0.2304e-10, -0.4564e-10, -0.8969e-10, -0.1735e-09, -0.3224e-09,           &
    -0.5933e-09, -0.9756e-09, -0.1428e-08, -0.1446e-08,  0.1156e-09,           &
     0.4499e-08,  0.8469e-08, -0.2720e-08, -0.4904e-07, -0.1401e-06,           &
    -0.2801e-06, -0.4681e-06, -0.6761e-06, -0.8387e-06, -0.8879e-06,           &
    -0.8040e-06/
   data ((o1(ip,iw,1),iw=1,21),ip=10,10)/                                      &
     0.7098e-05,  0.1415e-04,  0.2821e-04,  0.5615e-04,  0.1115e-03,           &
     0.2204e-03,  0.4325e-03,  0.8382e-03,  0.1593e-02,  0.2940e-02,           &
     0.5194e-02,  0.8666e-02,  0.1356e-01,  0.2006e-01,  0.2865e-01,           &
     0.4026e-01,  0.5631e-01,  0.7863e-01,  0.1093e+00,  0.1500e+00,           &
     0.2017e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=10,10)/                                      &
    -0.2661e-09, -0.5923e-09, -0.1426e-08, -0.3816e-08, -0.1159e-07,           &
    -0.3654e-07, -0.1143e-06, -0.3559e-06, -0.1074e-05, -0.3083e-05,           &
    -0.8159e-05, -0.1932e-04, -0.3998e-04, -0.7253e-04, -0.1194e-03,           &
    -0.1845e-03, -0.2718e-03, -0.3833e-03, -0.5176e-03, -0.6701e-03,           &
    -0.8354e-03/
   data ((o3(ip,iw,1),iw=1,21),ip=10,10)/                                      &
    -0.2232e-10, -0.4421e-10, -0.8695e-10, -0.1684e-09, -0.3141e-09,           &
    -0.5765e-09, -0.9606e-09, -0.1434e-08, -0.1551e-08, -0.2663e-09,           &
     0.3515e-08,  0.6549e-08, -0.5479e-08, -0.5312e-07, -0.1460e-06,           &
    -0.2883e-06, -0.4787e-06, -0.6863e-06, -0.8399e-06, -0.8703e-06,           &
    -0.7602e-06/
   data ((o1(ip,iw,1),iw=1,21),ip=11,11)/                                      &
     0.7088e-05,  0.1413e-04,  0.2817e-04,  0.5608e-04,  0.1114e-03,           &
     0.2203e-03,  0.4325e-03,  0.8390e-03,  0.1598e-02,  0.2955e-02,           &
     0.5242e-02,  0.8796e-02,  0.1386e-01,  0.2067e-01,  0.2978e-01,           &
     0.4224e-01,  0.5964e-01,  0.8406e-01,  0.1178e+00,  0.1627e+00,           &
     0.2197e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=11,11)/                                      &
    -0.4394e-09, -0.9330e-09, -0.2086e-08, -0.5054e-08, -0.1373e-07,           &
    -0.3971e-07, -0.1178e-06, -0.3546e-06, -0.1049e-05, -0.2976e-05,           &
    -0.7847e-05, -0.1860e-04, -0.3864e-04, -0.7038e-04, -0.1162e-03,           &
    -0.1798e-03, -0.2654e-03, -0.3754e-03, -0.5091e-03, -0.6621e-03,           &
    -0.8286e-03/
   data ((o3(ip,iw,1),iw=1,21),ip=11,11)/                                      &
    -0.2127e-10, -0.4216e-10, -0.8300e-10, -0.1611e-09, -0.3019e-09,           &
    -0.5597e-09, -0.9431e-09, -0.1450e-08, -0.1694e-08, -0.7913e-09,           &
     0.2144e-08,  0.3990e-08, -0.9282e-08, -0.5810e-07, -0.1525e-06,           &
    -0.2965e-06, -0.4869e-06, -0.6894e-06, -0.8281e-06, -0.8350e-06,           &
    -0.6956e-06/
   data ((o1(ip,iw,1),iw=1,21),ip=12,12)/                                      &
     0.7073e-05,  0.1411e-04,  0.2811e-04,  0.5598e-04,  0.1112e-03,           &
     0.2201e-03,  0.4324e-03,  0.8401e-03,  0.1604e-02,  0.2978e-02,           &
     0.5313e-02,  0.8987e-02,  0.1431e-01,  0.2158e-01,  0.3145e-01,           &
     0.4512e-01,  0.6440e-01,  0.9161e-01,  0.1293e+00,  0.1793e+00,           &
     0.2426e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=12,12)/                                      &
    -0.6829e-09, -0.1412e-08, -0.3014e-08, -0.6799e-08, -0.1675e-07,           &
    -0.4450e-07, -0.1235e-06, -0.3538e-06, -0.1014e-05, -0.2827e-05,           &
    -0.7407e-05, -0.1759e-04, -0.3676e-04, -0.6744e-04, -0.1120e-03,           &
    -0.1742e-03, -0.2585e-03, -0.3683e-03, -0.5034e-03, -0.6594e-03,           &
    -0.8290e-03/
   data ((o3(ip,iw,1),iw=1,21),ip=12,12)/                                      &
    -0.1985e-10, -0.3937e-10, -0.7761e-10, -0.1511e-09, -0.2855e-09,           &
    -0.5313e-09, -0.9251e-09, -0.1470e-08, -0.1898e-08, -0.1519e-08,           &
     0.2914e-09,  0.5675e-09, -0.1405e-07, -0.6359e-07, -0.1584e-06,           &
    -0.3020e-06, -0.4893e-06, -0.6821e-06, -0.8021e-06, -0.7834e-06,           &
    -0.6105e-06/
   data ((o1(ip,iw,1),iw=1,21),ip=13,13)/                                      &
     0.7053e-05,  0.1407e-04,  0.2804e-04,  0.5584e-04,  0.1110e-03,           &
     0.2198e-03,  0.4324e-03,  0.8420e-03,  0.1613e-02,  0.3011e-02,           &
     0.5416e-02,  0.9263e-02,  0.1495e-01,  0.2289e-01,  0.3384e-01,           &
     0.4918e-01,  0.7096e-01,  0.1018e+00,  0.1442e+00,  0.2004e+00,           &
     0.2708e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=13,13)/                                      &
    -0.1004e-08, -0.2043e-08, -0.4239e-08, -0.9104e-08, -0.2075e-07,           &
    -0.5096e-07, -0.1307e-06, -0.3520e-06, -0.9671e-06, -0.2630e-05,           &
    -0.6825e-05, -0.1624e-04, -0.3429e-04, -0.6369e-04, -0.1069e-03,           &
    -0.1680e-03, -0.2520e-03, -0.3635e-03, -0.5029e-03, -0.6647e-03,           &
    -0.8390e-03/
   data ((o3(ip,iw,1),iw=1,21),ip=13,13)/                                      &
    -0.1807e-10, -0.3587e-10, -0.7085e-10, -0.1385e-09, -0.2648e-09,           &
    -0.4958e-09, -0.8900e-09, -0.1473e-08, -0.2112e-08, -0.2399e-08,           &
    -0.2002e-08, -0.3646e-08, -0.1931e-07, -0.6852e-07, -0.1618e-06,           &
    -0.3021e-06, -0.4828e-06, -0.6634e-06, -0.7643e-06, -0.7177e-06,           &
    -0.5054e-06/
   data ((o1(ip,iw,1),iw=1,21),ip=14,14)/                                      &
     0.7029e-05,  0.1402e-04,  0.2795e-04,  0.5567e-04,  0.1107e-03,           &
     0.2195e-03,  0.4326e-03,  0.8447e-03,  0.1625e-02,  0.3056e-02,           &
     0.5554e-02,  0.9638e-02,  0.1584e-01,  0.2470e-01,  0.3713e-01,           &
     0.5470e-01,  0.7969e-01,  0.1149e+00,  0.1631e+00,  0.2263e+00,           &
     0.3045e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=14,14)/                                      &
    -0.1387e-08, -0.2798e-08, -0.5706e-08, -0.1187e-07, -0.2564e-07,           &
    -0.5866e-07, -0.1398e-06, -0.3516e-06, -0.9148e-06, -0.2398e-05,           &
    -0.6122e-05, -0.1459e-04, -0.3125e-04, -0.5923e-04, -0.1013e-03,           &
    -0.1620e-03, -0.2473e-03, -0.3631e-03, -0.5098e-03, -0.6800e-03,           &
    -0.8603e-03/
   data ((o3(ip,iw,1),iw=1,21),ip=14,14)/                                      &
    -0.1610e-10, -0.3200e-10, -0.6337e-10, -0.1245e-09, -0.2408e-09,           &
    -0.4533e-09, -0.8405e-09, -0.1464e-08, -0.2337e-08, -0.3341e-08,           &
    -0.4467e-08, -0.8154e-08, -0.2436e-07, -0.7128e-07, -0.1604e-06,           &
    -0.2945e-06, -0.4666e-06, -0.6357e-06, -0.7187e-06, -0.6419e-06,           &
    -0.3795e-06/
   data ((o1(ip,iw,1),iw=1,21),ip=15,15)/                                      &
     0.7003e-05,  0.1397e-04,  0.2785e-04,  0.5550e-04,  0.1104e-03,           &
     0.2192e-03,  0.4329e-03,  0.8481e-03,  0.1641e-02,  0.3112e-02,           &
     0.5729e-02,  0.1012e-01,  0.1698e-01,  0.2708e-01,  0.4146e-01,           &
     0.6189e-01,  0.9085e-01,  0.1313e+00,  0.1862e+00,  0.2571e+00,           &
     0.3433e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=15,15)/                                      &
    -0.1788e-08, -0.3588e-08, -0.7244e-08, -0.1479e-07, -0.3083e-07,           &
    -0.6671e-07, -0.1497e-06, -0.3519e-06, -0.8607e-06, -0.2154e-05,           &
    -0.5364e-05, -0.1276e-04, -0.2785e-04, -0.5435e-04, -0.9573e-04,           &
    -0.1570e-03, -0.2455e-03, -0.3682e-03, -0.5253e-03, -0.7065e-03,           &
    -0.8938e-03/
   data ((o3(ip,iw,1),iw=1,21),ip=15,15)/                                      &
    -0.1429e-10, -0.2843e-10, -0.5645e-10, -0.1115e-09, -0.2181e-09,           &
    -0.4200e-09, -0.7916e-09, -0.1460e-08, -0.2542e-08, -0.4168e-08,           &
    -0.6703e-08, -0.1215e-07, -0.2821e-07, -0.7073e-07, -0.1530e-06,           &
    -0.2791e-06, -0.4426e-06, -0.6027e-06, -0.6707e-06, -0.5591e-06,           &
    -0.2328e-06/
   data ((o1(ip,iw,1),iw=1,21),ip=16,16)/                                      &
     0.6981e-05,  0.1393e-04,  0.2777e-04,  0.5535e-04,  0.1102e-03,           &
     0.2190e-03,  0.4336e-03,  0.8527e-03,  0.1660e-02,  0.3177e-02,           &
     0.5930e-02,  0.1068e-01,  0.1836e-01,  0.3000e-01,  0.4686e-01,           &
     0.7083e-01,  0.1046e+00,  0.1511e+00,  0.2134e+00,  0.2924e+00,           &
     0.3861e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=16,16)/                                      &
    -0.2141e-08, -0.4286e-08, -0.8603e-08, -0.1737e-07, -0.3548e-07,           &
    -0.7410e-07, -0.1590e-06, -0.3537e-06, -0.8142e-06, -0.1935e-05,           &
    -0.4658e-05, -0.1099e-04, -0.2444e-04, -0.4948e-04, -0.9067e-04,           &
    -0.1538e-03, -0.2474e-03, -0.3793e-03, -0.5495e-03, -0.7439e-03,           &
    -0.9383e-03/
   data ((o3(ip,iw,1),iw=1,21),ip=16,16)/                                      &
    -0.1295e-10, -0.2581e-10, -0.5136e-10, -0.1019e-09, -0.2011e-09,           &
    -0.3916e-09, -0.7585e-09, -0.1439e-08, -0.2648e-08, -0.4747e-08,           &
    -0.8301e-08, -0.1499e-07, -0.3024e-07, -0.6702e-07, -0.1399e-06,           &
    -0.2564e-06, -0.4117e-06, -0.5669e-06, -0.6239e-06, -0.4748e-06,           &
    -0.7013e-07/
   data ((o1(ip,iw,1),iw=1,21),ip=17,17)/                                      &
     0.6967e-05,  0.1390e-04,  0.2772e-04,  0.5527e-04,  0.1101e-03,           &
     0.2191e-03,  0.4346e-03,  0.8579e-03,  0.1679e-02,  0.3244e-02,           &
     0.6139e-02,  0.1127e-01,  0.1986e-01,  0.3330e-01,  0.5315e-01,           &
     0.8139e-01,  0.1207e+00,  0.1741e+00,  0.2442e+00,  0.3311e+00,           &
     0.4312e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=17,17)/                                      &
    -0.2400e-08, -0.4796e-08, -0.9599e-08, -0.1927e-07, -0.3892e-07,           &
    -0.7954e-07, -0.1661e-06, -0.3540e-06, -0.7780e-06, -0.1763e-05,           &
    -0.4092e-05, -0.9512e-05, -0.2142e-04, -0.4502e-04, -0.8640e-04,           &
    -0.1525e-03, -0.2526e-03, -0.3955e-03, -0.5805e-03, -0.7897e-03,           &
    -0.9899e-03/
   data ((o3(ip,iw,1),iw=1,21),ip=17,17)/                                      &
    -0.1220e-10, -0.2432e-10, -0.4845e-10, -0.9640e-10, -0.1912e-09,           &
    -0.3771e-09, -0.7392e-09, -0.1420e-08, -0.2702e-08, -0.5049e-08,           &
    -0.9214e-08, -0.1659e-07, -0.3101e-07, -0.6162e-07, -0.1235e-06,           &
    -0.2287e-06, -0.3755e-06, -0.5274e-06, -0.5790e-06, -0.3947e-06,           &
     0.1003e-06/
   data ((o1(ip,iw,1),iw=1,21),ip=18,18)/                                      &
     0.6963e-05,  0.1389e-04,  0.2771e-04,  0.5526e-04,  0.1101e-03,           &
     0.2193e-03,  0.4359e-03,  0.8629e-03,  0.1698e-02,  0.3305e-02,           &
     0.6331e-02,  0.1183e-01,  0.2133e-01,  0.3671e-01,  0.5996e-01,           &
     0.9320e-01,  0.1389e+00,  0.1999e+00,  0.2778e+00,  0.3717e+00,           &
     0.4761e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=18,18)/                                      &
    -0.2557e-08, -0.5106e-08, -0.1020e-07, -0.2043e-07, -0.4103e-07,           &
    -0.8293e-07, -0.1697e-06, -0.3531e-06, -0.7531e-06, -0.1645e-05,           &
    -0.3690e-05, -0.8411e-05, -0.1902e-04, -0.4118e-04, -0.8276e-04,           &
    -0.1525e-03, -0.2601e-03, -0.4147e-03, -0.6149e-03, -0.8384e-03,           &
    -0.1042e-02/
   data ((o3(ip,iw,1),iw=1,21),ip=18,18)/                                      &
    -0.1189e-10, -0.2372e-10, -0.4729e-10, -0.9421e-10, -0.1873e-09,           &
    -0.3713e-09, -0.7317e-09, -0.1437e-08, -0.2764e-08, -0.5243e-08,           &
    -0.9691e-08, -0.1751e-07, -0.3122e-07, -0.5693e-07, -0.1076e-06,           &
    -0.1981e-06, -0.3324e-06, -0.4785e-06, -0.5280e-06, -0.3174e-06,           &
     0.2672e-06/
   data ((o1(ip,iw,1),iw=1,21),ip=19,19)/                                      &
     0.6963e-05,  0.1389e-04,  0.2771e-04,  0.5528e-04,  0.1102e-03,           &
     0.2196e-03,  0.4370e-03,  0.8672e-03,  0.1712e-02,  0.3355e-02,           &
     0.6488e-02,  0.1230e-01,  0.2262e-01,  0.3989e-01,  0.6677e-01,           &
     0.1056e+00,  0.1586e+00,  0.2276e+00,  0.3131e+00,  0.4124e+00,           &
     0.5188e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=19,19)/                                      &
    -0.2630e-08, -0.5249e-08, -0.1048e-07, -0.2096e-07, -0.4198e-07,           &
    -0.8440e-07, -0.1710e-06, -0.3513e-06, -0.7326e-06, -0.1562e-05,           &
    -0.3416e-05, -0.7637e-05, -0.1719e-04, -0.3795e-04, -0.7926e-04,           &
    -0.1524e-03, -0.2680e-03, -0.4344e-03, -0.6486e-03, -0.8838e-03,           &
    -0.1089e-02/
   data ((o3(ip,iw,1),iw=1,21),ip=19,19)/                                      &
    -0.1188e-10, -0.2369e-10, -0.4725e-10, -0.9417e-10, -0.1875e-09,           &
    -0.3725e-09, -0.7365e-09, -0.1445e-08, -0.2814e-08, -0.5384e-08,           &
    -0.1008e-07, -0.1816e-07, -0.3179e-07, -0.5453e-07, -0.9500e-07,           &
    -0.1679e-06, -0.2819e-06, -0.4109e-06, -0.4555e-06, -0.2283e-06,           &
     0.4283e-06/
   data ((o1(ip,iw,1),iw=1,21),ip=20,20)/                                      &
     0.6963e-05,  0.1389e-04,  0.2772e-04,  0.5529e-04,  0.1103e-03,           &
     0.2198e-03,  0.4377e-03,  0.8701e-03,  0.1723e-02,  0.3391e-02,           &
     0.6606e-02,  0.1266e-01,  0.2366e-01,  0.4262e-01,  0.7304e-01,           &
     0.1179e+00,  0.1789e+00,  0.2563e+00,  0.3487e+00,  0.4516e+00,           &
     0.5572e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=20,20)/                                      &
    -0.2651e-08, -0.5291e-08, -0.1056e-07, -0.2110e-07, -0.4221e-07,           &
    -0.8462e-07, -0.1705e-06, -0.3466e-06, -0.7155e-06, -0.1501e-05,           &
    -0.3223e-05, -0.7079e-05, -0.1581e-04, -0.3517e-04, -0.7553e-04,           &
    -0.1510e-03, -0.2746e-03, -0.4528e-03, -0.6789e-03, -0.9214e-03,           &
    -0.1124e-02/
   data ((o3(ip,iw,1),iw=1,21),ip=20,20)/                                      &
    -0.1193e-10, -0.2380e-10, -0.4748e-10, -0.9465e-10, -0.1886e-09,           &
    -0.3751e-09, -0.7436e-09, -0.1466e-08, -0.2872e-08, -0.5508e-08,           &
    -0.1038e-07, -0.1891e-07, -0.3279e-07, -0.5420e-07, -0.8711e-07,           &
    -0.1403e-06, -0.2248e-06, -0.3221e-06, -0.3459e-06, -0.1066e-06,           &
     0.5938e-06/
   data ((o1(ip,iw,1),iw=1,21),ip=21,21)/                                      &
     0.6964e-05,  0.1389e-04,  0.2772e-04,  0.5530e-04,  0.1103e-03,           &
     0.2199e-03,  0.4382e-03,  0.8719e-03,  0.1730e-02,  0.3416e-02,           &
     0.6690e-02,  0.1293e-01,  0.2445e-01,  0.4479e-01,  0.7837e-01,           &
     0.1291e+00,  0.1985e+00,  0.2846e+00,  0.3831e+00,  0.4875e+00,           &
     0.5902e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=21,21)/                                      &
    -0.2654e-08, -0.5296e-08, -0.1057e-07, -0.2111e-07, -0.4219e-07,           &
    -0.8445e-07, -0.1696e-06, -0.3428e-06, -0.7013e-06, -0.1458e-05,           &
    -0.3084e-05, -0.6678e-05, -0.1476e-04, -0.3284e-04, -0.7173e-04,           &
    -0.1481e-03, -0.2786e-03, -0.4688e-03, -0.7052e-03, -0.9506e-03,           &
    -0.1148e-02/
   data ((o3(ip,iw,1),iw=1,21),ip=21,21)/                                      &
    -0.1195e-10, -0.2384e-10, -0.4755e-10, -0.9482e-10, -0.1890e-09,           &
    -0.3761e-09, -0.7469e-09, -0.1476e-08, -0.2892e-08, -0.5603e-08,           &
    -0.1060e-07, -0.1942e-07, -0.3393e-07, -0.5508e-07, -0.8290e-07,           &
    -0.1182e-06, -0.1657e-06, -0.2170e-06, -0.1997e-06,  0.6227e-07,           &
     0.7847e-06/
   data ((o1(ip,iw,1),iw=1,21),ip=22,22)/                                      &
     0.6964e-05,  0.1389e-04,  0.2772e-04,  0.5530e-04,  0.1103e-03,           &
     0.2200e-03,  0.4385e-03,  0.8731e-03,  0.1735e-02,  0.3433e-02,           &
     0.6748e-02,  0.1311e-01,  0.2502e-01,  0.4642e-01,  0.8258e-01,           &
     0.1385e+00,  0.2160e+00,  0.3107e+00,  0.4146e+00,  0.5188e+00,           &
     0.6171e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=22,22)/                                      &
    -0.2653e-08, -0.5295e-08, -0.1057e-07, -0.2110e-07, -0.4215e-07,           &
    -0.8430e-07, -0.1690e-06, -0.3403e-06, -0.6919e-06, -0.1427e-05,           &
    -0.2991e-05, -0.6399e-05, -0.1398e-04, -0.3099e-04, -0.6824e-04,           &
    -0.1441e-03, -0.2795e-03, -0.4814e-03, -0.7282e-03, -0.9739e-03,           &
    -0.1163e-02/
   data ((o3(ip,iw,1),iw=1,21),ip=22,22)/                                      &
    -0.1195e-10, -0.2384e-10, -0.4756e-10, -0.9485e-10, -0.1891e-09,           &
    -0.3765e-09, -0.7483e-09, -0.1481e-08, -0.2908e-08, -0.5660e-08,           &
    -0.1075e-07, -0.1980e-07, -0.3472e-07, -0.5626e-07, -0.8149e-07,           &
    -0.1027e-06, -0.1136e-06, -0.1071e-06, -0.2991e-07,  0.2743e-06,           &
     0.1017e-05/
   data ((o1(ip,iw,1),iw=1,21),ip=23,23)/                                      &
     0.6964e-05,  0.1389e-04,  0.2772e-04,  0.5530e-04,  0.1103e-03,           &
     0.2200e-03,  0.4387e-03,  0.8738e-03,  0.1738e-02,  0.3445e-02,           &
     0.6786e-02,  0.1324e-01,  0.2541e-01,  0.4757e-01,  0.8567e-01,           &
     0.1459e+00,  0.2303e+00,  0.3331e+00,  0.4415e+00,  0.5444e+00,           &
     0.6377e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=23,23)/                                      &
    -0.2653e-08, -0.5294e-08, -0.1057e-07, -0.2109e-07, -0.4212e-07,           &
    -0.8420e-07, -0.1686e-06, -0.3388e-06, -0.6858e-06, -0.1406e-05,           &
    -0.2928e-05, -0.6206e-05, -0.1344e-04, -0.2961e-04, -0.6533e-04,           &
    -0.1399e-03, -0.2780e-03, -0.4904e-03, -0.7488e-03, -0.9953e-03,           &
    -0.1175e-02/
   data ((o3(ip,iw,1),iw=1,21),ip=23,23)/                                      &
    -0.1195e-10, -0.2384e-10, -0.4756e-10, -0.9485e-10, -0.1891e-09,           &
    -0.3767e-09, -0.7492e-09, -0.1485e-08, -0.2924e-08, -0.5671e-08,           &
    -0.1084e-07, -0.2009e-07, -0.3549e-07, -0.5773e-07, -0.8208e-07,           &
    -0.9394e-07, -0.7270e-07, -0.3947e-08,  0.1456e-06,  0.5083e-06,           &
     0.1270e-05/
   data ((o1(ip,iw,1),iw=1,21),ip=24,24)/                                      &
     0.6964e-05,  0.1389e-04,  0.2772e-04,  0.5531e-04,  0.1103e-03,           &
     0.2201e-03,  0.4388e-03,  0.8743e-03,  0.1740e-02,  0.3452e-02,           &
     0.6811e-02,  0.1332e-01,  0.2567e-01,  0.4835e-01,  0.8782e-01,           &
     0.1511e+00,  0.2412e+00,  0.3507e+00,  0.4627e+00,  0.5638e+00,           &
     0.6526e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=24,24)/                                      &
    -0.2653e-08, -0.5294e-08, -0.1056e-07, -0.2109e-07, -0.4210e-07,           &
    -0.8413e-07, -0.1684e-06, -0.3379e-06, -0.6820e-06, -0.1393e-05,           &
    -0.2889e-05, -0.6080e-05, -0.1307e-04, -0.2861e-04, -0.6310e-04,           &
    -0.1363e-03, -0.2758e-03, -0.4969e-03, -0.7681e-03, -0.1017e-02,           &
    -0.1186e-02/
   data ((o3(ip,iw,1),iw=1,21),ip=24,24)/                                      &
    -0.1195e-10, -0.2384e-10, -0.4756e-10, -0.9485e-10, -0.1891e-09,           &
    -0.3768e-09, -0.7497e-09, -0.1487e-08, -0.2933e-08, -0.5710e-08,           &
    -0.1089e-07, -0.2037e-07, -0.3616e-07, -0.5907e-07, -0.8351e-07,           &
    -0.8925e-07, -0.4122e-07,  0.8779e-07,  0.3143e-06,  0.7281e-06,           &
     0.1500e-05/
   data ((o1(ip,iw,1),iw=1,21),ip=25,25)/                                      &
     0.6964e-05,  0.1389e-04,  0.2772e-04,  0.5531e-04,  0.1103e-03,           &
     0.2201e-03,  0.4388e-03,  0.8745e-03,  0.1741e-02,  0.3456e-02,           &
     0.6827e-02,  0.1337e-01,  0.2584e-01,  0.4885e-01,  0.8924e-01,           &
     0.1547e+00,  0.2488e+00,  0.3632e+00,  0.4779e+00,  0.5771e+00,           &
     0.6627e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=25,25)/                                      &
    -0.2653e-08, -0.5293e-08, -0.1056e-07, -0.2108e-07, -0.4209e-07,           &
    -0.8409e-07, -0.1682e-06, -0.3373e-06, -0.6797e-06, -0.1383e-05,           &
    -0.2862e-05, -0.5993e-05, -0.1283e-04, -0.2795e-04, -0.6158e-04,           &
    -0.1338e-03, -0.2743e-03, -0.5030e-03, -0.7863e-03, -0.1038e-02,           &
    -0.1196e-02/
   data ((o3(ip,iw,1),iw=1,21),ip=25,25)/                                      &
    -0.1195e-10, -0.2383e-10, -0.4755e-10, -0.9484e-10, -0.1891e-09,           &
    -0.3768e-09, -0.7499e-09, -0.1489e-08, -0.2939e-08, -0.5741e-08,           &
    -0.1100e-07, -0.2066e-07, -0.3660e-07, -0.6002e-07, -0.8431e-07,           &
    -0.8556e-07, -0.1674e-07,  0.1638e-06,  0.4525e-06,  0.8949e-06,           &
     0.1669e-05/
   data ((o1(ip,iw,1),iw=1,21),ip=26,26)/                                      &
     0.6964e-05,  0.1389e-04,  0.2772e-04,  0.5531e-04,  0.1103e-03,           &
     0.2201e-03,  0.4389e-03,  0.8747e-03,  0.1741e-02,  0.3458e-02,           &
     0.6836e-02,  0.1340e-01,  0.2594e-01,  0.4916e-01,  0.9011e-01,           &
     0.1569e+00,  0.2536e+00,  0.3714e+00,  0.4877e+00,  0.5856e+00,           &
     0.6695e+00/
   data ((o2(ip,iw,1),iw=1,21),ip=26,26)/                                      &
    -0.2652e-08, -0.5292e-08, -0.1056e-07, -0.2108e-07, -0.4208e-07,           &
    -0.8406e-07, -0.1681e-06, -0.3369e-06, -0.6784e-06, -0.1378e-05,           &
    -0.2843e-05, -0.5944e-05, -0.1269e-04, -0.2759e-04, -0.6078e-04,           &
    -0.1326e-03, -0.2742e-03, -0.5088e-03, -0.8013e-03, -0.1054e-02,           &
    -0.1202e-02/
   data ((o3(ip,iw,1),iw=1,21),ip=26,26)/                                      &
    -0.1194e-10, -0.2383e-10, -0.4754e-10, -0.9482e-10, -0.1891e-09,           &
    -0.3768e-09, -0.7499e-09, -0.1489e-08, -0.2941e-08, -0.5752e-08,           &
    -0.1104e-07, -0.2069e-07, -0.3661e-07, -0.6012e-07, -0.8399e-07,           &
    -0.8183e-07,  0.1930e-08,  0.2167e-06,  0.5434e-06,  0.9990e-06,           &
     0.1787e-05/
   save c1,c2,c3,o1,o2,o3
   save h11,h12,h13,h21,h22,h23,h81,h82,h83
!
!-----tables co2 and h2o are only used with 'high' option
!
   if (high) then
!
     do iw = 1,nh
       do ip = 1,nx
         hh11(ip,iw,1)=1.0-h11(ip,iw,1)
         hh21(ip,iw,1)=1.0-h21(ip,iw,1)
         hh81(ip,iw,1)=1.0-h81(ip,iw,1)
       enddo
     enddo
!
     do iw = 1,nc
       do ip = 1,nx
         cc1(ip,iw,1)=1.0-c1(ip,iw,1)
       enddo
     enddo
!
   endif
!
!-----always use table look-up for ozone transmittance
!
   do iw = 1,no
     do ip = 1,nx
       oo1(ip,iw,1)=1.0-o1(ip,iw,1)
     enddo
   enddo
!
!-----compute layer pressure (pa) and layer temperature minus 250k (dt)
!
   do k = 1,np
     do i = 1,ipts   
       pa(i,k)=0.5*(prse(i,k)+prse(i,k+1))
       dt(i,k)=ta(i,k)-250.0
     enddo
   enddo
!
!-----compute layer absorber amount
!
!     dh2o : water vapor amount (g/cm**2)
!     dcont: scaled water vapor amount for continuum absorption (g/cm**2)
!     dco2 : co2 amount (cm-atm)stp
!     do3  : o3 amount (cm-atm)stp
!     dn2o : n2o amount (cm-atm)stp
!     dch4 : ch4 amount (cm-atm)stp
!     df11 : cfc11 amount (cm-atm)stp
!     df12 : cfc12 amount (cm-atm)stp
!     df22 : cfc22 amount (cm-atm)stp
!     the factor 1.02 is equal to 1000/980
!     factors 789 and 476 are for unit conversion
!       789=22.4*100*100/9.8/28.966 convert from l/l to (cm-atm)stp
!       476=100*100/9.8/(48/22.4)   convert from g/g to (cm-atm)stp
!     the factor 0.001618 is equal to 1.02/(.622*1013.25)
!     the factor 6.081 is equal to 1800/296
!
   do k = 1,np
     do i = 1,ipts   
       dp           = prse(i,k+1)-prse(i,k)
       dh2o (i,k) = 1.02*wa(i,k)*dp+1.e-10
       do3  (i,k) = 476.0*oa(i,k)*dp+1.e-10
       dco2 (i,k) = 789.*co2*dp+1.e-10
       dch4 (i,k) = 789.*ch4*dp+1.e-10
       dn2o (i,k) = 789.*n2o*dp+1.e-10
       df11 (i,k) = 789.*cfc11*dp+1.e-10
       df12 (i,k) = 789.*cfc12*dp+1.e-10
       df22 (i,k) = 789.*cfc22*dp+1.e-10
!
!-----compute scaled water vapor amount for h2o continuum absorption
!     following eq. (43).
!
       xx=pa(i,k)*0.001618*wa(i,k)*wa(i,k)*dp
       dcont(i,k) = xx*exp(1800./ta(i,k)-6.081)+1.e-10
!
!-----compute effective cloud-free fraction, clr, for each layer
!     the cloud diffuse transmittance is approximated by using 
!     a diffusivity factor of 1.66, cf. eq. (7)
!
       clr(i,k)=1.0-(fcld(i,k)*(1.0-exp(-1.66*taucl(i,k))))
     enddo
   enddo
!
!-----compute column-integrated h2o amoumt, h2o-weighted pressure
!     and temperature.  it follows eqs. (37) and (38).
!
   if (high) then
     call rad_lw_nasa_column(ipts,m,np,pa,dt,dh2o,sh2o,swpre,swtem)
   endif
!
!-----the surface (with an index np+1) is treated as a layer filled with
!     black clouds.
!
   do i = 1,ipts   
     clr(i,0)    = 1.0
     clr(i,np+1) = 0.0
     st4(i)      = 0.0
   enddo
!
!-----initialize fluxes
!
   do k = 1,np+1
     do i = 1,ipts   
       flx(i,k)  = 0.0
       flc(i,k)  = 0.0
       dfdts(i,k)= 0.0
       flxu(i,k) = 0.0
       flxd(i,k) = 0.0
       flxuc(i,k) = 0.0
       flxdc(i,k) = 0.0
       rflx(i,k) = 0.0
       rflc(i,k) = 0.0
     enddo
   enddo
!
!-----integration over spectral bands
!
   loop_1000: do ib = 1,10
!
!-----if h2otbl, compute h2o (line) transmittance using table look-up.
!     if conbnd, compute h2o (continuum) transmittance in bands 3-6.
!     if co2bnd, compute co2 transmittance in band 3.
!     if oznbnd, compute  o3 transmittance in band 5.
!     if n2obnd, compute n2o transmittance in band 3, 6 and 7.
!     if ch4bnd, compute ch4 transmittance in band 6 and 7.
!     if combnd, compute co2-minor transmittance in band 4 and 5.
!     if f11bnd, compute cfc11 transmittance in band 4 and 5.
!     if f12bnd, compute cfc12 transmittance in band 4 and 6.
!     if f22bnd, compute cfc22 transmittance in band 4 and 6.
!     if b10bnd, compute flux reduction due to n2o in band 10 only (540-620).
!
     h2otbl=high.and.(ib.eq.1.or.ib.eq.2.or.ib.eq.8)
     conbnd=ib.ge.3.and.ib.le.6
     co2bnd=ib.eq.3
     oznbnd=ib.eq.5
     n2obnd=ib.eq.6.or.ib.eq.7
     ch4bnd=ib.eq.6.or.ib.eq.7
     combnd=ib.eq.4.or.ib.eq.5
     f11bnd=ib.eq.4.or.ib.eq.5
     f12bnd=ib.eq.4.or.ib.eq.6
     f22bnd=ib.eq.4.or.ib.eq.6
     b10bnd=ib.eq.10
!
!-----blayer is the spectrally integrated planck flux of the mean layer
!     temperature derived from eq. (22)
!     the fitting for the planck flux is valid in the range 160-345 k.
!
     do k=1,np
       do i=1,ipts   
         blayer(i,k)=ta(i,k)*(ta(i,k)*(ta(i,k)                                 &
                       *(ta(i,k)*cb(5,ib)+cb(4,ib))+cb(3,ib))                  &
                       +cb(2,ib))+cb(1,ib)
       enddo
     enddo
!
!-----the earth's surface, with an index "np+1", is treated as a layer
!
     do i = 1,ipts   
       blayer(i,0)   = 0.0
       blayer(i,np+1)= ts(i)*(ts(i)*(ts(i)                                     &
                       *(ts(i)*cb(5,ib)+cb(4,ib))+cb(3,ib))                    &
                       +cb(2,ib))+cb(1,ib)
!
!-----dbs is the derivative of the surface planck flux with respect to
!     surface temperature (eq. 59).
!
       dbs(i)=ts(i)*(ts(i)*(ts(i)                                              &
              *4.*cb(5,ib)+3.*cb(4,ib))+2.*cb(3,ib))+cb(2,ib)
     enddo
!
!-----compute column-integrated absorber amoumt, absorber-weighted
!     pressure and temperature for co2 (band 3) and o3 (band 5).
!     it follows eqs. (37) and (38).
!
!-----this is in the band loop to save storage
!
     if (high .and. co2bnd) then
       call rad_lw_nasa_column(ipts,m,np,pa,dt,dco2,sco3,scopre,scotem)
     endif
!
     if (oznbnd) then
       call rad_lw_nasa_column(ipts,m,np,pa,dt,do3,sco3,scopre,scotem)
     endif
!
!-----compute the exponential terms (eq. 32) at each layer due to
!     water vapor line absorption when k-distribution is used
!
     if (.not.h2otbl .and. .not.b10bnd) then
       call rad_lw_nasa_h2oexps(ipts,ib,m,np,dh2o,pa,                          &
                                  dt,xkw,aw,bw,pm,mw,h2oexp)

     endif
!
!-----compute the exponential terms (eq. 46) at each layer due to
!     water vapor continuum absorption
!
     if (conbnd) then
       call rad_lw_nasa_conexps(ipts,ib,m,np,dcont,xke,ne,conexp)
     endif
!
!-----compute the exponential terms (eq. 32) at each layer due to
!     co2 absorption
!
     if (.not.high .and. co2bnd) then
       call rad_lw_nasa_co2exps(ipts,m,np,dco2,pa,dt,co2exp)
     endif
!
!***** for trace gases ****!
!
     if (trace) then
!
!-----compute the exponential terms at each layer due to n2o absorption
!
       if (n2obnd) then
         call rad_lw_nasa_n2oexps(ipts,ib,m,np,dn2o,pa,dt,n2oexp)
       endif
!
!-----compute the exponential terms at each layer due to ch4 absorption
!
       if (ch4bnd) then
         call rad_lw_nasa_ch4exps(ipts,ib,m,np,dch4,pa,dt,ch4exp)
       endif
!
!-----compute the exponential terms at each layer due to co2 minor absorption
!
       if (combnd) then
         call rad_lw_nasa_comexps(ipts,ib,m,np,dco2,dt,comexp)
       endif
!
!-----compute the exponential terms at each layer due to cfc11 absorption
!
       if (f11bnd) then
         a1  = 1.26610e-3
         b1  = 3.55940e-6
         fk1 = 1.57481e+1
         a2  = 8.19370e-4
         b2  = 4.67810e-6
         fk2 = 8.42338e+0
         call rad_lw_nasa_cfcexps(ipts,ib,m,np,a1,b1,fk1,                      &
                                     a2,b2,fk2,df11,dt,f11exp)

       endif
!
!-----compute the exponential terms at each layer due to cfc12 absorption
!
       if (f12bnd) then
         a1  = 8.77370e-4
         b1  =-5.88440e-6
         fk1 = 1.31266e+1
         a2  = 8.62000e-4
         b2  =-4.22500e-6
         fk2 = 3.07189e+1
         call rad_lw_nasa_cfcexps(ipts,ib,m,np,a1,b1,fk1,                      &
                                     a2,b2,fk2,df12,dt,f12exp)
       endif
!
!-----compute the exponential terms at each layer due to cfc22 absorption
!
       if (f22bnd) then
         a1  = 9.65130e-4
         b1  = 1.31280e-5
         fk1 = 5.13385e+0
         a2  =-3.00010e-5 
         b2  = 5.25010e-7
         fk2 = 2.72167e+1
         call rad_lw_nasa_cfcexps(ipts,ib,m,np,a1,b1,fk1,                      &
                                     a2,b2,fk2,df22,dt,f22exp)
       endif
!
!-----compute the exponential terms at each layer in band 10 due to
!     h2o line and continuum, co2, and n2o absorption
!
       if (b10bnd) then
         call rad_lw_nasa_b10exps(ipts,m,np,dh2o,dcont,dco2,dn2o,pa,dt         &
                                    ,h2oexp,conexp,co2exp,n2oexp)
       endif
!
     endif
!
!-----compute transmittances for regions between levels k1 and k2
!     and update the fluxes at the two levels.
!
     loop_2000: do k1 = 1,np
!
!-----initialize fclr, th2o, tcon, tco2, and tranal
!     fclr is the clear fraction of the line-of-sight
!     clrlw, clrmd, and clrhi are the clear-sky fractions of the 
!     low, middle, and high troposphere, respectively
!     tranal is the aerosol transmission function
!
       do i = 1,ipts   
         fclr(i)  = 1.0
         clrlw(i) = 1.0
         clrmd(i) = 1.0
         clrhi(i) = 1.0
         tranal(i)= 1.0
       enddo
!
!-----for h2o line transmission
!
       if (.not. h2otbl) then
         do ik = 1,6
           do i = 1,ipts   
             th2o(i,ik)=1.0
           enddo
         enddo
       endif
!
!-----for h2o continuum transmission
!
       if (conbnd) then
         do iq = 1,3
           do i = 1,ipts   
             tcon(i,iq)=1.0
           enddo
         enddo
       endif
!
!-----for co2 transmission using k-distribution method.
!     band 3 is divided into 3 sub-bands, but sub-bands 3a and 3c
!     are combined in computing the co2 transmittance.
!
       if (.not.high .and. co2bnd) then
         do isb = 1,2
           do ik = 1,6
             do i = 1,ipts   
               tco2(i,ik,isb)=1.0
             enddo
           enddo
         enddo
       endif
!
!***** for trace gases ****!
!
       if (trace) then
!
!-----for n2o transmission using k-distribution method.
!
         if (n2obnd) then
           do ik = 1,4
             do i = 1,ipts   
               tn2o(i,ik)=1.0
             enddo
           enddo
         endif
!
!-----for ch4 transmission using k-distribution method.
!
         if (ch4bnd) then
           do ik = 1,4
             do i = 1,ipts   
               tch4(i,ik)=1.0
             enddo
           enddo
         endif
!
!-----for co2-minor transmission using k-distribution method.
!
         if (combnd) then
           do ik = 1,2
             do i = 1,ipts   
               tcom(i,ik)=1.0
             enddo
           enddo
         endif
!
!-----for cfc-11 transmission using k-distribution method.
!
         if (f11bnd) then
           do i = 1,ipts   
             tf11(i)=1.0
           enddo
         endif
!
!-----for cfc-12 transmission using k-distribution method.
!
         if (f12bnd) then
           do i = 1,ipts   
             tf12(i)=1.0
           enddo
         endif
!
!-----for cfc-22 transmission when using k-distribution method.
!
         if (f22bnd) then
           do i = 1,ipts   
             tf22(i)=1.0
           enddo
         endif
!
!-----for the transmission in band 10 using k-distribution method.
!
         if (b10bnd) then
           do ik = 1,6
             do i = 1,ipts   
               th2o(i,ik)=1.0
               tco2(i,ik,1)=1.0
             enddo
           enddo
           do iq = 1,3
             do i = 1,ipts   
               tcon(i,iq)=1.0
             enddo
           enddo
           do ik = 1,4
             do i = 1,ipts   
               tn2o(i,ik)=1.0
             enddo
           enddo
         endif
!
       endif
!
!-----loop over the bottom level of the region (k2)
!
       loop_3000: do k2 = k1+1,np+1
!
!-----initialize total transmission function (trant)
!
         do i = 1,ipts   
           trant(i)=1.0
         enddo
!
         if (h2otbl) then
           w1=-8.0
           p1=-2.0
           dwe=0.3
           dpe=0.2
!
!-----compute water vapor transmittance using table look-up
!
           if (ib.eq.1) then
             call rad_lw_nasa_tablup(ipts,k1,k2,m,np,nx,nh,nt,                 &
                                          sh2o,swpre,swtem,                    &
                                          w1,p1,dwe,dpe,hh11,h12,h13,trant)

           endif
           if (ib.eq.2) then
             call rad_lw_nasa_tablup(ipts,k1,k2,m,np,nx,nh,nt,                 &
                                          sh2o,swpre,swtem,                    &
                                          w1,p1,dwe,dpe,hh21,h22,h23,trant)
           endif
           if (ib.eq.8) then
             call rad_lw_nasa_tablup(ipts,k1,k2,m,np,nx,nh,nt,                 &
                                          sh2o,swpre,swtem,                    &
                                          w1,p1,dwe,dpe,hh81,h82,h83,trant)
           endif
         else
!
!-----compute water vapor transmittance using k-distribution.
!
           if(ib.ne.10) then
             call rad_lw_nasa_h2okdis(ipts,ib,m,np,k2-1,                       &
                                           fkw,gkw,ne,h2oexp,conexp,           &
                                           th2o,tcon,trant)
           endif
         endif
!
         if (co2bnd) then
           if (high) then
!
!-----compute co2 transmittance using table look-up method
!
             w1=-4.0
             p1=-2.0
             dwe=0.3
             dpe=0.2
             call rad_lw_nasa_tablup(ipts,k1,k2,m,np,nx,nc,nt,                 &
                                          sco3,scopre,scotem,                  &
                                          w1,p1,dwe,dpe,cc1,c2,c3,trant)
           else
!
!-----compute co2 transmittance using k-distribution method
!
             call rad_lw_nasa_co2kdis(ipts,m,np,k2-1,co2exp,tco2,trant)
           endif
!
         endif
!
!-----compute o3 transmittance using table look-up
!
         if (oznbnd) then
           w1=-6.0
           p1=-2.0
           dwe=0.3
           dpe=0.2
           call rad_lw_nasa_tablup(ipts,k1,k2,m,np,nx,no,nt,                   &
                                       sco3,scopre,scotem,                     &
                                       w1,p1,dwe,dpe,oo1,o2,o3,trant)
         endif
!
!***** for trace gases ****!
!
         if (trace) then
!
!-----compute n2o transmittance using k-distribution method
!
           if (n2obnd) then
             call rad_lw_nasa_n2okdis(ipts,ib,m,np,k2-1,n2oexp,tn2o,trant)
           endif
!
!-----compute ch4 transmittance using k-distribution method
!
           if (ch4bnd) then
             call rad_lw_nasa_ch4kdis(ipts,ib,m,np,k2-1,ch4exp,tch4,trant)
           endif
!
!-----compute co2-minor transmittance using k-distribution method
!
           if (combnd) then
             call rad_lw_nasa_comkdis(ipts,ib,m,np,k2-1,comexp,tcom,trant)
           endif
!
!-----compute cfc11 transmittance using k-distribution method
!
           if (f11bnd) then
             call rad_lw_nasa_cfckdis(ipts,m,np,k2-1,f11exp,tf11,trant)
           endif
!
!-----compute cfc12 transmittance using k-distribution method
!
           if (f12bnd) then
             call rad_lw_nasa_cfckdis(ipts,m,np,k2-1,f12exp,tf12,trant)
           endif
!
!-----compute cfc22 transmittance using k-distribution method
!
           if (f22bnd) then
             call rad_lw_nasa_cfckdis(ipts,m,np,k2-1,f22exp,tf22,trant)
           endif
!
!-----compute transmittance in band 10 using k-distribution method
!     here, trant is the change in transmittance due to n2o absorption
!
           if (b10bnd) then
             call rad_lw_nasa_b10kdis(ipts,m,np,k2-1,                          &
                                           h2oexp,conexp,co2exp,n2oexp,        &
                                           th2o,tcon,tco2,tn2o,trant)
           endif
         endif
!
! ======= end of tranmisttance computation ==============
!
!-----include aerosol effect
!
         do i = 1,ipts   
           tranal(i)=tranal(i)*exp(-1.66*taual(i,k2-1,ib))
           trant (i)=trant(i) *tranal(i)
         enddo
         do i = 1,ipts   
           if( k2 .le. ict ) then
             clrhi(i)=min(clr(i,k2-1),clrhi(i))
           elseif( k2 .gt. ict .and. k2 .le. icb ) then
             clrmd(i)=min(clr(i,k2-1),clrmd(i))
           elseif( k2 .gt. icb ) then
             clrlw(i)=min(clr(i,k2-1),clrlw(i))
           endif
           fclr(i)=clrlw(i)*clrmd(i)*clrhi(i)
         enddo
!
!-----compute upward and downward fluxes (band 1-9).
!     add "boundary" terms to the net downward flux.
!     these are the first terms on the right-hand-side of
!     eqs. (56a) and (56b).  downward fluxes are positive.
!
         if (.not. b10bnd) then
           if (k2 .eq. k1+1) then
             do i = 1,ipts   
!
!-----compute net downward clear-sky fluxes
!
               flc(i,k1)=flc(i,k1)-blayer(i,k1)
               flc(i,k2)=flc(i,k2)+blayer(i,k1)
!
!-----compute upward and downward fluxes for all-sky situation
!
               flxu(i,k1)=flxu(i,k1)-blayer(i,k1)
               flxd(i,k2)=flxd(i,k2)+blayer(i,k1)
!
!-----compute upward and downward fluxes for clear-sky situation
!
               flxuc(i,k1)=flxuc(i,k1)-blayer(i,k1)
               flxdc(i,k2)=flxdc(i,k2)+blayer(i,k1)
             enddo
           endif
!
!-----add flux components involving the four layers above and below
!     the levels k1 and k2.  it follows eqs. (56a) and (56b).
!
           do i = 1,ipts   
             xx=trant(i)*(blayer(i,k2-1)-blayer(i,k2))
             flc(i,k1) =flc(i,k1)+xx
             flxu(i,k1) =flxu(i,k1)+xx*fclr(i)
             flxuc(i,k1) =flxuc(i,k1)+xx
             xx=trant(i)*(blayer(i,k1-1)-blayer(i,k1))
             flc(i,k2) =flc(i,k2)+xx
             flxd(i,k2) =flxd(i,k2)+xx*fclr(i)
             flxdc(i,k2) =flxdc(i,k2)+xx
           enddo
         else
!
!-----flux reduction due to n2o in band 10
!
           if (trace) then
             do i = 1,ipts   
               rflx(i,k1) = rflx(i,k1)+trant(i)*fclr(i)                        &
                        *(blayer(i,k2-1)-blayer(i,k2))
                     flxu(i,k1) =flxu(i,k1)+trant(i)*fclr(i)                   &
                        *(blayer(i,k2-1)-blayer(i,k2))
                     flxuc(i,k1) =flxuc(i,k1)+trant(i)                         &
                        *(blayer(i,k2-1)-blayer(i,k2))
                     rflx(i,k2) = rflx(i,k2)+trant(i)*fclr(i)                  &
                        *(blayer(i,k1-1)-blayer(i,k1))
                     flxd(i,k2) =flxd(i,k2)+trant(i)*fclr(i)                   &
                        *(blayer(i,k1-1)-blayer(i,k1))
                     flxdc(i,k2) =flxdc(i,k2)+trant(i)                         &
                        *(blayer(i,k1-1)-blayer(i,k1))
                     rflc(i,k1) = rflc(i,k1)+trant(i)                          &
                        *(blayer(i,k2-1)-blayer(i,k2))
                     rflc(i,k2) = rflc(i,k2)+trant(i)                          &
                        *(blayer(i,k1-1)-blayer(i,k1))
             enddo
           endif
!
         endif
!
       enddo loop_3000
!
!-----compute the partial derivative of fluxes with respect to
!     surface temperature (eq. 59).
!
       if (trace .or. (.not. b10bnd) ) then
         do i = 1,ipts   
           dfdts(i,k1) =dfdts(i,k1)-dbs(i)*trant(i)*fclr(i)
         enddo
       endif
!
     enddo loop_2000
!
!-----add contribution from the surface to the flux terms at the surface.
!
     if (.not. b10bnd) then
       do i = 1,ipts   
         flc(i,np+1)=flc(i,np+1)-blayer(i,np+1)
         flxu(i,np+1)=flxu(i,np+1)-blayer(i,np+1)
         flxuc(i,np+1)=flxuc(i,np+1)-blayer(i,np+1)
         st4(i)=st4(i)-blayer(i,np+1)
         dfdts(i,np+1)=dfdts(i,np+1)-dbs(i)
       enddo
     endif
!
   enddo loop_1000
!
!-----adjust fluxes due to n2o absorption in band 10
!
   do k = 1,np+1
     do i = 1,ipts   
       flx(i,k)=flxd(i,k)+flxu(i,k)
       flc(i,k)=flc(i,k)+rflc(i,k)
     enddo
   enddo
!
   return
   end subroutine rad_lw_nasa
#endif /* ~LWRMDC end */
