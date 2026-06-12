#include <define.h>
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [rad_main_solver]
!      |
!      |--- [rad_cloud_comp] *
!      |--- [rad_cld_property_ice] *
!      |--- [rad_cld_property_nasa] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   subroutine rad_cloud_comp(idimt,prslv,t,cldary,ibeg,ipts,xlatrd,            &
                     ktop,kbtm,nclds,cldlw,taucl,cfac,cldsw,                   &
                     taulw,cld0)
!-------------------------------------------------------------------------------
!
! subroutine: rad_cloud_comp
!
! program history log:
!   1993-02-00  y.h.         cloud radiative properties calculations
!                            after davis (1982)
!                            and harshvardhan et al. (1987)
!   1995-11-00  y.h.         modified to provide mixed cloud overlapping scheme
!                            for chous sw radiation scheme (1992)
!   2000-03-00  fanglin yang
!                            to use chou_uiuc's lw scheme, save
!                            taulw(i,l)   - cloud optical depth for each layer, 
!                              k=1-->top layer
!                            cld0 (i,l) - unweighted cloud fraction, 
!                              k=1-->top layer
!
!-------------------------------------------------------------------------------
!     input variables:
!        prslv(i,k)    - level pressure (mb) k=1 is toa
!        t (i,k)       - absolute temperature, k=1 is top layer (k)
!        cldary(i,k)   - cloud array contains compressed cloud
!                        fractions of 3 types (stratiform, conv
!                        and stratus), k=1 is the mdl sfc layer
!        ibeg,ipts     - indices for the beginnig no. and the
!                        total no. of array elements to be processed
! ---   modify xlatrd to general for regional and global (h.-m.h. juangk
! ***    xlatrd        - current latitude in radians (1st data pt)
!                         for models with diff lat at each pt, need to
!                         use the lat of each point....careful.....
!    output variables:
! --- clouds for lw rad. k=1 is the sfc, k=2 is the
!      lowest cloud layer, and so on
!        ktop,kbtm(i,k)- cloud top and bottom indeces, ktop and
!                        kbtm values from 1 to l model layers,
!                        with value of 1 being the top mdl layer
!        nclds(i)      - no. of separated cloud layers in a column
!        cldlw(i,k)    - cloud fractions for lw, emissivity adjusted
!        emis(i,k)     - cloud emissivity
!  ***   ityp(i,k)     - type of clouds, ityp=1, and 2 are for
!                        the rh, and conv types
! --- clouds for sw rad. k=1 is the top layer or level
!        taucl(i,k)    - cloud optical depth in every model layer
!                        k=1 is the top layer
!        cfac(i,k)     - feaction of clear sky view at the layer interface
!        cldsw(i,k)    - layer cloud fractions for sw
!
!-------------------------------------------------------------------------------
   use paramodel
   use rdparm
   use comcd1
!-------------------------------------------------------------------------------
   real      ::  prslv(imbx,lp1), cldary(idimt,l), xlatrd(idimt),              &
                 taucl(imbx,l),   t    (imbx,lp1), cldlw(imbx,lp1),            &
                 cldsw(imbx,l),    cfac(imbx,lp1), taulw(imbx,l),              &
                 cld0(imbx,l)
   integer   ::  ktop(imbx,lp1),  kbtm(imbx,lp1),                              &
                 ityp(imbx,lp1),  nclds(imax)
!
! ---  workspace ---
!
   real      ::  xamt (imax),tauc (imax),cl1 (imax),                           &
                 cl2  (imax),alfa (imax)
   integer   ::  mtyp (imax),kcld (imax),mbtm(imax)
!97  3, cl2  (imax),    ptopd(imax),    alfa (imax)
   logical   ::  bitx(imax),bitw(imax),bit1,bit2
!
!===>    begin here ................................................
!
   do i = 1,imax
     kcld(i) = 2
     mbtm(i) = 1
     mtyp(i) = 0
     xamt(i) = 0.0e0
     ityp(i,1) = 0
     ktop(i,1) = lp1
     kbtm(i,1) = lp1
     cldlw(i,1)= 1.0e0
     cfac(i,1) = 1.0e0
   enddo
   do k = 2,lp1
     do i = 1,imax
       ityp(i,k) = 0
       ktop(i,k) = 1
       kbtm(i,k) = 1
       cldlw(i,k)   = 0.0e0
!       cld0 (i,k)   = 0.0e0
       cldsw(i,k-1) = 0.0e0
       taucl(i,k-1) = 0.0e0
       taulw(i,k-1) = 0.0e0
       cld0 (i,k-1)   = 0.0e0
       cfac(i,k) = 1.0e0
     enddo
   enddo
!
! --- loop over mdl layers (bottom up)
!
   do 200 k = 2,l
!
     bit1 = .false.
     do i = 1,ipts
       ir = i + ibeg - 1
       bitx(i) = cldary(ir,k).gt.0.0e0
       bit1 = bit1 .or. bitx(i)
     enddo
     if (.not. bit1) go to 200
!
! --- decompress cloud array
!
     do i = 1,ipts
       cl1(i) = 0.0e0
       cl2(i) = 0.0e0
       bitw(i) = bitx(i)
     enddo
     do i = 1,ipts
       if (bitx(i)) then
         ir = i + ibeg - 1
         cl1(i) = amod(cldary(ir,k), 2.0e0)
         cltemp = amod(cldary(ir,k), 10.0e0)
         cl2(i) = 1.0e-4 * (cldary(ir,k) - cltemp)
!
! --- mtyp=1,2 for rh+stratus, and conv cloud types
!
         if (cl2(i) .gt. 0.0e0) then
           mtyp(i) = 2
         else
           mtyp(i) = 1
         endif
       endif
     enddo
     if(k.lt.l) then
       do i = 1,ipts
         ir = i + ibeg - 1
         if(bitw(i)) then
           bitw(i) = cldary(ir,k+1).le.0.0e0
         endif
       enddo
     endif
     bit2 = .false.
     do i = 1,ipts
       bit2 = bit2 .or. bitw(i)
       if (bitx(i)) then
         if(ityp(i,kcld(i)).eq.0) then
           ityp(i,kcld(i)) = mtyp(i)
           xamt(i) = cl1(i)
           if (mtyp(i) .eq. 2) xamt(i) = cl2(i)
           mbtm(i) = k
         else if(ityp(i,kcld(i)).ne.mtyp(i) .or.                               &
                (mtyp(i).eq.2 .and. xamt(i).ne.cl2(i)) ) then
           cldlw(i,kcld(i)) = xamt(i)
           ktop(i,kcld(i)) = lp1 - (k - 1)
           kbtm(i,kcld(i)) = lp1 - mbtm(i)
           ityp(i,kcld(i)+1) = mtyp(i)
           mbtm(i) = k
           xamt(i) = cl1(i)
           if (mtyp(i).eq.2) xamt(i) = cl2(i)
           kcld(i) = kcld(i) + 1
         else if(mtyp(i).eq.1) then
           xamt(i) = amax1(xamt(i), cl1(i))
         endif
       end if
     enddo
     if (.not. bit2) go to 200
     do i = 1,ipts
       if (bitw(i)) then
         cldlw(i,kcld(i)) = xamt(i)
         ktop(i,kcld(i)) = lp1 - k
         kbtm(i,kcld(i)) = lp1 - mbtm(i)
         kcld(i) = kcld(i) + 1
         mtyp(i) = 0
         mbtm(i) = 1
         xamt(i) = 0.0e0
       end if
     enddo
!
200 continue
!
! --- record num of cld lyrs and find max num of cld lyrs
!
   mclds = 0
   do i = 1,ipts
     nclds(i) = kcld(i) - 2
     mclds = max(mclds, nclds(i))
   enddo
!
! --- estimate cloud optical properties from t and q
!     (top down)
!
   do nncld = 1,mclds
     nc = mclds - nncld + 2
     do i = 1,ipts
       tauc(i) = 0.0e0
       bitx(i) = cldlw(i,nc) .gt. 0.0e0
       bitw(i) = bitx(i)
!
!conv - reduce conv cloud amount for sw rad
!
       if (ityp(i,nc) .eq. 2) then
         alfa(i) = amax1(0.25e0,                                               &
                   1.0e0-0.125e0*(kbtm(i,nc)-ktop(i,nc)))
       else
         alfa(i) = 1.0e0
       end if
     enddo
!
! --- find top pressure for mid cloud (3) domain=function of latitude
!
     minktp=1
     maxkbt=l
     if (nncld .eq. 1) kstrt = minktp
!
! --- calc cld thickness delp and mean temp (celsius)
!
     do kk = minktp,maxkbt
       do i = 1,ipts
         if (kk.ge.ktop(i,nc) .and. kk.le.kbtm(i,nc) .and.bitx(i)) then
           delp = prslv(i,kk+1) - prslv(i,kk)
           tcld = t(i,kk) - 273.16e0
!
! --- convective cloud
!
           if (ityp(i,nc) .eq. 2) then
             tau0 = delp * 0.06e0
!
! --- rh clouds
!
           else
             if (tcld .le. -10.0e0) then
               tau0 = delp                                                     &
                      * amax1(0.1e-3, 2.00e-6*(tcld+82.5e0)**2)
             else
               tau0 = delp*amin1(0.08e0, 6.949e-3*tcld+0.08e0)
             end if
           end if
           tauc(i) = tauc(i) + tau0
           cldsw(i,kk) = cldlw(i,nc)
           taucl(i,kk) = tau0 * alfa(i) * cldlw(i,nc)
           cld0 (i,kk) = cldlw(i,nc)
           taulw(i,kk) = tau0
           if (bitw(i)) then
             cfac(i,kk+1) = cfac(i,kk) * (1.0e0 - cldsw(i,kk))
             bitw(i) = .false.
           else
             cfac(i,kk+1) = cfac(i,kk)
           end if
         elseif (kk.gt.kbtm(i,nc) .and. bitx(i)) then
           cfac(i,kk+1) = cfac(i,kk)
         end if
       enddo
     enddo
     mkbtp1 = maxkbt + 1
     do k = mkbtp1,l
       do i = 1,ipts
         if (bitx(i)) cfac(i,k+1) = cfac(i,mkbtp1)
       enddo
     enddo
!
! --- calc cld emis
!
     do i = 1,ipts
       if (bitx(i))                                                            &
         cldlw(i,nc) = cldlw(i,nc)*(1.0e0-exp(-0.75e0*tauc(i)))
     enddo
!
   enddo
!
! --- cloud scaled for sw
!
   do kk = 1,l
     do i = 1,ipts
       if (cfac(i,lp1) .lt. 1.0e0)                                             &
         taucl(i,kk) = taucl(i,kk) / (1.0e0 - cfac(i,lp1))
     enddo
   enddo
   if (ipts .eq. imax) go to 565
   ipts1 = ipts + 1
   do i = ipts1,imax
     nclds(i) = nclds(ipts)
   enddo
   do k = 1,lp1
     do i = ipts1,imax
       cldlw(i,k) = cldlw(ipts,k)
       ktop(i,k) = ktop(ipts,k)
       kbtm(i,k) = kbtm(ipts,k)
       cfac(i,k) = cfac(ipts,k)
     enddo
   enddo
   do k = 1,l
     do i = ipts1,imax
       taulw(i,k) = taulw(ipts,k)
       cld0 (i,k) = cld0 (ipts,k)
       taucl(i,k) = taucl(ipts,k)
       cldsw(i,k) = cldsw(ipts,k)
     enddo
   enddo
565 continue
!
   return
   end subroutine rad_cloud_comp
!
!-------------------------------------------------------------------------------
   subroutine rad_cld_property_ice(idimt,prsi,prsl,prslv,                      &
                     t,cldtot,cldcnv,ibeg,ipts,lat,                            &
                     icwp,clwp,slmsk,xlat,ktop,kbtm,nclds,                     &
                     cldlw,taucl,cfac,cldsw,                                   &
                     cwp,cip,rew,rei,fice,                                     &
                     taulw,cld0)
!-------------------------------------------------------------------------------
!
! program history log:
!   1993-02-00  y.h.
!        cloud radiative properties calculations after davis (1982)
!        and harshvardhan et al. (1987).
!   1995-11-00  y.h.
!        modified to provide mixed cloud overlapping scheme for
!        chou's sw radiation scheme (1992)_.
!   1998-08-00  y.h.
!        modified to estimate effective radius of water/ice cloud
!        drop, fraction of ice water content, and cloud emissivity
!        for lw radiation calculation (kiehl et al. 1998,j.clim).
!   2000-03-00  fanglin yang
!        to use chou_uiuc's lw scheme, save
!        taulw(i,l)   - cloud optical depth for each layer, k=1-->top layer
!        cld0 (i,l) - unweighted cloud fraction, k=1-->top layer
!
!-------------------------------------------------------------------------------
!     input variables:
!        prsl(i,k)     - mdl sigma layer mean value, k=1 is sfc layer
!        prslv(i,k)    - level pressure (mb)         k=1 is toa
!        t (i,k)       - absolute temperature, k=1 is top layer (k)
!        cldtot(i,k)   - stratiform cloud      k=1 is mdl sfc layer
!        cldcnv(i,k)   - convective cloud      k=1 is mdl sfc layer
!        ibeg,ipts     - indices for the beginnig no. and the
!                        total no. of array elements to be processed
! ---   modify xlatrd to general for regional and global (h.-m.h. juangk
! ***    xlat        - current latitude in radians (1st data pt)
!                         for models with diff lat at each pt, need to
!                         use the lat of each point....careful.....
!        icwp          - flag indicates the method used for cloud
!                        properties, =0 use t-p; =1 use clwp.
!        clwp          - layer cloud water+ice path (g/m**2), k=1 is toa
!        slmsk         - land/sea/ice mask (0:sea.1:land,2:ice)
!    output variables:
! --- clouds for lw rad. k=1 is the sfc, k=2 is the
!      lowest cloud layer, and so on
!        ktop,kbtm(i,k)- cloud top and bottom indeces, ktop and
!                        kbtm values from 1 to l model layers,
!                        with value of 1 being the top mdl layer
!        nclds(i)      - no. of separated cloud layers in a column
!        cldlw(i,k)    - cloud fractions for lw, emissivity adjusted
!  ***   ityp(i,k)     - type of clouds, ityp=1, and 2 are for
!                        the rh, and conv types
! --- clouds for sw rad. k=1 is the top layer or level
!        taucl(i,k)    - cloud optical depth in every model layer
!                        k=1 is the top layer (used for icwp=0 only)
!        cfac(i,k)     - fraction of clear sky view at the layer interfa
!        cldsw(i,k)    - layer cloud fractions for sw
!  - - following are output variables for icwp=1
!        cwp           - layer cloud water path (g/m**2)
!        cip           - layer cloud ice path (g/m**2)
!        rew (i,k)     - effective layer cloud water drop size (micron)
!        rei (i,k)     - effective layer cloud ice drop size (micron)
!        fice(i,k)     - fraction of cloud ice content
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only : igrd12_,jgrd12_
#else
   use paramodel, only : latg2_,lonf2_
#endif
   use rdparm, only    : imax,imbx,lp1,l
   use comcd1
   real      ::  cldtot(idimt,l), cldcnv(idimt,l),slmsk(idimt), xlat(idimt),   &
                 prsi(idimt,l+1), prsl(idimt,l), prslv(imax,lp1), t(imax,lp1)
!
! - output array
!
   real      ::  taucl(imax,l),  cldlw(imax,lp1),   cldsw(imax,l),             &
                 cfac(imax,lp1),  taulw(imbx,l),     cld0(imbx,l),             &
                 rew (imax,l),     rei (imax,l),     fice(imax,l),             &
                 clwp(imax,l),     cwp (imax,l),     cip (imax,l)
   integer   ::  ktop(imax,lp1),  kbtm(imax,lp1),                              &
                 ityp(imax,lp1),  nclds(imax)
!
! ---  workspace ---
!
   real      ::  xamt (imax), tauc (imax), cl1  (imax),                        &
                 cl2  (imax), alfa (imax), wgt
   integer   ::  mtyp (imax), kcld (imax), mbtm (imax) 
   logical   ::  bitx(imax), bity(imax), bitw(imax), bit1, bit2
!
!===>    begin here ................................................
!
   do i = 1,imax
     kcld(i)    = 2
     mbtm(i)    = 1
     xamt(i)    = 0.0e0
     ityp(i,1)  = 0
     ktop(i,1)  = lp1
     kbtm(i,1)  = lp1
     cldlw(i,1) = 1.0e0
     cfac(i,1)  = 1.0e0
   enddo
!
   do k = 2,lp1
     do i = 1,imax
       ityp(i,k)  = 0
       ktop(i,k)  = 1
       kbtm(i,k)  = 1
       cldlw(i,k) = 0.0e0
       cfac(i,k)  = 1.0e0
     enddo
   enddo
!
   do k = 1,l
     do i = 1,imax
       cldsw(i,k) = 0.0e0
       taucl(i,k) = 0.0e0
       taulw(i,k) = 0.0e0
       cld0(i,k)  = 0.0e0
       cwp (i,k)  = 0.0e0
       cip (i,k)  = 0.0e0
       rew (i,k)  = 10.0e0
       rei (i,k)  = 10.0e0
       fice(i,k)  = 0.0e0
     enddo
   enddo
!
!0499 if (icwp .eq. 1) then ! 0898 - use cloud water/ice content
!
   do k = 1,l
     do i = 1,ipts
       ir   = i + ibeg -1
       delt = 263.16 - t(i,k)
       if (nint(slmsk(ir)) .eq. 1) then
         rew(i,k) = 5.0e0 + 5.0e0*amin1(1.0e0, amax1(0.0e0,                    &
                    delt*0.05e0 )) ! - effective radius for water
       end if
       wgt = amin1(1.0e0, amax1(0.0e0,                                      &
                  (prsl(ir,lp1-k)/prsi(ir,1) - 0.3e0) / 0.5e0 ))
       rei(i,k) = 50.0 - 40.0*wgt     ! - effective radius for ice
       fice(i,k) = amin1(1.0, amax1(0.0e0, delt/20.0e0 ))
                                   ! - fraction of ice in cloud
       cip (i,k) = clwp(i,k) * fice(i,k) 
       cwp (i,k) = clwp(i,k) - cip(i,k)
     enddo
   enddo
!
!0499 end if
!
! --- loop over mdl layers (bottom up)
!
   do k = 2,l
!
     do i = 1,ipts
       ir      = i + ibeg - 1
       cl1(i)  = cldtot(ir,k)
       cl2(i)  = cldcnv(ir,k)
     enddo
!
     bit1 = .false.
     do i = 1,ipts
       bitx(i) = cl1(i).gt.0.001e0 .or. cl2(i).gt.0.001e0
       bit1    = bit1 .or. bitx(i)
     enddo
     if (bit1) then
       do i = 1,ipts
!
! --- mtyp=1,2 for rh+stratus, and conv cloud types
!
         if (cl2(i) .gt. 0.001e0) then
           mtyp(i) = 2
         else if (cl1(i) .gt. 0.001e0) then
           mtyp(i) = 1
         else
           mtyp(i) = 0
         endif
       enddo
!
       if (k .lt. l) then
         do i = 1,ipts
           ir = i + ibeg - 1
!
! --- set bitw for clear gap above cloud
!
           bitw(i) = bitx(i) .and. cldtot(ir,k+1).le.0.001e0                   &
                             .and. cldcnv(ir,k+1).le.0.001e0
         enddo
       endif
       bit2 = .false.
       do i = 1,ipts
         bit2 = bit2 .or. bitw(i)
         if (bitx(i)) then
           kkcl = kcld(i)
           if(ityp(i,kkcl) .eq. 0) then
             ityp(i,kkcl) = mtyp(i)
             xamt(i) = cl1(i)
             if (mtyp(i) .eq. 2) xamt(i) = cl2(i)
             mbtm(i) = k
           else if(ityp(i,kkcl).ne.mtyp(i) .or.                                &
                  (mtyp(i).eq.2 .and. xamt(i).ne.cl2(i)) ) then
             cldlw(i,kkcl)  = xamt(i)
             ktop(i,kkcl)   = lp1 - (k - 1)
             kbtm(i,kkcl)   = lp1 - mbtm(i)
             ityp(i,kkcl+1) = mtyp(i)
             mbtm(i)        = k
             xamt(i)        = cl1(i)
             if (mtyp(i).eq.2) xamt(i) = cl2(i)
             kcld(i) = kkcl + 1
           else if(mtyp(i).eq.1) then
             xamt(i) = amax1(xamt(i), cl1(i))
           endif
         endif
       enddo
       if (bit2) then
         do i = 1,ipts
           if (bitw(i)) then
             kkcl = kcld(i)
             cldlw(i,kkcl) = xamt(i)
             ktop(i,kkcl)  = lp1 - k
             kbtm(i,kkcl)  = lp1 - mbtm(i)
             kcld(i)       = kkcl + 1
             mtyp(i)       = 0
             mbtm(i)       = 1
             xamt(i)       = 0.0e0
           endif
         enddo
       endif                             ! bit2 endif
     endif                               ! bit1 endif
   enddo                                 ! the k loop ends here!
!
! --- record num of cld lyrs and find max num of cld lyrs
!
   mclds = 0
   do i = 1,ipts
     nclds(i) = kcld(i) - 2
     mclds    = max(mclds, nclds(i))
   enddo
!
!     write(6,221) mclds, ibeg
!221  format(' in rad_cloud_comp: maxclds =',i4,' ibeg=',i4)
!     if (mclds .eq. 0) return
!
! --- estimate cloud optical properties from t and q  -- (top down)
!
   do nncld = 1,mclds
     nc = mclds - nncld + 2
     do i = 1,ipts
       bitx(i) = cldlw(i,nc) .ge. 0.001e0
!         bity(i) = bitx(i) .and. ityp(i,nc).eq.2
       bitw(i) = bitx(i)
     enddo
!
! --- find top pressure for mid cloud (3) domain=function of latitude
!
     minktp=1
     maxkbt=l
!
!       write(6,241) nc,minktp, maxkbt
!241    format(3x,'nc, minktp, maxkbt =',3i6/3x,'bitx, bitw :')
!
! --- find clear sky view at each levels
!
     do kk = minktp,maxkbt
       do i = 1,ipts
         if (kk.ge.ktop(i,nc) .and.                                            &
             kk.le.kbtm(i,nc) .and. bitx(i)) then
           cldsw(i,kk) = cldlw(i,nc)
           if (bitw(i)) then
             cfac(i,kk+1) = cfac(i,kk) * (1.0e0 - cldsw(i,kk))
             bitw(i) = .false.
           else
             cfac(i,kk+1) = cfac(i,kk)
           end if
         elseif (kk.gt.kbtm(i,nc) .and. bitx(i)) then
           cfac(i,kk+1) = cfac(i,kk)
         endif
       enddo
     enddo
!
     mkbtp1 = maxkbt + 1
     do k = mkbtp1,l
       do i = 1,ipts
         if (bitx(i)) cfac(i,k+1) = cfac(i,mkbtp1)
       enddo
     enddo
!
     do i = 1,ipts
       tauc(i) = 0.0e0
     enddo
!
     if (icwp .ne. 1) then
!
!conv - reduce conv cloud amount for sw rad
!
       do i = 1,ipts
         if (ityp(i,nc) .eq. 2) then
!0799         alfa(i) = cldlw(i,nc)
           alfa(i) =                                                           &
               amax1(0.25e0, 1.0e0-0.125e0*(kbtm(i,nc)-ktop(i,nc)))
         else
!             alfa(i) = sqrt(cldlw(i,nc))
           alfa(i) = 1.0
         endif
       enddo
!
! --- calc cld thickness delp and mean temp (celsius)
!
       do kk = minktp,maxkbt
         do i = 1,ipts
           if (kk.ge.ktop(i,nc) .and.                                          &
               kk.le.kbtm(i,nc) .and. bitx(i)) then
             delp = prslv(i,kk+1) - prslv(i,kk)
             tcld = t(i,kk) - 273.16e0
!
! --- convective cloud
!
             if (ityp(i,nc) .eq. 2) then
               tau0 = delp * 0.05e0
!
!0499 - if conv cld, set to water cloud only
!
               fice(i,kk) = 0.0e0
!
! --- rh clouds
!
             else
               if (tcld .le. -10.0e0) then
                 tau0 = delp                                                   &
                      * amax1(0.1e-3, 2.00e-6*(tcld+82.5e0)**2)
               else
                 tau0 = delp                                                   &
                      * amin1(0.08e0, 6.949e-3*tcld+0.08e0)
               endif
             endif
             taulw(i,kk) = tau0
             tauc(i)     = tauc(i) + tau0
             taucl(i,kk) = tau0 * alfa(i)
           endif
         enddo
       enddo
!
! --- calc cld emis and effective cloud cover for lw
!
       do i = 1,ipts
         if (bitx(i))                                                          &
           cldlw(i,nc) = cldlw(i,nc)*(1.0e0-exp(-0.75e0*tauc(i)))
       enddo
       do kk = minktp,maxkbt
         do i = 1,ipts
           if (bitx(i)) cld0 (i,kk) = cldlw(i,nc)
         enddo
       enddo
!
     else                         ! prognostic liquid water
!
! --- calc tauc and emis
!
       do kk = minktp,maxkbt
         do i = 1,ipts
           if (kk .ge. ktop(i,nc) .and.                                        &
               kk .le. kbtm(i,nc) .and. bitx(i)) then
!
!0499 - if conv cld, set to liquid content only
!
             if (ityp(i,nc) .eq. 2) then
               fice(i,kk) = 0.0e0
               cwp (i,kk) = cwp(i,kk) + cip(i,kk)
               cip (i,kk) = 0.0e0
             end if
             awlw    = 0.090361e0 * cwp(i,kk)
             ailw    = (0.005e0 + 1.0e0/rei(i,kk)) * cip(i,kk)
             tauc(i) = tauc(i) + awlw + ailw
             taulw(i,kk)=awlw+ailw
             taucl(i,kk)=awlw
           endif
         enddo
       enddo
       do i = 1,ipts
         if (bitx(i))                                                          &
!             cldlw(i,nc) = cldlw(i,nc) * amin1(1.0e0, amax1(0.2e0,
           cldlw(i,nc) = cldlw(i,nc) * amin1(1.0e0,                            &
                         1.0e0 - exp(-1.66e0*tauc(i)) )
       enddo
       do kk = minktp,maxkbt
         do i = 1,ipts
            if (bitx(i)) cld0 (i,kk) = cldlw(i,nc)
         enddo
       enddo
     endif
!
   enddo               ! end of nncld loop!
!
   if (icwp .eq. 1) then
     do k = 1,l
       do i = 1,ipts
         if (cldsw(i,k) .lt. 0.001e0) then
           taucl(i,k) = 0.0e0
           cldsw(i,k) = 0.0e0
           cwp  (i,k) = 0.0e0
           cip  (i,k) = 0.0e0
         endif
       enddo
     enddo
   endif
!
   if (ipts .lt. imax) then     ! --- fill up arrays
     ipts1 = ipts + 1
     do i = ipts1,imax
       nclds(i) = nclds(ipts)
     enddo
     do k = 1,lp1
       do i = ipts1,imax
         cldlw(i,k) = cldlw(ipts,k)
         ktop(i,k) = ktop(ipts,k)
         kbtm(i,k) = kbtm(ipts,k)
         cfac(i,k) = cfac(ipts,k)
       enddo
     enddo
     do k = 1,l
       do i = ipts1,imax
         taulw(i,k) = taulw(ipts,k)
         cld0 (i,k) = cld0 (ipts,k)
         taucl(i,k) = taucl(ipts,k)
         cldsw(i,k) = cldsw(ipts,k)
         cwp  (i,k) = cwp  (ipts,k)
         cip  (i,k) = cip  (ipts,k)
         rew  (i,k) = rew  (ipts,k)
         rei  (i,k) = rei  (ipts,k)
         fice (i,k) = fice (ipts,k)
       enddo
     enddo
   endif
!
#ifdef DBG
#ifndef MP
#ifdef RMP
   latp = jgrd12_
   ir = igrd12_/4   ! center point
#else
   latp = latg2_
   ir = lonf2_/4    ! (180,0) point
#endif
#ifdef SMP
   i = 1
   ir = 1
#else
   i   = ir - ibeg +1
   if(i.ge.1.and.i.le.ipts.and.lat.eq.latp) then
#endif
     write(6,501) ir,lat
     write(6,502)
501  format(//,'=== cloud properties in rad_cld_property_ice at (i,lat) ',     & 
            2i10,'   ===',/)
502  format('    k','      clwp','      fice','       cwp',                    &
            '       cip','       rew', '       rei','      ityp')
503  format(i5,6f10.5,i10)
     do k = 1,l
       write(6,503) k,clwp(i,k),fice(i,k),cwp(i,k),                            &
            cip(i,k),rew(i,k),rei(i,k),ityp(i,k)
     enddo
!
     write(6,505)
505  format(//,'    k','   tem (C)','     cldlw',' ktop',' kbtm',              &
            '     cldsw','      cfac','    cldtot','    cldcnv','     taulw')
506  format(i5,2f10.5,2i5,5f10.5)
     do k = 1,l
       write(6,506) k,t(i,k)-273.16,cldlw(i,k),ktop(i,k),kbtm(i,k),            &
         cldsw(i,k),cfac(i,k),cldtot(i,l-k+1),cldcnv(i,l-k+1),                 &
         taulw(i,l-k+1)
     enddo
#ifndef SMP
   endif
#endif
#endif
#endif
!
   return
   end subroutine rad_cld_property_ice
!
!-------------------------------------------------------------------------------
   subroutine rad_cld_property_nasa(idimt,prslv,t,cldtot,cldcnv,ibeg,ipts,lat, & 
                     icwp,clwp,slmsk,xlat,ktop,kbtm,nclds,                     &
                     cldlw,taucl,cfac,cldsw,                                   &
                     cwp,cip,rew,rei,fice,                                     &
                     taulw,cld0)
!-------------------------------------------------------------------------------
!
! program history log:
!   1993-02-00  y.h.
!        cloud radiative properties calculations after davis (1982)
!        and harshvardhan et al. (1987).
!   1995-11-00  y.h.
!        modified to provide mixed cloud overlapping scheme for
!        chou's sw radiation scheme (1992)_.
!   1998-08-00  y.h.
!        modified to estimate effective radius of water/ice cloud
!        drop, fraction of ice water content, and cloud emissivity
!        for lw radiation calculation (kiehl et al. 1998,j.clim).
!   2000-03-00  fanglin yang
!        to use chou_uiuc's lw scheme, save
!        taulw(i,l)   - cloud optical depth for each layer, k=1-->top layer
!        cld0 (i,l) - unweighted cloud fraction, k=1-->top layer
!   2000-08-00  y.h.
!        modified clear sky view 'cfac', use combined max/ran overlap
!        scheme on contiguous clouds. 'cldsw' will keep original
!        values. 'cldlw' will break up if difference between
!        two contiguous cloud is significant large enough.
!
!--------------------------------------------------------------------
!     input variables:
!        prslv(i,k)    - level pressure (mb)         k=1 is toa
!        t (i,k)       - absolute temperature, k=1 is top layer (k)
!        cldtot(i,k)   - stratiform cloud      k=1 is mdl sfc layer
!        cldcnv(i,k)   - convective cloud      k=1 is mdl sfc layer
!        ibeg,ipts     - indices for the beginnig no. and the
!                        total no. of array elements to be processed
! ---   modify xlatrd to general for regional and global (h.-m.h. juangk
! ***    xlat        - current latitude in radians (1st data pt)
!                         for models with diff lat at each pt, need to
!                         use the lat of each point....careful.....
!        icwp          - flag indicates the method used for cloud
!                        properties, =0 use t-p; =1 use clwp.
!        clwp          - layer cloud water+ice path (g/m**2), k=1 is toa
!        slmsk         - land/sea/ice mask (0:sea.1:land,2:ice)
!    output variables:
! --- clouds for lw rad. k=1 is the sfc, k=2 is the
!      lowest cloud layer, and so on
!        ktop,kbtm(i,k)- cloud top and bottom indeces, ktop and
!                        kbtm values from 1 to l model layers,
!                        with value of 1 being the top mdl layer
!        nclds(i)      - no. of separated cloud layers in a column
!        cldlw(i,k)    - cloud fractions for lw, emissivity adjusted
! --- clouds for sw rad. k=1 is the top layer or level
!        taucl(i,k)    - cloud optical depth in every model layer
!                        k=1 is the top layer (used for icwp=0 only)
!        cfac(i,k)     - fraction of clear sky view at the layer interfa
!        cldsw(i,k)    - layer cloud fractions for sw
!  ***   ityp(i,k)     - type of clouds, ityp=1, and 2 are for
!                        the rh, and conv types
!  - - following are output variables for icwp=1
!        cwp           - layer cloud water path (g/m**2)
!        cip           - layer cloud ice path (g/m**2)
!        rew (i,k)     - effective layer cloud water drop size (micron)
!        rei (i,k)     - effective layer cloud ice drop size (micron)
!        fice(i,k)     - fraction of cloud ice content
!
!-------------------------------------------------------------------------------
#ifdef DBG
#ifndef MP
#ifdef RMP
   use paramodel, only : igrd12_,jgrd12_
#else
   use paramodel, only : lonf2_,latg2_
#endif
#endif
#endif
   use constant, only  : pi_
   use rdparm99
   use comcd1
!-------------------------------------------------------------------------------
   real               ::  cldtot(idimt,l), cldcnv(idimt,l), slmsk(idimt),      &
                          prslv(imbx,lp1),    t (imbx,lp1),  xlat(idimt)
!
! - output array
!
   real               ::  taucl(imbx,l),  cldlw(imbx,lp1), cldsw(imbx,l),      &
                          cfac(imbx,lp1),   taulw(imbx,l),  cld0(imbx,l),      &
                            rew (imbx,l),    rei (imbx,l),  fice(imbx,l),      &
                            clwp(imbx,l),    cwp (imbx,l),  cip (imbx,l)     
   integer            ::  ktop(imbx,lp1),  kbtm(imbx,lp1),                     &
                          ityp(imbx,lp1),     nclds(imax)
!
! ---  workspace ---
!
   real               ::  xamt (imax),     tauc (imax),                        &
                          adlw (imbx,l), adsw (imbx,l),                        &
                          topl (imax),     topm (imax)
   integer            ::  mtyp (imax),     kcld (imax),      mbtm (imax)
   logical            ::  bitx(imax),       bitw(imax),       bit1,            &
                          bit2
!
!===>    begin here ................................................
!
   do i = 1,imax
     kcld(i)    = 2
     mbtm(i)    = 1
     xamt(i)    = 0.0e0
     ityp(i,1)  = 0
     mtyp(i)    = 0
     ktop(i,1)  = lp1
     kbtm(i,1)  = lp1
     cldlw(i,1) = 1.0e0
     cfac(i,1)  = 1.0e0
   enddo
!
   do k = 2,lp1
     do i = 1,imax
       ityp(i,k)  = 0
       ktop(i,k)  = 1
       kbtm(i,k)  = 1
       cldlw(i,k) = 0.0e0
       cfac(i,k)  = 1.0e0
     enddo
   enddo
!
   do k = 1,l
     do i = 1,imax
       cldsw(i,k) = 0.0e0
       taucl(i,k) = 0.0e0
       taulw(i,k) = 0.0e0
       cld0(i,k)  = 0.0e0
       cwp (i,k)  = 0.0e0
       cip (i,k)  = 0.0e0
       rew (i,k)  = 10.0e0
       rei (i,k)  = 10.0e0
       fice(i,k)  = 0.0e0
     enddo
   enddo
!
! --- find top pressure of low and mid cloud domains
!
   dptl = ptopc(2,2)-ptopc(2,1)
   dptm = ptopc(3,2)-ptopc(3,1)
   do i = 1,ipts
     ir = i + ibeg -1
     fac = max (0.0e0, 4.0e0*abs(xlat(ir))/pi_-1.0e0)
     topl(i) = ptopc(2,1) + dptl * fac
     topm(i) = ptopc(3,1) + dptm * fac
     bitw(i) = nint(slmsk(ir)) .ne. 1        ! sea/ice points
   enddo
!
! --- cloud overlapping criteria and adjustment over land/sea/ice
! ..... vertically contiguous clouds are treated with ran/max overlap
! ... crtcl- criterion of applying random or max overlaping (0.2-0.5)
! ... adlw,adsw - adjustment over land/sea for lw and sw:
! ...             0 (more ran) - 0.5 (more max)
!
   crtcl = 0.3
   do k = 1,l
     do i = 1,ipts
       ir = i + ibeg -1
       adlw(i,k) = 0.0
       adsw(i,k) = 0.0
       plyr = 0.5 * (prslv(i,k) + prslv(i,k+1))
       if (bitw(i)) then
         if (plyr .ge. topl(i)) then         ! low cloud domain
           adlw(i,k) = 0.2
           adsw(i,k) = 0.2
!           elseif (plyr .ge. topm(i)) then     ! mid cloud domain
!             adlw(i,k) = 0.0
!             adsw(i,k) = 0.0
!           else                                ! high cloud domain
!             adlw(i,k) = 0.0
!             adsw(i,k) = 0.0
         endif
       endif
     enddo
   enddo
!
!0499 if (icwp .eq. 1) then ! 0898 - use cloud water/ice content
!
   reimax = 80.0
   reimin = 15.0
   tem    = (reimax-reimin)/(60.0-20.0)
   do k = 1,l
     do i = 1,ipts
       ir   = i + ibeg -1
       delt = 273.16 - t(i,k)
       if (nint(slmsk(ir)) .eq. 1) then
         rew(i,k) = 5.0e0 + 5.0e0*amin1(1.0e0, amax1(0.0e0,                    &
                    delt*0.05e0 )) ! - effective radius for water
       endif
       rei(i,k) = reimin + (60 - delt) * tem
       rei(i,k) = max(reimin, min(reimax, rei(i,k)))
!                                     ! - effective radius for ice
       fice(i,k) = min(1.0, max(0.0e0, delt/20.0e0 ))
                                   ! - fraction of ice in cloud
       cip (i,k) = clwp(i,k) * fice(i,k) 
       cwp (i,k) = clwp(i,k) - cip(i,k)
     enddo
   enddo
!0499 end if
!
! --- sw cloud structure according to cloud types
!
   do k = 2,l
     kr = lp1 - k                 ! vertical indices for sw
!
     do i = 1,ipts
       ir = i + ibeg - 1
       cldsw(i,kr) = max(cldcnv(ir,k),cldtot(ir,k))
!
! --- mtyp=1,2 for rh+stratus, and conv cloud types
!
       if (cldtot(ir,k) .gt. 0.0e0) ityp(i,kr) = 1
       if (cldcnv(ir,k) .gt. 0.0e0) ityp(i,kr) = 2
     enddo                       ! end i loop
   enddo                         ! end k loop
!
! --- loop over mdl layers (bottom up)
!     lw cloud structure
!
   do k = 2,l
     kr = lp1 - k                 ! vertical indices for sw
!
     bit1 = .false.
     do i = 1,ipts
       bitx(i) = ityp(i,kr) .gt. 0
       bit1    = bit1 .or. bitx(i)
       bitw(i) = bitx(i)
     enddo
     if (bit1) then
       if (k .lt. l) then
         do i = 1,ipts
!
! --- set bitw for clear gap above cloud
!
           bitw(i) = bitx(i) .and. ityp(i,kr-1).eq.0
         enddo
       endif
       bit2 = .false.
       do i = 1,ipts
         bit2 = bit2 .or. bitw(i)
         if (bitx(i)) then
           kkcl = kcld(i)
           cl12 = cldsw(i,kr)
           if(mtyp(i) .eq. 0) then
             mtyp(i) = ityp(i,kr)
             xamt(i) = cl12
             mbtm(i) = k
           else
             sdif = 2.0*abs(xamt(i)-cl12)/(xamt(i)+cl12)-adlw(i,k)
             if(ityp(i,kr).ne.mtyp(i) .or. sdif.ge.crtcl) then
               cldlw(i,kkcl)  = xamt(i)
               ktop(i,kkcl)   = lp1 - (k - 1)
               kbtm(i,kkcl)   = lp1 - mbtm(i)
               mtyp(i)        = ityp(i,kr)
               mbtm(i)        = k
               xamt(i)        = cl12
               kcld(i) = kkcl + 1
             else
               xamt(i) = max(xamt(i), cl12)
             endif
           endif
         endif
       enddo
       if (bit2) then
         do i = 1,ipts
           if (bitw(i)) then
             kkcl = kcld(i)
             cldlw(i,kkcl) = xamt(i)
             ktop(i,kkcl)  = lp1 - k
             kbtm(i,kkcl)  = lp1 - mbtm(i)
             kcld(i)       = kkcl + 1
             mtyp(i)       = 0
             mbtm(i)       = 1
             xamt(i)       = 0.0e0
           endif
         enddo
       endif                             ! bit2 endif
     endif                               ! bit1 endif
   enddo                                 ! the k loop ends here!
!
! --- define fraction of clear sky view at each level for sw
!     (top down)
!
   do i = 1,ipts
     bitw(i) = .true.
   end do
!
   do k = 2,l
     do i = 1,ipts
       if (bitw(i)) then
         cr1 = cldsw(i,k)
         cr2 = cldsw(i,k-1)
         crf = 1.0 - cr1
         if (crf .ge. 1.0) then
           cfac(i,k+1) = cfac(i,k)
         else if (crf .gt. 0.0) then
           if (cr2 .le. 0.0) then
             cfac(i,k+1) = cfac(i,k) * crf
           else
             cl12 = min(cfac(i,k), crf)
             sdif = min(0.80, max(0.20,                                        &
                    2.0*abs(cr1-cr2)/(cr1+cr2)-adsw(i,k) ))
             cfac(i,k+1) = (1.0-sdif)*cl12 + sdif*cfac(i,k)*crf
           endif
         else
           do kk = k+1,lp1
             cfac(i,kk) = 0.0
           enddo
           bitw(i) = .false.
         endif
       endif
     enddo                          ! end of i loop
   enddo                          ! end of k loop
!
! --- estimate sw cloud optical depth (diagnostic method)
!
   if (icwp .ne. 1) then
     do k = 1,l
       do i = 1,ipts
         if (ityp(i,k) .gt. 0) then
!
! --- calc cld thickness delp and mean temp (celsius)
!
           delp = prslv(i,k+1) - prslv(i,k)
           tcld = t(i,k) - 273.16e0
!
! --- convective cloud
!
           if (ityp(i,k) .eq. 2) then
             tau0 = delp * 0.06e0
             fice(i,k) = 0.0e0
!
! --- rh clouds
!
           else
             if (tcld .le. -10.0e0) then
               tau0 = delp                                                     &
                    * amax1(0.1e-3, 2.00e-6*(tcld+82.5e0)**2)
             else
               tau0 = delp                                                     &
                    * amin1(0.08e0, 6.949e-3*tcld+0.08e0)
             endif
           endif
           taucl(i,k) = tau0
           taulw(i,k) = tau0
         endif
       enddo                       ! end of i loop
     enddo                       ! end of k loop
   endif
!
! --- record num of cld lyrs and find max num of cld lyrs
!
   mclds = 0
   do i = 1,ipts
     nclds(i) = kcld(i) - 2
     mclds = max(mclds, nclds(i))
   enddo
!
! --- estimate lw cloud optical properties from t and q
!
   do nncld = 1,mclds
     nc = mclds - nncld + 2
!
     do i = 1,ipts
       tauc(i) = 0.0e0
       bitx(i) = cldlw(i,nc) .gt. 0.0e0
     enddo
!
     minktp=1
     maxkbt=l
!
     if (icwp .ne. 1) then
       do kk = minktp,maxkbt
         do i = 1,ipts
           if (kk.ge.ktop(i,nc) .and.                                          &
               kk.le.kbtm(i,nc) .and. bitx(i)) then
             tauc(i) = tauc(i) + taucl(i,kk)
           endif
         enddo
       enddo
!
! --- calc cld emis and effective cloud cover for lw
!
       do i = 1,ipts
         if (bitx(i))                                                          &
           cldlw(i,nc) = cldlw(i,nc)*(1.0e0-exp(-0.75e0*tauc(i)))
       enddo
       do kk = minktp,maxkbt
         do i = 1,ipts
           if (bitx(i)) cld0 (i,kk) = cldlw(i,nc)
         enddo
       enddo
     else                         ! prognostic liquit water
!
! --- calc tauc and emis
!
       do kk = minktp,maxkbt
         do i = 1,ipts
           if (kk .ge. ktop(i,nc) .and.                                        &
               kk .le. kbtm(i,nc) .and. bitx(i)) then
!
             awlw    = 0.090361e0 * cwp(i,kk)
!               ailw    = (0.005e0 + 1.0e0/rei(i,kk)) * cip(i,kk)
             ailw    = (0.0029e0 + 1.0e0/rei(i,kk)) * cip(i,kk)
             tauc(i) = tauc(i) + awlw + ailw
             taulw(i,kk) = awlw + ailw
           endif
         enddo
       enddo
       do i = 1,ipts
         if (bitx(i))                                                          &
           cldlw(i,nc) = cldlw(i,nc) * min(1.0e0,                              &
                         1.0e0 - exp(-1.66e0*tauc(i)) )                        
       enddo
       do kk = minktp,maxkbt
         do i = 1,ipts
           if (bitx(i)) cld0 (i,kk) = cldlw(i,nc)
         enddo
       enddo
     endif
!
   enddo               ! end of nncld loop!
!     if (icwp .eq. 1) then
   do k = 1,l
     do i = 1,ipts
       if (cldsw(i,k) .lt. 0.001e0) then
         taucl(i,k) = 0.0e0
         cldsw(i,k) = 0.0e0
         cwp  (i,k) = 0.0e0
         cip  (i,k) = 0.0e0
       endif
     enddo
   enddo
!     end if
   if (ipts .lt. imax) then     ! --- fill up arrays
     ipts1 = ipts + 1
     do i = ipts1,imax
       nclds(i) = nclds(ipts)
     enddo
     do k = 1,lp1
       do i = ipts1,imax
         cldlw(i,k) = cldlw(ipts,k)
         ktop(i,k) = ktop(ipts,k)
         kbtm(i,k) = kbtm(ipts,k)
         cfac(i,k) = cfac(ipts,k)
       enddo
     enddo
     do k = 1,l
       do i = ipts1,imax
         taulw(i,k) = taulw(ipts,k)
         cld0 (i,k) = cld0 (ipts,k)
         taucl(i,k) = taucl(ipts,k)
         cldsw(i,k) = cldsw(ipts,k)
         cwp  (i,k) = cwp  (ipts,k)
         cip  (i,k) = cip  (ipts,k)
         rew  (i,k) = rew  (ipts,k)
         rei  (i,k) = rei  (ipts,k)
         fice (i,k) = fice (ipts,k)
       enddo
     enddo
   endif
!
#ifdef DBG
#ifndef MP
#ifdef RMP
   latp = jgrd12_
   ir = igrd12_/4   ! center point
#else
   latp = latg2_
   ir = lonf2_/4   ! (180,0) point
#endif
   i   = ir - ibeg +1
   if(i.ge.1.and.i.le.ipts.and.lat.eq.latp) then
     write(6,501) ir,lat
     write(6,502)
501  format(//,'=== cloud properties in rad_cld_property_ice at (i,lat) ',     &
       2i10,'   ===',/)
502  format('    k','      clwp','      fice','       cwp',                    &
   '       cip','       rew', '       rei','      ityp')
503  format(i5,6f10.5,i10)
     do k = 1,l
       write(6,503) k,clwp(i,k),fice(i,k),cwp(i,k),                            &
            cip(i,k),rew(i,k),rei(i,k),ityp(i,k)
     enddo
!
     write(6,505)
505  format(//,'    k','   tem (C)','     cldlw',' ktop',' kbtm',              &
   '     cldsw','      cfac','    cldtot','    cldcnv')
506  format(i5,2f10.5,2i5,4f10.5)
     do k = 1,l
       write(6,506) k,t(i,k)-273.16,cldlw(i,k),ktop(i,k),kbtm(i,k),            &
        cldsw(i,k),cfac(i,k),cldtot(ir,l-k+1),cldcnv(ir,l-k+1)
     enddo
   endif
#endif
#endif
!
   return
   end subroutine rad_cld_property_nasa
!
!-------------------------------------------------------------------------------
