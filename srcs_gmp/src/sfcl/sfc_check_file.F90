#include <define.h>
   subroutine sfc_check_file
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [sfc_check_file]
!           |
!           |-- [sfc_check_snowdepth] *
!           |-- [sfc_check_seaice] *
!           |-- [sfc_check_vicland] *
!           |-- [file_check_index] *
!           |-- [file_check_name] *
!           |-- [sfc_check_maxmin] *
!
!-------------------------------------------------------------------------------
   end subroutine sfc_check_file
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_check_snowdepth(snow,snwdph,ijdim)
!-------------------------------------------------------------------------------
   real                 ::  snow(ijdim),snwdph(ijdim)
!
   do ij = 1,ijdim
     if(snow(ij).gt.0.) then
       if(snwdph(ij).le.0.) then
         snwdph(ij)=snow(ij)*5.
       endif
     endif
   enddo
!
   return
   end subroutine sfc_check_snowdepth
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_check_seaice(ais,glacir,amxice,slmask,idim,jdim)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  idim,jdim
   real                 ::  ais(idim,jdim),glacir(idim,jdim)
   real                 ::  amxice(idim,jdim),slmask(idim,jdim)
!
   integer              ::  i,j,ij,kount,ip,im,jp,jm
   real                 ::  perr
!
!  check sea-ice cover mask against land-sea mask
!
   kount=1
   do j = 1,jdim
     do i = 1,idim
       if(slmask(i,j).eq.0..and.glacir(i,j).eq.1..and.ais(i,j).ne.1.) then
         ais(i,j)=1.
         kount=kount+1
       endif
       if(slmask(i,j).eq.1..and.ais(i,j).eq.1.) then
         ais(i,j)=0.
         kount=kount+1
       endif
       if(slmask(i,j).eq.0..and.amxice(i,j).eq.0..and. ais(i,j).eq.1.) then
         ais(i,j)=0.
         kount=kount+1
       endif
     enddo
   enddo
#ifndef MP
#undef DANGER
#ifdef DANGER
!
!  remove isolated open ocean surrounded by sea ice and/or land
!
   ij=0
   do j = 1,jdim
     jp=j+1
     jm=j-1
     if(jp.gt.jdim) jp=jdim-1
     if(jm.lt.1) jm=2
     do i = 1,idim
       ip=i+1
       im=i-1
       if(ip.gt.idim) ip=idim-1
       if(im.lt.1) im=2
       if(slmask(i,j).eq.0..and.ais(i,j).eq.0.) then
         if((slmask(ip,jp).eq.1..or.ais(ip,jp).eq.1.).and.                     &
            (slmask(i ,jp).eq.1..or.ais(i ,jp).eq.1.).and.                     &
            (slmask(im,jp).eq.1..or.ais(im,jp).eq.1.).and.                     &
            (slmask(ip,j ).eq.1..or.ais(ip,j ).eq.1.).and.                     &
            (slmask(im,j ).eq.1..or.ais(im,j ).eq.1.).and.                     &
            (slmask(ip,jm).eq.1..or.ais(ip,jm).eq.1.).and.                     &
            (slmask(i ,jm).eq.1..or.ais(i ,jm).eq.1.).and.                     &
            (slmask(im,jm).eq.1..or.ais(im,jm).eq.1.)) then
           ij=ij+1
           ais(i,j)=1.
         endif
       endif
     enddo
   enddo
!
   if (ij.ne.0)  then
     write(6,'(I8,A)') ij,' sea isolated point to sea-ice'
   endif
#endif
#endif
!
!  print results
!
#ifdef DBG
   perr=float(kount)/float(idim*jdim)*100.
   write(6,100) perr
#endif
  100 format("sfc_check_seaice:% qc points=",f4.1)
!
   return
   end subroutine sfc_check_seaice
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_check_vicland(idim,jdim)
!-------------------------------------------------------------------------------
!
! abstract: check the land surface with VIC parameters according to
!   non zero soil density
!
! program history log:
!   2004-01     ji chen ECPC/CRD/SIO/UCSD
!
!-------------------------------------------------------------------------------
   use varsfc
   use comsfc           !   sfcfcs
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  idim, jdim
#ifdef VICLSM1         /* to return */
   integer              ::  ijdim
   integer              ::  nvicland(idim*jdim), landnp(idim,jdim)
   integer              ::  landvic(idim,jdim),landsli(idim,jdim)

   integer              ::  m, n, np, nf, iv
   integer              ::  i, j, ij, k, kk
   integer              ::  ii1,jj1,ii21,ii22,jj21,jj22
!
!--- end of variable definition ---
!
#ifdef DBG
   open(1,file='check_vic_Ws.dat',status='unknown')
!
#endif
   ijdim = idim*jdim
   m = 0                                  ! total number of land grids
   n = 0
   do ij = 1,ijdim
     if(sfcfcs(ij,jsli).eq.1) then       ! land
       m = m + 1
       if(sfcfcs(ij,jsld).le.0.0) then  ! soil density
         n = n + 1
         nvicland(n) = ij             
       endif
     endif
   enddo
!
#ifdef DBG
   close(1)
#endif
   if(n.gt.0) then
     do ij = 1,ijdim
       i = mod(ij,idim)
       if(i.eq.0) i = idim
       j = (ij-i)/idim + 1
       if(sfcfcs(ij,jsld).gt.0.0) then
         landvic(i,j) = 1
       else
         landvic(i,j) = 0
       endif
       if(sfcfcs(ij,jsli).eq.1.) then
         landsli(i,j) = 1
       else
         landsli(i,j) = 0
       endif
     enddo
!
!         write(*, *) 'start to fill out VIC grids with sli land'
!
     do i = 1,idim
       do j =1,jdim
         k = 0
         if(landsli(i,j).eq.1.and.landvic(i,j).eq.0)then
           ii1 = i
           jj1 = j
           kk = 1
!
55         ii21 = ii1 + kk
           if(ii21.ge.1.and.ii21.le.idim) then
             if(landvic(ii21,jj1).eq.1) then
               k = (jj1 - 1)*idim + ii21
             endif
           endif
           if(k.eq.0) then
             ii22 = ii1 - kk
             if(ii22.ge.1.and.ii22.le.idim) then
               if(landvic(ii22,jj1).eq.1) then
                 k = (jj1 - 1)*idim + ii22
               endif
             endif
           endif
           if(k.eq.0) then
             jj21 = jj1 + kk
             if(jj21.ge.1.and.jj21.le.jdim) then
               if(landvic(ii1,jj21).eq.1) then
                 k = (jj21 - 1)*idim + ii1
               endif
             endif
           endif
           if(k.eq.0) then
             jj22 = jj1 - kk
             if(jj22.ge.1.and.jj22.le.jdim) then
               if(landvic(ii1,jj22).eq.1) then
                 k = (jj22 - 1)*idim + ii1
               endif
             endif
           endif
           if(k.eq.0) then
             if(ii21.ge.1.and.ii21.le.idim.and.                           &
                jj21.ge.1.and.jj21.le.jdim) then
               if(landvic(ii21,jj21).eq.1) then
                 k = (jj21 - 1)*idim + ii21
               endif
             endif
           endif
           if(k.eq.0) then
             if(ii22.ge.1.and.ii22.le.idim.and.                           &
                jj21.ge.1.and.jj21.le.jdim) then
               if(landvic(ii22,jj21).eq.1) then
                 k = (jj21 - 1)*idim + ii22
               endif
             endif
           endif
           if(k.eq.0) then
             if(ii21.ge.1.and.ii21.le.idim.and.                           &
                jj22.ge.1.and.jj22.le.jdim) then
               if(landvic(ii21,jj22).eq.1) then
                 k = (jj22 - 1)*idim + ii21
               endif
             endif
           endif
           if(k.eq.0) then
             if(ii22.ge.1.and.ii22.le.idim.and.                           &
                jj22.ge.1.and.jj22.le.jdim) then
               if(landvic(ii22,jj22).eq.1) then
                 k = (jj22 - 1)*idim + ii22
               endif
             endif
           endif
           if(k.eq.0) then
             if((ii21.ge.1.and.ii21.le.idim).or.                          &
                (ii22.ge.1.and.ii22.le.idim).or.                          &
                (jj21.ge.1.and.jj21.le.jdim).or.                          &
                (jj22.ge.1.and.jj22.le.jdim)) then
               kk = kk + 1
               !
               goto 55
               !
             else
               write(*,*) 'no vic land found for sli land'
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
           endif
         endif
         landnp(i,j) = k
       enddo
     enddo
!
#ifdef DBG
     open(77,file='VIC_land.dat',status='unknown')
     open(88,file='SLI_land.dat',status='unknown')
     open(99,file='VIC_SLI_land.dat',status='unknown')
     do j = 1,jdim
       write(77,177) (landvic(i,j),i=1,idim)
       write(88,177) (landsli(i,j),i=1,idim)
       write(99,199) (landnp(i,j),i=1,idim)
177    format(721i1)
199    format(721i7)
     enddo
     close(77)
     close(88)
     close(99)
     open(66,file='VIC_noland.dat',status='unknown')
     write(66,*) 'the number of land grids without VIC para ',n
#endif
!
     do ij = 1,n
       i = mod(nvicland(ij),idim)
       if(i.eq.0) i = idim
       j = (nvicland(ij)-i)/idim + 1
#ifdef DBG
       write(66,166) nvicland(ij), landnp(i,j)
166    format(2i15)
#endif
       np = nvicland(ij)
       nf = landnp(i,j)
!
       sfcfcs(np,jveg)=sfcfcs(nf,jveg)
       sfcfcs(np,jcpy)=sfcfcs(nf,jcpy)
       sfcfcs(np,jvet)=sfcfcs(nf,jvet)
       sfcfcs(np,jbif)=sfcfcs(nf,jbif)
       sfcfcs(np,jds )=sfcfcs(nf,jds )
       sfcfcs(np,jdsm)=sfcfcs(nf,jdsm)
       sfcfcs(np,jws )=sfcfcs(nf,jws )
       sfcfcs(np,jcef)=sfcfcs(nf,jcef)
       sfcfcs(np,jlai)=sfcfcs(nf,jlai)
       sfcfcs(np,jslz)=sfcfcs(nf,jslz)
       sfcfcs(np,jsnz)=sfcfcs(nf,jsnz)
!
       do m = 1,lsoil_
         sfcfcs(np,jrot+m-1)=sfcfcs(nf,jrot+m-1)
         sfcfcs(np,jexp+m-1)=sfcfcs(nf,jexp+m-1)
         sfcfcs(np,jkst+m-1)=sfcfcs(nf,jkst+m-1)
         sfcfcs(np,jdph+m-1)=sfcfcs(nf,jdph+m-1)
         sfcfcs(np,jbub+m-1)=sfcfcs(nf,jbub+m-1)
         sfcfcs(np,jqrt+m-1)=sfcfcs(nf,jqrt+m-1)
         sfcfcs(np,jbkd+m-1)=sfcfcs(nf,jbkd+m-1)
         sfcfcs(np,jsld+m-1)=sfcfcs(nf,jsld+m-1)
         sfcfcs(np,jwcr+m-1)=sfcfcs(nf,jwcr+m-1)
         sfcfcs(np,jwpw+m-1)=sfcfcs(nf,jwpw+m-1)
         sfcfcs(np,jsmr+m-1)=sfcfcs(nf,jsmr+m-1)
         sfcfcs(np,jsmx+m-1)=sfcfcs(nf,jsmx+m-1)
       enddo
       do n = 1,nsoil_
         sfcfcs(np,jdpn+n-1)=sfcfcs(nf,jdpn+n-1)
         sfcfcs(np,jsxn+n-1)=sfcfcs(nf,jsxn+n-1)
         sfcfcs(np,jepn+n-1)=sfcfcs(nf,jepn+n-1)
         sfcfcs(np,jbbn+n-1)=sfcfcs(nf,jbbn+n-1)
         sfcfcs(np,japn+n-1)=sfcfcs(nf,japn+n-1)
         sfcfcs(np,jbtn+n-1)=sfcfcs(nf,jbtn+n-1)
         sfcfcs(np,jgmn+n-1)=sfcfcs(nf,jgmn+n-1)
       enddo
     enddo
#ifdef DBG
     close(66)
#endif
   endif
#endif
!
   return
   end subroutine sfc_check_vicland
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine file_check_index(is2g)
!-------------------------------------------------------------------------------
   use vargrb
   use varsfc
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  is2g(numsfcs),n
!
   do n = 1,numsfcs
     is2g(n)=0
   enddo
!
   is2g(jtsf)=itsf  ! sst and land surface temperature
   is2g(jsno)=isno  ! snow depth
   is2g(jtg3)=itg3  ! deep soil temperature
   is2g(jzor)=izor  ! surface roughness
   is2g(jcv )=9999  ! convective cloud cover
   is2g(jcvb)=9999  ! convective cloud base
   is2g(jcvt)=9999  ! convective cloud top
   is2g(jsli)=iais  ! sea ice mask (to be converted to slimsk)
   is2g(jcpy)=9999  ! canopy water content
   is2g(jf10)=9999  ! 10m and lowest sigma level coversion factor
#ifndef OSULSM1
   is2g(jalf)=ialf  ! albedo fraction
   is2g(just)=9999  ! fricational velocity
   is2g(jffm)=9999  ! coefficient for momentun exchange
   is2g(jffh)=9999  ! coefficient for heat exchange
#endif
!
#ifdef OSULSM1
   is2g(jsmc)=iso2  ! soil moisture
   is2g(jstc)=ito2  ! soil temperature
   is2g(jalb)=iab1  ! albedo
   is2g(jplr)=iplr  ! stomata resistance
#endif
!
#ifdef OSULSM2
   is2g(jsmc)=iso2  ! soil moisture
   is2g(jstc)=ito2  ! soil temperature
   is2g(jalb)=iab4  ! albedo
   is2g(jveg)=iveg  ! vegetation cover
   is2g(jvet)=ivet  ! vegetation type
   is2g(jsot)=isot  ! soil type
#endif
!
#ifdef NOALSM1
   is2g(jsmc)=isn1  ! soil moisture
   is2g(jstc)=itn1  ! soil temperature
   is2g(jalb)=iab4  ! albedo
   is2g(jveg)=iveg  ! vegetation cover
   is2g(jvet)=ivet  ! vegetation type
   is2g(jsot)=isot  ! soil type
   is2g(jprc)=9999  ! precip for noah
   is2g(jsrf)=9999  ! sr flag for noah
   is2g(jsnd)=9999  ! snow depth for noah
   is2g(jslc)=9999  ! slc for noah
   is2g(jsmn)=ismn  ! min vegetation cover for noah
   is2g(jsmx)=ismx  ! max vegetation cover for noah
   is2g(jslo)=islo  ! slope type for noah
   is2g(jsna)=isna  ! snow albedo for noah
#endif
!
#ifdef VICLSM1
   is2g(jsmc)=isv1  ! soil moisture
   is2g(jstc)=itv1  ! soil temperature
   is2g(jalb)=iab4  ! albedo
   is2g(jveg)=ivgv  ! vegetation cover
   is2g(jvet)=ivtv  ! vegetation type
   is2g(jrot)=ivrt  ! vegetation root
   is2g(jprc)=9999  ! precip for noah
   is2g(jsrf)=9999  ! sr flag for noah
   is2g(jbif)=ibif  ! Variable infil curve parameter (N/A)
   is2g(jds )=ids   ! Fract of Dsm nonlinear baseflow begins
   is2g(jdsm)=idsm  ! Maximum velocity of baseflow (mm/day)
   is2g(jws )=iws   ! Fract maxi sm nonlinear baseflow occurs
   is2g(jcef)=icef  ! c
   is2g(jexp)=iexp  ! Para the vari of Ksat with sm
   is2g(jkst)=ikst  ! Saturated hydrologic conductivity (mm/day)
   is2g(jdph)=idph  ! soil layer thickness (m)
   is2g(jbub)=ibub  ! Bubbling pressure of soil layer (cm)
   is2g(jqrt)=iqrt  ! Quartz content of soil layer (fraction)
   is2g(jbkd)=ibkd  ! Bulk density of soil layer (kg/m3)
   is2g(jsld)=isld  ! Soil density of soil layer (kg/m3)
   is2g(jwcr)=iwcr  ! sm content at the critical point (mm)
   is2g(jwpw)=iwpw  ! sm content wilting point (mm)
   is2g(jsmr)=ismr  ! Soil moisture residual moisture (mm)
   is2g(jsmx)=9999  ! maximum soil moisture (mm)
   is2g(jdpn)=9999  ! maximum sm at soil node (m3/m3)
   is2g(jsxn)=9999  ! maximum sm at soil node (m3/m3)
   is2g(jepn)=9999  ! Para the vari of Ksat at soil node (N/A)
   is2g(jbbn)=9999  ! bubbling pressure at soil node (cm)
   is2g(japn)=9999  ! para alpha at soil node
   is2g(jbtn)=9999  ! para beta at soil node
   is2g(jgmn)=9999  ! para gamma at soil node
   is2g(jlai)=9999  ! lai
   is2g(jslz)=islz  ! soil roughness (m)
   is2g(jsnz)=isnz  ! snow roughness (m)
   is2g(jsic)=9999  ! soil ice content (mm)
   is2g(jcsn)=9999  ! canopy snow  (mm h2o)
   is2g(jrsn)=9999  ! snow density (kg/m^3)
   is2g(jtsn)=9999  ! snow surface temperature (K)
   is2g(jtpk)=9999  ! snow pack temperature (K)
   is2g(jsfw)=9999  ! surface snow water equivalent (mm h2o)
   is2g(jpkw)=9999  ! snow pack snow water equivalent (mm h2o)
   is2g(jlst)=9999  ! time step since last snow fall
#endif
!
#ifdef VICLSM2
   is2g(jsmc)=isv2  ! soil moisture
   is2g(jstc)=itv2  ! soil temperature
   is2g(jalb)=iab4  ! albedo
   is2g(jveg)=ivgv  ! vegetation cover
   is2g(jvet)=ivtv  ! vegetation type
   is2g(jrot)=ivrt  ! vegetation root
   is2g(jprc)=9999  ! precip for noah
   is2g(jsrf)=9999  ! sr flag for noah
   is2g(jbif)=ibif  ! Variable infil curve parameter (N/A)
   is2g(jds )=ids   ! Fract of Dsm nonlinear baseflow begins
   is2g(jdsm)=idsm  ! Maximum velocity of baseflow (mm/day)
   is2g(jws )=iws   ! Fract maxi sm nonlinear baseflow occurs
   is2g(jcef)=icef  ! c
   is2g(jexp)=iexp  ! Para the vari of Ksat with sm
   is2g(jkst)=ikst  ! Saturated hydrologic conductivity (mm/day)
   is2g(jdph)=idph  ! soil layer thickness (m)
   is2g(jbub)=ibub  ! Bubbling pressure of soil layer (cm)
   is2g(jqrt)=iqrt  ! Quartz content of soil layer (fraction)
   is2g(jbkd)=ibkd  ! Bulk density of soil layer (kg/m3)
   is2g(jsld)=isld  ! Soil density of soil layer (kg/m3)
   is2g(jwcr)=iwcr  ! sm content at the critical point (mm)
   is2g(jwpw)=iwpw  ! sm content wilting point (mm)
   is2g(jsmr)=ismr  ! Soil moisture residual moisture (mm)
   is2g(jsmx)=9999  ! maximum soil moisture (mm)
   is2g(jdpn)=9999  ! maximum sm at soil node (m3/m3)
   is2g(jsxn)=9999  ! maximum sm at soil node (m3/m3)
   is2g(jepn)=9999  ! Para the vari of Ksat at soil node (N/A)
   is2g(jbbn)=9999  ! bubbling pressure at soil node (cm)
   is2g(japn)=9999  ! para alpha at soil node
   is2g(jbtn)=9999  ! para beta at soil node
   is2g(jgmn)=9999  ! para gamma at soil node
   is2g(jlai)=9999  ! lai
   is2g(jslz)=islz  ! soil roughness (m)
   is2g(jsnz)=isnz  ! snow roughness (m)
   is2g(jsic)=9999  ! soil ice content (mm)
   is2g(jcsn)=9999  ! canopy snow  (mm h2o)
   is2g(jrsn)=9999  ! snow density (kg/m^3)
   is2g(jtsn)=9999  ! snow surface temperature (K)
   is2g(jtpk)=9999  ! snow pack temperature (K)
   is2g(jsfw)=9999  ! surface snow water equivalent (mm h2o)
   is2g(jpkw)=9999  ! snow pack snow water equivalent (mm h2o)
   is2g(jlst)=9999  ! time step since last snow fall
   is2g(jlai)=ilai  ! leaf area index
#endif
   is2g(joml)=ioml  ! ocean mixed layer depth
   is2g(jaer)=iaer  ! aerosol distribution
   is2g(jkpr)=ikpr  ! aerosol distribution
   is2g(jden)=iden  ! aerosol distribution
   is2g(jdxc)=idxc  ! aerosol distribution
   is2g(jmix)=imix  ! aerosol distribution
!
   return
   end subroutine file_check_index
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine file_check_name(condir,bindir,fn,numgrbs,                        &
                      fnmskg,fnorog,fnmask)
#include <abort.h>
!-------------------------------------------------------------------------------
!
!  checks sfc file names for consistency.
!
!  if condir and bin dir are given, add them in front
!  of fn, fntfc, fnmskg,fnorog and fnmask to
!  provide full directory structure.  
!
!  note that if filenames starts from '/', assume that
!  full directory is already given and condir, bindir are not added.
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  numgrbs
   character(len=128)   ::  condir,bindir
   character(len=128)   ::  fn(numgrbs)
   character(len=128)   ::  fnmskg,fnorog,fnmask
   integer              ::  n,nbin,ncon
!
   n=1
!
   do while (condir(n:n).ne.' '.and.n.le.128)
     n=n+1
   enddo
!
   ncon=n-1
   n=1
   do while (bindir(n:n).ne.' '.and.n.le.128)
     n=n+1
   enddo
!
   nbin=n-1
!
   do n = 1,numgrbs
     if(ncon.gt.0) then
       if(fn(n)(1:4).ne.'    '.and.fn(n)(1:1).ne.'/') then
         fn(n)=condir(1:ncon)//'/'//fn(n)
       endif
     endif
   enddo
!
   if(fnmskg(1:4).eq.'    ') then
     print *,'fnmskg empty'
     call MPABORT
   else
     if(ncon.gt.0.and.fnmskg(1:1).ne.'/') then
       fnmskg=condir(1:ncon)//'/'//fnmskg
     endif
   endif
!
   if(fnorog(1:4).eq.'    ') then
     print *,'fnorog empty'
     call MPABORT
   endif
!
   if(fnmask(1:4).eq.'    ') then
     print *,'fnmask empty'
     call MPABORT
   endif
!
   return
   end subroutine file_check_name
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_check_maxmin(fld,ijdim,loismsk,vmaxmin,svar)
!-------------------------------------------------------------------------------
!
!  loismsk=0  ... ocean without sea ice
!          1  ... land without snow
!          2  ... sea ice without snow
!          3  ... land with snow
!          4  ... sea ice with snow
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ijdim
   real                 ::  fld(ijdim)
   integer              ::  loismsk(ijdim)
   character(len=8)     ::  svar
!
!  first index for max and min
!  second index for loismsk
!
   real                 ::  vmaxmin(2,5)
   real                 ::  perr
   integer              ::  kount,ij,n
!
!  if criteria max .lt. min, do not check
!
   if(vmaxmin(1,1).lt.vmaxmin(2,1)) return
!
!  check against land-sea mask and ice cover mask
!
   kount=0
   do ij = 1,ijdim
     do n = 1,5
       if(loismsk(ij).eq.n-1) then
         if(fld(ij).gt.vmaxmin(1,n)) then
           fld(ij)=vmaxmin(1,n)
           kount=kount+1
         elseif(fld(ij).lt.vmaxmin(2,n)) then
           fld(ij)=vmaxmin(2,n)
           kount=kount+1
         endif
       endif
     enddo
   enddo
! 
!  print results
!
#ifdef DBG
   perr=float(kount)/float(ijdim)*100.
   write(6,100) svar,perr
100 format(a10,":% qc points=",f4.1)
#endif
!
   return
   end subroutine sfc_check_maxmin
!-------------------------------------------------------------------------------
