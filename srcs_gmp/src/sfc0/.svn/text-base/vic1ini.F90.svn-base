#include <define.h>
   subroutine vic1ini(ijdim,im)
!-----------------------------------------------------------------------
!  
!  subprogram: vic1ini
!
!  ABSTRACT: initialize the rest variables and parameters after merging
!            grib data, analysis data, and other land surface initial data
!
!   2003-02-06  ji chen                development
!   2008-03-09  kyeong hee seol        cvs verion setup
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!  The following parameters and variables should be initialized after 
!  calling sfc
!    
!    jtsf,jsmc,jsno,jstc,jtg3,jzor,jcv ,jcvb,jcvt,jalb,
!    jsli,jveg,jcpy,jf10,jvet,jrot,jalf,just,jffm,jffh,
!    jds ,jdsm,jws ,jcef,jexp,jkst,jbub,jqrt,jbkd,jsld,
!    jwcr,jwpw,jsmr,
!
!  We need initialize the following parameters and variables here:
!
!           Number        Comments
!  1. jsmx:  lsoil        maximum soil moisture (mm)
!  2. jdpn:  nsoil        thickness of soil node (m)
!  3. jsxn:  nsoil        maximum sm at soil node (m3/m3)
!  4. jepn:  nsoil        Para the vari of Ksat at soil node (N/A)
!  5. jbbn:  nsoil        bubbling pressure at soil node (cm)
!  6. japn:  nsoil        para alpha at soil node
!  7. jbtn:  nsoil        para beta at soil node
!  8. jgmn:  nsoil        para gamma at soil node
!  9. jlai:      1        leaf area index (monthly)
! 10. jsic:  lsoil        soil ice content (mm)
! 11. jcsn:      1        canopy snow  (m h2o)
! 12. jrsn:      1        snow density (kg/m^3)
! 13. jtsn:      1        snow surface temperature (K)
! 14. jtpk:      1        snow pack temperature (K)
! 15. jsfw:      1        surface snow water equivalent (m h2o)
! 16. jpkw:      1        snow pack snow water equivalent (m h2o)
! 17. jlst:      1        time step since last snow fall
!
!-----------------------------------------------------------------------
   use constant, only : t0c_
   use varsfc
   use comsfc         ! sfcfcs
#ifdef VICLSM1
   use vic_veglib
#endif
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
   integer              ::  ijdim,im
#ifdef VICLSM1
   integer, parameter   ::  lsoil=lsoil_, nsoil=nsoil_
   real                 ::  vdep(lsoil),vdzn(nsoil),sumd,sumdp,totdep
   real                 ::  smxlo(lsoil),bublo(lsoil),exptlo(lsoil)
   real                 ::  smmaxnd(nsoil),expnd(nsoil),bubnd(nsoil)
   real                 ::  alpha(nsoil),beta(nsoil),gamma(nsoil)
   real                 ::  tmp_T, tmp_smc, smmax, buble
   real                 ::  exple, sfc_max_unfwat
   integer              ::  ij, m, n
!-----------------------------------------------------------------------
!
!  in VICLSM1, some parameters are constant over the study regions,   
!  so we just give the fixed values of them here. In the next version 
!  VICLSM, we will read their grib files                              
!
   do ij = 1,ijdim
     if(sfcfcs(ij,jsli).eq.1.) then
       sfcfcs(ij,jlai) = veg_lai(im,sfcfcs(ij,jvet))
       sfcfcs(ij,jcsn) = 0.0
       sfcfcs(ij,jrsn) = 50.0
       if(sfcfcs(ij,jsno).gt.0.) then
         sfcfcs(ij,jlst) = 0.0
         sfcfcs(ij,jstc) = sfcfcs(ij,jstc+1)  ! 6/2004
       else
         sfcfcs(ij,jlst) = -1.0
       endif
       sfcfcs(ij,jtsn) = min(sfcfcs(ij,jtsf),t0c_)
       sfcfcs(ij,jtpk) = min(sfcfcs(ij,jtsf),t0c_)
       sfcfcs(ij,jsfw) = 0.0
       sfcfcs(ij,jpkw) = 0.0
!         
! the Wcr and Wpwp are the fraction of maximum moisture 05/25/2004
!
       do m = 1, lsoil
         if(sfcfcs(ij,jsld+m-1).gt.0) then
               sfcfcs(ij,jsmx+m-1) = (sfcfcs(ij,jsld+m-1)-                     &
                             sfcfcs(ij,jbkd+m-1))/sfcfcs(ij,jsld+m-1) ! fraction
         else
               sfcfcs(ij,jsmx+m-1) = 0.0
         end if
         if(sfcfcs(ij,jsmc+m-1).gt.sfcfcs(ij,jsmx+m-1)) then
            sfcfcs(ij,jsmc+m-1) = sfcfcs(ij,jsmx+m-1)
         endif
         sfcfcs(ij,jwcr+m-1)=sfcfcs(ij,jwcr+m-1)*                              &
                             sfcfcs(ij,jsmx+m-1)
         sfcfcs(ij,jwpw+m-1)=sfcfcs(ij,jwpw+m-1)*                              &
                             sfcfcs(ij,jsmx+m-1)
         sfcfcs(ij,jsmx+m-1)=sfcfcs(ij,jsmx+m-1)*                              &
                             sfcfcs(ij,jdph+m-1)*1000.0 ! maximum sm (mm)
       enddo
!
! check soil layer thickness effects on the VIC performancs 5/22/2004
!
!      do m = 1,lsoil_
!        if(sfcfcs(ij,jdph+m-1).gt.0)then
!          sfcfcs(ij,jwcr+m-1) = sfcfcs(ij,jwcr+m-1) /                         &
!                                sfcfcs(ij,jdph+m-1)
!          sfcfcs(ij,jwpw+m-1) = sfcfcs(ij,jwpw+m-1) /                         &
!                                sfcfcs(ij,jdph+m-1)
!          sfcfcs(ij,jsmr+m-1) = sfcfcs(ij,jsmr+m-1) /                         &
!                                sfcfcs(ij,jdph+m-1)
!          sfcfcs(ij,jsmx+m-1) = sfcfcs(ij,jsmx+m-1) /                         &
!                                sfcfcs(ij,jdph+m-1)
!        endif
!      enddo
!      sfcfcs(ij,jdph  ) = 0.1
!      sfcfcs(ij,jdph+1) = 0.3
!      sfcfcs(ij,jdph+2) = 1.0
!      do m = 1,lsoil_
!        sfcfcs(ij,jwcr+m-1) = sfcfcs(ij,jwcr+m-1) *                           &
!                              sfcfcs(ij,jdph+m-1)
!        sfcfcs(ij,jwpw+m-1) = sfcfcs(ij,jwpw+m-1) *                           &
!                              sfcfcs(ij,jdph+m-1)
!        sfcfcs(ij,jsmr+m-1) = sfcfcs(ij,jsmr+m-1) *                           &
!                              sfcfcs(ij,jdph+m-1)
!        sfcfcs(ij,jsmx+m-1) = sfcfcs(ij,jsmx+m-1) *                           &
!                              sfcfcs(ij,jdph+m-1)
!      enddo
!
! end check
!
       do m = 1,lsoil
         tmp_T= sfcfcs(ij,jstc+m)
         smmax= sfcfcs(ij,jsmx+m-1)
         buble= sfcfcs(ij,jbub+m-1)
         exple = sfcfcs(ij,jexp+m-1)
         if(tmp_T.lt.t0c_)then
           sfcfcs(ij,jsic+m-1)=sfcfcs(ij,jsmc+m-1)*                            &
                               sfcfcs(ij,jdph+m-1)*1000.0-                     &
                               sfc_max_unfwat(tmp_T,smmax,buble,exple)
           if(sfcfcs(ij,jsic+m-1).lt.0.0)sfcfcs(ij,jsic+m-1)=0.0
         else
           sfcfcs(ij,jsic+m-1) = 0.0
         endif
         tmp_smc=sfcfcs(ij,jsmc+m-1)*sfcfcs(ij,jdph+m-1)*1000.0
         if(sfcfcs(ij,jsic+m-1).gt.tmp_smc) then
            sfcfcs(ij,jsic+m-1) = tmp_smc
         endif
       enddo
     endif
   enddo
!   
!  define the depth of soil nodes, the number of soil nodes >= 3
!  currently, we are using a very simple way to find the soil temperture
!  profile for each soil node (Ji 2003)
!
   print *,'in vic1ini for computing node parameters'
!
   do ij = 1,ijdim
     if(sfcfcs(ij,jsli).eq.1.) then
       do m = 1,lsoil
         vdep(m) = sfcfcs(ij,jdph+m-1)
       enddo
       if(vdep(1) .gt.0.0) then
         vdzn(1) = 0.0
         vdzn(2) = vdep(1)
         sumd = 0
         do n = 3,nsoil
           vdzn(n) = (n-1)*0.1
           sumd = sumd + vdzn(n)
         enddo
         sumdp = 0
         do m = 2,lsoil
           sumdp = sumdp + vdep(m)
         enddo
         do n = 3,nsoil
           if(sumd.gt.0) then
             vdzn(n) = vdzn(n)/sumd*sumdp
           else
             vdzn(n) = 0
           end if
         enddo
       else
         do n = 1, nsoil
           vdzn(n) = 0.0
         enddo
       endif
       do n = 1, nsoil
         sfcfcs(ij,jdpn+n-1) = vdzn(n)
       enddo
     endif
   enddo
!
! compute thermal node soil parameters
!
   do ij = 1,ijdim
     if(sfcfcs(ij,jsli).eq.1.) then
       totdep = 0
       do m = 1,lsoil
         vdep(m) =  sfcfcs(ij,jdph+m-1)
         totdep  = totdep + vdep(m)
       enddo
       do n = 1,nsoil
         vdzn(n) = sfcfcs(ij,jdpn+n-1)
       enddo
       if(totdep.gt.0.0) then
         do m = 1, lsoil
           smxlo(m)  = sfcfcs(ij,jsmx+m-1)
           bublo(m)  = sfcfcs(ij,jbub+m-1)
           exptlo(m) = sfcfcs(ij,jexp+m-1)
         enddo
         call nodepara(lsoil, nsoil, vdep, vdzn,                               &
                       smxlo, bublo, exptlo, smmaxnd,                          &
                       expnd, bubnd, alpha, beta,                              &
                       gamma)
         do n = 1,nsoil
           sfcfcs(ij,jsxn+n-1) = smmaxnd(n)
           sfcfcs(ij,jepn+n-1) = expnd(n)
           sfcfcs(ij,jbbn+n-1) = bubnd(n)
           sfcfcs(ij,japn+n-1) = alpha(n)
           sfcfcs(ij,jbtn+n-1) = beta(n)
           sfcfcs(ij,jgmn+n-1) = gamma(n)
         enddo
       else
         do n = 1,nsoil
           sfcfcs(ij,jsxn+n-1) = 0.0
           sfcfcs(ij,jepn+n-1) = 0.0
           sfcfcs(ij,jbbn+n-1) = 0.0
           sfcfcs(ij,japn+n-1) = 0.0
           sfcfcs(ij,jbtn+n-1) = 0.0
           sfcfcs(ij,jgmn+n-1) = 0.0
         enddo
       endif
     endif
   enddo
!
#endif
   return
   end
!-----------------------------------------------------------------------
