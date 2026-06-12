#include <define.h>
   subroutine phys_flux_zero
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!    [phys_flux_zero]
!        |
!        |--- [dyn_tmax_zero_out] *
!        |--- [dyn_flux_zero_out] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!   2012-09-21  haiqin li              COUPLE_ROP option
!
!-------------------------------------------------------------------------------
   end subroutine phys_flux_zero
!
!-------------------------------------------------------------------------------
   subroutine dyn_tmax_zero_out
!-------------------------------------------------------------------------------
   use paramodel, only : LONF2S,LATG2S
   use comsfc
   use comfphys
!-------------------------------------------------------------------------------
   do j = 1,LATG2S
     do i = 1,LONF2S
       tmpmax(i,j) = 0.0
       tmpmin(i,j) = 1.0e10
     enddo
   enddo
!
   return                                                                    
   end subroutine dyn_tmax_zero_out  
!
!-------------------------------------------------------------------------------
   subroutine dyn_flux_zero_out(ifluxr)
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only : IGRD12S,JGRD12S
   use rscomf_rerun
#else
   use paramodel, only : IGRD12S=>LONF2S,JGRD12S=>LATG2S
   use comfphys
   use radiag
#endif
   use comsfc
#ifdef COUPLE_ROP
   integer             :: dhour,rtime,rtimsw,rtimlw
#endif
!-------------------------------------------------------------------------------
!
! regional spectral model developed by hann-ming henry juang
!       version with nested to global spectral model
!
!. this routine setup all the routine and constant for RMP.
!    and get input data from gmp or others include surface files.
!
!-------------------------------------------------------------------------------
#ifdef ORIGIN_THREAD
!$doacross share(dusfc,dvsfc,dtsfc,dqsfc,
!$&              raintot,raincps,ulwsfc,dlwsfc,
!$&              gflux,tmpmax,tmpmin,runoff,ep,
!$&              cldwrk,dugwd,dvgwd,psmean,cvavg,
!$&              snowmelt,snowevap,snowfall,qull,qvll,
!$&              alhtfl,evcnp,bgrun
!$&         local(i,lat)
#endif
!
#ifdef COUPLE_ROP
   dhour=dtpost/3600.
   if(dtpost.gt.0) then
     rtime=1./dtpost
   else
     rtime=0.
   endif
!
   secswr=max(real(dhour),dtswav) * 3600.
   seclwr=max(real(dhour),dtlwav) * 3600.
   if(secswr.gt.0.) then
     rtimsw=1./secswr
   else
     rtimsw=1.
   endif
   if(seclwr.gt.0.) then
     rtimlw=1./seclwr
   else
     rtimlw=1.
   endif
!
   if(ifluxr.eq.3) then
     do lat = 1,JGRD12S
       do i = 1,IGRD12S
         romsevap(i,lat)=romsevap(i,lat)+dqsfc(i,lat)
         romssens(i,lat)=romssens(i,lat)+dtsfc(i,lat)
         romsustr(i,lat)=romsustr(i,lat)+dusfc(i,lat)
         romsvstr(i,lat)=romsvstr(i,lat)+dvsfc(i,lat)
         romslwup(i,lat)=romslwup(i,lat)+ulwsfc(i,lat)
         romsswup(i,lat)=romsswup(i,lat)+fluxr(i,lat,3)
         romslwdn(i,lat)=romslwdn(i,lat)+dlwsfc(i,lat)
         romsswdn(i,lat)=romsswdn(i,lat)+fluxr(i,lat,4)
         romsprcp(i,lat)=romsprcp(i,lat)+raintot(i,lat)
       enddo
     enddo
     romsrtime =romsrtime + 1./rtime
     romsrtswup=romsrtswup+ 1./rtimsw
     romsrtswdn=romsrtswdn+ 1./rtimsw
   endif
#endif
   if(ifluxr.ne.2) then
#ifdef OPENMP
!$omp parallel do private(i,lat)
#endif
     do lat = 1,JGRD12S
       do i = 1,IGRD12S
         raintot(i,lat) = 0.0e0
         raincps(i,lat) = 0.0e0
         dusfc (i,lat) = 0.0e0
         dvsfc (i,lat) = 0.0e0
         dtsfc (i,lat) = 0.0e0
         dqsfc (i,lat) = 0.0e0
         ulwsfc(i,lat) = 0.0e0
         dlwsfc(i,lat) = 0.0e0
         gflux (i,lat) = 0.0e0
         tmpmax(i,lat) = 0.0e0
         tmpmin(i,lat) = 1.0e10
         runoff(i,lat) = 0.0e0
         ep    (i,lat) = 0.0e0
         cldwrk(i,lat) = 0.0e0
         dugwd (i,lat) = 0.0e0
         dvgwd (i,lat) = 0.0e0
#ifdef GWDC
         dugwdc(i,lat) = 0.0e0
         dvgwdc(i,lat) = 0.0e0
#endif
         cvavg (i,lat) = 0.0e0
         psmean(i,lat) = 0.0e0
         snowmelt(i,lat) = 0.0e0
         snowevap(i,lat) = 0.0e0
         snowfall(i,lat) = 0.0e0
         qull(i,lat) = 0.0e0
         qvll(i,lat) = 0.0e0
         alhtfl(i,lat) = 0.0e0
         evcnp(i,lat) = 0.0e0
         bgrun(i,lat) = 0.0e0
#ifdef NOAHYDRO
         raintot2(i,lat) = 0.0e0
#endif
#ifdef SMP_NUDGING_RAD
         scmsdn (i,lat) = 0.0e0
         scmsup (i,lat) = 0.0e0
#endif
#ifdef VICLSM1
         gheat (i,lat) = 0.0e0
#endif
       enddo
     enddo
     dtflux=0.0e0
   endif
!
   if((ifluxr.eq.0).or.(ifluxr.eq.3)) then
     do iv = 1,26
       do lat = 1,JGRD12S
         do i = 1,IGRD12S
           fluxr(i,lat,iv) = 0.e0
         enddo
       enddo
     enddo
     do lat = 1,JGRD12S
       do i = 1,IGRD12S
         dlwsfc(i,lat)=0.
         ulwsfc(i,lat)=0.
         cvavg (i,lat)=0.
       enddo
     enddo
   endif
!
#ifdef COUPLE_ROP
   if(ifluxr.ne.3) then
     do lat = 1,JGRD12S
       do i = 1,IGRD12S
         romsevap(i,lat)=0.
         romssens(i,lat)=0.
         romsustr(i,lat)=0.
         romsvstr(i,lat)=0.
         romslwup(i,lat)=0.
         romsswup(i,lat)=0.
         romslwdn(i,lat)=0.
         romsswdn(i,lat)=0.
         romsprcp(i,lat)=0.
       enddo
     enddo
     romsrtime=0.
     romsrtswup=0.
     romsrtswdn=0.
   endif
#endif
!
   return
   end subroutine dyn_flux_zero_out
!
!-------------------------------------------------------------------------------
