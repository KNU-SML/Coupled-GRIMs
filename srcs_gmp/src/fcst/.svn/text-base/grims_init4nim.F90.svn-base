#include "define.h"
!-------------------------------------------------------------------------------
   subroutine grims_init4nim(igms,igme,allowed_to_read)
!-------------------------------------------------------------------------------
!
! ::: structure :::
!
!   - grims_init      : initialize allocatable
!   - grims_variable  : treat variable from reading in grims_setup
!
! ::: history :::
!
!   - 2011-12-01  jung-eun kim
!
!-------------------------------------------------------------------------------
#ifdef NIM  /* NIM */
   use varsfc,    only   : lalbd_, numsfcs
   use paramodel, only   : para_init, nwater_, LONF2S, LATG2S
   use constant,  only   : rhoair0_, rhoh2o_, rhosnow_, cliq_, cvap_
   use comfspec_vr, only : comfspec_vr_init
   use comgpln,   only   : comgpln_init
   use comfcst,   only   : comfcst_init
   use co2dta,    only   : co2dta_init
   use comgda,    only   : comgda_init
   use comznl,    only   : comznl_init
   use comreduce, only   : comreduce_init
#ifndef SWRMDC
#ifndef ICECLOUD
   use rdparm,    only   : rdparm_init
   use comswaer,  only   : comswaer_init
#else
   use rdparm8,   only   : rdparm8_init
   use comswaer8, only   : comswaer8_init
#endif
#else
   use rdparm99,   only  : rdparm99_init
   use comswaer99, only  : comswaer99_init
#endif /* not SWRMDC */
   use lsm_init_module, only : dfkt_init, lsm_init_df, lsm_init_kt
   use comsfc            ! init_comsfc, sfcfcs
   use comfgrid
   use comfphys,   only  : comfphys_init, clstp
   use comfver
   use comio             ! iope
   use comgrad
   use radiag
#ifdef WSM3
   use module_mp_wsm3
#endif
#ifdef WSM5
   use module_mp_wsm5
#endif
#ifdef WSM6
   use module_mp_wsm6
#endif
#ifdef WDM5
   use module_mp_wdm5
#endif
#ifdef WDM6
   use module_mp_wdm6
#endif
!-------------------------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------------------------
   real*8               ::  swhr_int,lwhr_int
   integer, intent(in)  ::  igms,igme
   logical, intent(in)  ::  allowed_to_read
!-------------------------------------------------------------------------------
!
   if (allowed_to_read) then
     call para_init(igme)
     call rad_constant_init
     call comfspec_vr_init
     call comfgrid_init
     call comfphys_init
     call comfver_init
     call comgpln_init
     call comio_init
     call radiag_init
     call comfcst_init
     call co2dta_init
     call comgda_init
     call comznl_init
     call comgrad_init
     call comreduce_init
#ifndef SWRMDC
#ifndef ICECLOUD
     call rdparm_init
     call comswaer_init
#else
     call rdparm8_init
     call comswaer8_init
#endif
#else
     call rdparm99_init
     call comswaer99_init
#endif /* not SWRMDC */
     call init_comsfc(igme,igms)
     call funct_svp_init
     call funct_dew_point_temp_init
     call funct_pot_temp_init
     call funct_moist_adiabat_init
     call dfkt_init
     call lsm_init_df
     call lsm_init_kt
     call dyn_flux_zero_out(0)
!
#ifdef WSM3
     if (nwater_.eq.3) call wsm3init(rhoair0_, rhoh2o_, rhosnow_, cliq_, cvap_,&
                                    allowed_to_read)
#endif
#ifdef WSM5
     if (water_.eq.5) call wsm5init(rhoair0_, rhoh2o_, rhosnow_, cliq_, cvap_, &
                                    allowed_to_read)
#endif
#ifdef WSM6
     if (water_.eq.6) call wsm6init(rhoair0_, rhoh2o_, rhosnow_, cliq_, cvap_, &
                                    allowed_to_read)
#endif
#ifdef WDM5
     if (nwater_.eq.8) call wdm5init(rhoair0_, rhoh2o_, rhosnow_, cliq_, cvap_,&
                           allowed_to_read)
#endif
#ifdef WDM6
     if (nwater_.eq.9) call wdm6init(rhoair0_, rhoh2o_, rhosnow_, cliq_, cvap_,&
                           allowed_to_read)
#endif
   endif
!
   cowave=0.    ! osu index
   dtwave=0.    ! not used
#ifdef NIM_DBG
   if( iope ) then
     write(6,*) ' reading sfcfile (fixio) completed.'
   endif
#endif
!
! cmean,clstp control time averaging of convective clds in kuo
!
   clstp=99.
!
   return
#endif /* NIM */
   end subroutine grims_init4nim
#ifdef NIM  /* NIM */
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   subroutine grims_variable4nim(rhclnim,iernim,o3outnim,pstrnim,jerrnim)
!-------------------------------------------------------------------------------
   use comgrad
   use comreado3,  only :  o3out, pstr, jerr, jmout, loz
#ifdef NIM_DBGM
   use comio, only      :  iope,ncpu 
#endif
!-------------------------------------------------------------------------------
   real*8               :: rhclnim(nbin,nlon,nlat,mcld,nseal)
   real*8               :: o3outnim(jmout,loz),pstrnim(loz)
   integer              :: iernim,jerrnim
!-------------------------------------------------------------------------------
!
! for rad_cloud_read
!
!  rhcl(:,:,:,:,:)=rhclnim(:,:,:,:,:)
   do nsl = 1,nseal
      do ken = 1,mcld
         do l = 1,nlat
            do jl = 1,nlon
               do i = 1,nbin
                  rhcl(i,jl,l,ken,nsl) = rhclnim(i,jl,l,ken,nsl)
               enddo
            enddo
         enddo
      enddo
   enddo
   ier=iernim
#ifdef NIM_DBGM
   if(ncpu.eq.1) then  ! for 1cpu
     print*, 'grims_variable rhcl-1',rhcl(1,1,1,1,1)
     print*, 'grims_variable rhcl-1',rhclnim(1,1,1,1,1)
     print*, 'grims_variable rhcl-2',rhcl(2,2,2,2,2)
     print*, 'grims_variable rhcl-2',rhclnim(2,2,2,2,2)
     print*, 'grims_variable ier',ier
   else
     if(iope) then     ! for 2cpu
       print*, 'T,grims_variable rhcl-1',rhcl(1,1,1,1,1)
       print*, 'T,grims_variable rhcl-1',rhclnim(1,1,1,1,1)
       print*, 'T,grims_variable rhcl-2',rhcl(2,2,2,2,2)
       print*, 'T,grims_variable rhcl-2',rhclnim(2,2,2,2,2)
       print*, 'T,grims_variable ier',ier
     else
       print*, 'F,grims_variable rhcl-1',rhcl(1,1,1,1,1)
       print*, 'F,grims_variable rhcl-1',rhclnim(1,1,1,1,1)
       print*, 'F,grims_variable rhcl-2',rhcl(2,2,2,2,2)
       print*, 'F,grims_variable rhcl-2',rhclnim(2,2,2,2,2)
       print*, 'F,grims_variable ier',ier
     endif
   endif
#endif
!
! for rad_ozone_read
!
   o3out(:,:)=o3outnim(:,:)
   pstr(:)=pstrnim(:)
   jerr=jerrnim
!
#ifdef NIM_DBGM
   if(ncpu.eq.1) then  ! for 1cpu
     print*, 'grims_variable o3out-1',o3out(jmout,loz)
     print*, 'grims_variable o3out-1',o3out(jmout,loz)
     print*, 'grims_variable o3out-2',o3out(5,5)
     print*, 'grims_variable o3out-2',o3out(5,5)
     print*, 'grims_variable jerr',jerr
   else
     if(iope) then     ! for 2cpu
       print*, 'T,grims_variable o3out-1',o3out(jmout,loz)
       print*, 'T,grims_variable o3out-1',o3out(jmout,loz)
       print*, 'T,grims_variable o3out-2',o3out(5,5)
       print*, 'T,grims_variable o3out-2',o3out(5,5)
       print*, 'T,grims_variable jerr',jerr
     else
       print*, 'F,grims_variable o3out-1',o3out(jmout,loz)
       print*, 'F,grims_variable o3out-1',o3out(jmout,loz)
       print*, 'F,grims_variable o3out-2',o3out(5,5)
       print*, 'F,grims_variable o3out-2',o3out(5,5)
       print*, 'F,grims_variable jerr',jerr
     endif
   endif
#endif
   return
   end subroutine grims_variable4nim
#endif /* NIM */
