#include <define.h>   
   module diag_3d_module
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [diag_3d_module]
!      |
!      |--- [diag_3d_index] *
!      |--- [diag_3d_zero_out] *
!      |--- [diag_3d_get] *
!      |--- [diag_3d_arrange] *
!      |--- [diag_3d_archive] *
!      |--- [diag_3d_digital_filter] *
!
!-------------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine diag_3d_index
!-------------------------------------------------------------------------------
   use paramodel
   use comgda
!-------------------------------------------------------------------------------
!
!  set levels and indices for diagnostics
!
!  use kd as defined below when using subroutine diag_3d_archive
!
!  kd     cnmgda
!   1     dtlarg
!   2     dtconv
!   3     dqconv
!   4     dtshal
!   5     dqshal
!   6     dtvrdf
!   7     duvrdf
!   8     dvvrdf
!   9     dqvrdf
!   10    dthsw
!   11    dthlw
!   12    cloud
!   13    cvcld
!
!-------------------------------------------------------------------------------
!
!  set identification index
!
#ifdef DG3
!
!  large-scale heating
!
   cnmgda( 1)='dtlarg  '
   ipugda( 1)=241
   ibmgda( 1)=0
#endif
#ifdef DG3
!
!  deep convective heating
!
   cnmgda( 2)='dtconv  '
   ipugda( 2)=242
   ibmgda( 2)=0
#endif
#ifdef DG3
!
!  deep convective moistening
!
   cnmgda( 3)='dqconv  '
   ipugda( 3)=243
   ibmgda( 3)=0
#endif
#ifdef DG3
!
!  shallow convective heating
!
   cnmgda( 4)='dtshal  '
   ipugda( 4)=244
   ibmgda( 4)=0
#endif
#ifdef DG3
!
!  shallow convective moistening
!
   cnmgda( 5)='dqshal  '
   ipugda( 5)=245
   ibmgda( 5)=0
#endif
#ifdef DG3
!
!  vertical diffusion of temperature
!
   cnmgda( 6)='dtvrdf  '
   ipugda( 6)=246
   ibmgda( 6)=0
#endif
#ifdef DG3
!
!  vertical diffusion of zonal wind
!
   cnmgda( 7)='duvrdf  '
   ipugda( 7)=247
   ibmgda( 7)=0
#endif
#ifdef DG3
!
!  vertical diffusion of meridional wind
!
   cnmgda( 8)='dvvrdf  '
   ipugda( 8)=248
   ibmgda( 8)=0
#endif
#ifdef DG3
!
!  vertical diffusion of moisture
!
   cnmgda( 9)='dqvrdf  '
   ipugda( 9)=249
   ibmgda( 9)=0
#endif
#ifdef DG3
!
!  short wave radiation heating
!
   cnmgda(10)='dthsw   '
   ipugda(10)=250
   ibmgda(10)=0
#endif
#ifdef DG3
!
!  long wave radiation heating
!
   cnmgda(11)='dthlw   '
   ipugda(11)=251
   ibmgda(11)=0
#endif
#ifdef DG3
!
!  cloud amount
!
   cnmgda(12)='cloud   '
   ipugda(12)=71
   ibmgda(12)=0
#endif
#ifdef DG3
!
!  convective cloud amount
!
   cnmgda(13)='cvcld   '
   ipugda(13)=72
   ibmgda(13)=0
#endif
!
   return
   end subroutine diag_3d_index
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine diag_3d_zero_out
!-------------------------------------------------------------------------------
   use comgda
!-------------------------------------------------------------------------------
#ifdef DG3
   real                 ::  gda(nwgda*kdgda)
!
   do n = 1,nwgda*kdgda
     gda(n)=0.
   enddo
   do jr = 1,nrgda
     call diag_3d_arrange(jr,nwgda*kdgda,gda)
   enddo
!
#endif
   return
   end subroutine diag_3d_zero_out
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine diag_3d_get(j,l,a)
!-------------------------------------------------------------------------------
   use comgda
!-------------------------------------------------------------------------------
!
!  get diagnostics from disk or memory
!  call another diag_3d_get or diag_3d_arrange to finish i/o
!
!  j is record number
!  l is record length
!  a is array of length l to get
!
!-------------------------------------------------------------------------------
#ifdef DG3
   real                 ::  a(l)
!
   do n = j*l-l+1,j*l
     a(n-(j*l-l+1)+1)=gdd(n)
   enddo
!
#endif
   return
   end subroutine diag_3d_get
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine diag_3d_arrange(j,l,a)
!-------------------------------------------------------------------------------
   use comgda
!-------------------------------------------------------------------------------
!
!  put diagnostics onto disk or memory
!  call another diag_3d_get or diag_3d_arrange to finish i/o
!
!  j is record number
!  l is record length
!  a is array of length l to get
!
!-------------------------------------------------------------------------------
!
#ifdef DG3
   real                 ::  a(l)
   do n=j*l-l+1,j*l
     gdd(n)=a(n-(j*l-l+1)+1)
   enddo
!
#endif
   return
   end subroutine diag_3d_arrange
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine diag_3d_archive(lons2,lonf2,a,dt,kd,gda)
!-------------------------------------------------------------------------------
   use paramodel, only : levs_,LONF2S
   use comgda
!-------------------------------------------------------------------------------
!
!  accumulate diagnostics
!
!  nlonx is 1st dimension of a
!  a is array to accumulate
!  dt is factor by which to multiply before accumulating
!  kd is diagnostic number
!  gda is accumulation array
!
!-------------------------------------------------------------------------------
#ifdef DG3
   real                 ::  a(lonf2,levs_)
   real                 ::  gda(nwgda,kdgda)
!
   if(kd.gt.0.and.kd.le.kdgda) then
     do k = 1,levs_
       do i = 1,lons2
         ik=i+LONF2S*(k-1)
         gda(ik,kd)=gda(ik,kd)+a(i,k)*dt
       enddo
     enddo
   endif
!
#endif
   return
   end subroutine diag_3d_archive
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine diag_3d_digital_filter(fac)
!-------------------------------------------------------------------------------
#ifdef DG3
   use paramodel 
   use comgda 
!-------------------------------------------------------------------------------
   real  ::  gda(nwgda*kdgda)
!
   do jr = 1,nrgda
     call diag_3d_get(jr,nwgda*kdgda,gda)
     do n = 1,nwgda*kdgda
       gda(n)=fac*gda(n)
     enddo
     call diag_3d_arrange(jr,nwgda*kdgda,gda)
   enddo
!
#endif
   return
   end subroutine diag_3d_digital_filter
!-------------------------------------------------------------------------------
   end module diag_3d_module
