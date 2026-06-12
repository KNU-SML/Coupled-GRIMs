#include <define.h>
   subroutine sfc_usgs_rough(vegtyp,ijmax,z0)
!-------------------------------------------------------------------------------
!
! subprogram: sfc_usgs_rough
!
! abstract: compute roughness from vegetation type (13 types)
!
! usage:    call sfc_usgs_rough(vegtyp,ijmax,z0)
!   input argument list:
!     vegtyp   - real array of ijmax.  vegetation type.
!     ijmax    - dimension of vegtype
!   output argument list:
!     rough    - real array of ijmax.  surface roughness.
!
!-------------------------------------------------------------------------------
   real, dimension(ijmax)      ::  vegtyp, z0
   real, dimension(13)         ::  ri
!
   save
!
!  the clim values
!
!     data ri/0.05,0.5,0.2,0.5,0.5,0.1,0.5,2.7,1.1,0.8,                        &
!             0.85,0.6,0.001/
!
!  mm5 values
!     data RI/10.,0.05,0.12,0.13,0.2,0.1,0.1,0.5,0.5,0.4,                      &
!             0.5,0.5,0.005/
!
!   Fein Chen's values
!
   data RI/0.01,0.075,0.238,0.04,0.05,0.11,0.11,2.654,1.089                    &
          ,0.81,0.832,0.826,0.11/
!
   do ij = 1,ijmax
     iv=vegtyp(ij)+1
     if ( iv .gt. 13 ) then
       print *,'vegetation type greater than 13'
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
     z0(ij)=ri(iv)*100.
   enddo
!
   return
   end subroutine sfc_usgs_rough
!-------------------------------------------------------------------------------
