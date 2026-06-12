#include <define.h>
   subroutine dyn_transpose_coeff(a,kmax)
!------------------------------------------------------------------------------
!
! subprogram:    dyn_transpose_coeff  transposes spectral coefficients.
!
! abstract: transposes rhomboidal spectral coefficients.
!
! program history log:
!   1988-04-20  joseph sela
!   2000-01-01  song-you hong          usgs topo, kagwd options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call dyn_transpose_coeff (a, n)
!   input argument list:
!     a        - rhomboidal fields of spectral coefficients.
!     n        - number of  fields in array a.
!
!   output argument list:
!     a        - transposed fields of spectral coefficients.
!
!------------------------------------------------------------------------------
   use comchgr, only : mdim, indxnn, indxmm, b
!
   save
!
   real                 ::  a(mdim,kmax)
!------------------------------------------------------------------------------
!
   do k = 1,kmax
     do i = 1,mdim
       b(indxnn(i)) = a(i,k)
     enddo
     do i = 1,mdim
       a(i,k)=b(i)
     enddo
   enddo
!
   return
   end
!------------------------------------------------------------------------------
