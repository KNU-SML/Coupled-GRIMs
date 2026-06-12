#include <define.h>
   subroutine chgr_fourier_diff(a,b,n)
!-------------------------------------------------------------------------------
!
! subprogram: chgr_fourier_diff 
!
! abstract: forms sum and difference of 2 fouriers.
!
! program history log:
!   1988-04-08  joseph sela
!   2000-01-01  song-you hong          cvs version
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call symns (a, b, n)
!   input argument list:
!     a        - let xxx be array elements a(1,k) ... a(lonf_,k).
!                let yyy be array elements
!                a(lonf_+1,k) ... a(lonf2_,k).
!                on input xxx is northern hemisphere fouriers.
!                on input yyy is southern hemisphere fouriers.
!     n        - number of fields of arrays a and b.
!
!   output argument list:
!     a        - on output array elements xxx are the sum of
!                the northern and southern hemisphere fouriers,
!                i.e., the multipliers of symmetric pln.
!     b        - on output array b is the difference of
!                the northern and southern hemisphere fouriers,
!                i.e., the multipliers of antisymmetric pln.
!
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_
!-------------------------------------------------------------------------------
   save
   real  ::  a(lonf_,2,n), b(lonf_,n)
!-------------------------------------------------------------------------------
!
!    on input a corresp. to n. hemis. fouriers
!    on input b corresp. to s. hemis. fouriers
!    on output a corresp. to multipliers of sym. pln
!    on output b corresp. to multipliers of antisym. pln
!
   do k = 1,n
     do i = 1,lonf_
       b(i,k)   = a(i,1,k) - a(i,2,k)
     enddo
     do i = 1,lonf_
       a(i,1,k) = a(i,1,k) + a(i,2,k)
     enddo
   enddo
!
   return
   end
!-------------------------------------------------------------------------------
