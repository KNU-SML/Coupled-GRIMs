#include <define.h>
   subroutine chgr_hybrid_coeff(km,si,sl,ak,bk)
!-------------------------------------------------------------------------------
!
! subprogram: chgr_hybrid_coeff
!
! abstract: for sigma pressure re: hybrids implemented may 1st 2007 @ 12Z
!
!------------------------------------------------------------------------------
   real  ::  si(km+1),sl(km),ak(km+1),bk(km+1)
!------------------------------------------------------------------------------
!
! set ak and bk to related to centibar
!
   do k = 1,km
     ak(k) = si(k) / 1000.  ! Pa -> centibar 
     bk(k) = sl(k)
   enddo
!
   ak(km+1) = 0.0
   bk(km+1) = 0.0
!
   do k = 1,km+1
     print *,' k ak bk ',k,ak(k),bk(k)
   enddo
!
   return
   end subroutine chgr_hybrid_coeff
!------------------------------------------------------------------------------
