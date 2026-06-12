#include <define.h>
   subroutine mpn2nn(a,lntp,b,lntpp,nvar)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpn2nn
!            
! abstract:  transpose (lntp,nvar) to (lntpp,nvar)
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpn2nn(a,lntp,b,lntpp,nvar)
!
!    input argument lists:
!   a   - real (lntp,nvar) partial field 
!   lntp   - integer partial spectral grid 
!   lntpp   - integer sub partial spectral grid 
!   nvar   - integer total set of fields
!
!    output argument list:
!   b   - real (lntpp,nvar) sub partial field 
! 
! subprograms called:
!
!-------------------------------------------------------------------------------
!  use commpi, only : ncol, lntstr, lnpstr, mype, lntlen
   use commpi       ! ncol, lntstr, lnpstr, mype, lntlen
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer     ::  lntpp,lntp,nvar,n,m,ipe0,lntstr0
   real        ::  a(lntp,nvar),b(lntpp,nvar)
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do n = 1,nvar
       do m = 1,lntp
         b(m,n) = a(m,n)
       enddo
     enddo
     return
   endif
!
! otherwise
! cut the part and through away the rest.
!
   ipe0=int(mype/ncol)*ncol
   lntstr0=lntstr(mype)*2-lnpstr(ipe0)*2
   do n = 1,nvar
     do m = 1,lntlen(mype)*2
       b(m,n)=a(m+lntstr0,n)
     enddo
   enddo
!
   return
   end subroutine mpn2nn
!-------------------------------------------------------------------------------
