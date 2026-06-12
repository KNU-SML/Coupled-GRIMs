#include "define.h"
   subroutine mpsumall(dat,len,ityp,icom,nmem)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  len,ityp,icom,nmem
   real                 ::  dat(len)
   real(_mpi_real_)     ::  tmpsnd(len),tmprcv(len,nmem)
   integer              ::  i,k
!
   if (nmem.eq.1) then
     return
   endif
!
   tmpsnd(:)=dat(:)
!    
   CALL mpi_allgather(tmpsnd, len, ityp, tmprcv, len, ityp, icom)
!
   dat(:)=0.0
   do i = 1,nmem
     do k = 1,len
       dat(k)=dat(k)+tmprcv(k,i)
     enddo
   enddo
!
   return
   end subroutine mpsumall
!-------------------------------------------------------------------------------
