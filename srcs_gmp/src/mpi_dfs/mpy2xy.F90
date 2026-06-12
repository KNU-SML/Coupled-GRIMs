#include "define.h"
   subroutine mpy2xy(a,iba,b,ib,jbw,nvar)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! abstract:  yz2xy  ( global lon to global levs)
!
! usage:   call mpy2xy(a,iba,b,ib,jbw,nvar)
!
!    input argument lists:
!   a   - real (iba,jbw,nvar) partial field 
!   iba   - integer, global longitudal grid number
!   ib   - integer, partial longitudal grid number
!   jbw   - integer, partial latitudal grid number
!   nvar- integer total set of fields
!
!    output argument list:
!   b   - real (ib,jbw,nvar) sub partial field 
! 
! subprograms called:
!
!-------------------------------------------------------------------------------
   use commpi, only  :  ncol,mype,latlen,lonlen
   use dfsvar, only  :  igs,ige,ils,ile
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  iba,ib,jbw,nvar
   real     ::  a(igs:ige,jbw,nvar),b(ib,jbw,nvar)
!
! local
!
   integer  ::  is,i,j,n,ie,ip,jp
!
   ip=lonlen(mype)
   jp=latlen(mype)*2
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do n = 1,nvar
       do j = 1,jp
         do i = 1,iba
           b(igs+i-1,j,n) = a(i,j,n)
         enddo
       enddo
     enddo
     return
   endif
!
! otherwise
! cut the part and through away the rest.
!
   do n = 1,nvar
     do j = 1,jp
       b(1:ip,j,n)=a(ils:ile,j,n)
     enddo
   enddo
!
   return
   end subroutine mpy2xy
!-------------------------------------------------------------------------------
