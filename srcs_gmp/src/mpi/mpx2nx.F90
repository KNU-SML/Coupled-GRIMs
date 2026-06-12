#include <define.h>
   subroutine mpx2nx(a,lonf2,lota,b,lonf2p,lotb,latg2p,leva,levb,nlev)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:   mpx2nx
!
! abstract:   transpose (if,k,jp) to (ip,k,jp)
!
! program history log:
!    99-06-27  henry juang      finish entire test for gsm
!
! usage:        call mpx2nx(a,lonf2,lota,b,lonf2p,lotb,latg2p,
!        &                  leva,levb,nlev)
!
!    input argument lists:
!       a       - real (lonf2,lota,latg2p) partial field in j k
!       lonf2p  - integer partial longitude grid x 2
!       latg2p  - integer partial latitude grid / 2
!       lonf2   - integer total longitude grid x 2
!       leva    - integer starting vertical layer for a
!       levb    - integer starting vertical layer for b
!       nlev    - integer total set of fields
!
!    output argument list:
!       b       - real (lonf2p,lotb,latg2p) partial field in i j
!
! subprograms called:
!
!-------------------------------------------------------------------------------
!  use commpi, only : ncol, mype, latstr, latlen, lonstr, lonlen, latdef
   use commpi       ! ncol, mype, latstr, latlen, lonstr, lonlen, latdef
#ifdef REDUCE_GRID
   use comreduce, only : lonfd, lonfds, lonfdp
#else
   use paramodel, only : lonf_
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer     ::  lonf2,lonf2p,latg2p,lota,lotb,leva,levb,nlev,               &
                   i,j,n,ioff,noffa,noffb,lonf2d,lonfp,lonff
   real        ::  a(lonf2,lota,latg2p),b(lonf2p,lotb,latg2p)
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do j = 1,latg2p
       do n = 1,lota
         do i = 1,lonf2
           b(i,n,j) = a(i,n,j)
         enddo
       enddo
     enddo
     return
   endif
!
! otherwise
! cut the need and through away the remain.
!
#ifndef REDUCE_GRID
   lonff=lonf_
   lonfp=lonlen(mype)
   ioff=lonstr(mype)-1
#endif
   do n = 1,nlev
     noffa=leva-1+n
     noffb=levb-1+n
     do j = 1,latlen(mype)
#ifdef REDUCE_GRID
       lonff=lonfd(latdef(latstr(mype)+j-1))
       lonfp=lonfdp(j,mype)
       ioff=lonfds(j,mype)-1
#endif
       do i = 1,lonfp
         b(i      ,noffb,j)=a(i+ioff      ,noffa,j)
         b(i+lonfp,noffb,j)=a(i+ioff+lonff,noffa,j)
       enddo
     enddo
   enddo
!
   return
   end subroutine mpx2nx
!-------------------------------------------------------------------------------
