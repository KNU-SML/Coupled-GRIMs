#include <define.h>
   subroutine mpsp2f(a,lntpp,b,lnt,ntotal)
!-------------------------------------------------------------------------------
! subprogram documentation block
!
! subprogram:   mpsp2f
!
! abstract: transpose (lnpp,kf) to (lnf,kf)
!
! program history log:
!    99-06-27  henry juang      finish entire test for gsm
!
! usage:        call mpsf2p(a,lnt,b,lntpp,ntotal)
!
!    input argument lists:
!       a       - real (lntpp,ntotal) sub partial field
!       lnt     - integer total spectral grid
!       lntpp   - integer sub partial spectral grid
!       ntotal  - integer total set of fields
!
!    output argument list:
!       b       - real (lnt,ntotal) total field
!
! subprograms called:
!   spcshflo    - shafle spectral for output
!   mpi_gatherv  - gather all pe in the same comm
!
!-------------------------------------------------------------------------------
   use paramodel, only : jcap_
   use commpi          ! lntlen, lntstr, lnpstr, mype, master, npes,           &
                       ! mpi_comm_world, MPI_REAL, lwvdef
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer            ::  lnt,lntpp,ntotal,                                    &
                          n,m,k,mk,ierr
   real               ::  a(lntpp,ntotal),b(lnt,ntotal)
!
   real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
   integer,allocatable::  len(:),loc(:)
!
   allocate(tmpsnd(lntpp*ntotal))
   allocate(tmprcv(lnt*ntotal))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then
     mk=0
     do n = 0,npes-1
       loc(n)=mk
       do m = 1,lntlen(n)*2
         do k = 1,ntotal
           mk=mk+1
         enddo
       enddo
       len(n)=ntotal*lntlen(n)*2
     enddo
   else
     len(mype)=ntotal*lntlen(mype)*2
   endif
!
   mk=0
   do m = 1,lntlen(mype)*2
     do k = 1,ntotal
       mk=mk+1
       tmpsnd(mk)=a(m,k)
     enddo
   enddo
!
!soojin_couple
!   call mpi_gatherv(tmpsnd,len(mype),MPI_REAL,                                 &
!                    tmprcv,len(0),loc(0),MPI_REAL,0,mpi_comm_world,ierr)
   call mpi_gatherv(tmpsnd,len(mype),MPI_REAL,                                 &
                    tmprcv,len(0),loc(0),MPI_REAL,0,mpi_comm_private,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   if( mype.eq.0 ) then
     mk=0
     do n = 0,npes-1
       do m = 1,lntlen(n)*2
         do k = 1,ntotal
           mk=mk+1
           b(m+lntstr(n)*2,k)=tmprcv(mk)
         enddo
       enddo
     enddo
   endif
!
   if( mype.eq.master ) call spcshflo(b,lnt,ntotal,jcap_,lwvdef)
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
   return
   end subroutine mpsp2f
!-------------------------------------------------------------------------------
