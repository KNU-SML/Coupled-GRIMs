#include <define.h>
   subroutine mpsf2p(a,lnt,b,lntpp,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpsf2p
!            
! abstract: transpose (lnf,kf) to (lnpp,kf)
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpsf2p(a,lnt,b,lntpp,ntotal)
!
!    input argument lists:
!   a   - real (lnt,ntotal) total field
!   lnt   - integer total spectral grid
!   lntpp   - integer sub partial spectral grid
!   ntotal   - integer total set of fields
!
!    output argument list:
!   b   - real (lntpp,ntotal) sub partial field
! 
! subprograms called:
!   spcshfli      - shafle spectral for balance from input
!   spcshflo      - shafle spectral for output
!   mpi_scatterv  - send to all pe in the same comm
!
!-------------------------------------------------------------------------------
   use paramodel, only : jcap_
   use commpi       ! npes, lntlen, lntstr, mype, master,                      &
                    ! mpi_comm_world, MPI_REAL, lwvdef
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer            ::  lnt,lntpp,ntotal,                                    &
                          n,m,k,mk,ierr
   real               ::  a(lnt,ntotal),b(lntpp,ntotal)
!
   real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
   integer,allocatable::  len(:),loc(:)
!
   allocate(tmpsnd(lnt*ntotal))
   allocate(tmprcv(lntpp*ntotal))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then
      call spcshfli(a,lnt,ntotal,jcap_,lwvdef)
      mk=0
      do n=0,npes-1
         loc(n)=mk
         do m=1,lntlen(n)*2
            do k=1,ntotal
               mk=mk+1
               tmpsnd(mk)=a(m+lntstr(n)*2,k)
            enddo
         enddo
         len(n)=ntotal*lntlen(n)*2
      enddo
   else
      len(mype)=ntotal*lntlen(mype)*2
   endif
!
!soojin_couple
!   call mpi_scatterv(tmpsnd,len(0),loc(0),MPI_REAL,                            &
!        tmprcv,len(mype),MPI_REAL,0,mpi_comm_world,ierr)
   call mpi_scatterv(tmpsnd,len(0),loc(0),MPI_REAL,                            &
        tmprcv,len(mype),MPI_REAL,0,mpi_comm_private,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   mk=0
   do m=1,lntlen(mype)*2
      do k=1,ntotal
         mk=mk+1
         b(m,k)=tmprcv(mk)
      enddo
   enddo
!
   if( mype.eq.0 ) call spcshflo(a,lnt,ntotal,jcap_,lwvdef)
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
   return
   end
