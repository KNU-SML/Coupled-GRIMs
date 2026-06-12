#include <define.h>
   subroutine mplatall(a,latgp,b,latg,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mplatall
!            
! abstract: transpose (jp,ntotal) to (jf,ntotal)
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mplatall(a,latgp,b,latg,ntotal)
!
!    input argument lists:
!   a   - real (latgp,ntotal) partial field for each pe 
!   latgp   - integer partial grid in latitude divide 2
!   latg   - integer total grid in latitude divide 2
!   ntotal   - integer total set of fields
!
!    output argument list:
!   b   - real (latg,ntotal) full field 
! 
! subprograms called:
!   mpi_scatterv   - scatter message from master pe to all pe
!   mpi_bcast      - boardcast to all
!
!-------------------------------------------------------------------------------
!  use commpi, only : nrow, myrow, mype, latlen, latstr, comm_column,          &
   use commpi       ! nrow, myrow, mype, latlen, latstr, comm_column,          &
                    ! MPI_REAL, ncol, latdef
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer            ::  latgp,latg,ntotal,                                   &
                          jj,ji,n,j,k,mk,nk,lens,ip,ierr
   real               ::  a(latgp,ntotal),b(latg,ntotal)
!
   real(_mpi_real_),allocatable::  tmpsnd(:),tmprcv(:)
   integer,allocatable::  len(:),loc(:)
!
   allocate(tmpsnd(latgp*ntotal))
   allocate(tmprcv(latg *ntotal))
   allocate(len(0:nrow-1))
   allocate(loc(0:nrow-1))
!
   mk=0
   do k = 1,ntotal
     do j = 1,latlen(mype)
       mk=mk+1
       tmpsnd(mk)=a(j,k)
     enddo
   enddo
!
   if( myrow.eq.0 ) then
     mk=0
     do n = 0,nrow-1
       loc(n)=mk
       ip=mype+ncol*n
       len(n)=ntotal*latlen(ip)
       mk=mk+len(n)
     enddo
   endif
   lens=ntotal*latlen(mype)
!
! --- not working for allgatherv, dont know why ----
!     call mpi_allgatherv(tmpsnd,lens,MPI_REAL,
!    1     tmprcv,len(0),loc(0),MPI_REAL,comm_column,ierr)
!
! --- works for gatherv then bcast ----
!
   call mpi_gatherv(tmpsnd,lens,MPI_REAL,                                      &
        tmprcv,len(0),loc(0),MPI_REAL,0,comm_column,ierr)
   lens=ntotal*latg
   call mpi_bcast(tmprcv,lens,MPI_REAL,0,comm_column,ierr)
!
   mk=0
   do n = 0,nrow-1
     ip=ncol*n
     do k = 1,ntotal
       do j = 1,latlen(ip)
         jj=j+latstr(ip)-1
         ji=latdef(jj)
         mk=mk+1
         b(ji,k)=tmprcv(mk)
       enddo
     enddo
   enddo
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
   return
   end subroutine mplatall
!-------------------------------------------------------------------------------
