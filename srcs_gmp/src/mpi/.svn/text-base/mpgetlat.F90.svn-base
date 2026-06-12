#include <define.h>
   subroutine mpgetlat(a,latgp,b,latg,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpgetlat
!            
! abstract: transpose (ntotal,jp) to (ntotal,jf)
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpgetlat(a,latgp,b,latg,ntotal)
!
!    input argument lists:
!   a   - real (ntotal,latgp) any array
!   latgp   - integer partial number of latitudes
!   ntotal   - integer group of latitude
!
!    output argument list:
!   b   - real (ntotal,latg) any array
!   latg   - integer total number of latitudes
! 
! subprograms called:
!   mpi_gatherv   - as gather message with variated length
!
!-------------------------------------------------------------------------------
   use commpi          ! npes, ncol, mycol, nrow, mype, latstr, latlen,        &
                       ! MPI_REAL, comm_column
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer            ::  latgp,latg,ntotal,                                   &
                          jj,n,j,k,mk,nk,ip,ierr
   real               ::  a(ntotal,latgp),b(ntotal,latg)
   real(_mpi_real_),allocatable::  tmpsnd(:),tmprcv(:)
   integer,allocatable::  len(:),loc(:)
!
   allocate(tmpsnd(latgp*ntotal))
   allocate(tmprcv(latg *ntotal))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mycol.eq.0 ) then
!
     mk=0
     do k = 1,ntotal
       do j = 1,latlen(mype)
         mk=mk+1
         tmpsnd(mk)=a(k,j)
         mk=mk+1
         tmpsnd(mk)=a(k,latgp+1-j)
       enddo
     enddo
!
     if( mype.eq.0 ) then
       mk=0
       do n = 0,nrow-1
         loc(n)=mk
         ip=mype+ncol*n
         len(n)=ntotal*latlen(ip)*2
         mk=mk+len(n)
       enddo
     else
       len(mype)=ntotal*latlen(mype)*2
     endif
!
     call mpi_gatherv(tmpsnd,len(mype),MPI_REAL,                               &
        tmprcv,len(0),loc(0),MPI_REAL,0,comm_column,ierr)
!
     if( mype.eq.0 ) then
       mk=0
       do n = 0,nrow-1
         ip=mype+ncol*n
         do k = 1,ntotal
           do j = 1,latlen(ip)
             jj=j+latstr(ip)-1
             mk=mk+1
             b(k,jj)=tmprcv(mk)
             jj=latg+1-jj
             mk=mk+1
             b(k,jj)=tmprcv(mk)
           enddo
         enddo
       enddo
     endif
!
   endif
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
   return
   end subroutine mpgetlat
!-------------------------------------------------------------------------------
