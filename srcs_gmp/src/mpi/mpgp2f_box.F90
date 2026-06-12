#include <define.h>
   subroutine mpgp2f_box(a,lonfp,latgp,b,lonf,latg,ntotal)

#ifdef AOMG
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpgp2f
!            
! abstract:  transpose (ip,jp,kf) to (if,jf,kf)
!            return decomposed domain in box shape to entire domain
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpgp2f(a,lonf2p,latg2p,b,lonf2,latg2,ntotal)
!
!    input argument lists:
!   a   - real (lonf2p,latg2p,ntotal) partial field for each pe 
!   lonf2p   - integer partial grid in longitude time 2
!   latg2p   - integer partial grid in latitude divide 2
!   ntotal   - integer total set of fields
!
!    output argument list:
!   b   - real (lonf2,latg2,ntotal) full field 
!   lonf2   - integer total grid in longitude time 2
!   latg2   - integer total grid in latitude divide 2
! 
! subprograms called:
!   mpi_gatherv   - gather message from all pe to master
!
!-------------------------------------------------------------------------------
   use commpi       ! npes, mype, latlen, latstr, master, lonlen, lonstr,      &
                    ! MPI_REAL, mpi_comm_world, latdef
   use paramodel, only : lonf_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer            ::  lonfp,latgp,lonf,latg,ntotal,                    &
                          ii,jj,n,i,j,k,mk,ierr,lonlend,lonstrd,lonff
   real               ::  a(lonfp,latgp,ntotal),b(lonf,latg,ntotal)
   real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
   integer,allocatable::  len(:),loc(:)
!
#ifdef REDUCE_GRID
   print*, "coupled mode cannot run with REDUCE_GRID"
   call abort
#endif /* REDUCE_GRID */

   allocate(tmpsnd(lonfp*latgp*ntotal))
   allocate(tmprcv(lonf *latg *ntotal))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then

!- check oasis
   print*, "=== in mpgp2f_box ==="
!- end check

     mk=0
     do n = 0,npes-1
       loc(n)=mk
       len(n)=ntotal*lonlen(n)*2*latlen(n)
       mk=loc(n)+len(n)
     enddo
   endif
!
   lonlend=lonlen(mype)
   mk=0
   do k = 1,ntotal
     do j = 1,latlen(mype)*2
       do i = 1,lonlend
         mk=mk+1
         tmpsnd(mk)=a(i,j,k)
!         mk=mk+1
!         tmpsnd(mk)=a(i+lonlend,j,k)
       enddo
     enddo
   enddo
   len(mype)=mk
!
   call mpi_gatherv(tmpsnd,len(mype),MPI_REAL,                                 &
        tmprcv,len(0),loc(0),MPI_REAL,0,mpi_comm_private,ierr)
!
   if( mype.eq.0 ) then
     mk=0
     do n = 0,npes-1
       lonlend=lonlen(n)
       lonstrd=lonstr(n)-1
       lonff=lonf_

       do k = 1,ntotal
         do j = 1,latlen(n)*2
           jj=j+(latstr(n)-1)*2
           do i = 1,lonlend
             ii=i+lonstrd
             mk=mk+1
             b(ii,jj,k)=tmprcv(mk)
           enddo
         enddo
       enddo

     enddo
   endif
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
#endif /* AOMG */
   return
   end subroutine mpgp2f_box
