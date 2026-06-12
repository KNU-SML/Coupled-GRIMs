#include <define.h>
   subroutine mpnx2x(a,lonf2p,lota,b,lonf2,lotb,latg2p,leva,levb,nvar)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpnx2x
!            
! abstract:   transpose (ip,k,jp) to (if,k,jp)
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpnx2x(a,lonf2p,lota,b,lonf2,lotb,latg2p,
!        &                  leva,levb,nvar)
!
!    input argument lists:
!   a   - real (lonf2p,lota,latg2p) partial field in i j
!   lonf2p   - integer partial longitude grid x 2
!   latg2p   - integer partial latitude grid / 2
!   lonf2   - integer total longitude grid x 2
!   leva   - integer starting vertical layer for a
!   levb   - integer starting vertical layer for b
!   nvar   - integer total set of fields
!
!    output argument list:
!   b   - real (lonf2,lotb,latg2p) partial field in j k
! 
! subprograms called:
!   mpi_alltoallv  - send and receive from all pe in the same comm
!
!-------------------------------------------------------------------------------
   use commpi       ! ncol, mype, latstr, latlen, lonstr, lonlen,              &
                    ! MPI_REAL, comm_row, latdef
#ifdef REDUCE_GRID
   use comreduce    ! lonfdp, lonfd, lonfds
#else
   use paramodel, only : lonf_
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer            ::  lonf2p,lonf2,latg2p,lota,lotb,nvar,leva,levb,status, &
                          i,j,k,ii,ipe,ipe0,ioff,noffa,noffb,n,mn,len,ierr,    &
                          lonf2d,lonfp,lonff,latlon
   real               ::  a(lonf2p,lota,latg2p),b(lonf2,lotb,latg2p)
!
   real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
   integer,allocatable::lensnd(:),lenrcv(:)
   integer,allocatable::locsnd(:),locrcv(:)
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do j = 1,latg2p
       do k = 1,lota
         do i = 1,lonf2
           b(i,k,j) = a(i,k,j)
         enddo
       enddo
     enddo
     return
   endif
!
! otherwise
!
   allocate(tmpsnd(lonf2p*latg2p*nvar*ncol))
   allocate(tmprcv(lonf2p*latg2p*nvar*ncol))
   allocate(lensnd(ncol))
   allocate(lenrcv(ncol))
   allocate(locsnd(ncol))
   allocate(locrcv(ncol))
!
   ipe0 = int(mype/ncol)*ncol
!
   mn = 0
#ifndef REDUCE_GRID
   lonfp = lonlen(mype)
#endif
   do ii = 1,ncol
     locsnd(ii) = mn
     do n = 1,nvar
       noffa = leva-1+n
       do j = 1,latlen(mype)
#ifdef REDUCE_GRID
         lonfp = lonfdp(j,mype)
#endif
         do i = 1,lonfp
           mn = mn+1
           tmpsnd(mn) = a(i      ,noffa,j)
           mn = mn+1
           tmpsnd(mn) = a(i+lonfp,noffa,j)
         enddo
       enddo
     enddo
     lensnd(ii) = mn-locsnd(ii)
   enddo
!
   mn = 0
   do ii = 1,ncol
     locrcv(ii) = mn
     ipe = (ii-1)+ipe0
#ifdef REDUCE_GRID
     latlon = 0
     do j = 1,latlen(mype)
       latlon = latlon+lonfdp(j,ipe)
     enddo
#else
     latlon = latlen(mype)*lonlen(ipe)
#endif
     lenrcv(ii) = latlon*2*nvar
     mn = mn+lenrcv(ii)
   enddo
!
   call mpi_alltoallv(tmpsnd,lensnd,locsnd,MPI_REAL,                           &
                      tmprcv,lenrcv,locrcv,MPI_REAL,comm_row,ierr)
!
   mn = 0
   do ii = 1,ncol
     ipe = (ii-1)+ipe0
#ifndef REDUCE_GRID
     lonff = lonf_
     lonfp = lonlen(ipe)
     ioff = lonstr(ipe)-1
#endif
     do n = 1,nvar
       noffb = levb-1+n
       do j = 1,latlen(mype)
#ifdef REDUCE_GRID
         lonff = lonfd(latdef(latstr(mype)+j-1))
         lonfp = lonfdp(j,ipe)
         ioff = lonfds(j,ipe)-1
#endif
         do i = 1,lonfp
           mn = mn+1
           b(i+ioff      ,noffb,j) = tmprcv(mn)
           mn = mn+1
           b(i+ioff+lonff,noffb,j) = tmprcv(mn)
         enddo
       enddo
     enddo
   enddo
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(lensnd)
   deallocate(lenrcv)
   deallocate(locsnd)
   deallocate(locrcv)
!
   return
   end subroutine mpnx2x
!-------------------------------------------------------------------------------
