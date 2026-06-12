#include <define.h>
   subroutine mpnx2nk(a,lonf2p,lota,b,lonf2,lotb,latg2p,levs,levsp,            &
                      leva,levb,nvar)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:   mpnx2nk
!
! abstract:  transpose (ip,kf) to (if,kp)
!
! program history log:
!    99-06-27  henry juang      finish entire test for gsm
!
! usage:    call mpnx2nk(a,lonf2p,lota,b,lonf2,lotb,latg2p,levs,levsp,
!          &                   leva,levb,nvar)
!
!    input argument lists:
!       a       - real (lonf2p,lota,latg2p) partial field in j k
!       lonf2   - integer total longitude grid x 2
!       lota    - integer total field of nvar*levsp
!       lonf2p  - integer partial longitude grid x 2
!       lotb    - integer total field of nvar*levs
!       latg2p  - integer partial latitude grid / 2
!       levsp   - integer partial vertical grid
!       levs    - integer total vertical grid
!       leva    - integer start vertical grid for a
!       levb    - integer start vertical grid for b
!       nvar    - integer total set of fields
!
!    output argument list:
!       b       - real (lonf2,lotb,latg2p) partial field in i j
!
! subprograms called:
!   mpi_alltoallv  - send and receive from all pe in the same comm
!
!-------------------------------------------------------------------------------
!  use commpi, only : ncol, mype, comm_row, levstr, levlen, latstr, latlen,    &
   use commpi       ! ncol, mype, comm_row, levstr, levlen, latstr, latlen,    &
                    ! MPI_REAL, latdef
#ifdef REDUCE_GRID
   use comreduce    ! lonfd, lonfds, lonfdp
#else
   use paramodel, only : lonf_
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer            ::  lonf2,levsp,latg2p,lonf2p,levs,nvar,status,          &
                          lota,lotb,leva,levb,ioff,koff,                       &
                          i,j,k,ii,ipe,ipe0,n,mn,len,kp,kf,ierr,               &
                          lonf2d,lonff,lonfp,latlon
   real               ::  a(lonf2p,lota,latg2p),b(lonf2,lotb,latg2p)
!
   real(_mpi_real_),allocatable::  tmpsnd(:),tmprcv(:)
   integer,allocatable::  lensnd(:),lenrcv(:)
   integer,allocatable::  locsnd(:),locrcv(:)
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
   allocate(tmpsnd(lonf2p*latg2p*levsp*nvar*ncol))
   allocate(tmprcv(lonf2p*latg2p*levsp*nvar*ncol))
   allocate(lensnd(ncol))
   allocate(lenrcv(ncol))
   allocate(locsnd(ncol))
   allocate(locrcv(ncol))
!
   ipe0=int(mype/ncol)*ncol
!
! cut in k
!
   mn=0
#ifndef REDUCE_GRID
   lonfp=lonlen(mype)
#endif
   do ii = 1,ncol
     locsnd(ii)=mn
     ipe=(ii-1)+ipe0
     koff=levstr(ipe)-1
     do n = 1,nvar
       kf=(n-1)*levs+leva-1+koff
       do j = 1,latlen(mype)
#ifdef REDUCE_GRID
         lonfp=lonfdp(j,mype)
#endif
         do k = 1,levlen(ipe)
           do i = 1,lonfp
             mn=mn+1
             tmpsnd(mn) = a(i      ,k+kf,j)
             mn=mn+1
             tmpsnd(mn) = a(i+lonfp,k+kf,j)
           enddo
         enddo
       enddo
     enddo
     lensnd(ii)=mn-locsnd(ii)
   enddo
!
   mn=0
   do ii = 1,ncol
     locrcv(ii)=mn
     ipe=(ii-1)+ipe0
#ifdef REDUCE_GRID
     latlon=0
     do j = 1,latlen(mype)
       latlon=latlon+lonfdp(j,ipe)
     enddo
#else
     latlon=latlen(mype)*lonlen(ipe)
#endif
     lenrcv(ii)=levlen(mype)*latlon*2*nvar
     mn=mn+lenrcv(ii)
   enddo
!
   call mpi_alltoallv(tmpsnd,lensnd,locsnd,MPI_REAL,                           &
                      tmprcv,lenrcv,locrcv,MPI_REAL,comm_row,ierr)
!
! put to x
!
   mn=0
   do ii = 1,ncol
     ipe=(ii-1)+ipe0
#ifndef REDUCE_GRID
     lonff=lonf_
     lonfp=lonlen(ipe)
     ioff=lonstr(ipe)-1
#endif
     do n = 1,nvar
       kp=(n-1)*levsp+levb-1
       do j = 1,latlen(mype)
#ifdef REDUCE_GRID
         lonff=lonfd(latdef(latstr(mype)+j-1))
         lonfp=lonfdp(j,ipe)
         ioff=lonfds(j,ipe)-1
#endif
         do k = 1,levlen(mype)
           do i = 1,lonfp
             mn=mn+1
             b(i+ioff      ,k+kp,j) = tmprcv(mn)
             mn=mn+1
             b(i+ioff+lonff,k+kp,j) = tmprcv(mn)
           enddo
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
   end subroutine mpnx2nk
!-------------------------------------------------------------------------------
