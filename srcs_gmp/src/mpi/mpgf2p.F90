#include <define.h>
   subroutine mpgf2p(a,lonf2,latg2,b,lonf2p,latg2p,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpgf2p
!            
! abstract: transpose (if,jf,kf) to (ip,jp,kf)
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpgf2p(a,lonf2,latg2,b,lonf2p,latg2p,ntotal)
!
!    input argument lists:
!   a   - real (lonf2,latg2,ntotal) full field 
!   lonf2   - integer total grid in longitude time 2
!   latg2   - integer total grid in latitude divide 2
!   ntotal   - integer total set of fields
!
!    output argument list:
!   b   - real (lonf2p,latg2p,ntotal) partial field for each pe 
!   lonf2p   - integer partial grid in longitude time 2
!   latg2p   - integer partial grid in latitude divide 2
! 
! subprograms called:
!   mpi_scatterv   - scatter message from master pe to all pe
!
!-------------------------------------------------------------------------------
   use commpi       ! npes, mype, master, lonstr, lonlen, latstr, latlen,      &
                    ! MPI_REAL, latdef, mpi_comm_world
#ifdef REDUCE_GRID
   use comreduce    ! lonfds, lonfdp, lonfd
#else
   use paramodel, only : lonf_
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer            ::  lonf2,latg2,lonf2p,latg2p,ntotal,                    &
                          ii,jj,n,mk,i,j,k,ierr,lonlend,lonstrd,lonff
   real               ::  a(lonf2,latg2,ntotal),b(lonf2p,latg2p,ntotal)
   real(_mpi_real_)   ::  tmpsnd(lonf2 *latg2 *ntotal)
   real(_mpi_real_)   ::  tmprcv(lonf2p*latg2p*ntotal)
   integer            ::  len(0:npes-1)
   integer            ::  loc(0:npes-1)
!
   if( mype.eq.master ) then
     mk=0
     do n = 0,npes-1
       loc(n)=mk
#ifndef REDUCE_GRID
       lonlend=lonlen(n)
       lonstrd=lonstr(n)-1
       lonff=lonf_
#endif
       do k = 1,ntotal
         do j = 1,latlen(n)
           jj=j+latstr(n)-1
           jj=latdef(jj)
#ifdef REDUCE_GRID
           lonlend=lonfdp(j,n)
           lonstrd=lonfds(j,n)-1
           lonff=lonfd(jj)
#endif
           do i = 1,lonlend
             ii=i+lonstrd
             mk=mk+1
             tmpsnd(mk)=a(ii,jj,k)
             ii=ii+lonff
             mk=mk+1
             tmpsnd(mk)=a(ii,jj,k)
           enddo
         enddo
       enddo
       len(n)=mk-loc(n)
     enddo
   else
#ifdef REDUCE_GRID
     mk=0
     do j = 1,latlen(mype)
       mk=mk+lonfdp(j,mype)
     enddo
     len(mype)=ntotal*mk*2
#else
     len(mype)=ntotal*lonlen(mype)*2*latlen(mype)
#endif
   endif
!
!soojin_couple
!   call mpi_scatterv(tmpsnd,len(0),loc(0),MPI_REAL,                            &
!                        tmprcv,len(mype),MPI_REAL,0,mpi_comm_world,ierr)
   call mpi_scatterv(tmpsnd,len(0),loc(0),mpi_real8,                           &
                        tmprcv,len(mype),mpi_real8,0,mpi_comm_private,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   mk=0
#ifndef REDUCE_GRID
   lonlend=lonlen(mype)
#endif
   do k = 1,ntotal
     do j = 1,latlen(mype)
#ifdef REDUCE_GRID
       lonlend=lonfdp(j,mype)
#endif
       do i = 1,lonlend
         mk=mk+1
         b(i        ,j,k)=tmprcv(mk)
         mk=mk+1
         b(i+lonlend,j,k)=tmprcv(mk)
       enddo
       do i = lonlend*2+1,lonf2p
         b(i,j,k)=tmprcv(mk)
       enddo
     enddo
     do j = latlen(mype)+1,latg2p
       do i = 1,lonf2p
         b(i,j,k)=tmprcv(mk)
       enddo
     enddo
   enddo
!
   return
   end subroutine mpgf2p
!-------------------------------------------------------------------------------
