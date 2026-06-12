#include <define.h>
   subroutine mpgfpk2fk(a,idim,nlen,nstr,nend,b,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpgfpk2fk
!            
! abstract:  gather pgb records to master processor
!
! program history log:
!    04-08-27  H.Kanamaru and M.Kanamitsu  the first version
!
! usage:   call mpgfpk2fk(a,idim,nlen,nstr,nend,b,ntotal)
!
!    input argument lists:
!   a   - 
!   idim   - 
!   nlen   - 
!   nstr   - 
!   nend   - 
!   ntotal   - 
!
!    output argument list:
!   b   - 
!   nlen   - 
! 
! subprograms called:
!   mpi_gatherv   - gather message from all pe to master
!   mpi_allgatherv - gather message from all pe to all pe
!
!-------------------------------------------------------------------------------
   use paramodel, only : npes_
   use commpi       ! mype, mpi_integer, mpi_comm_world, mpi_character, npes_
#ifdef REDUCE_GRID
   use comreduce    ! lonfd, lonfdp, lonfds
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer               ::  idim,ntotal,                                      &
                             n,i,j,k,mk,kk,ierr,nlocac,lenac
   character             ::  a(idim,ntotal),b(idim,ntotal)
   integer               ::  nstr(0:npes_-1),nend(0:npes_-1),nstrend(0:npes_-1)
   integer               ::  nlen(ntotal)
   character, allocatable::  tmpsnd(:)
   character, allocatable::  tmprcv(:)
   integer,allocatable   ::  len(:),loc(:),mlen(:),nloc(:)
!
   npes=npes_
   allocate(tmpsnd(idim*(nend(mype)-nstr(mype)+1)))
   allocate(tmprcv(idim*ntotal))
   allocate(len(npes_))
   allocate(loc(0:npes_-1))
   allocate(mlen(ntotal))
   allocate(nloc(0:npes_-1))
!
! make nstrend from nstr and nend, number of record fileds in each pe
! also make nloc, accumulated number of record files from pe=0 
!
   nlocac=0 
   do i = 0,npes_-1
     nloc(i)=nlocac
     nstrend(i)=nend(i)-nstr(i)+1
     nlocac=nlocac+nstrend(i)
   enddo
!
! gather and bcst nlen into mlen, array count for each record
!
!soojin_couple
!   call mpi_allgatherv(nlen,nstrend(mype),mpi_integer,                         &
!       mlen,nstrend(0),nloc(0),mpi_integer,mpi_comm_world,ierr)
   call mpi_allgatherv(nlen,nstrend(mype),mpi_integer,                         &
       mlen,nstrend(0),nloc(0),mpi_integer,mpi_comm_private,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! rearrange a into tmpsnd
!
   mk=0
   do k = nstr(mype),nend(mype)
     do i = 1,nlen(k-nstr(mype)+1)
       mk=mk+1
       tmpsnd(mk)=a(i,k-nstr(mype)+1)
     enddo
   enddo
!
! make len from mlen, array count for each pe
! rearrange len into loc, for use in gatherv of tmpsnd
!
   do i = 1,npes
     len(i)=0
   enddo
   lenac=0
   do j = 0,npes-1
     do i = nloc(j)+1,nloc(j)+nstrend(j)
       len(j+1)=len(j+1)+mlen(i)
     enddo
     loc(j)=lenac
     lenac=lenac+len(j+1)
   enddo
!
! gather to master, tmpsnd to tmprcv
!
!soojin_couple
!   call mpi_gatherv(tmpsnd,len(mype+1),mpi_character,                          &
!        tmprcv,len,loc(0),mpi_character,0,mpi_comm_world,ierr)
   call mpi_gatherv(tmpsnd,len(mype+1),mpi_character,                          &
        tmprcv,len,loc(0),mpi_character,0,mpi_comm_private,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! rearrange tmprcv into b
!
   if( mype.eq.0 ) then
     mk=0
     do k = 1,ntotal
       do i = 1,mlen(k)
         mk=mk+1
         b(i,k)=tmprcv(mk)
       enddo
     enddo
   endif
!
! pass mlen to nlen
!
   do i = 1,ntotal
     nlen(i)=mlen(i)
   enddo
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
   deallocate(mlen)
   deallocate(nloc)
!
   return
   end subroutine mpgfpk2fk
!-------------------------------------------------------------------------------
