#include <define.h>
   subroutine mpgp2f(a,lonf2p,latg2p,b,lonf2,latg2,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpgp2f
!            
! abstract:  transpose (ip,jp,kf) to (if,jf,kf)
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpgp2f(a,lonf2p,latg2p,b,lonf2,latg2,ntotal)
!
!    input argument lists:
!   b   - real (lonf2p,latg2p,ntotal) partial field for each pe 
!   lonf2p   - integer partial grid in longitude time 2
!   latg2p   - integer partial grid in latitude divide 2
!   ntotal   - integer total set of fields
!
!    output argument list:
!   a   - real (lonf2,latg2,ntotal) full field 
!   lonf2   - integer total grid in longitude time 2
!   latg2   - integer total grid in latitude divide 2
! 
! subprograms called:
!   mpi_gatherv   - gather message from all pe to master
!
!-------------------------------------------------------------------------------
   use commpi, only : real_type,comm_world,lonlen,lonstr,latlen,latstr,        &
                      npes,mype,master
   use dfsvar, only : ib,iba,latdef,jgs
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer            ::  lonf2p,latg2p,lonf2,latg2,ntotal
   integer            ::  ii,jj,n,i,j,k,mk,ierr,ilen,ist,iend
   real               ::  a(lonf2p,latg2p,ntotal),b(lonf2,latg2,ntotal)
!
   real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
   integer,allocatable::  len(:),loc(:)
!
   allocate(tmpsnd(lonf2p*latg2p*ntotal))
   allocate(tmprcv(lonf2 *latg2 *ntotal))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then
     mk=0
     do n = 0,npes-1
       loc(n)=mk
       len(n)=ntotal*lonlen(n)*latlen(n)*2
       mk=loc(n)+len(n)
     enddo
   endif
!
   mk=0
   do k = 1,ntotal
     do j = 1,latlen(mype)
       tmpsnd(mk+1:mk+ib)=a(1:ib,j,k)
       mk=mk+ib
       tmpsnd(mk+1:mk+ib)=a(ib+1:2*ib,j,k)
       mk=mk+ib
     enddo
   enddo
   len(mype)=mk
!
   call mpi_gatherv(tmpsnd,len(mype),real_type,                                &
                    tmprcv,len(0),loc(0),real_type,0,comm_world,ierr)
!
   if( mype.eq.0 ) then
     mk=0
     do n = 0,npes-1
       ilen=lonlen(n)
       ist=lonstr(n)
       do k = 1,ntotal
         do j = 1,latlen(n)
           jj=latdef(j+latstr(n)-jgs)
           b(ist:ist+ilen-1,jj,k)=tmprcv(mk+1:mk+ilen)
           mk=mk+ilen
           b(iba+ist:iba+ist+ilen-1,jj,k)=tmprcv(mk+1:mk+ilen)
           mk=mk+ilen
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
   return
   end subroutine mpgp2f
!-------------------------------------------------------------------------------
