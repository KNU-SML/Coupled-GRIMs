#include "define.h"
   subroutine mpgf2pd(tot,iba,jbwa,par,ib,jbw,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpgf2p
!            
! abstract: transpose (if,jf,kf) to (ip,jp,kf)
!
! usage:   call mpgf2p(a,iba,jbwa,b,ib,jbw,ntotal)
!
!    input argument lists:
!   a   - real (iba,jbwa,ntotal) full field 
!   iba   - integer total grid in longitude time 2
!   jbwa   - integer total grid in latitude divide 2
!   ntotal   - integer total set of fields
!
!    output argument list:
!   b   - real (ib,jbw,ntotal) partial field for each pe 
!   ib   - integer partial grid in longitude time 2
!   jbw   - integer partial grid in latitude divide 2
! 
! subprograms called:
!   mpi_scatterv   - scatter message from master pe to all pe
!
!-------------------------------------------------------------------------------
   use commpi, only : npes,mype,master,comm_world,real_type,lonstr,lonlen,     &
                      latstr,latlen,ncol
   use dfsvar, only : igs,ige,jgs,latdef
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  iba,jbwa,ib,jbw,ntotal
   real  ::  tot(igs:ige,jbwa,ntotal),par(ib,jbw,ntotal)
!
! local
!
   integer  ::  ilen,jlen,ip,jp
   integer  ::  ii,jj,n,mk,i,j,k,ierr,ist
   real(_mpi_real_),allocatable  ::  tmpsnd(:),tmprcv(:)
   integer,allocatable  ::  len(:),loc(:)
!
   ip=lonlen(mype)
   jp=latlen(mype)
   allocate(tmpsnd(iba *jbwa *ntotal))
   allocate(tmprcv(ip*jp*2*ntotal))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then
     mk=0
     do n = 0,npes-1
       loc(n)=mk
       ilen=lonlen(n)
       ist=lonstr(n)-1
       do k = 1,ntotal
         do j = 1,latlen(n)
           jj=latdef(latstr(n)+j-jgs)
           tmpsnd(mk+1:mk+ilen)=tot(ist+1:ist+ilen,jj,k)
           mk=mk+ilen
           tmpsnd(mk+1:mk+ilen)=tot(ist+1:ist+ilen,jbwa-jj+1,k)
           mk=mk+ilen
         enddo
       enddo
       len(n)=mk-loc(n)
     enddo
   else
     len(mype)=ntotal*lonlen(mype)*latlen(mype)*2
   endif
!
   call mpi_scatterv(tmpsnd,len(0),loc(0),real_type,tmprcv,len(mype),          &
                     real_type,0,comm_world,ierr)
!
   mk=0
   do k = 1,ntotal
     do j = 1,jp
       par(1:ip,j,k)=tmprcv(mk+1:mk+ip)
       mk=mk+ip
       par(1:ip,jbw-j+1,k)=tmprcv(mk+1:mk+ip)
       mk=mk+ip
     enddo
   enddo
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
   return
   end subroutine mpgf2pd
!-------------------------------------------------------------------------------
