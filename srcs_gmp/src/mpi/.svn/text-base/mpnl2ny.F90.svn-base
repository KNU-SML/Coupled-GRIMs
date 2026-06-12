#include <define.h>
   subroutine mpnl2ny(a,lonfp,latg,b,lonf,latgp,ntotal,kstr,klen)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpnl2ny
!            
! abstract:   transpose (ip,kp,jf) to (if,kp,jp)
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call mpnl2ny(a,lonfp,latg,b,lonf,latgp,ntotal,kstr,klen)
!
!    input argument lists:
!   a   - real (lonfp,ntotal,latg) partial field in i k
!   lonfp   - integer partial longitude grid x 2
!   latg   - integer total latitude grid / 2
!   lonf   - integer total longitude grid x 2
!   latgp   - integer partial latitude grid / 2
!   kstr   - integer starting set of fields
!   klen   - integer partial vertical layer
!   ntotal   - integer total set of fields
!
!    output argument list:
!   b   - real (lonf,ntotal,latgp) partial field in j k
! 
! subprograms called:
!   mpi_alltoallv  - send and receive from all pe in the same comm
!
!-------------------------------------------------------------------------------
   use paramodel, only : jcapp_, lonf_
   use commpi       ! ncol, mype, latlen, latstr, lwvlen, lwvstr,              &
                    ! latdef, MPI_REAL, comm_column, nrow, lwvdef
#ifdef REDUCE_GRID
   use comreduce    ! lonfd, lonfds, lonfdp, lcapdp
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer            ::  lonf,latgp,ntotal,lonfp,latg,status,lcaps,           &
                          i,i1,i2,j1,j,k,jj,jpe,jpe0,n,mn,len,ierr,            &
                          kstr,klen,ij,ji,llens                                
   real               ::  a(lonfp,ntotal,latg),b(lonf,ntotal,latgp)
!
   real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
   integer,allocatable::  lensnd(:),lenrcv(:)
   integer,allocatable::  locsnd(:),locrcv(:)
!
   allocate(tmpsnd(lonfp*latgp*klen*nrow))
   allocate(tmprcv(lonfp*latgp*klen*nrow))
   allocate(lensnd(nrow))
   allocate(lenrcv(nrow))
   allocate(locsnd(nrow))
   allocate(locrcv(nrow))
!
   jpe0=mod(mype,ncol)
   lcaps=(jcapp_+1)*2
!
! cut in y
!
   mn=0
#ifndef REDUCE_GRID
   llens=lwvlen(mype)*2
#endif
   do jj = 1,nrow
     locsnd(jj)=mn
     jpe=jpe0+(jj-1)*ncol
     do n = kstr,kstr+klen-1
       do j = 1,latlen(jpe)
         j1=j+latstr(jpe)-1
         ji=latdef(j1)
#ifdef REDUCE_GRID
         llens=lcapdp(ji,mype)*2
#endif
         do i = 1,llens
           mn=mn+1
           tmpsnd(mn)=a(i      ,n,ji)
           mn=mn+1
           tmpsnd(mn)=a(i+lcaps,n,ji)
         enddo
       enddo
     enddo
     lensnd(jj)=mn-locsnd(jj)
   enddo
!
   mn=0
   do jj = 1,nrow
     locrcv(jj)=mn
     jpe=jpe0+(jj-1)*ncol
#ifdef REDUCE_GRID
     llens=0
     do j = 1,latlen(mype)
       j1=j+latstr(mype)-1
       ji=latdef(j1)
       llens=llens+lcapdp(ji,jpe)
     enddo
     lenrcv(jj)=llens*4*klen
#else
     lenrcv(jj)=latlen(mype)*lwvlen(jpe)*4*klen
#endif
     mn=mn+lenrcv(jj)
   enddo
!
   call mpi_alltoallv(tmpsnd,lensnd,locsnd,MPI_REAL,                           &
                      tmprcv,lenrcv,locrcv,MPI_REAL,comm_column,ierr)
!
! restore l
!
   mn=0
   do jj = 1,nrow
     jpe=jpe0+(jj-1)*ncol
#ifndef REDUCE_GRID
     llens=lwvlen(jpe)
#endif
     do n = kstr,kstr+klen-1
       do j = 1,latlen(mype)
         i1=lwvstr(jpe)
#ifdef REDUCE_GRID
         j1=latdef(latstr(mype)+j-1)
         llens=lcapdp(j1,jpe)
#endif
         do i = 1,llens
           ij=lwvdef(i+i1)*2
           mn=mn+1
           b(ij+1      ,n,j)=tmprcv(mn)
           mn=mn+1
           b(ij+1+lonf_,n,j)=tmprcv(mn)
           mn=mn+1
           b(ij+2      ,n,j)=tmprcv(mn)
           mn=mn+1
           b(ij+2+lonf_,n,j)=tmprcv(mn)
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
   end subroutine mpnl2ny
!-------------------------------------------------------------------------------
