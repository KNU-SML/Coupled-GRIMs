#include "define.h"
   subroutine mpmz2z(a,mt,lats,latd,b,mtg,latgp,ntotal,kstr,klen)
!------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpnl2ny
!            
! abstract:   transpose (ip,kp,jf) to (if,kp,jp)
!
! usage:   call mpnl2ny(a,mt,lats,latd,b,mtg,latgp,ntotal,kstr,klen)
!
!    input argument lists:
!   a   - real (mt,ntotal,latg) partial field in i k
!   mt   - integer partial longitude grid x 2
!   latg   - integer total latitude grid / 2
!   mtg   - integer total longitude grid x 2
!   latgp   - integer partial latitude grid / 2
!   kstr   - integer starting set of fields
!   klen   - integer partial vertical layer
!   ntotal   - integer total set of fields
!
!    output argument list:
!   b   - real (mtg,ntotal,latgp) partial field in j k
! 
! subprograms called:
!   mpi_alltoallv  - send and receive from all pe in the same comm
!
!------------------------------------------------------------------------------
   use commpi, only  :  mype,real_type,comm_col,nrow,ncol,npes,mycol,lmstr,    &
                        lmlen,latstr
   use dfsvar, only  :  jgs,jge,mgs,mge
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer  ::  mtg,latgp,ntotal,mt,lats,latd,status,                          &
                i,i1,i2,j,k,jj,jpe,n,mn,len,ierr,                              &
                kstr,klen,ij,ji,llens,latl
   real  ::  a(mt,latd,ntotal),b(mgs:mge,latgp,ntotal)
!
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer  ::  lensnd,locsnd
   integer, allocatable  ::  locrcv(:),lenrcv(:)
!
   latl=latd-lats
!
   allocate(tmpsnd(mt*latl*klen))
   allocate(tmprcv(mtg*latl*klen))
   allocate(lenrcv(nrow))
   allocate(locrcv(nrow))
!
! cut in y
!
   mn=0
   do n = kstr,kstr+klen-1
     do j = lats+1,latd
       tmpsnd(mn+1:mn+mt)=a(1:mt,j,n)
       mn=mn+mt
     enddo
   enddo
   lensnd=mn
!
   mn=0
   do jj = 1,nrow
     locrcv(jj)=mn
     jpe=(jj-1)*ncol
     lenrcv(jj)=latl*lmlen(jpe)*klen
     mn=mn+lenrcv(jj)
   enddo
!
   call mpi_gatherv(tmpsnd,lensnd,real_type,                                   &
                    tmprcv,lenrcv,locrcv,real_type,0,comm_col,ierr)
!
! restore l
!
   if (mycol.eq.ncol-1) then
     mn=0
     do jj = 1,nrow
       jpe=(jj-1)*ncol
       i1=lmstr(jpe)
       llens=lmlen(jpe)
       do n = kstr,kstr+klen-1
         do j = 1,latl
           b(i1:i1+llens-1,j,n)=tmprcv(mn+1:mn+llens)
           mn=mn+llens
         enddo
       enddo
     enddo
   endif
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(lenrcv)
   deallocate(locrcv)
!
   return
   end subroutine mpmz2z
!-------------------------------------------------------------------------------
