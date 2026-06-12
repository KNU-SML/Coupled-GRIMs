#include "define.h"
   subroutine mpm2f(par,lmp,tot,lmg,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:   mpm2f
!
! usage:        call mpm2f(par,lmp,tot,lmg,ntotal)
!
!    input argument lists:
!       par     - real (lmp,ntotal) sub partial field
!       lmp     - integer total spectral grid
!       ntotal  - integer total set of fields
!
!    output argument list:
!       tot       - real (lmg,ntotal) total field
!
! subprograms called:
!   mpi_gatherv  - gather all pe in the same comm
!
!-------------------------------------------------------------------------------
   use commpi, only  :  comm_col,real_type,master,mype,myrow,                  &
                        lmlen,lmstr,ncol,nrow
   use dfsvar, only  :  mgs,mge
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  lmg,lmp,ntotal
   real  ::  par(lmp,ntotal),tot(mgs:mge,ntotal)
!
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer  ::  ms,me,jpe0,jpe,msm,mem,lmm,n,k,mk,ierr
   integer  ::  nr,lm
   integer,allocatable  ::  len(:),loc(:)
!
   allocate(tmpsnd(lmp*ntotal))
   allocate(tmprcv(lmg*ntotal))
   allocate(len(0:nrow-1))
   allocate(loc(0:nrow-1))
!
   jpe0=mod(mype,ncol)
   if (myrow.eq.0) then
     do nr = 0,nrow-1
       jpe=jpe0+nr*ncol
       len(nr)=ntotal*lmlen(jpe)
     enddo
     loc(0)=0
     do n = 1,nrow-1
       loc(n)=loc(n-1)+len(n-1)
     enddo
   else
     len(myrow)=ntotal*lmp
   endif
!
   mk=0
   do k = 1,ntotal
     tmpsnd(mk+1:mk+lmp)=par(1:lmp,k)
     mk=mk+lmp
   enddo
!
   call mpi_gatherv(tmpsnd,len(myrow),real_type,                               &
        tmprcv,len,loc,real_type,0,comm_col,ierr)
   if (ierr.ne.0) then
     if (mype.eq.master) write(6,*)'mpsp2f.F:mpi_gatherv failed'
     call mpabort
   endif
!
   if( myrow.eq.0) then
     mk=0
     do nr = 0,nrow-1
       jpe=jpe0+nr*ncol
       ms=lmstr(jpe)
       lmm=lmlen(jpe)/2
       lm=lmm
       if (ms.eq.0) lm=lm+1
       me=ms+lm-1
       msm=-me
       mem=msm+lmm-1
       do k = 1,ntotal
         tot(msm:mem,k)=tmprcv(mk+1:mk+lmm)
         mk=mk+lmm
         tot(ms :me ,k)=tmprcv(mk+1:mk+lm )
         mk=mk+lm
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
   end subroutine mpm2f
!-------------------------------------------------------------------------------
