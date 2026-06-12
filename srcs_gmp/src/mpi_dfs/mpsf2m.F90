#include "define.h"
   subroutine mpsf2m(tot,lmg,par,lmp,ntotal)
!-------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!
!!!!!not yet completed
!!!!!!!!!!!!!!!!!!!!!!!
!
! subprogram documentation block
!
! subprogram:   mpsf2m
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
                        lmlen,lmstr,ncol,nrow,comm_row
   use dfsvar, only  :  mgs,mge,jlg
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  lmg,lmp,ntotal
   real     ::  tot(mgs:mge,ntotal),par(lmp,ntotal)
!
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer  ::  ms,me,jpe0,jpe,msm,mem,lmm,n,k,mk,ierr
   integer  ::  nr,lm
   integer,allocatable  ::  len(:),loc(:)
!
   allocate(tmpsnd(lmg*ntotal))
   allocate(tmprcv(lmp*ntotal))
   allocate(len(0:nrow-1))
   allocate(loc(0:nrow-1))
!
   call mpi_bcast(tot,lmg*ntotal,real_type,0,comm_row,ierr)
!
   jpe0=mod(mype,ncol)
   if(myrow.eq.0) then
     mk=0
     do nr = 0,nrow-1
       loc(nr)=mk
       jpe=jpe0+nr*ncol
       ms=lmstr(jpe)
       lmm=lmlen(jpe)/2
       lm=lmm
       if(ms.eq.0) lm=lm+1
       me=ms+lm-1
       msm=-me
       mem=msm+lmm-1
       do k = 1,ntotal
         tmpsnd(mk+1:mk+lmm)=tot(msm:mem,k)
         mk=mk+lmm
         tmpsnd(mk+1:mk+lm) =tot(ms :me,k)
         mk=mk+lm
       enddo
       len(nr)=(lm+lmm)*ntotal
     enddo
   else
     len(myrow)=lmp*ntotal
   endif
!
   call mpi_scatterv(tmpsnd,len,loc,real_type,                                 &
                     tmprcv,len(myrow),real_type,0,comm_col,ierr)
!
   if (ierr.ne.0) then
     if (mype.eq.master) write(6,*)'mpsf2m.F:mpi_scatterv failed'
     call mpabort
   endif
!
   mk=0
   lm=lmlen(myrow)*ntotal
   par(1:lm,1)=tmprcv(mk+1:mk+lm)
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
   return
   end subroutine mpsf2m
!-------------------------------------------------------------------------------
