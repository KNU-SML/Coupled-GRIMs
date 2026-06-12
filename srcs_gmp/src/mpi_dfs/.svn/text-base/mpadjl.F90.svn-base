#include "define.h"
   subroutine mpadjl(par,lmp,lnp,tot,lmg,lng,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! add aditional meridional wave components to tot after JLG index
!
! usage:        call mpadjl(par,lmp,tot,lmg,lng,ntotal)
!
!    input argument lists:
!       par     - real (lmp,lnp,ntotal) sub partial field
!       lmp     - integer total spectral grid
!       lnp     - integer total spectral grid
!       ntotal  - integer total set of fields
!
!    output argument list:
!       tot       - real (lmg,lng,ntotal) total field
!
! subprograms called:
!   shufmg      - shuffle zonal wave after gather
!   mpi_gatherv  - gather all pe in the same comm
!
!-------------------------------------------------------------------------------
   use commpi, only  :  nrow,ncol,myrow,mype,comm_col,real_type,master,        &
                        lmlen,lmstr
   use dfsvar, only  :  mgs,mge,jlg
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  lmg,lmp,lng,lnp,ntotal
   real     ::  par(lmp,lnp,ntotal),tot(mgs:mge,lng,ntotal)
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer  ::  np,ms,me,lm,msm,mem,lmm,n,k,mk,ierr,nn
   integer, allocatable  ::  len(:),loc(:)
!
   allocate(tmpsnd(lmp*lnp*ntotal))
   allocate(tmprcv(lmg*lnp*ntotal))
   allocate(len(0:nrow-1))
   allocate(loc(0:nrow-1))
!
   if( myrow.eq.0 ) then
     do nn = 0,nrow-1
       np=nn*ncol
       len(nn)=ntotal*lmlen(np)*lnp
     enddo
     loc(0)=0
     do nn = 1,nrow-1
       loc(nn)=loc(nn-1)+len(nn-1)
     enddo
   else
     len(myrow)=ntotal*lmp*lnp
   endif
!
   mk=0
   do k = 1,ntotal
     do nn = 1,lnp
       tmpsnd(mk+1:mk+lmp)=par(1:lmp,nn,k)
       mk=mk+lmp
     enddo
   enddo
!
   call mpi_gatherv(tmpsnd,len(myrow),real_type,tmprcv,len,loc,real_type,0,    &
                    comm_col,ierr)
   if (ierr.ne.0) then
     if (mype.eq.master) write(96,*)'mpadjl.F:mpi_gatherv failed'
     call mpabort
   endif
!
   if (mype.eq.master) then
     mk=0
     do nn = 0,nrow-1
       np=nn*ncol
       ms=lmstr(np)
       lmm=lmlen(np)/2
       lm=lmm
       if (ms.eq.0) lm=lm+1
       me=ms+lm-1
       msm=-me
       mem=msm+lmm-1
       do k = 1,ntotal
         do n = 1,lnp
           np=jlg+n
           tot(msm:mem,np,k)=tmprcv(mk+1:mk+lmm)
           mk=mk+lmm
           tot(ms :me ,np,k)=tmprcv(mk+1:mk+lm )
           mk=mk+lm
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
   end subroutine mpadjl
!-------------------------------------------------------------------------------
