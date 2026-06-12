#include "define.h"
   subroutine mpsp2f(par,lmp,lnp,tot,lmg,lng,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:   mpsp2f
!
! usage:        call mpsp2f(par,lmp,lnp,tor,lmg,lng,ntotal)
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
   use commpi, only  :  npes,mype,comm_world,real_type,master,lnlen,lnstr,     &
                        lmlen,lmstr,ncol
   use dfsvar, only  :  mgs,mge,ngs,jlg
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  lmg,lmp,lng,lnp,ntotal
   real  ::  par(lmp,lnp,ntotal),tot(mgs:mge,lng,ntotal)
!
   real(_mpi_real_),allocatable  ::  tmpsnd(:),tmprcv(:)
   integer  ::  np,ms,me,lm,msm,mem,lmm,n,k,mk,ierr,nn
   integer,allocatable  ::  len(:),loc(:)
!
   allocate(tmpsnd(lmp*lnp*ntotal))
   allocate(tmprcv(lmg*jlg*ntotal))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then
     do n = 0,npes-1
       len(n)=ntotal*lmlen(n)*lnlen(n)
     enddo
     loc(0)=0
     do n = 1,npes-1
       loc(n)=loc(n-1)+len(n-1)
     enddo
   else
     len(mype)=ntotal*lmp*lnp
   endif
!
   mk=0
   do k = 1,ntotal
     do n = 1,lnp
       tmpsnd(mk+1:mk+lmp)=par(1:lmp,n,k)
       mk=mk+lmp
     enddo
   enddo
!
   call mpi_gatherv(tmpsnd,len(mype),real_type,tmprcv,len,loc,real_type,0,     &
                    comm_world,ierr)
   if (ierr.ne.0) then
     if (mype.eq.master) write(96,*)'mpsp2f.F:mpi_gatherv failed'
     call mpabort
   endif
!
   if( mype.eq.master) then
     mk=0
     do np = 0,npes-1
       ms=lmstr(np)
       lmm=lmlen(np)/2
       lm=lmm
       if (ms.eq.0) lm=lm+1
       me=ms+lm-1
       msm=-me
       mem=msm+lmm-1
       do k = 1,ntotal
         do n = 1,lnlen(np)
           nn=lnstr(np)-ngs+n
           tot(msm:mem,nn,k)=tmprcv(mk+1:mk+lmm)
           mk=mk+lmm
           tot(ms :me ,nn,k)=tmprcv(mk+1:mk+lm )
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
   end subroutine mpsp2f
!-------------------------------------------------------------------------------
