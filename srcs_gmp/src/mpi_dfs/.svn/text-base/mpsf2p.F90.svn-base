#include "define.h"
   subroutine mpsf2p(tot,lmg,lng,par,lmp,lnp,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! usage:   call mpsf2p(tot,lmg,lng,par,lmp,lnp,ntotal)
!
!    input argument lists:
!   tot   - real (lmg,lng,ntotal) total field
!   lmg   - global zonal wave number
!   lng   - global meridional wave number
!   ntotal   - integer total set of fields
!
!    output argument list:
!   par  - real (lmp,lnp,ntotal) sub partial field
! 
! subprograms called:
!   mpi_scatterv  - send to all pe in the same comm
!
!-------------------------------------------------------------------------------
   use commpi, only  :  npes,mype,comm_world,master,real_type,lnlen,lnstr,     &
                        lmlen,lmstr,ncol
   use dfsvar, only  :  mgs,mge,ngs,nge
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  lmg,lng,lmp,lnp,ntotal
   real  ::  tot(mgs:mge,ngs:nge,ntotal),par(lmp,lnp,ntotal)
!
   integer  ::  nend,np,n,m,k,mk,ms,me,msm,mem,lmm,ierr,lm
   real(_mpi_real_),allocatable  ::  tmpsnd(:),tmprcv(:)
   integer,allocatable  ::  len(:),loc(:)
!
   allocate(tmpsnd(lmg*lng*ntotal))
   allocate(tmprcv(lmp*lnp*ntotal))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then
     mk=0
     do np = 0,npes-1
       loc(np)=mk
       ms=lmstr(np)
       lmm=lmlen(np)/2
       lm=lmm
       if (ms.eq.0) lm=lm+1
       me=ms+lm-1
       msm=-me
       mem=msm+lmm-1
       nend=lnstr(np)+lnlen(np)-1
       do k = 1,ntotal
         do n = lnstr(np),nend
           tmpsnd(mk+1:mk+lmm)=tot(msm:mem,n,k)
           mk=mk+lmm
           tmpsnd(mk+1:mk+lm )=tot(ms :me ,n,k)
           mk=mk+lm
         enddo
       enddo
       len(np)=ntotal*(lm+lmm)*lnlen(np)
     enddo
   else
     len(mype)=ntotal*lmp*lnp
   endif
!
   call mpi_scatterv(tmpsnd,len,loc,real_type,                                 &
                     tmprcv,len(mype),real_type,0,comm_world,ierr)
!
   mk=0
   do k = 1,ntotal
     do n = 1,lnp
       par(1:lmp,n,k)=tmprcv(mk+1:mk+lmp)
       mk=mk+lmp
     enddo
   enddo
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
   return
   end subroutine mpsf2p
!-------------------------------------------------------------------------------
