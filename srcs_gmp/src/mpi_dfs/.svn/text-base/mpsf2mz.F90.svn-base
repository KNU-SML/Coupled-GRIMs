#include "define.h"
   subroutine mpsf2mz(tot,mtg,jlg,levs,lota,par,mt,jlm,lotb)
!------------------------------------------------------------------------------
!
! Distribute total fileds (mtg,jlg,lota) to (mt,jlg,lotb) on each PE !
!
! usage:   call mpsf2mz(tot,mtg,jlg,levs,lota,par,mt,jlm,lotb)
!
!    input argument lists:
!   tot   - real (mtg,jlg,lota) total field
!   mtg   - global zonal wave number
!   jlg   - global meridional wave number
!   lotb,lota   - integer total set of fields
!
!    output argument list:
!   par  - real (mt,jlg,levsp,nvar) sub partial field
! 
! subprograms called:
!   mpi_scatterv  - send to all pe in the same comm
!
!------------------------------------------------------------------------------
   use commpi, only  :  npes,mype,comm_world,master,real_type,                 &
                        lmlen,lmstr,ncol,levstr,levlen
   use dfsvar, only  :  mgs,mge
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer  ::  mtg,jlg,jlm,levs,lota,mt,lotb
   real  ::  tot(mgs:mge,jlg,lota),par(mt,jlm,lotb)
!
! local
!
   integer  ::  mk,n,k,kk,ierr,lm,lmm,ms,me,msm,mem,np,nv
   integer  ::  nsin,nvar,ksing
   real(_mpi_real_),allocatable  ::  tmpsnd(:),tmprcv(:)
   integer,allocatable  ::  len(:),loc(:)
!
   allocate(tmpsnd(mtg*jlg*lota))
   allocate(tmprcv(mt*jlg*lotb))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   nvar=lota/levs
   nsin=lota-nvar*levs
   ksing=levs*nvar+1
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
       do nv = 1,nvar
         do k = levstr(np),levstr(np)+levlen(np)-1
           kk=(nv-1)*levs+k
           do n = 1,jlg
             tmpsnd(mk+1:mk+lmm)=tot(msm:mem,n,kk)
             mk=mk+lmm
             tmpsnd(mk+1:mk+lm )=tot(ms :me ,n,kk)
             mk=mk+lm
           enddo
         enddo
       enddo
       do nv = 1,nsin
         kk=ksing+nv-1
         do n = 1,jlg
           tmpsnd(mk+1:mk+lmm)=tot(msm:mem,n,kk)
           mk=mk+lmm
           tmpsnd(mk+1:mk+lm )=tot(ms :me ,n,kk)
           mk=mk+lm
         enddo
       enddo
       len(np)=((lm+lmm)*jlg)*(nsin+nvar*levlen(np))
     enddo
   else
     len(mype)=lotb*mt*jlg
   endif
!
   call mpi_scatterv(tmpsnd,len,loc,real_type,                                 &
                     tmprcv,len(mype),real_type,0,comm_world,ierr)
!
   mk=0
   do k = 1,lotb
     lm=mt*jlg
     par(1:lm,1,k)=tmprcv(mk+1:mk+lm)
     mk=mk+lm
   enddo
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
   return
   end subroutine mpsf2mz
!-------------------------------------------------------------------------------
