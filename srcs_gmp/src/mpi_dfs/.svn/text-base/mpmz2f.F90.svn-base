#include "define.h"
   subroutine mpmz2f(par,mt,jlm,lota,tot,mtg,jlg,levs,lotb)
!------------------------------------------------------------------------------
!
! Collect (mt,jlg,levsp) data from all PE to (mtg,jlg,levs)
!
! usage:        call mpmz2f(par,mt,jlm,lota,tot,mtg,jlg,levs,lotb)
!
!    input argument lists:
!       par     - real (mt,jlm,lota) sub partial field
!       mt      - integer local spectral grid
!       jlm     - integer second dimension of par
!       lota    - integer total set of par fields
!       mtg     - integer total zonal spectral grid
!       jlg     - integer total meridional wave
!       levs   - integer total number of level
!       lota    - integer total set of tot fields
!
!    output argument list:
!       tot     - real (mtg,jlg,lotb) total field
!
!------------------------------------------------------------------------------
   use commpi, only  :  npes,mype,comm_world,real_type,master,                 &
                        lmlen,lmstr,ncol,levstr,levlen
   use dfsvar, only  :  mgs,mge
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer  ::  mt,jlm,jlg,mtg,levs,lota,lotb
   real  ::   par(mt,jlm,lota),tot(mgs:mge,jlg,lotb)
!
   real(_mpi_real_),allocatable  ::  tmpsnd(:),tmprcv(:)
   integer  ::  np,ms,me,lm,msm,mem,lmm,n,mk,ierr,ks,ke,k
   integer,allocatable  ::  len(:),loc(:)
   integer  ::  nvar,nsin,ksing,nv
!
   allocate(tmpsnd(mt*jlg*lota))
   allocate(tmprcv(mtg*jlg*lotb))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   nvar=lotb/levs
   nsin=lotb-nvar*levs
   ksing=levs*nvar+1

   if( mype.eq.master ) then
     do n = 0,npes-1
       len(n)=lota*lmlen(n)*jlg
     enddo
     loc(0)=0
     do n = 1,npes-1
       loc(n)=loc(n-1)+len(n-1)
     enddo
   else
     len(mype)=lota*mt*jlg
   endif
!
   mk=0
   do k = 1,lota
     lm=jlg*mt
     tmpsnd(mk+1:mk+lm)=par(1:lm,1,k)
     mk=mk+lm
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
       do nv = 1,nvar
         ks=(nv-1)*levs+levstr(np)
         ke=ks+levlen(np)-1
         do k = ks,ke
           do n = 1,jlg
             tot(msm:mem,n,k)=tmprcv(mk+1:mk+lmm)
             mk=mk+lmm
             tot(ms :me ,n,k)=tmprcv(mk+1:mk+lm )
             mk=mk+lm
           enddo
         enddo
       enddo
       do nv = 1,nsin
         k=nvar*levs+nv
         do n = 1,jlg
           tot(msm:mem,n,k)=tmprcv(mk+1:mk+lmm)
           mk=mk+lmm
           tot(ms :me ,n,k)=tmprcv(mk+1:mk+lm )
           mk=mk+lm
         enddo
       enddo
     enddo   ! np
   endif      ! master
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
   return
   end subroutine mpmz2f
!-------------------------------------------------------------------------------
