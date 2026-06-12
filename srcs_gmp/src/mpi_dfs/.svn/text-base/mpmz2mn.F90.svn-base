#include "define.h"
   subroutine mpmz2mn(a,jld,levsp,b,jl,mt,levs,nvar)
!------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:   mpnn2nk
!
! usage:        call mpnn2nk(a,jld,levsp,b,jl,mt,levs,nvar)
!
!    input argument lists:
!       a       - real (mt,jld,levsp,nvar) partial field
!       levsp   - integer partial vertical grid
!       jld      - second dimension of a
!       mt,jl   - integer partial spectral grid
!       levs    - integer full vertical grid
!       nvar    - integer total set of fields
!
!    output argument list:
!       b       - real (mt,jl,levs,nvar) partial field
!
! subprograms called:
!   mpi_alltoallv  - send and receive from all pe in the same comm
!
!------------------------------------------------------------------------------
   use commpi, only  :  npes,mype,master,ncol,nrow,real_type,                  &
                        comm_row,levlen,levstr,lnstr,lnlen
   use dfsvar, only  :  ngs,jlg
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer  ::  jl,levs,jld,mt,levsp,nvar
   real  ::  a(mt,jld,levsp,nvar),b(mt,jl,levs,nvar)
!
   integer  ::  np,ii,ipe,ipe0,n,k,m,mn,len,ierr,j,ll,kp
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer, allocatable  ::  lensnd(:),lenrcv(:)
   integer, allocatable  ::  locsnd(:),locrcv(:)
!
   kp=levlen(mype)
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do n = 1,nvar
       do k = 1,levs
         b(1:mt*jl,1,k,n)=a(1:mt*jl,1,k,n)
       enddo
     enddo
     return
   endif
!
! otherwise
!
   allocate(tmpsnd(mt*jlg*kp*nvar))
   allocate(tmprcv(mt*jl*levs*nvar))
   allocate(lensnd(ncol))
   allocate(lenrcv(ncol))
   allocate(locsnd(ncol))
   allocate(locrcv(ncol))
!
   ipe0=int(mype/ncol)*ncol
!
! cuting  k
!
   mn=0
   do ii = 1,ncol
     locsnd(ii)=mn
     ipe=(ii-1)+ipe0
     np=lnstr(ipe)-ngs
     do n = 1,nvar
       do k = 1,levlen(mype)
         do j = 1,lnlen(ipe)
           tmpsnd(mn+1:mn+mt)=a(1:mt,j+np,k,n)
           mn=mn+mt
         enddo
       enddo
     enddo
     lensnd(ii)=mn-locsnd(ii)
   enddo
!
   mn=0
   do ii = 1,ncol
     locrcv(ii)=mn
     ipe=(ii-1)+ipe0
     lenrcv(ii)=mt*jl*nvar*levlen(ipe)
     mn=mn+lenrcv(ii)
   enddo
!
   call mpi_alltoallv(tmpsnd,lensnd,locsnd,real_type,tmprcv,lenrcv,locrcv,     &
                      real_type,comm_row,ierr)
   if(ierr.ne.0) then
     print *,' error after alltoallv in mpmz2mn.F:',ierr
     call mpabort
   endif
!
! restore  n
!
   mn=0
   do ii = 1,ncol
     ipe=(ii-1)+ipe0
     do n = 1,nvar
       do k = 1,levlen(ipe)
         do j = 1,jl
           b(1:mt,j,k+levstr(ipe)-1,n)=tmprcv(mn+1:mn+mt)
           mn=mn+mt
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
   end subroutine mpmz2mn
!-------------------------------------------------------------------------------
