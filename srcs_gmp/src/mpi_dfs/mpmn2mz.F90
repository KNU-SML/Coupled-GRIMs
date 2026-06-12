#include "define.h"
   subroutine mpmn2mz(a,jt,levs,b,jlg,mt,levsp,nvar)
!------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:   mpnn2nk
!
! abstract:  transpose (jt,levs,nvar) to (jlg,mt,levsp,nvar)
!
! usage:        call mpnn2nk(a,jt,levs,b,jlg,mt,levsp,nvar)
!
!    input argument lists:
!       a       - real (jt,levs,nvar) partial field
!       jt    - integer sub partial spectral grid
!       levsp   - integer partial vertical grid
!       jlg,mt   - integer partial spectral grid
!       levs    - integer full vertical grid
!       nvar    - integer total set of fields
!
!    output argument list:
!       b       - real (jlg,mt,levsp,nvar) partial field
!
! subprograms called:
!   mpi_alltoallv  - send and receive from all pe in the same comm
!
!------------------------------------------------------------------------------
   use commpi, only  :  mype,ncol,real_type,comm_row,mycol,levlen,levstr,      &
                        lnlen,lnstr
   use dfsvar, only : ngs
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer  ::  jt,levs,jlg,mt,levsp,nvar
   real  ::  a(mt,jt,levs,nvar),b(mt,jlg,levsp,nvar)
!
   integer  ::  ii,ipe,ipe0,n,k,m,mn,len,ierr,j,ll,ns,ne,kp,nj
   integer  ::  jl
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer, allocatable  ::  lensnd(:),lenrcv(:)
   integer, allocatable  ::  locsnd(:),locrcv(:)
!
   kp=levlen(mype)
   jl=lnlen(mype)
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do n = 1,nvar
       do k = 1,levs
         do j = 1,jlg
           b(1:mt,j,k,n)=a(1:mt,j,k,n)
         enddo
       enddo
     enddo
     return
   endif
!
! otherwise
!
   allocate(tmpsnd(mt*jl*levs*nvar))
   allocate(tmprcv(mt*jlg*kp*nvar))
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
     do n = 1,nvar
       do k = 1,levlen(ipe)
         do j = 1,jl
           tmpsnd(mn+1:mn+mt)=a(1:mt,j,k+levstr(ipe)-1,n)
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
     lenrcv(ii)=mt*lnlen(ipe)*nvar*kp
     mn=mn+lenrcv(ii)
   enddo
!
   call mpi_alltoallv(tmpsnd,lensnd,locsnd,real_type,tmprcv,lenrcv,locrcv,     &
                      real_type,comm_row,ierr)
   if(ierr.ne.0) then
     print *,' error after alltoallv in mpmn2mz.F:',ierr
     call mpabort
   endif
!
! restore  n
!
   mn=0
   do ii = 1,ncol
     ipe=(ii-1)+ipe0
     ns=lnstr(ipe)+(1-ngs)
     ne=ns+lnlen(ipe)-1
     do n = 1,nvar
       do k = 1,kp
         do j = ns,ne
           b(1:mt,j,k,n)=tmprcv(mn+1:mn+mt)
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
   end subroutine mpmn2mz
!-------------------------------------------------------------------------------
