#include "define.h"
   subroutine mpm2mz(a,jl,levs,b,jtg,mt,levsp,nvar)
!------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:   mpnn2nk
!
! abstract:  transpose (jl,levs,nvar) to (jla,mt,levsp,nvar)
!
! usage:        call mpnn2nk(a,jl,levs,b,jla,mt,levsp,nvar)
!
!    input argument lists:
!       a       - real (jl,levs,nvar) partial field
!       jl    - integer sub partial spectral grid
!       levsp   - integer partial vertical grid
!       jla,mt   - integer partial spectral grid
!       levs    - integer full vertical grid
!       nvar    - integer total set of fields
!
!    output argument list:
!       b       - real (jla,mt,levsp,nvar) partial field
!
! subprograms called:
!   mpi_alltoallv  - send and receive from all pe in the same comm
!
!------------------------------------------------------------------------------
   use commpi, only  :  mype,master,ncol,nrow,real_type,comm_col,              &
                        levlen,levstr,lnstr,lnlen
   use dfsvar, only  :  ngs
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer  ::  jl,levs,jtg,mt,levsp,nvar,np,                                  & 
                ii,ipe,ipe0,n,k,m,mn,len,ierr,j,ll,kp                              
   real  ::  a(mt,jl,levs,nvar),b(mt,ngs:ngs+jtg-1,levsp,nvar)
!
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer, allocatable  ::  lensnd(:),locsnd(:)
   integer  ::  lenrcv
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do n = 1,nvar
       do k = 1,levs
         do m = 1,mt
           b(m,jtg+ngs-1,k,n)=a(m,1,k,n)
         enddo
       enddo
     enddo
     return
   endif
!
! otherwise
!
   kp=levlen(mype)
   allocate(tmpsnd(mt*levs*nvar))
   allocate(tmprcv(mt*kp*nvar))
   allocate(lensnd(ncol))
   allocate(locsnd(ncol))
!
   ipe0=int(mype/ncol)*ncol
!
! cuting  k
!
   mn=0
   ll=mt*jl
   do ii = 1,ncol
     locsnd(ii)=mn
     ipe=(ii-1)+ipe0
     do n = 1,nvar
       do k = 1,levlen(ipe)
         tmpsnd(mn+1:mn+ll)=a(1:ll,1,k+levstr(ipe)-1,n)
         mn=mn+ll
       enddo
     enddo
     lensnd(ii)=mn-locsnd(ii)
   enddo
!
   lenrcv=mt*nvar*levlen(mype)
!
   call mpi_scatterv(tmpsnd,lensnd,locsnd,real_type,                           &
                     tmprcv,lenrcv,       real_type,                           &
                     0,comm_col,ierr)
   if(ierr.ne.0) then
     print *,' error after mpi_scatterv in mpm2mz.F:',ierr
     call mpabort
   endif
!
! restore  n
!
   mn=0
   do n = 1,nvar
     do k = 1,kp
       b(1:mt,ngs+jtg-1,k,n)=tmprcv(mn+1:mn+mt)
       mn=mn+mt
     enddo
   enddo
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(lensnd)
   deallocate(locsnd)
!
   return
   end subroutine mpm2mz
!-------------------------------------------------------------------------------
