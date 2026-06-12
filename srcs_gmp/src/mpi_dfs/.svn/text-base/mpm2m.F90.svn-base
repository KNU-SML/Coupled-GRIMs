#include "define.h"
   subroutine mpm2m(a,jl,b,jtg,mt,nlev)
!-------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------
   use commpi, only  :  mycol,master,ncol,nrow,real_type,comm_col,             &
                        levlen,levstr,lnstr,lnlen
   use dfsvar, only  :  ngs
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  jl,nlev,jtg,mt,                                                &
                n,k,m,mn,len,ierr,j
   real  ::  a(mt,jl,nlev),b(mt,ngs:ngs+jtg-1,nlev)
!
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer, allocatable  ::  lensnd(:),locsnd(:)
   integer  ::  lenrcv
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do k = 1,nlev
       do m = 1,mt
         b(m,jtg+ngs-1,k)=a(m,1,k)
       enddo
     enddo
     return
   endif
!
! otherwise
!
   allocate(tmprcv(mt*nlev))
   if (mycol.eq.0) tmprcv(:)=a(1:mt*nlev,1,1)
!
   call mpi_bcast(tmprcv,mt*nlev,real_type,0,comm_col,ierr)

   if(ierr.ne.0) then
     print *,' error after mpi_scatterv in mpm2mz.F:',ierr
     call mpabort
   endif
!
! restore  n
!
   mn=0
   do n = 1,nlev
     b(1:mt,ngs+jtg-1,n)=tmprcv(mn+1:mn+mt)
     mn=mn+mt
   enddo
!
   deallocate(tmprcv)
!
   return
   end subroutine mpm2m
!-------------------------------------------------------------------------------
