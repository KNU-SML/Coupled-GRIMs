#include "define.h"
   subroutine mpmn2m(a,jl,b,jlg,mt,nvar)
!------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpmn2m
!            
! usage:   call mpmn2m(a,jl,b,jlg,mt,nvar)
!
!    input argument lists:
!   a   - real (mt,jl,nvar) partial field in m,n
!   mt   - integer partial zonal wave
!   jl   - integer partial meridional wave
!   jlg   - integer total meridional wave
!   nvar- integer total number of variables
!
!    output argument list:
!   b   - real (mt,jlg,nvar) partial field in zonal
! 
! subprograms called:
!   mpi_allgatherv
!
!------------------------------------------------------------------------------
   use dfsvar, only  :  ngs
   use commpi, only  :  mype,ncol,real_type,comm_row,lnlen,lnstr
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer  ::  jl,jlg,mt,nvar
   real  ::  a(mt,jl,nvar),b(mt,jlg,nvar)
!
! local
!
   integer  ::  i,j,k,ii,ipe,ipe0,n,mn,len,ierr,ln,ne,ns
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer  ::  lensnd
   integer, allocatable  ::  lenrcvs(:),locrcvs(:)
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do k = 1,nvar
       do j = 1,jlg
         b(1:mt,j,k) = a(1:mt,j,k)
       enddo
     enddo
     return
   endif
!
! otherwise
!
   allocate(tmpsnd(jl*mt*nvar))
   allocate(tmprcv(jlg*mt*nvar))
   allocate(lenrcvs(ncol))
   allocate(locrcvs(ncol))
!
   lensnd=0
   do n = 1,nvar
     do j = 1,jl
       tmpsnd(lensnd+1:lensnd+mt)=a(1:mt,j,n)
         lensnd=lensnd+mt
     enddo
   enddo
!
   ipe0=(mype/ncol)*ncol
   mn=0
   do ii = 1,ncol
     locrcvs(ii)=mn
     ipe=(ii-1)+ipe0
     ln=mt*lnlen(ipe)
     lenrcvs(ii)=ln*nvar
     mn=mn+lenrcvs(ii)
   enddo
!
   call mpi_allgatherv(tmpsnd, lensnd, real_type,tmprcv,lenrcvs, locrcvs,      &
                       real_type,comm_row, ierr)
!
   mn=0
   do ii = 1,ncol
     ipe=(ii-1)+ipe0
     ns=lnstr(ipe)+(1-ngs)
     ne=ns+lnlen(ipe)-1
     do k = 1,nvar
       do n = ns,ne
         b(1:mt,n,k)=tmprcv(mn+1:mn+mt)
         mn=mn+mt
       enddo
     enddo
   enddo
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(lenrcvs)
   deallocate(locrcvs)
!
   return
   end subroutine mpmn2m
!-------------------------------------------------------------------------------
