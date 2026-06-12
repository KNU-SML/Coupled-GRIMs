#include "define.h"
   subroutine mpmz2m(a,levsp,b,levs,mt,jlg)
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
   use dfsvar, only  :  kgs
   use commpi, only  :  mype,ncol,real_type,comm_row,levlen,levstr
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer       ::  levsp,levs,mt,jlg
   real          ::  a(mt,jlg,levsp),b(mt,jlg,levs)
!
! local
!
   integer  ::  i,j,k,ii,ipe,ipe0,n,mn,len,ierr,ln,ke,ks
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer  ::  lensnd
   integer, allocatable  ::  lenrcvs(:),locrcvs(:)
!
   forall(i=1:mt,j=1:jlg,k=1:levsp) b(i,j,k)=0.
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do k = 1,levsp
       do j = 1,jlg
         do i = 1,mt
           b(i,j,k) = a(i,j,k)
         enddo
       enddo
     enddo
     return
   endif
!
! otherwise
!
   allocate(tmpsnd(levsp*jlg*mt))
   allocate(tmprcv(levs*jlg*mt))
   allocate(lenrcvs(ncol))
   allocate(locrcvs(ncol))
!
   lensnd=0
   do n = 1,levsp
     do j = 1,jlg
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
     lenrcvs(ii)=mt*jlg*levlen(ipe)
     mn=mn+lenrcvs(ii)
   enddo
!
   call mpi_allgatherv(tmpsnd, lensnd, real_type,tmprcv,lenrcvs, locrcvs,      &
                       real_type,comm_row, ierr)
!
   mn=0
   do ii = 1,ncol
     ipe=(ii-1)+ipe0
     ks=levstr(ipe)
     ke=ks+levlen(ipe)-1
     do k = ks,ke
       do j = 1,jlg
         b(1:mt,j,k)=tmprcv(mn+1:mn+mt)
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
   end subroutine mpmz2m
!-------------------------------------------------------------------------------
