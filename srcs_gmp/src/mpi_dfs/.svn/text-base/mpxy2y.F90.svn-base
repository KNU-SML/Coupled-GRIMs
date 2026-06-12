#include "define.h"
   subroutine mpxy2y(a,lonp,b,long,latp,nvar)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpxy2y
!            
! usage:   call mpxy2y(a,lonp,b,long,latp,nvar)
!
!    input argument lists:
!   a   - real (lonp,latp,nvar) partial field in i j
!   lonp   - integer partial longitude grid
!   latp   - integer partial latitude grid
!   long   - integer total longitude grid
!   nvar   - integer total set of fields
!
!    output argument list:
!   b   - real (long,latp,nvar) partial field in j
! 
! subprograms called:
!   mpi_alltoallv  - send and receive from all pe in the same comm
!
!-------------------------------------------------------------------------------
   use dfsvar, only  :  igs,ige
   use commpi, only  :  mype,master,ncol,nrow,real_type,comm_row,lonlen,       &
                        latlen,levlen,levstr,lonstr
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                ::  lonp,long,latp,nvar
   real                   ::  a(lonp,latp,nvar),b(igs:ige,latp,nvar)
!
! local
!
   integer                ::  i,j,k,ii,ipe,ipe0,ioff,n,mn,len,ierr,            &
                              lonfp,latlon,ip,jp
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer               ::  lensnd
   integer, allocatable  ::  lenrcvs(:),locrcvs(:)
!
   ip=lonlen(mype)
   jp=latlen(mype)*2
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do k = 1,nvar
       do j = 1,jp
         do i = 1,long
           b(i+igs-1,j,k) = a(i,j,k)
         enddo
       enddo
     enddo
     return
   endif
!
! otherwise
!
   allocate(tmpsnd(ip*jp*nvar))
   allocate(tmprcv(long*jp*nvar))
   allocate(lenrcvs(ncol))
   allocate(locrcvs(ncol))
!
   lensnd=0
   do n = 1,nvar
     do j = 1,jp
       tmpsnd(lensnd+1:lensnd+ip)=a(1:ip,j,n)
       lensnd=lensnd+ip
     enddo
   enddo
!
   ipe0=(mype/ncol)*ncol
   mn=0
   do ii = 1,ncol
     locrcvs(ii)=mn
     ipe=(ii-1)+ipe0
     latlon=jp*lonlen(ipe)
     lenrcvs(ii)=latlon*nvar
     mn=mn+lenrcvs(ii)
   enddo
!
   call mpi_allgatherv(tmpsnd, lensnd, real_type,tmprcv,lenrcvs, locrcvs,      &
                       real_type, comm_row, ierr)
!
   mn=0
   do ii = 1,ncol
     ipe=(ii-1)+ipe0
     lonfp=lonlen(ipe)
     ioff=lonstr(ipe)-1
     do n = 1,nvar
       do j = 1,jp
         b(1+ioff:lonfp+ioff,j,n)=tmprcv(mn+1:mn+lonfp)
         mn=mn+lonfp
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
   end subroutine mpxy2y
!-------------------------------------------------------------------------------
