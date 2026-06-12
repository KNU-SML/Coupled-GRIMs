#include "define.h"
   subroutine mpyz2xy(a,lonf2,levsp,b,lonf2p,levs,latg2p,nvar)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpnk2nx
!            
! abstract:  transpose (if,kp) to (ip,kf)
!
! usage:   call mpyz2xy(a,lonf2,levsp,b,lonf2p,levs,latg2p,nvar)
!
!    input argument lists:
!   a   - real (lonf2,latg2p,levsp,nvar) partial field in j k
!   lonf2   - integer total longitude grid
!   lonf2p   - integer partial longitude grid
!   latg2p   - integer partial latitude grid
!   levsp   - integer partial vertical grid 
!   levs   - integer total vertical grid 
!   levsp   - integer partial vertical grid 
!   nvar   - integer total set of fields
!
!    output argument list:
!   b   - real (lonf2p,latg2p,levs,nvar) partial field in i j
! 
! subprograms called:
!   mpi_alltoallv  - send and receive from all pe in the same comm
!
!-------------------------------------------------------------------------------
   use commpi, only  :  mype,master,ncol,real_type,comm_row,                   &
                        latlen,lonlen,levlen,levstr,lonstr
   use dfsvar, only  :  igs,ige
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  lonf2,levsp,latg2p,lonf2p,levs,nvar
   real     ::  a(igs:ige,latg2p,levsp,nvar),b(lonf2p,latg2p,levs,nvar)
!
   integer  ::  i,j,k,ii,ipe,ipe0,n,m,mn,len,ierr,ks,ke,ioff,lonfp,ip,jp,kp
!
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer, allocatable  ::  lensnd(:),lenrcv(:)
   integer, allocatable  ::  locsnd(:),locrcv(:)
!
   ip=lonlen(mype)
   jp=latlen(mype)*2
   kp=levlen(mype)
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do n = 1,nvar
       do k = 1,levs
         do j = 1,jp
           do i = 1,ip
             b(i,j,k,n) = a(i+igs-1,j,k,n)
           enddo
         enddo
       enddo
     enddo
     return
   endif
!
! otherwise
!
   allocate(tmpsnd(lonf2*jp*kp*nvar))
   allocate(tmprcv(ip*jp*levs*nvar))
   allocate(lensnd(ncol))
   allocate(lenrcv(ncol))
   allocate(locsnd(ncol))
   allocate(locrcv(ncol))
!
   ipe0=int(mype/ncol)*ncol
!
! cut in x
!
   mn=0
   do ii = 1,ncol
     locsnd(ii)=mn
     ipe=(ii-1)+ipe0
     lonfp=lonlen(ipe)
     ioff=lonstr(ipe)-1
     do n = 1,nvar
       do k = 1,kp
         do j = 1,jp
           tmpsnd(mn+1:mn+lonfp)=a(1+ioff:lonfp+ioff,j,k,n)
           mn=mn+lonfp
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
     lenrcv(ii)=levlen(ipe)*ip*jp*nvar
     mn=mn+lenrcv(ii)
   enddo
!
   call mpi_alltoallv(tmpsnd,lensnd,locsnd,real_type,tmprcv,lenrcv,locrcv,     &
                      real_type,comm_row,ierr)
! restore to k
!
   mn=0
   do ii = 1,ncol
     ipe=(ii-1)+ipe0
     ks=levstr(ipe)
     ke=ks+levlen(ipe)-1
     do n = 1,nvar
       do k = ks,ke
         do j = 1,jp
           b(1:ip,j,k,n)=tmprcv(mn+1:mn+ip)
           mn=mn+ip
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
   end subroutine mpyz2xy
!-------------------------------------------------------------------------------
