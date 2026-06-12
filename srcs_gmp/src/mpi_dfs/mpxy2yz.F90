#include "define.h"
   subroutine mpxy2yz(a,lonp,lota,b,long,lotb,latp,levs,levsp,nvar)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! usage:    call mpnx2nk(a,lonp,lota,b,long,lotb,latp,levs,levsp,nvar)
!
!    input argument lists:
!       a       - real (lonp,latp,lota) partial field in j k
!       long   - integer total longitude grid x 2
!       lota    - integer total field of nvar*levsp
!       lonp  - integer partial longitude grid x 2
!       lotb    - integer total field of nvar*levs
!       latp  - integer partial latitude grid / 2
!       levsp   - integer partial vertical grid
!       levs    - integer total vertical grid
!       nvar    - integer total set of fields
!
!    output argument list:
!       b       - real (long,atp,lotb) partial field in i j
!
! subprograms called:
!   mpi_alltoallv  - send and receive from all pe in the same comm
!
!-------------------------------------------------------------------------------
   use commpi, only  :  mype,master,ncol,nrow,real_type,comm_row,lonlen,       &
                        latlen,levlen,levstr,lonstr
   use dfsvar, only  :  igs,ige
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  long,levsp,latp,lonp,levs,nvar,lota,lotb
   real     ::  a(lonp,latp,lota),b(igs:ige,latp,lotb)
!
   integer  ::  ioff,koff,lonfp,latlon,kl,kf,ip,jp,kp,                         &
                i,j,k,ii,ipe,ipe0,n,mn,len,ierr
!
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer, allocatable  ::  lensnd(:),lenrcv(:)
   integer, allocatable  ::  locsnd(:),locrcv(:)
!
   ip=lonlen(mype)
   jp=latlen(mype)*2
   kl=levlen(mype)
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do k = 1,nvar*levs
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
   allocate(tmpsnd(ip*jp*levs*nvar))
   allocate(tmprcv(long*jp*kl*nvar))
   allocate(lensnd(ncol))
   allocate(lenrcv(ncol))
   allocate(locsnd(ncol))
   allocate(locrcv(ncol))
!
   ipe0=int(mype/ncol)*ncol
!
! cut in k
!
   mn=0
   lonfp=lonlen(mype)
   do ii = 1,ncol
     locsnd(ii)=mn
     ipe=(ii-1)+ipe0
     koff=levstr(ipe)-1
     do n = 1,nvar
       kf=(n-1)*levs+koff
       do k = 1,levlen(ipe)
         do j = 1,jp
           tmpsnd(mn+1:mn+ip)=a(1:ip,j,k+kf)
           mn=mn+ip
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
     latlon=latlen(mype)*2*lonlen(ipe)
     lenrcv(ii)=levlen(mype)*latlon*nvar
     mn=mn+lenrcv(ii)
   enddo
!
   call mpi_alltoallv(tmpsnd,lensnd,locsnd,real_type,tmprcv,lenrcv,locrcv,     &
                      real_type,comm_row,ierr)
!
! put to x
!
   mn=0
   do ii = 1,ncol
     ipe=(ii-1)+ipe0
     lonfp=lonlen(ipe)
     ioff=lonstr(ipe)-1
     do n = 1,nvar
       kp=(n-1)*levsp
       do k = 1,kl
         do j = 1,jp
           b(1+ioff:lonfp+ioff,j,k+kp)=tmprcv(mn+1:mn+lonfp)
           mn=mn+lonfp
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
   end subroutine mpxy2yz
!-------------------------------------------------------------------------------
