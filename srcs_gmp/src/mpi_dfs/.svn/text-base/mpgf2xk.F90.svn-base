#include "define.h"
   subroutine mpgf2xk(a,iba,jbwa,levs,nvar,b,ib,levsp)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram: 	mpgf2xk
!            
! usage:	call mpgf2xk(a,iba,jbwa,levs,nvar,b,ib,levsp)
!
!    input argument lists:
!	a	- real (iba,jbwa,levs,nvar) full field 
!	iba	- integer total grid in longitude time 2
!	jbwa	- integer total grid in latitude divide 2
!	levs	- integer total number of model layers
!	nvar	- integer total set of fields
!
!    output argument list:
!	b	- real (ib,jbwa,levsp,nvar) partial field for each pe 
!	ib	- integer partial grid in longitude time 2
!	levsp    - integer partial number of model layer(divided by row PEs)
! 
! subprograms called:
!   mpi_scatterv	- scatter message from master pe to all pe
!
!-------------------------------------------------------------------------------
   use commpi, only : mype,master,comm_world,npes,real_type,                   &
                      lonlen,lerstr,lerlen,lonstr
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                       ::  iba,jbwa,ib,levs,levsp,nvar
   integer                       ::  nv,n,mk,i,j,k,ierr,lonlend,lonstrd
   integer                       ::  is,ie,ks
   real                          ::  a(iba,jbwa,levs,nvar)
   real                          ::  b(ib,jbwa,levsp,nvar)
   real(_mpi_real_),allocatable  ::  tmpsnd(:),tmprcv(:)
   integer,allocatable           ::  len(:),loc(:)
!
   allocate(tmpsnd(iba *jbwa *levs*nvar))
   allocate(tmprcv(ib*jbwa*levsp*nvar))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then
     mk=0
     do n = 0,npes-1
       loc(n)=mk
       lonlend=lonlen(n)
       lonstrd=lonstr(n)
       is=lonstrd
       ie=is+lonlend-1
#ifdef DBG
       write(6,'(A,5I8)')'n,is,ie=',n,is,ie
       call flush(6)
#endif
       do nv = 1,nvar
         do k = 1,lerlen(n)
           ks=lerstr(n)+k-1
#ifdef DBG
           write(6,'(A,3I8)')'ks,lerlen,nv=',ks,lerlen(n),nv
           call flush(6)
#endif
           do j = 1,jbwa
             tmpsnd(mk+1:mk+lonlend)=a(is:ie,j,ks,nv)
             mk=mk+lonlend
           enddo
         enddo
       enddo
       len(n)=mk-loc(n)
     enddo
   else
     len(mype)=nvar*lerlen(mype)*lonlen(mype)*jbwa
   endif
!
   call mpi_scatterv(tmpsnd,len(0),loc(0),real_type,                           &
        tmprcv,len(mype),real_type,0,comm_world,ierr)
!
   mk=0
   lonlend=lonlen(mype)
   do nv = 1,nvar
     do k = 1,lerlen(mype)
       do j = 1,jbwa
         b(1:lonlend,j,k,nv)=tmprcv(mk+1:mk+lonlend)
         mk=mk+lonlend
       enddo
     enddo
   enddo
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
   return
   end subroutine mpgf2xk
!-------------------------------------------------------------------------------
