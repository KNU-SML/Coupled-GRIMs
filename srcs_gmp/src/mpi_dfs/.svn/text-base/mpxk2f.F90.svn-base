#include "define.h"
   subroutine mpxk2f(a,ib,jbwa,levsp,nvar,b,iba,levs)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram: 	mpxk2f
!            
! abstract:  transpose (ip,jf,kp) to (if,jf,kf)
!
! usage:	call mpxk2fa,ib,jbwa,levsp,nvar,b,iba,levs)
!
!    input argument lists:
!	a	- real (ib,jbwa,levsp,nvar) partial field for each PE
!	ib	- integer partial grid in longitude time 2
!	jbwa 	- integer total grid in latitude divide 2
!	levsp	   - integer partial number of model layer
!	nvar	   - integer total set of fields
!
!    output argument list:
!	b	- real (iba,jbwa,ntotal) full field 
!	iba	- integer total grid in longitude time 2
!	jbwa	- integer total grid in latitude divide 2
! 
! subprograms called:
!   mpi_gatherv	- gather message from all pe to master
!
!-------------------------------------------------------------------------------
   use commpi, only : mype,master,comm_world,npes,real_type,                   &
                      lonlen,lonstr,lerlen,lerstr
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ib,iba,jbwa,nvar,levs,levsp
   integer              ::  nv,n,i,j,k,mk,ierr,lonlend,lonstrd,is,ie,ks
   real                 ::  a(ib,jbwa,levsp,nvar),b(iba,jbwa,levs,nvar)
!
   real(_mpi_real_),allocatable  ::  tmpsnd(:),tmprcv(:)
   integer,allocatable           ::  len(:),loc(:)
!
   allocate(tmpsnd(ib*jbwa*levsp*nvar))
   allocate(tmprcv(iba *jbwa *levs*nvar))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then
     mk=0
     do n = 0,npes-1
       loc(n)=mk
       len(n)=lerlen(n)*nvar*lonlen(n)*jbwa
       mk=loc(n)+len(n)
     enddo
   endif
!
   lonlend=lonlen(mype)
   mk=0
   do nv = 1,nvar
     do k = 1,lerlen(mype)
       do j = 1,jbwa
         tmpsnd(mk+1:mk+lonlend)=a(1:lonlend,j,k,nv)
         mk=mk+lonlend
       enddo
     enddo
   enddo
   len(mype)=mk
!
   call mpi_gatherv(tmpsnd,len(mype),real_type,                                &
                    tmprcv,len(0),loc(0),real_type,0,comm_world,ierr)
!
   if( mype.eq.0 ) then
     mk=0
     do n = 0,npes-1
       lonlend=lonlen(n)
       lonstrd=lonstr(n)
       is=lonstrd
       ie=is+lonlend-1
       do nv = 1,nvar
         do k = 1,lerlen(n)
           ks=lerstr(n)+k-1
           do j = 1,jbwa
             b(is:ie,j,ks,nv)=tmprcv(mk+1:mk+lonlend)
             mk=mk+lonlend
           enddo
         enddo
       enddo
     enddo
   endif
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(len)
   deallocate(loc)
!
   return
   end subroutine mpxk2f
!-------------------------------------------------------------------------------
