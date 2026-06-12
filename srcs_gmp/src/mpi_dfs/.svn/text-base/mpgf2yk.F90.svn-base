#include "define.h"
   subroutine mpgf2yk(a,iba,jbwa,levs,nvar,b,jbw,levsp)
!-------------------------------------------------------------------------------
!
!  subprogram documentation block
!
! subprogram: 	mpgf2yk
!            
! usage:	call mpgf2yk(a,iba,jbwa,levs,nvar,b,jbw,levsp)
!
!    input argument lists:
!	a	   - real (iba,jbwa,levs,nvar) full field 
!	iba	- integer total grid in longitude time 2
!	jbwa	- integer total grid in latitude divide 2
!	levs	- integer total number of model layers
!	nvar	- integer total set of fields
!
!    output argument list:
!	b	   - real (jbw,jbwa,levsp,nvar) partial field for each pe 
!	jbw   - integer partial grid in latitude time 2
!	levsp - integer partial number of model layer(divided by row PEs)
! 
! subprograms called:
!   mpi_scatterv	- scatter message from master pe to all pe
!
!-------------------------------------------------------------------------------
   use commpi, only : mype,master,comm_world,npes,real_type,                   &
                      latlen,latstr,levlen,levstr
   use dfsvar, only : jgs,latdef
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                       ::  iba,jbwa,jbw,levs,levsp,nvar
   integer                       ::  nv,n,mk,i,j,k,ierr,jj,ks
   real                          ::  a(iba,jbwa,levs,nvar)
   real                          ::  b(iba,jbw,levsp,nvar)
   real(_mpi_real_),allocatable  ::  tmpsnd(:),tmprcv(:)
   integer,allocatable           ::  len(:),loc(:)
!
   allocate(tmpsnd(iba*jbwa*levs *nvar))
   allocate(tmprcv(iba*jbw *levsp*nvar))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then
     mk=0
     do n = 0,npes-1
       loc(n)=mk
       do nv = 1,nvar
         do k = 1,levlen(n)
           ks=levstr(n)+k-1
#ifdef DBG
           write(6,'(A,3I8)')'ks,levlen,nv=',ks,levlen(n),nv
           call flush(6)
#endif
           do j = 1,latlen(n)
             jj=latdef(latstr(n)-jgs+j)
             tmpsnd(mk+1:mk+iba)=a(1:iba,jj,ks,nv)
             mk=mk+iba
             tmpsnd(mk+1:mk+iba)=a(1:iba,jbwa-jj+1,ks,nv)
             mk=mk+iba
           enddo
         enddo
       enddo
       len(n)=mk-loc(n)
     enddo
   else
     len(mype)=nvar*levlen(mype)*latlen(mype)*2*iba
   endif
!
   call mpi_scatterv(tmpsnd,len(0),loc(0),real_type,                           &
        tmprcv,len(mype),real_type,0,comm_world,ierr)
!
   mk=0
   do nv = 1,nvar
     do k = 1,levlen(mype)
       do j = 1,latlen(mype)
         jj=jbw-j+1
         b(1:iba,j,k,nv)=tmprcv(mk+1:mk+iba)
         mk=mk+iba
         b(1:iba,jj,k,nv)=tmprcv(mk+1:mk+iba)
         mk=mk+iba
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
   end subroutine mpgf2yk
!-------------------------------------------------------------------------------
