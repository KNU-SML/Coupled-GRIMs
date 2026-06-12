#include "define.h"
   subroutine mpyk2f(a,iba,jbw,levsp,nvar,b,jbwa,levs)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram: 	mpyk2f
!            
! abstract:  transpose (ip,jf,kp) to (if,jf,kf)
!
! usage:	call mpyk2f(a,iba,jbw,levsp,nvar,b,jbwa,levs)
!
!    input argument lists:
!	a	- real (iba,jbw,levsp,nvar) partial field for each PE
!	iba  	- integer total grid in longitude time 2
!	jbw 	- integer partial grid in latitude divide 2
!	levsp	   - integer partial number of model layer
!	nvar	   - integer total set of fields
!
!    output argument list:
!	b	- real (iba,jbwa,ntotal) full field 
!	jbwa	- integer total grid in latitude divide 2
! 
! subprograms called:
!   mpi_gatherv	- gather message from all pe to master
!
!-------------------------------------------------------------------------------
   use commpi, only : mype,master,comm_world,npes,real_type,                   &
                      latlen,latstr,levlen,levstr
   use dfsvar, only : jgs,latdef
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  iba,jbw,jbwa,nvar,levs,levsp
   integer              ::  nv,n,i,j,k,mk,ierr,jj,ks
   real                 ::  a(iba,jbw,levsp,nvar)
   real                 ::  b(iba,jbwa,levs,nvar)
!
   real(_mpi_real_),allocatable  ::  tmpsnd(:),tmprcv(:)
   integer,allocatable           ::  len(:),loc(:)
!
   allocate(tmpsnd(iba*jbw *levsp*nvar))
   allocate(tmprcv(iba*jbwa*levs *nvar))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then
     mk=0
     do n = 0,npes-1
       loc(n)=mk
       len(n)=levlen(n)*nvar*iba*latlen(n)*2
       mk=loc(n)+len(n)
     enddo
   endif
!
   mk=0
   do nv = 1,nvar
     do k = 1,levlen(mype)
       do j = 1,latlen(mype)
         jj=jbw-j+1
         tmpsnd(mk+1:mk+iba)=a(1:iba,j ,k,nv)
         mk=mk+iba
         tmpsnd(mk+1:mk+iba)=a(1:iba,jj,k,nv)
         mk=mk+iba
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
       do nv = 1,nvar
         do k = 1,levlen(n)
           ks=levstr(n)+k-1
           do j = 1,latlen(n)
             jj=latdef(latstr(n)-jgs+j)
             b(1:iba,jj,ks,nv)=tmprcv(mk+1:mk+iba)
             mk=mk+iba
             b(1:iba,jbwa-jj+1,ks,nv)=tmprcv(mk+1:mk+iba)
             mk=mk+iba
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
   end subroutine mpyk2f
!-------------------------------------------------------------------------------
