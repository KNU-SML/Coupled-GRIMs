#include <define.h>
   subroutine mpgf2xk(a,lonf,latg,levs,nvar,b,lonfp,levsp)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram: 	mpgf2xk
!            
! program history log:
!    07-09-15  hoon park
!
! usage:	call mpgf2xk(a,lonf2,latg2,levs,nvar,b,lonf2p,levsp)
!
!    input argument lists:
!	a	- real (lonf2,latg2,levs,nvar) full field 
!	lonf2	- integer total grid in longitude time 2
!	latg2	- integer total grid in latitude divide 2
!	levs	- integer total number of model layers
!	nvar	- integer total set of fields
!
!    output argument list:
!	b	- real (lonf2p,latg2,levsp,nvar) partial field for each pe 
!	lonf2p	- integer partial grid in longitude time 2
!	levsp    - integer partial number of model layer(divided by row PEs)
! 
! subprograms called:
!   mpi_scatterv	- scatter message from master pe to all pe
!
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_
   use commpi          ! npes, mype, lerlen, master, lonstr, lonlen,           &
                       ! MPI_REAL, mpi_comm_world, lerstr
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                       ::  lonf,latg,lonfp,levs,levsp,nvar
   integer                       ::  nv,n,mk,i,j,k,ierr,lonlend,lonstrd,lonff
   integer                       ::  is,ie,ks
   real                          ::  a(lonf,latg,levs,nvar)
   real                          ::  b(lonfp,latg,levsp,nvar)
   real(_mpi_real_),allocatable  ::  tmpsnd(:),tmprcv(:)
   integer,allocatable           ::  len(:),loc(:)
!
   allocate(tmpsnd(lonf *latg *levs*nvar))
   allocate(tmprcv(lonfp*latg*levsp*nvar))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if( mype.eq.master ) then
     mk=0
     do n = 0,npes-1
       loc(n)=mk
       lonlend=lonlen(n)
       lonstrd=lonstr(n)
       lonff=lonf_
       is=lonstrd
       ie=is+lonlend-1
       do nv = 1,nvar
         do k = 1,lerlen(n)
           ks=lerstr(n)+k-1
           do j = 1,latg
             tmpsnd(mk+1:mk+lonlend)=a(is:ie,j,ks,nv)
             mk=mk+lonlend
           enddo
         enddo
       enddo
       len(n)=mk-loc(n)
     enddo
   else
     len(mype)=nvar*lerlen(mype)*lonlen(mype)*latg
   endif
!
!soojin_couple
!   call mpi_scatterv(tmpsnd,len(0),loc(0),MPI_REAL,                            &
!        tmprcv,len(mype),MPI_REAL,0,mpi_comm_world,ierr)
   call mpi_scatterv(tmpsnd,len(0),loc(0),MPI_REAL,                            &
        tmprcv,len(mype),MPI_REAL,0,mpi_comm_private,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   mk=0
   lonlend=lonlen(mype)
   do nv = 1,nvar
     do k = 1,lerlen(mype)
       do j = 1,latg
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
