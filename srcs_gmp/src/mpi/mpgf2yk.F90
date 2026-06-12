#include <define.h>
   subroutine mpgf2yk(a,lonf,latg,levs,nvar,b,latgp,levsp)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram: 	mpgf2yk
!            
! program history log:
!    07-09-15  hoon park
!
! usage:	call mpgf2yk(a,lonf2,latg2,levs,nvar,b,latg2p,levsp)
!
!    input argument lists:
!	a	- real (lonf2,latg2,levs,nvar) full field 
!	lonf2	- integer total grid in longitude time 2
!	latg2	- integer total grid in latitude divide 2
!	levs	- integer total number of model layers
!	nvar	- integer total set of fields
!
!    output argument list:
!	b	- real (latg2p,latg2,levsp,nvar) partial field for each pe 
!	latg2p	- integer partial grid in latitude time 2
!	levsp    - integer partial number of model layer(divided by row PEs)
! 
! subprograms called:
!   mpi_scatterv	- scatter message from master pe to all pe
!
!-------------------------------------------------------------------------------
   use commpi       ! npes, mype, master, levlen, levstr, latstr, latlen,      &
                    ! MPI_REAL, mpi_comm_world
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                       ::  lonf,latg,latgp,levs,levsp,nvar
   integer                       ::  nv,n,mk,i,j,k,ierr
   integer                       ::  ks,jj
   real                          ::  a(lonf,latg,levs,nvar)
   real                          ::  b(lonf,latgp,levsp,nvar)
   real(_mpi_real_),allocatable  ::  tmpsnd(:),tmprcv(:)
   integer,allocatable           ::  len(:),loc(:)
!
   allocate(tmpsnd(lonf *latg *levs*nvar))
   allocate(tmprcv(lonf*latgp*levsp*nvar))
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
           do j = latstr(n),latstr(n)+latlen(n)-1
             tmpsnd(mk+1:mk+lonf)=a(1:lonf,j,ks,nv)
             mk=mk+lonf
             jj=latg-j+1
             tmpsnd(mk+1:mk+lonf)=a(1:lonf,jj,ks,nv)
             mk=mk+lonf
           enddo
         enddo
       enddo
       len(n)=mk-loc(n)
     enddo
   else
     len(mype)=nvar*levlen(mype)*latlen(mype)*lonf*2
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
   do nv = 1,nvar
     do k = 1,levlen(mype)
       do j = 1,latlen(mype)
         b(1:lonf,j,k,nv)=tmprcv(mk+1:mk+lonf)
         mk=mk+lonf
         jj=latgp-j+1
         b(1:lonf,jj,k,nv)=tmprcv(mk+1:mk+lonf)
         mk=mk+lonf
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
