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
!   mpi_isend
!   mpi_irecv
!   mpi_barrier
!   mpi_wait
!
!------------------------------------------------------------------------------
   use commpi, only : mype,ncol,real_type,levlen,levstr,comm_world,mycol     ,&
                      istat_size
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer                           ::  levsp,levs,mt,jlg
   real   , dimension(mt,jlg,levsp)  ::  a
   real   , dimension(mt,jlg,levs )  ::  b
!
! local
!
   integer, parameter             ::  tag=17
   integer                        ::  i,j,k,nc,ipe,ipe0,n,mn,size,ierr,ke,ks,kk
   integer                        ::  ilen,send_req,recv_req,status(istat_size)
   integer                        ::  dest,source,offset1,offset2
   real(_mpi_real_), allocatable  ::  tmpdata(:)
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
   allocate(tmpdata(mt*jlg*levs))
   tmpdata=0.
!
! input
!
   size=mt*jlg*levsp
   ilen=mycol*size
   do k = 1,levsp
     do j = 1,jlg
       tmpdata(ilen+1:ilen+mt)=a(1:mt,j,k)
       ilen=ilen+mt
     enddo
   enddo
!
! index
!
   ipe0=(mype/ncol)*ncol
   dest=ipe0+mod(mycol+1,ncol)
   source=ipe0+mod(mycol-1+ncol,ncol)
   offset1=mod(mycol+ncol,ncol)*size+1
   offset2=mod(mycol+ncol-1,ncol)*size+1
!
! non-blocking send-receive
!
   call mpi_barrier(comm_world,ierr)
   do nc = 0,ncol-2
     call mpi_isend(tmpdata(offset1),size,real_type,dest                      ,&
                                     tag,comm_world,send_req,ierr)
     call mpi_irecv(tmpdata(offset2),size,real_type,source                    ,&
                                     tag,comm_world,recv_req,ierr)
     offset1=mod(mycol-(nc+1)+ncol,ncol)*size+1
     offset2=mod(mycol-(nc+1)+ncol-1,ncol)*size+1
     call mpi_wait(send_req,status,ierr)
     call mpi_wait(recv_req,status,ierr)
   enddo
!
! output
!
   do nc = 0,ncol-1
     ipe=ipe0+nc
     ks=levstr(ipe)
     ke=ks+levlen(ipe)-1
     mn=nc*size
     do k = ks,ke
       do j = 1,jlg
         b(1:mt,j,k)=tmpdata(mn+1:mn+mt)
         mn=mn+mt
       enddo
     enddo
   enddo
!
   deallocate(tmpdata)
!
   return
   end subroutine mpmz2m
!-------------------------------------------------------------------------------
