!
   subroutine mpnbnd(a,mt,jl,levs,bl,br,n)
!------------------------------------------------------------------------------
!
! Purpose :  exchange meridional wave coeff. at boundary with
!            adjacent PE.
! Output    : bl  : from left PE
!             br  : from right PE
! Input      
!         n : number of boundary
!
!------------------------------------------------------------------------------
   use commpi, only : ncol,mycol,comm_row,real_type,istat_size
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer  ::  mt,jl,levs,n,jlb
   real  ::  a(mt,jl,levs),bl(mt,n,levs),br(mt,n,levs)
!
!local variables
!
   integer  ::  nc,k,ie,is,icnt,ndes,ierr
   integer  ::  istat(istat_size),isreqf,irreql,isreql,irreqf,ifirst,ilast
   real, allocatable, dimension(:)  ::  sndf,rcvf,sndl,rcvl
!
   if (mycol.ne.0) then
!
! send first n
!
     allocate(sndf(mt*levs*n))
     !allocate(rcvl(mt*levs*n))
     icnt=1
     do k = 1,levs
       do nc = 1,n
         sndf(icnt:icnt+mt-1)=a(1:mt,nc,k)
         icnt=icnt+mt
       enddo
     enddo
     icnt=icnt-1
     call mpi_isend(sndf,icnt,real_type,mycol-1,1,comm_row,isreqf,ierr)
     if (ierr.ne.0) call mpabort
     !call mpi_irecv(rcvl,icnt,real_type,mycol-1,ilast,
     call mpi_irecv(bl,icnt,real_type,mycol-1,1,comm_row,irreql,ierr)
     if (ierr.ne.0) call mpabort
     ndes=jl+n+1
   else
     ndes=jl+n
   endif
!
   if (mycol.ne.ncol-1) then
!
! send last n
!
     allocate(sndl(mt*levs*n))
     !allocate(rcvf(mt*levs*n))
     icnt=1
     do k = 1,levs
       do nc = n,1,-1
         sndl(icnt:icnt+mt-1)=a(1:mt,jl-nc+1,k)
         icnt=icnt+mt
       enddo
     enddo
     icnt=icnt-1
     call mpi_isend(sndl,icnt,real_type,mycol+1,1,comm_row,isreql,ierr)
     if (ierr.ne.0) call mpabort
     !call mpi_irecv(rcvf,icnt,real_type,mycol+1,
     call mpi_irecv(br,icnt,real_type,mycol+1,1,comm_row,irreqf,ierr)
     if (ierr.ne.0) call mpabort
   endif
!
   if (mycol.ne.0) then
     call MPI_WAIT(irreql,istat,ierr)
     if (ierr.ne.0) call mpabort
     !is=1
     !do k = 1,levs
     !  do nc = n,1,-1
     !      ie=is+mt-1
     !      bl(1:mt,nc,k)=rcvl(is:ie)
     !      is=ie+1
     !  enddo
     !enddo
     !deallocate(rcvl)
     call MPI_WAIT(isreqf,istat,ierr)
     if (ierr.ne.0) call mpabort
     deallocate(sndf)
   endif
!
   if (mycol.ne.ncol-1) then
     call MPI_WAIT(irreqf,istat,ierr)
     if (ierr.ne.0) call mpabort
     !is=1
     !do k = 1,levs
     !  do nc = 1,n
     !    br(1:mt,nc,k)=rcvf(is:ie)
     !    is=ie+1
     !  enddo
     !enddo
     !deallocate(rcvf)
     call MPI_WAIT(isreql,istat,ierr)
     if (ierr.ne.0) call mpabort
     deallocate(sndl)
   endif
!
   return
   end subroutine mpnbnd
!-------------------------------------------------------------------------------
