!
  subroutine mpinit(strwtime)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpinit
!            
! abstract: initializing the mpi by each pe.
!
! usage:   call mpinit
!
! subprograms called:
!   mpi_comm_size   - to get total number of pe for the comm
!   mpi_comm_rank   - to get my own pe number 
!   mpi_comm_split   - to split entire comm into sub comm
!
!-------------------------------------------------------------------------------
   use commpi, only : npes,mype,master,ncol,nrow,comm_row,comm_col,istat_size, &
                     mpivarini,comm_world,real_type,int_type,mycol,myrow
!-------------------------------------------------------------------------------
   include 'mpif.h'
!
   real*8               ::  strwtime
!
   integer              ::  ncoli
   data ncoli/2/
   namelist /nammpi/ncoli
!
   call mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
   call mpi_comm_rank(MPI_COMM_WORLD,mype,jerr)
!
   if( ierr.ne.0 .or. jerr.ne.0) then
     print *,'PE',mype,':********* Error stop in mpinit ********* '
     print *,'PE',mype,':error code from mpi_comm_size = ',ierr
     print *,'PE',mype,':error code from mpi_comm_rank = ',jerr
     print *,'PE',mype,':******* End of output for mpinit ******* '
     call mpabort
   endif
!
   master=0
   comm_world=MPI_COMM_WORLD
#if (_mpi_real_ == 4)
   real_type=MPI_REAL
#else
   real_type=MPI_DOUBLE_PRECISION
#endif
   int_type=MPI_INTEGER
   istat_size=MPI_STATUS_SIZE
   call mpivarini(npes)
!
   if (mype.eq.master) read(95,nammpi)
   call mpbcasti(ncoli,1)
!
   strwtime=mpi_wtime()
   ncol=ncoli
   nrow=npes/ncol
   msgtag=0
   myrow=mype/ncol
   mycol=mod(mype,ncol)
   call mpi_comm_split(MPI_COMM_WORLD,myrow,mycol,comm_row,ierr)
   call mpi_comm_split(MPI_COMM_WORLD,mycol,myrow,comm_col,jerr)
!
   if( ierr.ne.0 .or. jerr.ne.0 ) then
     print *,'PE',mype,':********* Error stop in mpinit ********* '
     print *,'PE',mype,':error code for doing comm_row = ',ierr
     print *,'PE',mype,':error code for doing comm_col = ',jerr
     print *,'PE',mype,':******* End of output for mpinit ******* '
     call mpabort
   endif
!
   return
   end subroutine mpinit
!-------------------------------------------------------------------------------
