#include "define.h"
   subroutine mpgp2fd(par,ib,jbw,tot,iba,jbwa,ntotal)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpgp2f
!            
! usage:   call mpgp2f(a,ib,jbw,b,iba,jbwa,ntotal)
!
!    input argument lists:
!   b   - real (ib,jbw,ntotal) partial field for each pe 
!   ib   - integer partial grid in longitude time 2
!   jbw   - integer partial grid in latitude divide 2
!   ntotal   - integer total set of fields
!
!    output argument list:
!   a   - real (iba,jbwa,ntotal) full field 
!   iba   - integer total grid in longitude time 2
!   jbwa   - integer total grid in latitude divide 2
! 
! subprograms called:
!   mpi_gatherv   - gather message from all pe to master
!
!-------------------------------------------------------------------------------
   use commpi, only  :  npes,mype,master,comm_world,real_type,lonstr,lonlen,   &
                        latstr,latlen,ncol
   use dfsvar, only  :  igs,ige,jgs,jge,latdef
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  ib,jbw,iba,jbwa,ntotal
   integer  ::  jj,n,i,j,k,mk,ierr,ilen,ist,iend,nr,ip,jp,js
   real  ::  par(ib,jbw,ntotal),tot(igs:ige,jbwa,ntotal)
!
   real(_mpi_real_),allocatable  ::  tmpsnd(:),tmprcv(:)
   integer,allocatable  ::  len(:),loc(:)
!
   ip=lonlen(mype)
   jp=latlen(mype)
   allocate(tmpsnd(ip*jp*2*ntotal))
   allocate(tmprcv(iba*jbwa*ntotal))
   allocate(len(0:npes-1))
   allocate(loc(0:npes-1))
!
   if (ige-igs+1.ne.iba.or.jge-jgs+1.ne.jbwa) then
     write(96,'(A,4I5)') 'dimension error : ige,igs,iba,jge,jgs,jbwa=',        &
                          ige,igs,iba,jge,jgs,jbwa
     call mpabort
   endif

   mk=0
   do n = 0,npes-1
     loc(n)=mk
     len(n)=ntotal*lonlen(n)*latlen(n)*2
     mk=loc(n)+len(n)
   enddo

   mk=0
   do k = 1,ntotal
     do j = 1,jp
       tmpsnd(mk+1:mk+ip)=par(1:ip,j,k)
       mk=mk+ip
       tmpsnd(mk+1:mk+ip)=par(1:ip,jbw-j+1,k)
       mk=mk+ip
     enddo
   enddo
!
   call mpi_gatherv(tmpsnd,len(mype),real_type,tmprcv,len,loc,real_type,       &
                    master,comm_world,ierr)
   if (ierr /= 0) then
     print*,'mpi_gatherv error:',ierr
     call mpabort
   endif
!
   if( mype.eq.0 ) then
     mk=0
     do n = 0,npes-1
       nr=n/ncol
       ilen=lonlen(n)
       ist=lonstr(n)-1
       do k = 1,ntotal
         do j = 1,latlen(n)
           jj=latdef(latstr(n)+j-jgs)
           js=jbwa-jj+1
           iend=ist+ilen
           tot(ist+1:iend,jj,k)=tmprcv(mk+1:mk+ilen)
           mk=mk+ilen
           tot(ist+1:iend,js,k)=tmprcv(mk+1:mk+ilen)
           mk=mk+ilen
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
   end subroutine mpgp2fd
!-------------------------------------------------------------------------------
