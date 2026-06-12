#include <define.h>
   subroutine preduceg
!-------------------------------------------------------------------------------
!
! subroutine preduceg      programmer: hann-ming juang
!
! purpose: partial reduce grid routine for mpi to compute lcapd,lonfd
!          lcapdp and lonfdp for each node.
!
! input:
!   lcapd   (latg2) maximal wave number needed for all latitude
!   lonfd   (latg2) maximal grid number for all latitudes
!   latg2   latitude dimension
! output: all into commpi
!   lonfd   (latg2) full lonf for reduced grid as the same as lonfd
!   lcapdp   (latg2 ,npes) partial wave resolution for reduced grid
!   lonfdp   (latg2p,npes) partial lonf length for reduced grid
!   lonfds   (latg2p,npes) partial lonf start point for reduced grid
!
!-------------------------------------------------------------------------------
   use paramodel, only : latg2_
!  use commpi, only : npes, ncol, nrow, lwvdef, lwvstr, lwvlen, latstr,        &
   use commpi       ! npes, ncol, nrow, lwvdef, lwvstr, lwvlen, latstr,        &
                    ! latlen, latdef
#ifdef REDUCE_GRID
   use comreduce    ! lcapd, lcapdp, lonfds, lonfd, lonfdp
!-------------------------------------------------------------------------------
   integer, allocatable :: lentmp(:)
   allocate(lentmp(ncol))
!
   latg2 = latg2_
   do n = 0,npes-1
     llstr = lwvstr(n)
     llens = lwvlen(n)
     do lat = 1,latg2
       ls = 0
       do ll = 1,llens
         l = lwvdef(llstr+ll)
         if( l.lt.lcapd(lat) ) ls = ls+1
       enddo
       lcapdp(lat,n) = ls
     enddo
   enddo
!
   do n = 0,nrow-1
     jpe = n*ncol
     jjstr = latstr(jpe)
     jjens = latlen(jpe)
     do lat = 1,jjens
       j = latdef(jjstr+lat-1)
       call equdiv(lonfd(j),ncol,lentmp)
       mk = 1
       do m = 0,ncol-1
         np = n*ncol+m
         lonfds(lat,np) = mk
         lonfdp(lat,np) = lentmp(m+1)
         mk = mk+lonfdp(lat,np)
       enddo
     enddo
   enddo
!
#ifdef DBG
   if(mype.eq.master) then
     do j = 1,latg2_
       print *,' lat lcapd lonfd ',j,lcapd(j),lonfd(j)
     enddo
     do n = 0,npes-1
       print *,' ----- n= ',n,' ---------'
       do j = 1,latg2
         print *,' lat lcapdp ',j,lcapdp(j,n)
       enddo
       do j = 1,latlen(n)
         jj = latstr(n)+j-1
         print *,' lat lonfds lonfdp ',jj,lonfds(j,n),lonfdp(j,n)
       enddo
     enddo
   endif
#endif
#endif
!
   return
   end subroutine preduceg
!-------------------------------------------------------------------------------
