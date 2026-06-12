#include <define.h>
   SUBROUTINE dfs_sincos_lat
!-------------------------------------------------------------------------------
   use constant, only : pi=>pi_
   use dfsvar, only   : jbwA,jbA,jgs,jge,SINLAT,COSLAT,TANLAT   ! update
!-------------------------------------------------------------------------------
!                                                                    
! abstract :     SINE OF LATITUDE, COSINE OF LATITUDE                           
!                                                                    
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real     ::  dy,bb
   integer  ::  j
!
!**********************************************************!
!       MERIDIONAL EIGENFUNCTION ------SIN and COS         !
!**********************************************************!
!
   if (.not.allocated(coslat)) allocate(coslat(jgs:jge,3))
   if (.not.allocated(sinlat)) allocate(sinlat(jgs:jge))
   if (.not.allocated(tanlat)) allocate(tanlat(jgs:jge))
   DY = PI/JBWA
!
   do j = 0,JBWA-1
     bb= 0.5*pi-DY*(j+0.5)
#ifdef IBMSP
     SINLAT(j+jgs  )= SIN(dble(bb))      ! SIN(LAT)
     COSLAT(j+jgs,1)= COS(dble(bb))      ! COS(LAT)
#else
     SINLAT(j+jgs  )= DSIN(dble(bb))      ! SIN(LAT)
     COSLAT(j+jgs,1)= DCOS(dble(bb))      ! COS(LAT)
#endif
     COSLAT(j+jgs,2)= 1./COSLAT(j+jgs,1)
     COSLAT(j+jgs,3)= COSLAT(j+jgs,2)*COSLAT(j+jgs,2)
     TANLAT(j+jgs  )= SINLAT(j+jgs  )*COSLAT(j+jgs,2)
   end do
!
   RETURN
   END SUBROUTINE dfs_sincos_lat
