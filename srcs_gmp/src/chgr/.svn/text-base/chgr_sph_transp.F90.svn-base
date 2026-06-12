#include <define.h>
   subroutine chgr_sph_transp(tri, rho, iswtch)
!-------------------------------------------------------------------------------
!
! subprogram:    chgr_sph_transp  
!
! abstract: scalar triangle <--> rhomboid shuffle.
!
! program history log:
!   1988-04-21  joseph sela
!   2000-01-01  song-you hong          usgs topo, kagwd options
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call chgr_sph_transp (tri, rho, iswitch)
!   input argument list:
!     tri      - triangle.
!                when iswitch=-1, triangle tri is stored in
!                the lower triangle of rhomboid rho, and
!                the upper triangle of rhomboid rho is set to zeros.
!     rho      - rhomboid.
!                when iswitch=1, the lower triangle of rhomboid rho
!                is stored in triangle tri.
!     iswitch  - integer switch.
!                shufle must be called first with iswitch=0
!                to initialize bit vector bitv.
!
!   output argument list:
!     tri      - triangle.
!                when iswitch=1, the lower triangle of rhomboid rho
!                is stored in triangle tri.
!     rho      - rhomboid.
!                when iswitch=-1, triangle tri is stored in
!                the lower triangle of rhomboid rho, and
!                the upper triangle of rhomboid rho is set to zeros.
!
!-------------------------------------------------------------------------------
   use paramter
   use comchgr, only : mdim, mwvp2, mwavep
!-------------------------------------------------------------------------------
   save
!
   real                 ::  tri( mdim ), rho( mwvp2 , mwavep )
!-------------------------------------------------------------------------------
!
!     when iswtch=1, the lower triangle of rhomboid rho
!     is stored in triangle tri.
!
!     when iswtch=-1, triangle tri is stored in
!     the lower triangle of rhomboid rho, and
!     the upper triangle of rhomboid rho is set to zeros.
!
   if (iswtch) 200, 400, 600
!
   200 continue
!
   do k = 1, mwavep
     do i = 1, mwvp2
       rho(i,k) = 0.0e0
     enddo
   enddo
   iplus = 0
   len =  mwvp2
   do k = 1, mwavep
     do i = 1,len
       rho(i,k) = tri(i+iplus)
     enddo
     iplus = iplus + len
     len = len - 2
   enddo
   return
!
   400 continue
   return
!
   600 continue
!
   iplus = 0
   len =  mwvp2
   do k = 1, mwavep
     do i = 1,len
       tri(i+iplus) = rho(i,k)
     enddo
     iplus = iplus + len
     len = len - 2
   enddo
!
   return
   end
!-------------------------------------------------------------------------------
