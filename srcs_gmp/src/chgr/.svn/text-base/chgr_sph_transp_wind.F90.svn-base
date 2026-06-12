#include <define.h>
   subroutine chgr_sph_transp_wind(tri, rho, iswtch)
!-------------------------------------------------------------------------------
!
! subprogram:    chgr_sph_transp_wind  
!
! abstract: vector triangle <--> rhomboid shuffle.
!
! program history log:
!   1988-01-01  sela
!   2000-01-01  song-you hong          usgs topo, kagwd options
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call chgr_sph_transp_wind (tri, rho, iswitch)
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
!-------------------------------------------------------------------------------
   use varsfc, only : lsoil_
   use paramter
   use comchgr, only : mdimv, mwvp2, mwave2
   save
!-------------------------------------------------------------------------------
   real                 :: tri( mdimv ), rho( mwvp2 , mwave2 )
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
   do k = 1, mwave2
     do i = 1, mwvp2
       rho(i,k) = 0.0
     enddo
   enddo
   do i = 1, mwvp2
     rho(i,1) = tri(i)
   enddo
   iplus =  mwvp2
   len   =  mwvp2
   do k = 2, mwave2
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
   do i = 1, mwvp2
     tri(i) = rho(i,1)
   enddo
   iplus =  mwvp2
   len   =  mwvp2
   do k = 2, mwave2
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
