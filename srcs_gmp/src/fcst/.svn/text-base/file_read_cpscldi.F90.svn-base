#include "define.h"
   subroutine file_read_cpscldi(qcicpsi,qrscpsi,taucldi,cldwpi,cldipi)
#if defined(EXPLICIT_CLOUDINESS) && defined(CPS_QINI)
!-------------------------------------------------------------------------------
!
! subroutine: file_read_cpscldi
!
! program history log:
!   2012-09-25  suryun ham         read initial cloud properties 
!
! usage: call file_read_cpscldi(n,qcicpsi,qrscpsi,taucldi,cldwpi,cldipi)
!   output argument list:
!     qcicpsi   - cloud/ice water from CPS
!     qrscpsi   - rain/snow water from CPS
!     taucldi   - cloud optical depth for each layer
!     cldwpi    - cloud water path
!     cldipi    - cloud ice path
!
!-------------------------------------------------------------------------------
#ifdef MP
   use comio,only        : iope
   use paramodel, only   : lonf2p_,latg2p_
#endif
   use paramodel,only    : LONF2S,LATG2S,levs_,lonf2_,latg2_
!
   integer, parameter    :: nvar=5
   real                  :: qcicpsi(LONF2S,levs_,LATG2S)
   real                  :: qrscpsi(LONF2S,levs_,LATG2S)
   real                  :: taucldi(LONF2S,levs_,LATG2S)
   real                  :: cldwpi(LONF2S,levs_,LATG2S)
   real                  :: cldipi(LONF2S,levs_,LATG2S)
!
   real                  :: work1(lonf2_,latg2_,levs_)
   real                  :: work2(lonf2_,latg2_,levs_)
   real                  :: work3(lonf2_,latg2_,levs_)
   real                  :: work4(lonf2_,latg2_,levs_)
   real                  :: work5(lonf2_,latg2_,levs_)
#ifdef MP
   real                  :: work1p(lonf2p_,latg2p_,levs_)
   real                  :: work2p(lonf2p_,latg2p_,levs_)
   real                  :: work3p(lonf2p_,latg2p_,levs_)
   real                  :: work4p(lonf2p_,latg2p_,levs_)
   real                  :: work5p(lonf2p_,latg2p_,levs_)
#endif
!   
   work1 = 0. ; work2 = 0. ; work3 = 0. ; work4 = 0. ; work5 = 0.
#ifdef MP
   work1p = 0. ; work2p = 0. ; work3p = 0. ; work4p = 0. ; work5p = 0.
#endif
!
#ifdef MP
   if( iope ) then
#endif
     read(44) (((work1(i,j,k),i=1,lonf2_),j=1,latg2_),k=1,levs_)
     read(44) (((work2(i,j,k),i=1,lonf2_),j=1,latg2_),k=1,levs_)
     read(44) (((work3(i,j,k),i=1,lonf2_),j=1,latg2_),k=1,levs_)
     read(44) (((work4(i,j,k),i=1,lonf2_),j=1,latg2_),k=1,levs_)
     read(44) (((work5(i,j,k),i=1,lonf2_),j=1,latg2_),k=1,levs_)
     print*,'cpscldi read finished'
     close(44)
#ifdef MP
   endif  !/* iope */
#endif
!
#ifdef MP
   call MPGF2P(work1,lonf2_,latg2_,work1p,lonf2p_,latg2p_,levs_)
   call MPGF2P(work2,lonf2_,latg2_,work2p,lonf2p_,latg2p_,levs_)
   call MPGF2P(work3,lonf2_,latg2_,work3p,lonf2p_,latg2p_,levs_)
   call MPGF2P(work4,lonf2_,latg2_,work4p,lonf2p_,latg2p_,levs_)
   call MPGF2P(work5,lonf2_,latg2_,work5p,lonf2p_,latg2p_,levs_)
   call mpbcastr(work1p,levs_*LATG2S*LONF2S)
   call mpbcastr(work2p,levs_*LATG2S*LONF2S)
   call mpbcastr(work3p,levs_*LATG2S*LONF2S)
   call mpbcastr(work4p,levs_*LATG2S*LONF2S)
   call mpbcastr(work5p,levs_*LATG2S*LONF2S)
#endif
!
   do k = 1,levs_
     do j = 1,LATG2S
       do i = 1,LONF2S
#ifdef MP
         qcicpsi(i,k,j) = work1p(i,j,k)
         qrscpsi(i,k,j) = work2p(i,j,k)
         taucldi(i,k,j) = work3p(i,j,k)
         cldwpi(i,k,j)  = work4p(i,j,k)
         cldipi(i,k,j)  = work5p(i,j,k)
#else
         qcicpsi(i,k,j) = work1(i,j,k)
         qrscpsi(i,k,j) = work2(i,j,k)
         taucldi(i,k,j) = work3(i,j,k)
         cldwpi(i,k,j)  = work4(i,j,k)
         cldipi(i,k,j)  = work5(i,j,k)
#endif
       enddo
     enddo
   enddo
!
#endif
   end 
