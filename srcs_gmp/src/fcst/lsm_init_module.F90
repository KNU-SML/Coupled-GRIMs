#include <define.h>
   module lsm_init_module
!-------------------------------------------------------------------------------
   use paramodel, only   :  ntype=>nstype_
   use comfcst, only     :  ngrid, dfk, ktk, b, satpsi, satkt, tsat, dfkt
!-------------------------------------------------------------------------------
!
!  abstract ::
!    lsm_init_df sets up moisture diffusivity and hydrolic conductivity
!    for all soil types
!    grddfs sets up thermal diffusivity for all soil types
!
!  ::: structure :::
!    
!    [gmp_start_setup] --- [lsm_init_module]
!                                |
!                                |---[dfkt_init] *
!                                |---[lsm_init_df] *
!                                |---[lsm_init_kt] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine dfkt_init
!-------------------------------------------------------------------------------
!
!  the nine soil types are:
!    1  ... loamy sand (coarse)
!    2  ... silty clay loam (medium)
!    3  ... light clay (fine)
!    4  ... sandy loam (coarse-medium)
!    5  ... sandy clay (coarse-fine)
!    6  ... clay loam  (medium-fine)
!    7  ... sandy clay loam (coarse-med-fine)
!    8  ... loam  (organic)
!    9  ... ice (use loamy sand property)
!
!    nstype = 16 from stagsgo soil data
!    soil types    statsgo (miller ??, 199?)  cosby et al (1984)
!             1          sand                  sand
!             2          loamy sand            loamy sand
!             3          sandy loam            sandy loam
!             4          silt loam             silty loam
!             5          silt                  silty loam
!             6          loam                  loam
!             7          sandy clay loam       sandy clay loam
!             8          silty clay loam       silty clay loam
!             9          clay loam             clay loam
!            10          sandy clay            sandy clay
!            11          silty clay            silty clay
!            12          clay                  light clay
!            13          organic materials     loam
!            14          water
!            15          bedrock
!            16          other (land-ice)
!                        the value of this class is the same as in classe-2
!
!-------------------------------------------------------------------------------
#ifndef OSULSM1
#ifdef STATSGO_SOIL
   b=(/2.79,  4.26,  4.74,  5.33,  5.33,  5.25,  6.66,                         &
       8.72,  8.17,  10.73, 10.39, 11.55, 5.25,  0.0, 2.79,  4.26/)
   satpsi=(/ 0.069, 0.036, 0.141, 0.759, 0.759, 0.355, 0.135,                  &
             0.617, 0.263, 0.098, 0.324, 0.468, 0.355, 0.0, 0.069, 0.036/)
   satkt=(/ 1.07e-6, 1.41e-5, 5.23e-6, 2.81e-6, 2.81e-6,                       &
            3.38e-6, 4.45e-6, 2.04e-6, 2.45e-6, 7.22e-6, 1.34e-6,              &
            9.74e-7, 3.38e-6, 0.0, 1.41e-4, 1.41e-5/)
   tsat=(/0.339, 0.421, 0.434, 0.476, 0.476, 0.439, 0.404,                     &
          0.464, 0.465, 0.406, 0.468, 0.468, 0.439, 1.0, 0.20, 0.421/)
#else
   b=(/ 4.26,8.72,11.55,4.74,10.73,8.17,6.77,5.25,4.26 /)
   satpsi=(/.04,.62,.47,.14,.10,.26,.14,.36,.04/)
   satkt=(/1.41e-5,.20e-5,.10e-5,.52e-5,.72e-5,                                &
              .25e-5,.45e-5,.34e-5,1.41e-5/)
   tsat=(/.421,.464,.468,.434,.406,.465,.404,.439,.421/)
#endif
#else
   b=(/4.05,4.38,4.9,5.3,5.39,7.12,7.75,8.52,10.4/)
!    &       ,10.4,11.4/)
   satpsi=(/.121,.09,.218,.786,.478,.299,.356,.63,.153/)
!    &            ,.49,.405/)
   satkt=(/1.76e-4,1.5633e-4,3.467e-5,7.2e-6,6.95e-6                           &
              ,6.3e-6,1.7e-6,2.45e-6,2.167e-6/)
!    &           ,1.033e-6,1.283e-6/)
   tsat=(/.395,.41,.435,.485,.451,.42,.477,.476,.426/)
!    &          ,.492,.482/)
#endif
   end subroutine dfkt_init
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine lsm_init_df
!-------------------------------------------------------------------------------
   do k = 1,ntype
     dynw = tsat(k) * .05
     f1 = b(k) * satkt(k) * satpsi(k) / tsat(k) ** (b(k) + 3.)
     f2 = satkt(k) / tsat(k) ** (b(k) * 2. + 3.)
!
!  convert from m/s to kg m-2 s-1 unit
!
     f1 = f1 * 1000.
     f2 = f2 * 1000.
     do i = 1, ngrid
       theta = float(i-1) * dynw
       theta = min(tsat(k),theta)
       dfk(i,k) = f1 * theta ** (b(k) + 2.)
       ktk(i,k) = f2 * theta ** (b(k) * 2. + 3.)
     enddo
   enddo
!
   return
   end subroutine lsm_init_df
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine lsm_init_kt
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer             ::  k,i
   real                ::  dynw,theta,pf,f1
!
   do k = 1,ntype
     if(satpsi(k).gt.0.) then
       dynw = tsat(k) * .05
       f1 = log10(satpsi(k)) + b(k) * log10(tsat(k)) + 2.
       do i = 1,ngrid
         theta = float(i-1) * dynw
         theta = min(tsat(k),theta)
         if(theta.gt.0.) then
           pf = f1 - b(k) * log10(theta)
         else
           pf = 5.2
         endif
         if(pf.le.5.1) then
           dfkt(i,k) = exp(-(2.7+pf)) * 420.
         else
           dfkt(i,k) = .1744
         endif
       enddo
     endif
   enddo
!
   return
   end subroutine lsm_init_kt
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   end module lsm_init_module
