!
   module co2dta
!-------------------------------------------------------------------------------
   use rdparm
!-------------------------------------------------------------------------------
!
!   the following common blocks contain pretabulated co2 transmission
!       functions, evaluated using the methods of fels and
!       schwarzkopf (1981) and schwarzkopf and fels (1985),
!   common co2bd3 contains co2 transmission functions and temperature
!   and pressure derivatives for the 560-800 cm-1 band. also included
!   are the standard temperatures and the weighting function. these
!   data are in block data bd3:
!         co251    =  transmission fctns for t0 (std. profile)
!                       with p(sfc)=1013.25 mb
!         co258    =  transmission fctns. for t0 (std. profile)
!                       with p(sfc)= ^810 mb
!         cdt51    =  first temperature derivative of co251
!         cdt58    =  first temperature derivative of co258
!         c2d51    =  second temperature derivative of co251
!         c2d58    =  second temperature derivative of co251
!         co2m51   =  transmission fctns for t0 for adjacent pressure
!                        levels, with no pressure quadrature. used for
!                        nearby layer computations. p(sfc)=1013.25 mb
!         co2m58   =  same as co2m51,with p(sfc)= ^810 mb
!         cdtm51   =  first temperature derivative of co2m51
!         cdtm58   =  first temperature derivative of co2m58
!         c2dm51   =  second temperature derivative of co2m51
!         c2dm58   =  second temperature derivative of co2m58
!         stemp    =  standard temperatures for model pressure level
!                        structure with p(sfc)=1013.25 mb
!         gtemp    =  weighting function for model pressure level
!                        structure with p(sfc)=1013.25 mb.
!         b0       =  temp. coefficient used for co2 trans. fctn.
!                        correction for t(k). (see ref. 4 and bd3)
!         b1       =  temp. coefficient, used along with b0
!         b2       =  temp. coefficient, used along with b0
!         b3       =  temp. coefficient, used along with b0
!
! --- common co2bd3 ---
!
   real               ::  b0,b1,b2,b3
   real, allocatable  ::  co251(:,:),co258(:,:),cdt51(:,:),                    &
                          cdt58(:,:),c2d51(:,:),c2d58(:,:),co2m51(:),          &
                          co2m58(:),cdtm51(:),cdtm58(:),c2dm51(:),c2dm58(:),   &
                          stemp(:),gtemp(:)
!
!   common co2bd2 contains co2 transmission functions and temperature
!   and pressure derivatives for the 560-670 cm-1 part of the 15 um
!   co2 band.  these data are in block data bd2.
!         co231    =  transmission fctns for t0 (std. profile)
!                       with p(sfc)=1013.25 mb
!         co238    =  transmission fctns. for t0 (std. profile)
!                       with p(sfc)= ^810 mb
!         cdt31    =  first temperature derivative of co231
!         cdt38    =  first temperature derivative of co238
!         c2d31    =  second temperature derivative of co231
!         c2d38    =  second temperature derivative of co231
!
! --- common co2bd2 ---
!
   real, allocatable  ::  co231(:),co238(:),cdt31(:),                          &
                          cdt38(:),c2d31(:),c2d38(:)
!
!   common co2bd4 contains co2 transmission functions and temperature
!   and pressure derivatives for the 670-800 cm-1 part of the 15 um
!   co2 band.  these data are in block data bd4.
!         co271    =  transmission fctns for t0 (std. profile)
!                       with p(sfc)=1013.25 mb
!         co278    =  transmission fctns. for t0 (std. profile)
!                       with p(sfc)= ^810 mb
!         cdt71    =  first temperature derivative of co271
!         cdt78    =  first temperature derivative of co278
!         c2d71    =  second temperature derivative of co271
!         c2d78    =  second temperature derivative of co271
!
! --- common co2bd4 ---
!
   real, allocatable  ::  co271(:),co278(:),cdt71(:),                          &
                          cdt78(:),c2d71(:),c2d78(:)
!
!   common co2bd5 contains co2 transmission functions for the 2270-
!   2380 part of the 4.3 um co2 band. these data are in block data bd5.
!         co211    =  transmission fctns for t0 (std. profile)
!                       with p(sfc)=1013.25 mb
!         co218    =  transmission fctns. for t0 (std. profile)
!                       with p(sfc)= ^810 mb
!
!  --- common co2bd5 ---
!
   real, allocatable  ::  co211(:),co218(:)
!
   contains
!-------------------------------------------------------------------------------
   subroutine co2dta_init
   call rdparm_init
   allocate(   co251(lp1,lp1),co258(lp1,lp1),cdt51(lp1,lp1),                   &
               cdt58(lp1,lp1),c2d51(lp1,lp1),c2d58(lp1,lp1),co2m51(l),         &
               co2m58(l),cdtm51(l),cdtm58(l),c2dm51(l),c2dm58(l),              &
               stemp(lp1),gtemp(lp1),                                          &
               co231(lp1),co238(lp1),cdt31(lp1),                               &
               cdt38(lp1),c2d31(lp1),c2d38(lp1),                               &
               co271(lp1),co278(lp1),cdt71(lp1),                               &
               cdt78(lp1),c2d71(lp1),c2d78(lp1),                               &
               co211(lp1),co218(lp1)    )
!
   end subroutine co2dta_init
!-------------------------------------------------------------------------------
   end module co2dta
