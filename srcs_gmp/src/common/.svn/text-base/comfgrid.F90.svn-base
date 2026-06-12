#include "define.h"
   module comfgrid
!-------------------------------------------------------------------------------
   use paramodel, only   :  LONF2S,LATG2S,latg2_,jcap_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   private              ::  LONF2S,LATG2S,latg2_,jcap_
!
   real   , allocatable, dimension(:)    ::  colrad                           ,&
                                             sinlat                           ,&
                                             rcs2                             ,&
                                             colrab
   real   , allocatable, dimension(:,:)  ::  sinlab                           ,&
                                             coslab
#ifdef DFS
   real   , allocatable, dimension(:)    ::  rbs2
#else
#ifndef MP
   integer, allocatable, dimension(:)    ::  lwvdef                           ,&
                                             latdef
#endif
   real   , allocatable, dimension(:)    ::  wgt                              ,&
                                             wgtcs                            ,&
                                             wgb                              ,&
                                             wgbcs                            ,&
                                             rbs2
#endif
!
   contains
!-------------------------------------------------------------------------------
   subroutine comfgrid_init
!-------------------------------------------------------------------------------
   allocate(   colrad   (latg2_)                                              ,&
               sinlat   (latg2_)                                              ,&
               rcs2     (latg2_)                                              ,&
               colrab   (LATG2S)           )
   allocate(   sinlab   (LONF2S,LATG2S)                                       ,&
               coslab   (LONF2S,LATG2S)    )
#ifdef DFS
   allocate(   rbs2     (LATG2S))
#else
#ifndef MP
   allocate(   lwvdef   (jcap_+1)                                             ,&
               latdef   (latg2_)           )
#endif
   allocate(   wgt      (latg2_)                                              ,&
               wgtcs    (latg2_)                                              ,&
               wgb      (LATG2S)                                              ,&
               wgbcs    (LATG2S)                                              ,&
               rbs2     (LATG2S)           )
#endif
   end subroutine comfgrid_init
!-------------------------------------------------------------------------------
   end module comfgrid
