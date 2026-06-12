#include <define.h>   
   module rscomf
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only : lnwav_,levs_,lngrd_,igrd12_,jgrd12_
   use varsfc, only    : lsoil_
   use comio
   use rscomver
   use rscomrad
   use rscomloc
   use rscompln
!-------------------------------------------------------------------------------
   private            :: lnwav_,levs_,lngrd_,igrd12_,jgrd12_
   private            :: lsoil_
!
!  begin two rloop (comf)
!     in-code nesting regional to global
!
! fillowing common has to be in order with the call of sums in
!     rmp_dynamics_driver, rmp_physics_main_driver, rloopr, rloopz
!
   integer,              dimension(4)      ::  idate
   real   , allocatable, dimension(:)      ::  z                              ,&
                                               q                              ,&
                                               qm                             ,&
                                               gz                             ,&
                                               dpdlam                         ,&
                                               dpdphi
   real   , allocatable, dimension(:,:)    ::  y                              ,&
                                               rt                             ,&
                                               di                             ,&
                                               te                             ,&
                                               rq                             ,&
                                               tem                            ,&
                                               rm                             ,&
                                               x                              ,&
                                               uu                             ,&
                                               uum                            ,&
                                               w                              ,&
                                               vv                             ,&
                                               vvm
   real   , allocatable, dimension(:)      ::  flat                           ,&
                                               flon                           ,&
                                               fm2                            ,&
                                               fm2x                           ,&
                                               fm2y
   logical                                 ::  lastep
   real                                    ::  sdec,cdec,slag,solhr,clstp
   real   , allocatable, dimension(:,:)    ::  slmsk                          ,&
                                               hprime                         ,&
                                               sfcnsw                         ,&
                                               sfcdlw                         ,&
#ifdef VIC
                                               sfcusw                         ,&
                                               sfcdsw                         ,&
#endif
                                               sinlar                         ,&
                                               coslar                         ,&
                                               albed                          ,&
                                               coszer                         ,&
                                               cv                             ,&
                                               cvt                            ,&
                                               cvb                            ,&
                                               tsflw
   real   , allocatable, dimension(:,:,:)  ::  swh                            ,&
                                               hlw
   real                                    ::  dtflux
   real   , allocatable, dimension(:,:)    ::  dusfc                          ,&
                                               dvsfc                          ,&
                                               dtsfc                          ,&
                                               dqsfc                          ,&
                                               dlwsfc                         ,&
                                               ulwsfc                         ,&
                                               raintot                        ,&
                                               tsea                           ,&
#ifdef NOAHYDRO
                                               raintot2                       ,&
#endif
                                               f10m                           ,&
                                               u10m                           ,&
                                               v10m                           ,&
                                               t2m                            ,&
                                               q2m                            ,&
                                               dugwd                          ,&
                                               dvgwd                          ,&
                                               psurf                          ,&
                                               psmean
   real   , allocatable, dimension(:,:)    ::  tg3                            ,&
                                               z0cm                           ,&
                                               plantr                         ,&
                                               snoweq                         ,&
                                               raincps                        ,&
                                               gflux                          ,&
                                               slrad                          ,&
                                               canopy                         ,&
                                               runoff                         ,&
                                               tmpmax                         ,&
                                               tmpmin                         ,&
                                               ep                             ,&
                                               cldwrk                         ,&
                                               hpbl                           ,&
                                               pwat
   real   , allocatable, dimension(:,:,:)  ::  smc                            ,&
                                               stc
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscomf_init
!-------------------------------------------------------------------------------
   allocate(   z        (lnwav_)                                              ,&
               q        (lnwav_)                                              ,&
               qm       (lnwav_)                                              ,&
               gz       (lnwav_)                                              ,&
               dpdlam   (lnwav_)                                              ,&
               dpdphi   (lnwav_)                                               )
   allocate(   y        (lnwav_,levs_)                                        ,&
               rt       (lnwav_,levs_)                                        ,&
               di       (lnwav_,levs_)                                        ,&
               te       (lnwav_,levs_)                                        ,&
               rq       (lnwav_,levs_)                                        ,&
               tem      (lnwav_,levs_)                                        ,&
               rm       (lnwav_,levs_)                                        ,&
               x        (lnwav_,levs_)                                        ,&
               uu       (lnwav_,levs_)                                        ,&
               uum      (lnwav_,levs_)                                        ,&
               w        (lnwav_,levs_)                                        ,&
               vv       (lnwav_,levs_)                                        ,&
               vvm      (lnwav_,levs_)                                         )
   allocate(   flat     (lngrd_)                                              ,&
               flon     (lngrd_)                                              ,&
               fm2      (lngrd_)                                              ,&
               fm2x     (lngrd_)                                              ,&
               fm2y     (lngrd_)                                               )
   allocate(   slmsk    (igrd12_,jgrd12_)                                     ,&
               hprime   (igrd12_,jgrd12_)                                     ,&
               sfcnsw   (igrd12_,jgrd12_)                                     ,&
               sfcdlw   (igrd12_,jgrd12_)                                     ,&
#ifdef VIC
               sfcusw   (igrd12_,jgrd12_)                                     ,&
               sfcdsw   (igrd12_,jgrd12_)                                     ,&
#endif
               sinlar   (igrd12_,jgrd12_)                                     ,&
               coslar   (igrd12_,jgrd12_)                                     ,&
               albed    (igrd12_,jgrd12_)                                     ,&
               coszer   (igrd12_,jgrd12_)                                     ,&
               cv       (igrd12_,jgrd12_)                                     ,&
               cvt      (igrd12_,jgrd12_)                                     ,&
               cvb      (igrd12_,jgrd12_)                                     ,&
               tsflw    (igrd12_,jgrd12_)                                      )
#ifdef RMPVECTORIZE
   allocate(   swh      (igrd12p_,jgrd12p_,levs_)                             ,&
               hlw      (igrd12p_,jgrd12p_,levs_)                              )
#else
   allocate(   swh      (igrd12_,levs_,jgrd12_)                               ,&
               hlw      (igrd12_,levs_,jgrd12_)                                )
#endif
   allocate(   dusfc    (igrd12_,jgrd12_)                                     ,&
               dvsfc    (igrd12_,jgrd12_)                                     ,&
               dtsfc    (igrd12_,jgrd12_)                                     ,&
               dqsfc    (igrd12_,jgrd12_)                                     ,&
               dlwsfc   (igrd12_,jgrd12_)                                     ,&
               ulwsfc   (igrd12_,jgrd12_)                                     ,&
               raintot  (igrd12_,jgrd12_)                                     ,&
               tsea     (igrd12_,jgrd12_)                                     ,&
#ifdef NOAHYDRO
               raintot2 (igrd12_,jgrd12_)                                     ,&
#endif
               f10m     (igrd12_,jgrd12_)                                     ,&
               u10m     (igrd12_,jgrd12_)                                     ,&
               v10m     (igrd12_,jgrd12_)                                     ,&
               t2m      (igrd12_,jgrd12_)                                     ,&
               q2m      (igrd12_,jgrd12_)                                     ,&
               dugwd    (igrd12_,jgrd12_)                                     ,&
               dvgwd    (igrd12_,jgrd12_)                                     ,&
               psurf    (igrd12_,jgrd12_)                                     ,&
               psmean   (igrd12_,jgrd12_)                                      )
   allocate(   tg3      (igrd12_,jgrd12_)                                     ,&
               z0cm     (igrd12_,jgrd12_)                                     ,&
               plantr   (igrd12_,jgrd12_)                                     ,&
               snoweq   (igrd12_,jgrd12_)                                     ,&
               raincps  (igrd12_,jgrd12_)                                     ,&
               gflux    (igrd12_,jgrd12_)                                     ,&
               slrad    (igrd12_,jgrd12_)                                     ,&
               canopy   (igrd12_,jgrd12_)                                     ,&
               runoff   (igrd12_,jgrd12_)                                     ,&
               tmpmax   (igrd12_,jgrd12_)                                     ,&
               tmpmin   (igrd12_,jgrd12_)                                     ,&
               ep       (igrd12_,jgrd12_)                                     ,&
               cldwrk   (igrd12_,jgrd12_)                                     ,&
               hpbl     (igrd12_,jgrd12_)                                     ,&
               pwat     (igrd12_,jgrd12_)                                      )
   allocate(   smc      (igrd12_,jgrd12_,lsoil_)                              ,&
               stc      (igrd12_,jgrd12_,lsoil_)                               )
!
   end subroutine rscomf_init
!-------------------------------------------------------------------------------
#endif /* RMP end */
   end module rscomf
