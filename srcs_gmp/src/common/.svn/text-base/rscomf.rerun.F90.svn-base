#include <define.h>
   module rscomf_rerun
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only  :  IGRD12S,JGRD12S,LNWAVS,LNGRDS,                      &
                           levs_,levh_,mtnvar_,levp1_
   use comio
   use rscomver
   use rscompln
   use rscomrad
   use rscomloc
#ifdef MP
   use paramodel, only  :  llwavp_,levsp_,levhp_,levp1p_
!-------------------------------------------------------------------------------
   private             ::  llwavp_,levsp_,levhp_,levp1p_
#else
!-------------------------------------------------------------------------------
#endif
   private             :: IGRD12S,JGRD12S,LNWAVS,LNGRDS,levs_,levh_,mtnvar_
!
!  following common has to be in order with the call of sums in
!    rmp_dynamics_driver, rmp_physics_main_driver, rloopr, rloopz
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
#ifdef MP
   real   , allocatable, dimension(:)      ::  za                             ,&
                                               qa                             ,&
                                               qma                            ,&
                                               gza                            ,&
                                               dpdlama                        ,&
                                               dpdphia
   real   , allocatable, dimension(:,:)    ::  ya                             ,&
                                               rta                            ,&
                                               dia                            ,&
                                               tea                            ,&
                                               rqa                            ,&
                                               tema                           ,&
                                               rma                            ,&
                                               xa                             ,&
                                               uua                            ,&
                                               uuma                           ,&
                                               wa                             ,&
                                               vva                            ,&
                                               vvma
#endif
#ifdef NONHYD
   real   , allocatable, dimension(:,:)    ::  p                              ,&
                                               t                              ,&
                                               o                              ,&
                                               pn                             ,&
                                               tn                             ,&
                                               on                             ,&
                                               pnm                            ,&
                                               tnm                            ,&
                                               onm
#ifdef MP
   real   , allocatable, dimension(:,:)    ::  pa                             ,&
                                               ta                             ,&
                                               oa                             ,&
                                               pna                            ,&
                                               tna                            ,&
                                               ona                            ,&
                                               pnma                           ,&
                                               tnma                           ,&
                                               onma
#endif /* MP end */
#endif /* NONHYD end */
   real   , allocatable, dimension(:)      ::  flat                           ,&
                                               flon                           ,&
                                               fm2                            ,&
                                               fm2x                           ,&
                                               fm2y
!
! put into module of comsfc
!#include <comsfc.h>
!
   logical                                 ::  lastep
   real                                    ::  sdec,cdec,slag,solhr,clstp
   real   , allocatable, dimension(:,:)    ::  sfcnsw                         ,&
                                               sfcdlw                         ,&
#ifdef VIC
                                               sfcusw                         ,&
                                               sfcdsw                         ,&
#endif
                                               coszer                         ,&
                                               coszdg                         ,&
                                               sinlar                         ,&
                                               coslar                         ,&
                                               albed                          ,&
                                               tsflw
   real   , allocatable, dimension(:,:,:)  ::  swh                            ,&
                                               hlw                            ,&
                                               ozon
   real                                    ::  dtflux
   real   , allocatable, dimension(:,:)    ::  dusfc                          ,&
                                               dvsfc                          ,&
                                               dtsfc                          ,&
                                               dqsfc                          ,&
                                               dlwsfc                         ,&
                                               ulwsfc                         ,&
                                               raintot                        ,&
#ifdef NOAHYDRO
                                               raintot2                       ,&
#endif
                                               u10m                           ,&
                                               v10m                           ,&
                                               t2m                            ,&
                                               q2m                            ,&
                                               dugwd                          ,&
                                               dvgwd                          ,&
#ifdef GWDC
                                               dugwdc                         ,&
                                               dvgwdc                         ,&
#endif
                                               psurf                          ,&
                                               psmean
   real   , allocatable, dimension(:,:)    ::  raincps                        ,&
                                               gflux                          ,&
                                               slrad                          ,&
                                               runoff                         ,&
                                               tmpmax                         ,&
                                               tmpmin                         ,&
                                               ep                             ,&
                                               cldwrk                         ,&
                                               hpbl                           ,&
                                               pwat                           ,&
                                               snowmelt                       ,&
                                               snowfall                       ,&
                                               snowevap                       ,&
                                               qull                           ,&
                                               qvll                           ,&
                                               alhtfl                         ,&
                                               evcnp                          ,&
                                               bgrun
   real   , allocatable, dimension(:,:,:)  ::  hprime
#ifdef COUPLE_ROP
   real                                    ::  romsrtime, romsrtswup, romsrtswdn
   real   , allocatable, dimension(:,:)    ::  romsevap, romssens             ,&
                                               romsustr, romsvstr             ,&
                                               romslwup, romsswup             ,&
                                               romslwdn, romsswdn             ,&
                                               romsprcp, romssgz
#endif
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscomf_rerun_init
!-------------------------------------------------------------------------------
   allocate(   z       (LNWAVS)                                               ,&
               q       (LNWAVS)                                               ,&
               qm      (LNWAVS)                                               ,&
               gz      (LNWAVS)                                               ,&
               dpdlam  (LNWAVS)                                               ,&
               dpdphi  (LNWAVS)                                                )
   allocate(   y       (LNWAVS,levs_)                                         ,&
               rt      (LNWAVS,levh_)                                         ,&
               di      (LNWAVS,levs_)                                         ,&
               te      (LNWAVS,levs_)                                         ,&
               rq      (LNWAVS,levh_)                                         ,&
               tem     (LNWAVS,levs_)                                         ,&
               rm      (LNWAVS,levh_)                                         ,&
               x       (LNWAVS,levs_)                                         ,&
               uu      (LNWAVS,levs_)                                         ,&
               uum     (LNWAVS,levs_)                                         ,&
               w       (LNWAVS,levs_)                                         ,&
               vv      (LNWAVS,levs_)                                         ,&
               vvm     (LNWAVS,levs_)                                          )
#ifdef MP
   allocate(   za      (llwavp_)                                              ,&
               qa      (llwavp_)                                              ,&
               qma     (llwavp_)                                              ,&
               gza     (llwavp_)                                              ,&
               dpdlama (llwavp_)                                              ,&
               dpdphia (llwavp_)                                               )
   allocate(   ya      (llwavp_,levsp_)                                       ,&
               rta     (llwavp_,levhp_)                                       ,&
               dia     (llwavp_,levsp_)                                       ,&
               tea     (llwavp_,levsp_)                                       ,&
               rqa     (llwavp_,levhp_)                                       ,&
               tema    (llwavp_,levsp_)                                       ,&
               rma     (llwavp_,levhp_)                                       ,&
               xa      (llwavp_,levsp_)                                       ,&
               uua     (llwavp_,levsp_)                                       ,&
               uuma    (llwavp_,levsp_)                                       ,&
               wa      (llwavp_,levsp_)                                       ,&
               vva     (llwavp_,levsp_)                                       ,&
               vvma    (llwavp_,levsp_)                                        )
#endif
#ifdef NONHYD
   allocate(   p       (LNWAVS,levs_)                                         ,&
               t       (LNWAVS,levs_)                                         ,&
               o       (LNWAVS,levp1_)                                        ,&
               pn      (LNWAVS,levs_)                                         ,&
               tn      (LNWAVS,levs_)                                         ,&
               on      (LNWAVS,levp1_)                                        ,&
               pnm     (LNWAVS,levs_)                                         ,&
               tnm     (LNWAVS,levs_)                                         ,&
               onm     (LNWAVS,levp1_)                                         )
#ifdef MP
   allocate(   pa      (llwavp_,levsp_)                                       ,&
               ta      (llwavp_,levsp_)                                       ,&
               oa      (llwavp_,levp1p_)                                      ,&
               pna     (llwavp_,levsp_)                                       ,&
               tna     (llwavp_,levsp_)                                       ,&
               ona     (llwavp_,levp1p_)                                      ,&
               pnma    (llwavp_,levsp_)                                       ,&
               tnma    (llwavp_,levsp_)                                       ,&
               onma    (llwavp_,levp1p_)                                       )
#endif /* MP end */
#endif /* NONHYD end */
   allocate(   flat    (LNGRDS)                                               ,&
               flon    (LNGRDS)                                               ,&
               fm2     (LNGRDS)                                               ,&
               fm2x    (LNGRDS)                                               ,&
               fm2y    (LNGRDS)                                                )
   allocate(   sfcnsw  (IGRD12S,JGRD12S)                                      ,&
               sfcdlw  (IGRD12S,JGRD12S)                                      ,&
#ifdef VIC
               sfcusw  (IGRD12S,JGRD12S)                                      ,&
               sfcdsw  (IGRD12S,JGRD12S)                                      ,&
#endif
               coszer  (IGRD12S,JGRD12S)                                      ,&
               coszdg  (IGRD12S,JGRD12S)                                      ,&
               sinlar  (IGRD12S,JGRD12S)                                      ,&
               coslar  (IGRD12S,JGRD12S)                                      ,&
               albed   (IGRD12S,JGRD12S)                                      ,&
               tsflw   (IGRD12S,JGRD12S)                                       )
#ifdef RMPVECTORIZE
   allocate(   swh     (igrd12p_,jgrd12p_,levs_)                              ,&
               hlw     (igrd12p_,jgrd12p_,levs_)                              ,&
               ozon    (igrd12p_,jgrd12p_,levs_)                               )
#else
   allocate(   swh     (IGRD12S,levs_,JGRD12S)                                ,&
               hlw     (IGRD12S,levs_,JGRD12S)                                ,&
               ozon    (IGRD12S,levs_,JGRD12S)                                 )
#endif
   allocate(   dusfc   (IGRD12S,JGRD12S)                                      ,&
               dvsfc   (IGRD12S,JGRD12S)                                      ,&
               dtsfc   (IGRD12S,JGRD12S)                                      ,&
               dqsfc   (IGRD12S,JGRD12S)                                      ,&
               dlwsfc  (IGRD12S,JGRD12S)                                      ,&
               ulwsfc  (IGRD12S,JGRD12S)                                      ,&
               raintot (IGRD12S,JGRD12S)                                      ,&
#ifdef NOAHYDRO
               raintot2(IGRD12S,JGRD12S)                                      ,&
#endif
               u10m    (IGRD12S,JGRD12S)                                      ,&
               v10m    (IGRD12S,JGRD12S)                                      ,&
               t2m     (IGRD12S,JGRD12S)                                      ,&
               q2m     (IGRD12S,JGRD12S)                                      ,&
               dugwd   (IGRD12S,JGRD12S)                                      ,&
               dvgwd   (IGRD12S,JGRD12S)                                      ,&
#ifdef GWDC
               dugwdc  (IGRD12S,JGRD12S)                                      ,&
               dvgwdc  (IGRD12S,JGRD12S)                                      ,&
#endif
               psurf   (IGRD12S,JGRD12S)                                      ,&
               psmean  (IGRD12S,JGRD12S)                                       )
   allocate(   raincps (IGRD12S,JGRD12S)                                      ,&
               gflux   (IGRD12S,JGRD12S)                                      ,&
               slrad   (IGRD12S,JGRD12S)                                      ,&
               runoff  (IGRD12S,JGRD12S)                                      ,&
               tmpmax  (IGRD12S,JGRD12S)                                      ,&
               tmpmin  (IGRD12S,JGRD12S)                                      ,&
               ep      (IGRD12S,JGRD12S)                                      ,&
               cldwrk  (IGRD12S,JGRD12S)                                      ,&
               hpbl    (IGRD12S,JGRD12S)                                      ,&
               pwat    (IGRD12S,JGRD12S)                                      ,&
               snowmelt(IGRD12S,JGRD12S)                                      ,&
               snowfall(IGRD12S,JGRD12S)                                      ,&
               snowevap(IGRD12S,JGRD12S)                                      ,&
               qull    (IGRD12S,JGRD12S)                                      ,&
               qvll    (IGRD12S,JGRD12S)                                      ,&
               alhtfl  (IGRD12S,JGRD12S)                                      ,&
               evcnp   (IGRD12S,JGRD12S)                                      ,&
               bgrun   (IGRD12S,JGRD12S)                                       )
   allocate(   hprime  (IGRD12S,JGRD12S,mtnvar_)                               )
#ifdef COUPLE_ROP
   allocate(   romsevap(IGRD12S,JGRD12S)                                      ,&
               romssens(IGRD12S,JGRD12S)                                      ,&
               romsustr(IGRD12S,JGRD12S)                                      ,&
               romsvstr(IGRD12S,JGRD12S)                                      ,&
               romslwup(IGRD12S,JGRD12S)                                      ,&
               romsswup(IGRD12S,JGRD12S)                                      ,&
               romslwdn(IGRD12S,JGRD12S)                                      ,&
               romsswdn(IGRD12S,JGRD12S)                                      ,&
               romsprcp(IGRD12S,JGRD12S)                                      ,&
               romssgz (IGRD12S,JGRD12S)                                       )
#endif
!
   end subroutine rscomf_rerun_init
!-------------------------------------------------------------------------------
#endif /* RMP end*/
   end module rscomf_rerun
