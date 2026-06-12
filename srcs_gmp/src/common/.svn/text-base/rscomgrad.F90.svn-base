#include <define.h>
   module rscomgrad
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only : IGRD12S,JGRD12S,LONF2S,LATG2S,levs_,jgrd1_
   use aerparm
#ifdef RMPVECTORIZE
   use paramodel, only : igrd12p_,jgrd12p_
!-------------------------------------------------------------------------------
   private            :: igrd12p_,jgrd12p_
#else
!-------------------------------------------------------------------------------
#endif
   private            :: IGRD12S,JGRD12S,LONF2S,LATG2S,levs_,jgrd1_
!
   integer, parameter                      ::  mcld  = 3                      ,&
                                               nseal = 2                      ,&
                                               nbin  = 100                    ,&
                                               nlon  = 2                      ,&
                                               nlat  = 4
   integer                 ::  ier
#ifndef SWRMDC
   integer                                 ::  nalaer
#else
   integer                                 ::  nfalb,nfaer
#endif
   integer                                 ::  idtln,idtls,istrat,jo3         ,&
                                               kalb,jsno,itimsw,itimlw        ,&
                                               icfc,icwp,iswsrc(5),ibnd,ko3
   real                                    ::  solc,rsin1,rcos1,rcos2         ,&
                                               tsmin,tsmax,shmin,shmax        ,&
                                               raddt,fjd,r1,alf               ,&
                                               dlt,dlon,dtlw,rrs2
#if defined (RRTMGSW) || defined (RRTMGLW)
   real                                    ::  jd
#endif
   real                                    ::  rhcl(nbin,nlon,nlat,mcld,nseal)
   real   , allocatable, dimension(:,:)    ::  cvr                            ,&
                                               cvtr                           ,&
                                               cvbr                           ,&
                                               alvbr                          ,&
                                               alnbr                          ,&
                                               alvdr                          ,&
                                               alndr                          ,&
                                               avecld                         ,&
                                               cldl                           ,&
                                               avecv                          ,&
                                               zonht                          ,&
                                               cldsig                         ,&
                                               facsf                          ,&
                                               facwf       
   real   , allocatable, dimension(:,:,:)  ::  alvsf                          ,&
                                               alnsf                          ,&
                                               alvwf                          ,&
                                               alnwf                          ,&
                                               paerr
#ifndef SWRMDC
   real   , allocatable, dimension(:,:,:)  ::  paerf
#else
   real   , allocatable, dimension(:,:,:)  ::  idxcg                          ,&
                                               cmixg                          ,&
                                               denng                          ,&
   real   , allocatable, dimension(:,:)    ::  kprfg
#endif
   logical                                 ::  runrad
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscomgrad_init
!-------------------------------------------------------------------------------
   allocate(   cvr     (IGRD12S,JGRD12S)                                      ,&
               cvtr    (IGRD12S,JGRD12S)                                      ,&
               cvbr    (IGRD12S,JGRD12S)                                      ,&
               alvbr   (IGRD12S,JGRD12S)                                      ,&
               alnbr   (IGRD12S,JGRD12S)                                      ,&
               alvdr   (IGRD12S,JGRD12S)                                      ,&
               alndr   (IGRD12S,JGRD12S)                                      ,&
               avecld  (levs_,jgrd1_)                                         ,&
               cldl    (4,jgrd1_)                                             ,&
               avecv   (3,jgrd1_)                                             ,&
               zonht   (levs_,jgrd1_)                                         ,&
               cldsig  (levs_,18)                                             ,&
               facsf   (LONF2S,LATG2S)                                        ,&
               facwf   (LONF2S,LATG2S)                                         )
   allocate(   alvsf   (LONF2S,LATG2S,4)                                      ,&
               alnsf   (LONF2S,LATG2S,4)                                      ,&
               alvwf   (LONF2S,LATG2S,4)                                      ,&
               alnwf   (LONF2S,LATG2S,4)                                      ,&
#ifdef RMPVECTORIZE
               paerr   (igrd12p_*jgrd12p_,5,1)                                 )
#else
               paerr   (IGRD12S,5,JGRD12S)                                     )
#endif
#ifndef SWRMDC
   allocate(   paerf   (LONF2S,LATG2S,5)                                       )
#else
   allocate(   idxcg   (nxc,IGRD12S,JGRD12S)                                  ,&
               cmixg   (nxc,IGRD12S,JGRD12S)                                  ,&
               denng   (ndn,IGRD12S,JGRD12S)                                   )
   allocate(   kprfg   (IGRD12S,JGRD12S)                                       )
#endif
!
   end subroutine rscomgrad_init
!-------------------------------------------------------------------------------
#endif /* RMP end */
   end module rscomgrad
