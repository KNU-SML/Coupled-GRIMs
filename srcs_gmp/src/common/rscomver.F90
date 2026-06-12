#include <define.h>
   module rscomver
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only  :  levs_,levp1_,levm1_
   use rscomcpu
!-------------------------------------------------------------------------------
   private             ::  levs_,levp1_,levm1_
!
   integer               ::  kdt,inistp,limlow,numsum,nummax,ncldb1           ,&
                             nout0,nout,nrsmi0,nrsmi1,nrsmi2                  ,&
                             nrsmo1,nrsmo2,nrsmop,nrflxi,nrflxf,nrflxp,nrsfcp ,&
                             nrflip,nrflop,nrinit,nr2dda,nrpken
   real                  ::   thour,deltim,dt2,dthr,hdthr                     ,&
                              dtltb,dtcvav,dtswav,dtlwav                      ,&
                              rcl,sl1,fhour,shour,dtpost,filta,filtb          ,&
                              dk,tk,percut,bndrlx,cowave,dtwave
   real   , allocatable, dimension(:,:)    ::  am                             ,&
                                               hm                             ,&
                                               tm                             ,&
                                               bm                             ,&
                                               cm
   real   , allocatable, dimension(:)      ::  spdmax                         ,&
                                               si                             ,&
                                               sl                             ,&
                                               del                            ,&
                                               rdel2                          ,&
                                               rmsdot                         ,&
                                               ci                             ,&
                                               cl                             ,&
                                               tov                            ,&
                                               sv                             ,&
                                               rpi                            ,&
                                               p1                             ,&
                                               p2                             ,&
                                               h1                             ,&
                                               h2                             ,&
                                               rpirec
#ifdef NONHYD
   real   , allocatable, dimension(:,:)    ::  cmn                            ,&
                                               dmn                            ,&
                                               emn                            ,&
                                               fmn                            ,&
                                               gmn                            ,&
                                               hmn
#endif /* NONHYD */
#ifdef HYBRID
   real   , allocatable, dimension(:)      ::  ak5                            ,&
                                               bk5
#endif
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscomver_init
!-------------------------------------------------------------------------------
   allocate(    am     (levs_,levs_)                                          ,&
                hm     (levs_,levs_)                                          ,&
                tm     (levs_,levs_)                                          ,&
                bm     (levs_,levs_)                                          ,&
                cm     (levs_,levs_)                                          ,&
                spdmax (levs_)                                                ,&
                si     (levp1_)                                               ,&
                sl     (levs_)                                                ,&
                del    (levs_)                                                ,&
                rdel2  (levp1_)                                               ,&
                rmsdot (levm1_)                                               ,&
                ci     (levp1_)                                               ,&
                cl     (levs_)                                                ,&
                tov    (levs_)                                                ,&
                sv     (levs_ )                                               ,&
                rpi    (levm1_)                                               ,&
                p1     (levs_ )                                               ,&
                p2     (levs_)                                                ,&
                h1     (levs_)                                                ,&
                h2     (levs_ )                                               ,&
                rpirec (levm1_)                                                )
#ifdef NONHYD
   allocate(    cmn    (levp1_,levs_)                                         ,&
                dmn    (levp1_,levs_)                                         ,&
                emn    (levs_,levp1_)                                         ,&
                fmn    (levs_ ,levs_)                                         ,&
                gmn    (levs_,levp1_)                                         ,&
                hmn    (levs_ ,levs_)                                          )
#endif /* NONHYD */
#ifdef HYBRID
   allocate(    ak5    (levp1_)                                               ,&
                bk5    (levp1_)                                                )
#endif
!
   end subroutine rscomver_init
!-------------------------------------------------------------------------------
#endif /* RMP end */
   end module rscomver
