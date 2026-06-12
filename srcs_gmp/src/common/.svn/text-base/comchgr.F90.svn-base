#include <define.h>
   module comchgr
!-------------------------------------------------------------------------------
   use paramodel, only : jcap1_,jcap2_, lnt2_, lnut2_, twoj1_
   use paramter, only  : idim ,jdim ,kdim ,mwave ,kdimq
   use parmchgr, only  : idimi,jdimi,kdimi,mwavei,kdimqi
   use padchgr, only   : padchgr_init, nscrch, ngwri
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  idate(4)
   integer              ::  idimt
   integer              ::  mdim
   integer              ::  mdim2
   integer              ::  mdimi
#ifdef DFS
   integer              ::  idimtd
   integer              ::  mdimd
#ifdef ALIASED
   integer              ::  mw
#endif
#endif
   integer              ::  jdimhf                                            ,&
                            ijdim                                             ,&
                            mdimv                                             ,&
                            kdimp                                             ,&
                            kdimm                     
   integer              ::  jdimhfi                                           ,&
                            ijdimi                                            ,&
                            mdimvi                                            ,&
                            kdimpi                                            ,&
                            kdimmi                     
   integer              ::  mwvp2                                             ,&
                            mdimhf                                            ,&
                            mdimvh                                            ,&
                            kdimt                                             ,&
                            idimti                                            ,&
                            mwvp2i                                            ,&
                            mdimhfi                                           ,&
                            mdimvhi
   integer              ::  npad
   integer              ::  mwavep, mwave2
!
!  common for input arrays
!
   real   , allocatable, dimension(:)      ::  qi                             ,&
                                               gzi
   real   , allocatable, dimension(:,:)    ::  tei                            ,&
                                               dii                            ,&
                                               zei                            ,&
                                               rqi
!
!  common for modified horizontal resolution
!
   real   , allocatable, dimension(:)      ::  q                              ,&
                                               gz
   real   , allocatable, dimension(:,:)    ::  te                             ,&
                                               di                             ,&
                                               ze                             ,&
                                               rq
   real   , allocatable, dimension(:)      ::  ps
   real   , allocatable, dimension(:,:)    ::  tf                             ,&
                                               rqf                            ,&
                                               uf                             ,&
                                               vf
!
!  common for arrays of modified vertical resolution
!
   real   , allocatable, dimension(:)      ::  qo                             ,&
                                               gzo
   real   , allocatable, dimension(:,:)    ::  teo                            ,&
                                               dio                            ,&
                                               zeo                            ,&
                                               rqo
   real   , allocatable, dimension(:)      ::  pso
   real   , allocatable, dimension(:,:)    ::  tfo                            ,&
                                               rqfo                           ,&
                                               ufo                            ,&
                                               vfo
!
!  common label
!
   character(len=4)                        ::  lab(8)
!
!  common cchgs
!
#ifndef DFS
   real   , allocatable, dimension(:)      ::  zss
#endif
   real   , allocatable, dimension(:)      ::  pss
   real   , allocatable, dimension(:,:)    ::  tts                            ,&
                                               qqs                            ,&
                                               uus                            ,&
                                               vvs
!
!  common cchga
!
#ifndef DFS
   real   , allocatable, dimension(:)      ::  zsa
#endif
   real   , allocatable, dimension(:)      ::  psa
   real   , allocatable, dimension(:,:)    ::  tta                            ,&
                                               qqa                            ,&
                                               uua                            ,&
                                               vva
!
!  common inpver
!
   real   , allocatable, dimension(:)      ::  siin                           ,&
                                               slin                           ,&
                                               delin                          ,&
                                               ciin                           ,&
                                               clin                           ,&
                                               rpiin                          ,&
                                               ak                             ,&
                                               bk
!
!  common comver
!
   real   , allocatable, dimension(:)      ::  si                             ,&
                                               sl                             ,&
                                               del                            ,&
                                               ci                             ,&
                                               cl                             ,&
                                               rpi                            ,&
                                               ak5                            ,&
                                               bk5
!
!  common plncom
!
   real   , allocatable, dimension(:)      ::  eps                            ,&
                                               colrad                         ,&
                                               wgt                            ,&
                                               wgtcs                          ,&
                                               rcs2
!
!  common comind
!
   integer, allocatable, dimension(:)      ::  indxnn                         ,&
                                               indxmm
!
!  common scrtch
!
   real   , allocatable, dimension(:,:,:)  ::  d                              ,&
                                               z                              ,&
                                               u                              ,&
                                               v
!
!  common gozcom
!
   real   , allocatable, dimension(:)      ::  dxa                            ,&
                                               dxb                            ,&
                                               dxc                            ,&
                                               dxd
!
!  common scratch
!
   real   , allocatable, dimension(:)      ::  b                              ,&
                                               c                              ,&
                                               pad
!
!  common chgr_sph_functc
!
   real                                    ::  jfir
   real   , allocatable, dimension(:)      ::  deps                           ,&
                                               rdeps                          ,&
                                               dx                             ,&
                                               y
!
   contains
!-------------------------------------------------------------------------------
   subroutine comchgr_init
!-------------------------------------------------------------------------------
   call padchgr_init
!
   idimt=idim*2
#ifdef DFS
   idimtd=idim*jdim
#ifdef ALIASED
   mw=(3*mwave+1)/2
#ifdef ALIASED2
   mdimd=(2*mwave+1)*((jdim-2)+1)
#else
   mdimd=(2*mwave+1)*(mw-mod(mw,2)+1)
#endif
#else
   mdimd=(2*mwave+1)*(mwave+1)
#endif /* ALIASED end */
#endif /* DFS end */
   mdim=(mwave+1)*(mwave+2)
   mdimi=(mwavei+1)*(mwavei+2)
   mdim2=(mwave+1)*(mwave+1)*2
!
   jdimhf=jdim/2
   ijdim=idim*jdim
   mdimv=(mwave+1)*(mwave+4)
   kdimp=kdim+1
   kdimm=kdim-1
   jdimhfi=jdimi/2
   ijdimi=idimi*jdimi
   mdimvi=(mwavei+1)*(mwavei+4)
   kdimpi=kdimi+1
   kdimmi=kdimi-1
!
   mwvp2=(mwave+1)*2
   mdimhf=mdim/2
   mdimvh=mdimv/2
   kdimt=kdim*2
   idimti=idimi*2
   mwvp2i=(mwavei+1)*2
   mdimhfi=mdimi/2
   mdimvhi=mdimvi/2
!
   mwavep=mwave+1
   mwave2=mwave+2
!
   npad=nscrch-ngwri
!
   allocate(   qi   (mdimi)                                                   ,&
               gzi  (mdimi)                                                   ,&
               tei  (mdimi,kdimi)                                             ,&
               dii  (mdimi,kdimi)                                             ,&
               zei  (mdimi,kdimi)                                             ,&
               rqi  (mdimi,kdimi)   )
   allocate(   q    (mdim)                                                    ,&
               gz   (mdim)                                                    ,&
               te   (mdim,kdimi)                                              ,&
               di   (mdim,kdimi)                                              ,&
               ze   (mdim,kdimi)                                              ,&
               rq   (mdim,kdimi)    )
   allocate(   ps   (idimt)                                                   ,&
               tf   (idimt,kdimi)                                             ,&
               rqf  (idimt,kdimi)                                             ,&
               uf   (idimt,kdimi)                                             ,&
               vf   (idimt,kdimi)   )
#ifdef DFS
#define MDIM mdimd
#else
#define MDIM mdim
#endif
   allocate(   qo    (MDIM)                                                   ,&
               gzo   (MDIM)                                                   ,&
               dio   (MDIM,kdim)                                              ,&
               zeo   (MDIM,kdim)                                              ,&
               teo   (MDIM,kdim)                                              ,&
               rqo   (MDIM,kdimq)   )
   allocate(   pso   (idimt)                                                  ,&
               tfo   (idimt,kdim)                                             ,&
               rqfo  (idimt,kdim)                                             ,&
               ufo   (idimt,kdim)                                             ,&
               vfo   (idimt,kdim)   )
#ifndef DFS
   allocate(   zss   (idimt)        )
#endif
   allocate(   pss   (idimt)                                                  ,&
               tts   (idimt,kdim)                                             ,&
               qqs   (idimt,kdimq)                                            ,&
               uus   (idimt,kdim)                                             ,&
               vvs   (idimt,kdim)   )
#ifndef DFS
   allocate(   zsa   (idimt)        )
#endif
   allocate(   psa   (idimt)                                                  ,&
               tta   (idimt,kdim)                                             ,&
               qqa   (idimt,kdimq)                                            ,&
               uua   (idimt,kdim)                                             ,&
               vva   (idimt,kdim)   )
   allocate(   siin  (kdimpi)                                                 ,&
               slin  (kdimi)                                                  ,&
               delin (kdimi)                                                  ,&
               ciin  (kdimpi)                                                 ,&
               clin  (kdimi)                                                  ,&
               rpiin (kdimmi)                                                 ,&
               ak    (kdimpi)                                                 ,&
               bk    (kdimpi)       )
   allocate(   si    (kdimp)                                                  ,&
               sl    (kdim )                                                  ,&
               del   (kdim )                                                  ,&
               ci    (kdimp)                                                  ,&
               cl    (kdim )                                                  ,&
               rpi   (kdimm)                                                  ,&
               ak5   (kdimp)                                                  ,&
               bk5   (kdimp)        )
   allocate(   eps   (mdimv )                                                 ,&
               colrad(jdimhf)                                                 ,&
               wgt   (jdimhf )                                                ,&
               wgtcs (jdimhf)                                                 ,&
               rcs2  (jdimhf)       )
   allocate(   indxnn(MDIM)                                                   ,&
               indxmm(MDIM)         )
   allocate(   d     (2,jcap1_,jcap1_)                                        ,&
               z     (2,jcap1_,jcap1_)                                        ,&
               u     (2,jcap1_,jcap2_)                                        ,&
               v     (2,jcap1_,jcap2_)  )
   allocate(   dxa   (lnt2_)                                                  ,&
               dxb   (lnt2_)                                                  ,&
               dxc   (lnt2_)                                                  ,&
               dxd   (lnt2_)        )
   allocate(   b     (MDIM)                                                   ,&
               c     (MDIM)                                                   ,&
               pad   (npad)         )
   allocate(   deps  (lnut2_)                                                 ,&
               rdeps (lnut2_)                                                 ,&
               dx    (twoj1_)                                                 ,&
               y     (jcap1_)       )
#undef MDIM
   end subroutine comchgr_init
!-------------------------------------------------------------------------------
   end module comchgr
