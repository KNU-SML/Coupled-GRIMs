#include <define.h>
   module comfcst
!-------------------------------------------------------------------------------
   use paramodel, only : lnt_, lnt2_, lnt22_, lnut2_                          ,&
                         lonf_, latg_, levs_, levh_                           ,&
                         LATG2S, LONF2S, LONF22S                              ,&
                         jcap_, jcap1_, twoj1_                                ,&
#ifdef RMP
                         IGRD12S, JGRD12S, igrd1_, igrd2_                     ,&
#endif
#ifdef MP
                         lonf2p_,latg2p_                                      ,&
#else
                         lonf2_,latg2_                                        ,&
#endif
#ifdef RIVER
                         io2_, jo2_                                           ,&
#endif
                         nstype_
   use rdparm
!-------------------------------------------------------------------------------
   integer, parameter  ::  ngrid = 22     ,&
                           nl    = 81     ,&
                           nlp1  = nl+1   ,&
#ifndef O3CHEM
                           jo3   = 18     ,&
                           ko3   = 46     ,&
#endif
#ifdef KF2
                           kfnt  = 250    ,&
                           kfnp  = 220    ,&
#endif
                           loz   = 17     ,&
                           mct   = 15     ,&  ! mct is a number of cloud types
                           kdum2 = 21     ,&
                           kens  = 2
   integer             ::  kdum
#ifndef DFS
   real                                    ::  dt
   real   , allocatable, dimension(:)      ::  gvdt                           ,&
                                               svdt
   real   , allocatable, dimension(:,:)    ::  amdt                           ,&
                                               bmdt
   real   , allocatable, dimension(:,:,:)  ::  dm
   real   , allocatable, dimension(:,:)    ::  spdlat
#ifdef CLM_CWF
#ifdef RMP
   real   , allocatable, dimension(:,:)    ::  cgs
#else
   real   , allocatable, dimension(:,:,:)  ::  cgs
#endif
#endif /* CLM_CWF end */
!
#define DEFAULT
#ifdef DCRFT
#undef DEFAULT
#ifndef RMP
   real   , allocatable, dimension(:)      ::  crscale                        ,&
                                               rcscale
#endif
#endif
#ifdef FFT_SGIMATH
#undef DEFAULT
   real   , allocatable, dimension(:)      ::  scale
   real   , allocatable, dimension(:,:)    ::  trig
#endif
#ifdef FFTW
#undef DEFAULT
   integer, allocatable, dimension(:)      ::  iplan_c_to_r                   ,&
                                               iplan_r_to_c
   real   , allocatable, dimension(:)      ::  scale
#endif
#ifdef RMP
#ifdef ASLES_OLD
#undef ASLES
#undef DEFAILT
   integer,              dimension(100)    ::  ifax_r
   real   , allocatable, dimension(:,:)    ::  trigs_r
#endif
#ifdef ASLES
#undef DEFAILT
   integer,              dimension(100)    ::  ifaxc
   real   , allocatable, dimension(:)      ::  trigsc
   integer,              dimension(100)    ::  ifaxs
   real   , allocatable, dimension(:)      ::  trigss
#endif
#endif /* RMP end */
#ifdef DEFAULT
#ifdef RMP
   integer,              dimension(100)    ::  ifax_r
   real   , allocatable, dimension(:,:)    ::  trigs_r
#endif
   integer, allocatable, dimension(:,:)    ::  ifax
   real   , allocatable, dimension(:,:,:)  ::  trigs
#endif
#endif /* ~DFS end */
!
#ifdef SMP
   real                                    ::  dtbdy                          ,&
                                               curtime                        ,&
                                               hour1                          ,&
                                               hour2                          ,&
                                               hours                          ,&
                                               houre                          ,&
                                               ftim1                          ,&
                                               ftim2
   real   , allocatable, dimension(:,:)    ::  vvel
#ifdef CLM_CWF
   real   , allocatable, dimension(:,:)    ::  wdiv                           ,&
                                               hadq
#endif
#ifdef SMP_NUDGING_FLX
   real   , allocatable, dimension(:)      ::  scmsh                          ,&
                                               scmlh
#endif
#ifdef SMP_NUDGING_RAD
   real   , allocatable, dimension(:)      ::  psfc                           ,&
                                               tsfc                           ,&
                                               qsfc                           ,&
                                               prec                           ,&
                                               scmlwup                        ,&
                                               scmlwdn                        ,&
                                               scmswup                        ,&
                                               scmswdn
#endif
   real   , allocatable, dimension(:)      ::  ps1, ps2
   real   , allocatable, dimension(:,:)    ::  uo1, uo2                       ,&
                                               vo1, vo2                       ,&
                                               to1, to2                       ,&
                                               qo1, qo2                       ,&
                                               wo1, wo2                       ,&
#ifdef REV_FRC
                                               ac1, ac2                       ,&
#endif
#ifdef RLX_FRC
                                               tscl1, tscl2                   ,&
#endif
                                               at1, at2                       ,&
                                               aq1, aq2                       ,&
                                               dv1, dv2
#endif  /* SMP end */
!
#ifdef EXPLICIT_CLOUDINESS
   real   , allocatable, dimension(:,:,:)  ::  qcicps                         ,&
                                               qrscps                         ,&
                                               taucld                         ,&
                                               cldwp                          ,&
                                               cldip
#endif
   real   , allocatable, dimension(:,:)    ::  dfk
   real   , allocatable, dimension(:,:)    ::  ktk
   real   , allocatable, dimension(:,:)    ::  dfkt
   real   , allocatable, dimension(:)      ::  b                              ,&
                                               satpsi                         ,&
                                               satkt                          ,&
                                               tsat
#ifdef SAS_DIAG
   real   , allocatable, dimension(:,:,:)  ::  dcu, dcv, dct, dcq, dch        ,&
                                               fcu, fcd, dlt, dlq, dlh
   real   , allocatable, dimension(:,:)    ::  deltb, delqb, delhb, cbmf
#else
#ifdef GWDC
   real   , allocatable, dimension(:,:,:)  ::  dct, fcu, fcd
#endif
#endif
!
#ifdef RASV2
   integer                                 ::  krmin                          ,&
                                               krmax                          ,&
                                               kfmax                          ,&
                                               ncrnd
   integer, allocatable, dimension(:)      ::  kctop
   real   , allocatable, dimension(:)      ::  sig                            ,&
                                               sgb                            ,&
                                               prj                            ,&
                                               rasal
   real   , allocatable, dimension(:,:)    ::  dsfc
#endif
#ifdef RIVER
   integer, allocatable,dimension(:)       ::  iindx1                         ,&
                                               iindx2                         ,&
                                               jindx1                         ,&
                                               jindx2
   real   , allocatable, dimension(:)      ::  ddx                            ,&
                                               ddy
#endif
!
!   *************************************************************!
!   seasonal climatologies of o3 (obtained on user vertical coord)
!   defined as 5 deg lat means n.p.->s.p.
!   *************************************************************!
!                                                   winter  spring
!                                                   summer  fall
!
   real   , allocatable, dimension(:,:), target ::  dduo3n, ddo3n2            ,&
                                                    ddo3n3, ddo3n4
!
!   *************************************************************!
!   seasonal climatologies of o3 on the detailed gfdl coordinate...
!   interpolation to each point pressure profile done in rad_ozone_gfdl
!   (see gloopr)
!   *************************************************************!
!                                                   winter  spring
!                                                   summer  fall
!
   real   , allocatable, dimension(:,:), target ::  xduo3n, xdo3n2            ,&
                                                    xdo3n3, xdo3n4
   real   , allocatable, dimension(:)           ::  prgfdl
!
   real                                    ::  degrad                         ,&
                                               hsigma                         ,&
                                               daysec                         ,&
                                               rco2
!
   integer                                 ::  kmaxvm
   integer, allocatable, dimension(:)      ::  ind                            ,&
                                               indx2                          ,&
                                               kmaxv                          ,&
                                               idummy2
   real   , allocatable, dimension(:,:)    ::  em1                            ,&
                                               em1wde                         ,&
                                               table1                         ,&
                                               table2                         ,&
                                               table3                         ,&
                                               em3                            ,&
                                               source                         ,&
                                               dsrce
!
!  downward sw fluxes from sw sib rad
!
   real   , allocatable, dimension(:,:)    ::  dfvbr                          ,&
                                               dfnbr                          ,&
                                               dfvdr                          ,&
                                               dfndr
   integer, allocatable, dimension(:)      ::  indxnn                         ,&
                                               indxmm
   real   , allocatable, dimension(:)      ::  dxa                            ,&
                                               dxb
   integer, allocatable, dimension(:)      ::  indxmv
   real   , allocatable, dimension(:)      ::  deps                           ,&
                                               rdeps                          ,&
                                               dx                             ,&
                                               y
   real   ,              dimension(loz)    ::  psnasa
   real   ,              dimension(37,loz) ::  o3nasa
#ifndef O3CHEM
   real   , allocatable, dimension(:,:)    ::  prdin                          ,&
                                               disin
#endif
#ifdef KF2
   real                                    ::  rdpr                           ,&
                                               rdthk                          ,&
                                               plutop
   real   , allocatable, dimension(:)      ::  the0k                          ,&
                                               alu
   real   , allocatable, dimension(:,:)    ::  ttab                           ,&
                                               qstab
#endif
!
   contains
!-------------------------------------------------------------------------------
   subroutine comfcst_init
!-------------------------------------------------------------------------------
   call rdparm_init
#ifndef HYBRID
   kdum=201-levs_-1-levs_
#else
   kdum=201-levs_-2-levs_
#endif
!
#ifndef DFS
   allocate(  gvdt   (levs_)                                                  ,&
              svdt   (levs_)                                                  ,&
              amdt   (levs_,levs_)                                            ,&
              bmdt   (levs_,levs_)                                            ,&
              dm     (levs_,levs_,0:jcap_)     )
   allocate(  spdlat (levs_,LATG2S)            )
#ifdef CLM_CWF
#ifdef RMP
   allocate(  cgs    (IGRD12S*JGRD12S,levs_)   )
#else
   allocate(  cgs    (LONF2S,LATG2S,levs_)     )
#endif
#endif
!
#ifdef REDUCE_GRID
#define LATS latg_/2
#else
#define LATS 1
#endif
!
#ifdef DCRFT
   allocate(  crscale (LATS)                                                  ,&
              rcscale (LATS)                   )
#endif
#ifdef FFT_SGIMATH
   allocate(  scale   (LATS)                                                  ,&
              trig    (lonf_+15,LATS)          )
#endif
#ifdef FFTW
   allocate(  scale   (LATS)                                                  ,&
              iplan_c_to_r(LATS)                                              ,&
              iplan_r_to_c(LATS)               )
#endif
#ifdef RMP
#ifdef ASLES_OLD
   allocate(   trigs_r (igrd2_,2)              )
#endif         
#ifdef ASLES   
   allocate(   trigsc  (igrd1_*3)              )
   allocate(   trigss  ((igrd1_*7)/2+1)        )
#endif         
#endif /* RMP end */
#ifdef DEFAULT
#ifdef RMP
   allocate(  trigs_r (igrd2_,2)               )
#endif
   allocate(  ifax    (100,LATS)                                              ,&
              trigs   (lonf_,2,LATS)           )
#endif
#endif /* not DFS end */
!
#ifdef SMP
   allocate(  vvel    (lnt22_,levs_)           )
#ifdef CLM_CWF
   allocate(  wdiv    (lnt22_,levs_)                                          ,&
              hadq    (lnt22_,levs_)           )
#endif
#ifdef SMP_NUDGING_FLX
   allocate(  scmsh   (lnt22_)                                                ,&
              scmlh   (lnt22_)                 )
#endif
#ifdef SMP_NUDGING_RAD
   allocate(  psfc    (lnt22_)                                                ,&
              tsfc    (lnt22_)                                                ,&
              qsfc    (lnt22_)                                                ,&
              prec    (lnt22_)                                                ,&
              scmlwup (lnt22_)                                                ,&
              scmlwdn (lnt22_)                                                ,&
              scmswup (lnt22_)                                                ,&
              scmswdn (lnt22_)                 )
#endif
   allocate(  ps1(lnt22_)      , ps2(lnt22_)         )
   allocate(  uo1(lnt22_,levs_), uo2(lnt22_,levs_)                            ,&
              vo1(lnt22_,levs_), vo2(lnt22_,levs_)                            ,&
              to1(lnt22_,levs_), to2(lnt22_,levs_)                            ,&
              qo1(lnt22_,levh_), qo2(lnt22_,levh_)                            ,&
              wo1(lnt22_,levs_), wo2(lnt22_,levs_)                            ,&
#ifdef REV_FRC
              ac1(lnt22_,levs_), ac2(lnt22_,levs_)                            ,&
#endif
#ifdef RLX_FRC
              tscl1(lnt22_,levs_), tscl2(lnt22_,levs_)                        ,&
#endif
              at1(lnt22_,levs_), at2(lnt22_,levs_)                            ,&
              aq1(lnt22_,levh_), aq2(lnt22_,levh_)                            ,&
              dv1(lnt22_,levs_), dv2(lnt22_,levs_)   ) 
#endif  /* SMP end */
!
#ifdef EXPLICIT_CLOUDINESS
   allocate(  qcicps  (LONF2S,levs_,LATG2S)                                   ,&
              qrscps  (LONF2S,levs_,LATG2S)                                   ,&
              taucld  (LONF2S,levs_,LATG2S)                                   ,&
              cldwp   (LONF2S,levs_,LATG2S)                                   ,&
              cldip   (LONF2S,levs_,LATG2S)          )
#endif
   allocate(  dfk     (ngrid,nstype_)                )
   allocate(  ktk     (ngrid,nstype_)                )
   allocate(  dfkt    (ngrid,nstype_)                )
   allocate(  b       (nstype_)                                               ,&
              satpsi  (nstype_)                                               ,&
              satkt   (nstype_)                                               ,&
              tsat    (nstype_)                      )
#ifdef SAS_DIAG
   allocate(  dcu     (LONF2S,levs_,LATG2S)                                   ,&
              dcv     (LONF2S,levs_,LATG2S)                                   ,&
              dct     (LONF2S,levs_,LATG2S)                                   ,&
              dcq     (LONF2S,levs_,LATG2S)                                   ,&
              dch     (LONF2S,levs_,LATG2S)                                   ,&
              fcu     (LONF2S,levs_,LATG2S)                                   ,&
              fcd     (LONF2S,levs_,LATG2S)                                   ,&
              dlt     (LONF2S,levs_,LATG2S)                                   ,&
              dlq     (LONF2S,levs_,LATG2S)                                   ,&
              dlh     (LONF2S,levs_,LATG2S)                                   ,&
              deltb   (LONF2S,LATG2S)                                         ,&
              delqb   (LONF2S,LATG2S)                                         ,&
              delhb   (LONF2S,LATG2S)                                         ,&
              cbmf    (LONF2S,LATG2S)                )
#else
#ifdef GWDC
   allocate(  dct     (LONF2S,levs_,LATG2S)                                   ,&
              fcu     (LONF2S,levs_,LATG2S)                                   ,&
              fcd     (LONF2S,levs_,LATG2S)          )
#endif
#endif
#ifdef RASV2
   allocate(  kctop   (mct+1)                        )
   allocate(  sig     (levs_+1)                                               ,&
              sgb     (levs_)                                                 ,&
              prj     (levs_+1)                                               ,&
              rasal   (levs_)                                                 ,&
              dsfc    (LONF2S,LATG2S)                )
#endif
#ifdef RIVER
   allocate(  iindx1  (io2_)                                                  ,&
              iindx2  (io2_)                                                  ,&
              jindx1  (jo2_)                                                  ,&
              jindx2  (jo2_)                                                  ,&
              ddx     (io2_)                                                  ,&
              ddy     (jo2_)                         )
#endif
   allocate(  dduo3n  (37,l)                                                  ,&
              ddo3n2  (37,l)                                                  ,&
              ddo3n3  (37,l)                                                  ,&
              ddo3n4  (37,l)                         )
   allocate(  xduo3n  (37,nl)                                                 ,&
              xdo3n2  (37,nl)                                                 ,&
              xdo3n3  (37,nl)                                                 ,&
              xdo3n4  (37,nl)                                                 ,&
              prgfdl  (nl)                           )
   allocate(  ind     (imax)                                                  ,&
              indx2   (lp1v)                                                  ,&
              kmaxv   (lp1)                                                   ,&
              idummy2 (imax+lp1v+lp1+1)                                       ,&
              em1     (28,180)                                                ,&
              em1wde  (28,180)                                                ,&
              table1  (28,180)                                                ,&
              table2  (28,180)                                                ,&
              table3  (28,180)                                                ,&
              em3     (28,180)                                                ,&
              source  (28,nbly)                                               ,&
              dsrce   (28,nbly)                      )
#ifdef DG3
#ifdef RMPVECTORIZE
   
   allocate(  dfvbr   (LONF2S,1)                                              ,&
              dfnbr   (LONF2S,1)                                              ,&
              dfvdr   (LONF2S,1)                                              ,&
              dfndr   (LONF2S,1)                     )
#else
   allocate(  dfvbr   (LONF2S,LATG2S)                                         ,&
              dfnbr   (LONF2S,LATG2S)                                         ,&
              dfvdr   (LONF2S,LATG2S)                                         ,&
              dfndr   (LONF2S,LATG2S)                )
#endif
#endif
   allocate(  indxnn (lnt2_)                                                  ,&
              indxmm (lnt2_)                   )
   allocate(  dxa    (lnt2_)                                                  ,&
              dxb    (lnt2_)                   )
   allocate(  deps   (lnut2_)                                                 ,&
              rdeps  (lnut2_)                                                 ,&
              dx     (twoj1_)                                                 ,&
              y      (jcap1_)                                                 ,&
              indxmv (lnut2_)                  )
#ifndef O3CHEM
   allocate(  prdin  (jo3,ko3)                                                ,&
              disin  (jo3,ko3)                 )
#endif
#ifdef KF2
   allocate(  the0k  (1:kfnp)                                                 ,&
              alu    (1:200)                   )
   allocate(  ttab   (1:kfnt,1:kfnp)                                          ,&
              qstab  (1:kfnt,1:kfnp)           )
#endif
   end subroutine comfcst_init
!-------------------------------------------------------------------------------
   end module comfcst
