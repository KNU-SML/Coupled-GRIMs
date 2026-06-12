#include <define.h>
   module rscomrsml
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_, levs_, levh_, levp1_                          ,&
                         lonf2_, latg2_                                       ,&
                         LONF2S                                               ,&
                         lnt2_, lnt22_, lnut2_                                ,&
                         jcap1_ , jcap2_, twoj1_                              ,&
                         igrd1_ , igrd2_                                      ,&
                                  jgrd12_                                     ,&
                         IGRD12S, JGRD12S                                     ,&
                         lngrd_ , LNGRDS                                      ,&
                         lnwav_ , LNWAVS
#ifdef MP
   use paramodel, only : levsp_, levp1p_, llwavp_
#endif
   use rsparltb, only  : rsparltb_init, lngrdb
!-------------------------------------------------------------------------------
   integer, allocatable, dimension(:)      ::  nc2b00                         ,&
                                               nc2b10                         ,&
                                               nc2b11                         ,&
                                               nc2b01
   real   , allocatable, dimension(:)      ::  dc2b00                         ,&
                                               dc2b10                         ,&
                                               dc2b11                         ,&
                                               dc2b01
   integer, allocatable, dimension(:)      ::  ncb11
   real   , allocatable, dimension(:)      ::  xcb0                           ,&
                                               xcb1                           ,&
                                               xcb2                           ,&
                                               xcb3                           ,&
                                               ycb0                           ,&
                                               ycb1                           ,&
                                               ycb2                           ,&
                                               ycb3
   integer, allocatable, dimension(:)      ::  nc2r00                         ,&
                                               nc2r10                         ,&
                                               nc2r11                         ,&
                                               nc2r01
   real   , allocatable, dimension(:)      ::  dc2r00                         ,&
                                               dc2r10                         ,&
                                               dc2r11                         ,&
                                               dc2r01
   integer, allocatable, dimension(:)      ::  ncr11
   real   , allocatable, dimension(:)      ::  xcr0                           ,&
                                               xcr1                           ,&
                                               xcr2                           ,&
                                               xcr3                           ,&
                                               ycr0                           ,&
                                               ycr1                           ,&
                                               ycr2                           ,&
                                               ycr3
   integer, allocatable, dimension(:)      ::  ng2b00                         ,&
                                               ng2b10                         ,&
                                               ng2b11                         ,&
                                               ng2b01
   real   , allocatable, dimension(:)      ::  dg2b00                         ,&
                                               dg2b10                         ,&
                                               dg2b11                         ,&
                                               dg2b01
   integer, allocatable, dimension(:)      ::  ngb11
   real   , allocatable, dimension(:)      ::  xgb0                           ,&
                                               xgb1                           ,&
                                               xgb2                           ,&
                                               xgb3                           ,&
                                               ygb0                           ,&
                                               ygb1                           ,&
                                               ygb2                           ,&
                                               ygb3
   integer, allocatable, dimension(:)      ::  ng2r00                         ,&
                                               ng2r10                         ,&
                                               ng2r11                         ,&
                                               ng2r01
   real   , allocatable, dimension(:)      ::  dg2r00                         ,&
                                               dg2r10                         ,&
                                               dg2r11                         ,&
                                               dg2r01
   integer, allocatable, dimension(:)      ::  ngr11
   real   , allocatable, dimension(:)      ::  xgr0                           ,&
                                               xgr1                           ,&
                                               xgr2                           ,&
                                               xgr3                           ,&
                                               ygr0                           ,&
                                               ygr1                           ,&
                                               ygr2                           ,&
                                               ygr3
   real   , allocatable, dimension(:)      ::  gm                             ,&
                                               xm                             ,&
                                               csln                           ,&
                                               snln                           ,&
                                               colrad                         ,&
                                               wgt                            ,&
                                               rcslb                          ,&
                                               rsnlb                          ,&
                                               rlatb                          ,&
                                               rlonb
   real   , allocatable, dimension(:,:)    ::  epsi
   integer, allocatable, dimension(:)      ::  indxmv
   real   , allocatable, dimension(:)      ::  deps                           ,&
                                               rdeps                          ,&
                                               dx                             ,&
                                               y
   real   , allocatable, dimension(:)      ::  dxa                            ,&
                                               dxb
   real   , allocatable, dimension(:)      ::  sm
   real   , allocatable, dimension(:,:)    ::  spdlat
   real                                    ::  dthour                         ,&
                                               dshour                         ,&
                                               totsum                         ,&
                                               dsolhr
   real   , allocatable, dimension(:)      ::  qs
   real   , allocatable, dimension(:,:)    ::  tes                            ,&
                                               rqs                            ,&
                                               uus                            ,&
                                               vvs
#ifdef NONHYD
   real                                    ::  dtx2
   real   , allocatable, dimension(:)      ::  tor
   real   , allocatable, dimension(:,:)    ::  wcm                            ,&
                                               x3minv                         ,&
                                               x4m
   real   , allocatable, dimension(:,:,:)  ::  dcm
#endif
#ifdef SAS_DIAG
   real   , allocatable, dimension(:,:)    ::  deltb                          ,&
                                               delqb                          ,&
                                               delhb                          ,&
                                               cbmf
   real   , allocatable, dimension(:,:,:)  ::  dcu                            ,&
                                               dcv                            ,&
                                               dct                            ,&
                                               dcq                            ,&
                                               dch                            ,&
                                               fcu                            ,&
                                               fcd                            ,&
                                               dlt                            ,&
                                               dlq                            ,&
                                               dlh
#else
#ifdef GWDC 
   real   , allocatable, dimension(:,:,:)  ::  dct                            ,&
                                               fcu                            ,&
                                               fcd
#endif /* GWDC end */
#endif /* SAS_DIAG end */
!
#ifdef NONHYD
   real   , allocatable, dimension(:,:)    ::  pnx                            ,&
                                               tnx                            ,&
                                               onx                            ,&
                                               pny                            ,&
                                               tny                            ,&
                                               ony
#ifdef MP
   real   , allocatable, dimension(:,:)    ::  pnxk                           ,&
                                               tnxk                           ,&
                                               onxk                           ,&
                                               pnyk                           ,&
                                               tnyk                           ,&
                                               onyk
#endif
#endif /* NONHYD end */
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscomrsml_init
!-------------------------------------------------------------------------------
   call rsparltb_init
   allocate(   nc2b00  (lngrdb)                                               ,&
               nc2b10  (lngrdb)                                               ,&
               nc2b11  (lngrdb)                                               ,&
               nc2b01  (lngrdb)                                               ,&
               dc2b00  (lngrdb)                                               ,&
               dc2b10  (lngrdb)                                               ,&
               dc2b11  (lngrdb)                                               ,&
               dc2b01  (lngrdb)                                                )
   allocate(   ncb11   (lngrdb)                                               ,&
               xcb0    (lngrdb)                                               ,&
               xcb1    (lngrdb)                                               ,&
               xcb2    (lngrdb)                                               ,&
               xcb3    (lngrdb)                                               ,&
               ycb0    (lngrdb)                                               ,&
               ycb1    (lngrdb)                                               ,&
               ycb2    (lngrdb)                                               ,&
               ycb3    (lngrdb)                                                )
   allocate(   nc2r00  (lngrd_)                                               ,&
               nc2r10  (lngrd_)                                               ,&
               nc2r11  (lngrd_)                                               ,&
               nc2r01  (lngrd_)                                               ,&
               dc2r00  (lngrd_)                                               ,&
               dc2r10  (lngrd_)                                               ,&
               dc2r11  (lngrd_)                                               ,&
               dc2r01  (lngrd_)                                                )
   allocate(   ncr11   (lngrd_)                                               ,&
               xcr0    (lngrd_)                                               ,&
               xcr1    (lngrd_)                                               ,&
               xcr2    (lngrd_)                                               ,&
               xcr3    (lngrd_)                                               ,&
               ycr0    (lngrd_)                                               ,&
               ycr1    (lngrd_)                                               ,&
               ycr2    (lngrd_)                                               ,&
               ycr3    (lngrd_)                                                )
   allocate(   ng2b00  (lngrdb)                                               ,&
               ng2b10  (lngrdb)                                               ,&
               ng2b11  (lngrdb)                                               ,&
               ng2b01  (lngrdb)                                               ,&
               dg2b00  (lngrdb)                                               ,&
               dg2b10  (lngrdb)                                               ,&
               dg2b11  (lngrdb)                                               ,&
               dg2b01  (lngrdb)                                                )
   allocate(   ngb11   (lngrdb)                                               ,&
               xgb0    (lngrdb)                                               ,&
               xgb1    (lngrdb)                                               ,&
               xgb2    (lngrdb)                                               ,&
               xgb3    (lngrdb)                                               ,&
               ygb0    (lngrdb)                                               ,&
               ygb1    (lngrdb)                                               ,&
               ygb2    (lngrdb)                                               ,&
               ygb3    (lngrdb)                                                )
   allocate(   ng2r00  (LNGRDS)                                               ,&
               ng2r10  (LNGRDS)                                               ,&
               ng2r11  (LNGRDS)                                               ,&
               ng2r01  (LNGRDS)                                               ,&
               dg2r00  (LNGRDS)                                               ,&
               dg2r10  (LNGRDS)                                               ,&
               dg2r11  (LNGRDS)                                               ,&
               dg2r01  (LNGRDS)                                                )
   allocate(   ngr11   (lngrd_)                                               ,&
               xgr0    (lngrd_)                                               ,&
               xgr1    (lngrd_)                                               ,&
               xgr2    (lngrd_)                                               ,&
               xgr3    (lngrd_)                                               ,&
               ygr0    (lngrd_)                                               ,&
               ygr1    (lngrd_)                                               ,&
               ygr2    (lngrd_)                                               ,&
               ygr3    (lngrd_)                                                )
   allocate(   gm      (lnt22_)                                               ,&
               xm      (lngrdb)                                               ,&
               csln    (lonf2_)                                               ,&
               snln    (lonf2_)                                               ,&
               colrad  (latg2_)                                               ,&
               wgt     (latg2_)                                               ,&
               rcslb   (lngrdb)                                               ,&
               rsnlb   (lngrdb)                                               ,&
               rlatb   (lngrdb)                                               ,&
               rlonb   (lngrdb)                                               ,&
               epsi    (jcap2_,jcap1_)                                         )
   allocate(   indxmv  (lnut2_)                                               ,&
               deps    (lnut2_)                                               ,&
               rdeps   (lnut2_)                                               ,&
               dx      (twoj1_)                                               ,&
               y       (jcap1_)                                                )
   allocate(   dxa     (lnt2_)                                                ,&
               dxb     (lnt2_)                                                 )
   allocate(   sm      (LNWAVS)                                                )
#ifdef RMPVECTORIZE
   allocate(   spdlat  (levs_,1)                                               )
#else
   allocate(   spdlat  (levs_,JGRD12S)                                         )
#endif
   allocate(   qs      (lnwav_)                                                )
   allocate(   tes     (lnwav_,levs_)                                         ,&
               rqs     (lnwav_,levh_)                                         ,&
               uus     (lnwav_,levs_)                                         ,&
               vvs     (lnwav_,levs_)                                          )
#ifdef NONHYD
   allocate(   tor     (levs_)                                                 )
   allocate(   wcm     (levs_,levp1_)                                         ,& 
               x3minv  (levp1_,levp1_)                                        ,&
               x4m     (levp1_,levs_ )                                         )
   allocate(   dcm     (LNWAVS,levs_,levs_)                                    )
#endif /* NONHYD end */
#ifdef SAS_DIAG
   allocate(   deltb   (LONF2S,jgrd12_)                                       ,&
               delqb   (LONF2S,jgrd12_)                                       ,&
               delhb   (LONF2S,jgrd12_)                                       ,&
               cbmf    (LONF2S,jgrd12_)                                        )
   allocate(   dcu     (LONF2S,levs_,jgrd12_)                                 ,&
               dcv     (LONF2S,levs_,jgrd12_)                                 ,&
               dct     (LONF2S,levs_,jgrd12_)                                 ,&
               dcq     (LONF2S,levs_,jgrd12_)                                 ,&
               dch     (LONF2S,levs_,jgrd12_)                                 ,&
               fcu     (LONF2S,levs_,jgrd12_)                                 ,&
               fcd     (LONF2S,levs_,jgrd12_)                                 ,&
               dlt     (LONF2S,levs_,jgrd12_)                                 ,&
               dlq     (LONF2S,levs_,jgrd12_)                                 ,&
               dlh     (LONF2S,levs_,jgrd12_)                                  )
#else
#ifdef GWDC 
   allocate(   dct     (LONF2S,levs_,jgrd12_)                                 ,&
               fcu     (LONF2S,levs_,jgrd12_)                                 ,&
               fcd     (LONF2S,levs_,jgrd12_)                                  )
#endif /* GWDC end */
#endif /* SAS_DIAG end */
!
#ifdef NONHYD
   allocate(   pnx     (LNWAVS,levs_)                                         ,&
               tnx     (LNWAVS,levs_)                                         ,&
               onx     (LNWAVS,levp1_)                                        ,&
               pny     (LNWAVS,levs_)                                         ,&
               tny     (LNWAVS,levs_)                                         ,&
               ony     (LNWAVS,levp1_)                                         )
#ifdef MP
   allocate(   pnxk    (llwavp_,levsp_)                                       ,&
               tnxk    (llwavp_,levsp_)                                       ,&
               onxk    (llwavp_,levp1p_)                                      ,&
               pnyk    (llwavp_,levsp_)                                       ,&
               tnyk    (llwavp_,levsp_)                                       ,&
               onyk    (llwavp_,levp1p_)                                       )
#endif
#endif /* NONHYD end */
! 
   end subroutine rscomrsml_init
!-------------------------------------------------------------------------------
   end module rscomrsml
