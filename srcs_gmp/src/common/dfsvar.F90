#include "define.h"
   module dfsvar
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   logical                                              ::  iope
   integer                                              ::  ntrac,ncldg,ntotal
! ------------------------------------------------------------------------------
! global index
! ------------------------------------------------------------------------------
   integer                                              ::  levs,levh         ,&
                                                            tlevs             ,&
                                                            tld  ,tlv         ,&
                                                            tlt  ,tlq         ,&
                                                            mta  ,jla         ,&
                                                            mtg  ,jlg         ,&
                                                            mtha ,jlha        ,&
                                                            mt1a              ,&
                                                            lnt2a,lnt22a      ,&
                                                            iba  ,jba  ,jbwa  ,&
                                                            ib1a ,jb1a ,jbw1a ,&
                                                            kgs  ,kge         ,&
                                                            mgs  ,mge         ,&
                                                            ngs  ,nge         ,&
                                                            igs  ,ige         ,&
                                                            jgs  ,jge
! ------------------------------------------------------------------------------
! local index
! ------------------------------------------------------------------------------
   integer                                              ::  levsp ,levhp      ,&
#ifdef MP
                                                            kler              ,&
#endif
                                                            tlevsp            ,&
                                                            tldp,tlvp         ,&
                                                            tltp,tlqp         ,&
                                                            mt   ,jl          ,&
                                                            mth  ,jlh ,jlhp   ,&
                                                            mt1               ,&
                                                            lnt2 ,lnt22       ,&
                                                            ib   ,jb  ,jbw    ,&
                                                            ib1  ,jb1 ,jbw1   ,&
                                                            kls  ,kle         ,&
                                                            mls  ,mle         ,&
                                                            nls  ,nle         ,&
                                                            ils  ,ile         ,&
                                                            jls  ,jle         ,&
                                                            mlsm,mlem,mlm,ml
!-------------------------------------------------------------------------------
   integer, allocatable         , dimension(:)          ::  latdef
#ifdef MP
   integer, allocatable, target , dimension(:,:)        ::  midxr
   integer             , pointer, dimension(:)          ::  midx
#else
   integer, allocatable         , dimension(:)          ::  midx
#endif
   integer                                              ::  nod,nev
   integer, allocatable         , dimension(:)          ::  mevs,mods,mop
#ifdef ALIASED
   integer, allocatable         , dimension(:)          ::  nmxa
#endif
!
   character(len=20)                                    ::  VERTCORD
   integer                                              ::  lap_dim
   real                                                 ::  PI2,PI180         ,&
                                                            aOmega,aOmega2    ,&
                                                            DISSI,GENER       ,&
                                                            cooling,press     ,&
                                                            damp1,damp2  
   real                         , dimension(3)          ::  gamma_v           ,&
                                                            gamma_d           ,&
                                                            gamma_t           ,&
                                                            gamma_q
   integer, allocatable         , dimension(:,:)        ::  JCOL2JS
   real   , allocatable         , dimension(:)          ::  aMSQUAR           ,&
                                                            SINLAT,TANLAT
   real   , allocatable         , dimension(:,:)        ::  COSLAT
   real   , allocatable         , dimension(:,:,:,:)    ::  AMATm,AMATm_      ,&
                                                            DMATm,DMATm_
   real   , allocatable         , dimension(:,:,:,:,:)  ::  DMATmi_
   real   , allocatable, target , dimension(:,:,:,:,:)  ::  ZDMATmi_a_RE_v    ,&
                                                            ZDMATmi_a_RE_d    ,&
                                                            ZDMATmi_a_RE_t    ,&
                                                            ZDMATmi_a_RE_q    ,&
                                                            ZDMATmi_a_IM_v    ,&
                                                            ZDMATmi_a_IM_d    ,&
                                                            ZDMATmi_a_IM_t    ,&
                                                            ZDMATmi_a_IM_q
!-------------------------------------------------------------------------------
! for 3d-model 
!-------------------------------------------------------------------------------
   real   , allocatable         , dimension(:)          ::  sigma             ,&
                                                            sigmafull         ,&
                                                            tsigma            ,&
                                                            delsig            ,&
                                                            valpha, beta      ,&
                                                            gama1 , gama2     ,&
                                                            TAVEXY, FAVEXY
   real   , allocatable         , dimension(:)          ::  GV
   real   , allocatable         , dimension(:)          ::  PHS
   real   , allocatable         , dimension(:,:)        ::  TRS,GRS,GRSINV
   real   , allocatable         , dimension(:,:)        ::  PHSDM
   real   , allocatable         , dimension(:,:,:)      ::  TRSDM
   real   , allocatable         , dimension(:,:)        ::  DIVN
   real   , allocatable         , dimension(:,:,:)      ::  sigdot,APSK       ,&
                                                            wR    ,dwR        ,&
                                                            ADT1
#ifdef HDSZ
   real   , allocatable         , dimension(:)          ::  Z_VIS
   real   , allocatable         , dimension(:,:)        ::  T_VIS
   real   , allocatable         , dimension(:,:,:)      ::  T_REF
#endif
   real   , allocatable         , dimension(:)          ::  EIGVAL
   real   , allocatable         , dimension(:,:)        ::  VERMAT,EIGVEC     ,&
                                                            AEIGVEC
!
! spectral variables
!
   real   , allocatable, target , dimension(:,:,:)      ::  V1  ,V2  ,V3  ,VJ
   real                , pointer, dimension(:,:,:)      ::  d1  ,d2  ,d3  ,dj ,&
#ifndef NISLQ
                                                            q1  ,q2  ,q3  ,qj ,&
#endif
                                                            t1  ,t2  ,t3  ,tj
   real                , pointer, dimension(:,:)        ::  PSL1,SF1          ,&
                                                            PSL2, XPS2, LAP2  ,&
                                                            PSL3, XPS3, LAP3  ,&
                                                            PSJ
#ifdef MP
   real   , allocatable         , dimension(:,:)        ::  SF1l
#endif
#ifdef PERT
!
! for pertubation
!
   real   , allocatable, target , dimension(:,:,:)      ::  V0
   real                , pointer, dimension(:,:,:)      ::  D0                ,&
                                                            T0                ,&
#ifndef NISLQ
                                                            Q0                ,&
#endif
                                                            PSL0
#endif
!
! grid variables
!
   real   , allocatable         , dimension(:,:)        ::  xlon,xlat
   real   , allocatable, target , dimension(:,:,:)      ::  VOR
   real                , pointer, dimension(:,:,:)      ::  dio               ,&
#ifndef NISLQ
                                                            qai               ,&
#endif
                                                            tai
   real                , pointer, dimension(:,:)        ::  prs               ,&
                                                            xpsl              ,&
                                                            ypsl              ,&
                                                            plap
!
! uv wind on grid space
!
   real   , allocatable, target , dimension(:,:,:)      ::  VXC
   real                , pointer, dimension(:,:,:)      ::  VYC
   real   , allocatable         , dimension(:)          ::  spdmax
#ifdef CLM_CWF
   real   , allocatable         , dimension(:,:,:)      ::  cgs
#endif
#ifdef BARO_TEST
   real   , allocatable         , dimension(:,:)        ::  SFAI
#endif 
   real   , allocatable         , dimension(:,:)        ::  apsn
#ifdef NISLQ
   real   , allocatable         , dimension(:,:,:)      ::  q1,q2,q3,qj,qai
#ifdef PERT
   real   , allocatable         , dimension(:,:,:)      ::  q0
#endif
#endif /* NISLQ end */
!
   contains
!-------------------------------------------------------------------------------
   subroutine get_dfs_dim(jcap,ilevs,ntr,ncld,isiz,jsiz)
!-------------------------------------------------------------------------------
#ifdef MP
   use commpi       ! mype,master,lmlen,lmstr,lnlen,lnstr,                     &
                    ! lonlen,lonstr,latlen,latstr,levstr,levlen,               &
                    ! ncol,nrow,myrow,lerlen,lerstr
!-------------------------------------------------------------------------------
#endif
   implicit none
!-------------------------------------------------------------------------------
   integer                                            ::  jcap,ilevs,ntr      ,&
                                                          ncld,isiz,jsiz      ,&
                                                          i,ip0,nr,mx,mh      ,&
                                                          me,ms,jj,j,m
! -----------------------------------------------------------------------------
! set ntotal=ntrac+ncld
! -----------------------------------------------------------------------------
   ntrac=ntr
   ncldg=ncld
   ntotal=ntr+ncld
! -----------------------------------------------------------------------------
! vertical layer     
! -----------------------------------------------------------------------------
   levs=ilevs
   levh=levs*ntotal
#ifdef NISLQ
   tlevs=3*levs        ! total layers for prognostic variables
#else
   tlevs=3*levs+levh   ! total layers for prognostic variables
#endif
   tld=1
   tlv=tld+levs
   tlt=tlv+levs
#ifndef NISLQ
   tlq=tlt+levs
#endif
! -----------------------------------------------------------------------------
! East-West wave size
! -----------------------------------------------------------------------------
   mta=jcap
   mtg=2*mta+1            ! -MTA:MTA
   mt1a=mta+1
   mtha=mta/2
! -----------------------------------------------------------------------------
! South-North wave size
! -----------------------------------------------------------------------------
#ifdef ALIASED
#ifdef ALIASED2
   jla=jsiz-2
#else
   jla=3*(mta+1)/2
   jla=jla-mod(jla,2)
#endif
#else
   jla=mta
#endif /* ALIASED end */
   jlg=jla+1            ! 0:jla
   jlha=jla/2
! ------------------------------------------------------------------------------
! total wave number
! ------------------------------------------------------------------------------
   lnt2a=mtg*jlg
   lnt22a=lnt2a
! ------------------------------------------------------------------------------
! I size
! ------------------------------------------------------------------------------
   iba=isiz
   ib1a=iba-1
! ------------------------------------------------------------------------------
! J size
! ------------------------------------------------------------------------------
   jbwa=jsiz
   jbw1a=jbwa-1
   jba=jbwa/2
   jb1a=jba-1
! ------------------------------------------------------------------------------
! set global start & end index
! ------------------------------------------------------------------------------
   mgs=-mta
   mge=mta
   ngs=0
   nge=jla
   kgs=1
   kge=levs
   igs=1
   ige=iba
   jgs=1
   jge=jbwa
! ------------------------------------------------------------------------------
! set local domain size
! ------------------------------------------------------------------------------
   allocate(latdef(jba))
   do j = jgs,jgs+jba-1
      latdef(j)=j
   enddo
!
#ifdef MP
   call mpdimset(mta,jla,levs,iba,jbwa,igs,jgs,0,ngs,latdef)
   if (mype.eq.master) then
     iope=.true.
   else
     iope=.false.
   endif
!
! vertical
!
   levsp=levlen(mype)    ! levs/ncol
   kler=lerlen(mype)    ! levs/nrow
   kls=levstr(mype)
   kle=kls+levsp-1
!
! east-west
!
   mt=lmlen(mype) 
   mth=mt/2
   mt1=mt-1
   mls=lmstr(mype)
   if (mls.eq.0) then
     ml=mth+1    ! odd
   else
     ml=mth
   endif
   mle=mls+ml-1
   mlsm=-mle        ! native start
   mlm=mth
   mlem=mlsm+mlm-1  ! native end
!
! south-north
!
   jl=lnlen(mype)
   nls=lnstr(mype)
   nle=nls+jl-1
!
! wave number indice
!
   allocate(midxr(maxval(lmlen),0:nrow-1))
   do nr = 0,nrow-1
     ip0=nr*ncol
     mx=lmlen(ip0)
     mh=mx/2
     ms=lmstr(ip0)
     me=ms+mh-1
     if (ms.eq.0) me=me+1
     do i = 1,mh
       midxr(i,nr)=-(me-i+1)
     enddo
     do i = mh+1,mx
       midxr(i,nr)=ms+i-mh-1
     enddo
#ifdef DBG
     if (mype.eq.master) then
       write(6,*)'nr=',nr
       write(6,*)'midxr=',midxr(:,nr)
     endif
#endif
   enddo
!
   midx=>midxr(1:mt,myrow)
!
! grid size
!
   ib=lonlen(mype)
   ib1=ib-1
   ils=lonstr(mype)
   ile=ils+ib-1
   jbw=latlen(mype)*2
   jbw1=jbw-1
   jb=jbw/2
   jb1=jb-1
   jls=latstr(mype)
   jle=jls+jb1
#else /* ~MP */
   iope=.true.
!
! vertical
!
   levsp=levs
   kls=kgs
   kle=kge
!
! east-west
!
   mt=mtg           ! -MT:MT
   mth=mtha
   mt1=mt-1
   mls=0
   mle=mta
   ml=mta+1
   mlsm=-mta
   mlem=-1
   mlm=mta
!
! south-north
!
   jl=jla+1         ! 0:JLA
   nls=ngs
   nle=nge
!
! wave number indice
!
   allocate(midx(mt))
   do i = 1,mt
     midx(i)=-mta+i-1
   enddo
!
! grid size
!
   ib=iba
   ib1=ib1a
   ils=igs
   ile=ige
   jbw=jbwa
   jbw1=jbw1a
   jb=jba
   jb1=jb1a
   jls=jgs
   jle=jls+jb1
#endif
   levhp=levsp*ntotal
#ifdef NISLQ
   tlevsp=3*levsp         ! total layers for prognostic variables
#else
   tlevsp=3*levsp+levhp   ! total layers for prognostic variables
#endif
   tldp=1
   tlvp=tldp+levsp
   tltp=tlvp+levsp
#ifndef NISLQ
   tlqp=tltp+levsp
#endif
!
   lnt2=mt*jlg
   lnt22=lnt2
   jlh=(nle-nls)/2+1
   jlhp=jlh+1
! ------------------------------------------------------------------------------
! index for even and odd wavenumber
! ------------------------------------------------------------------------------
   nev=0; nod=0
   allocate(mop(mt))
   allocate(mods(mt/2+1))
   allocate(mevs(mt/2+1))
   do i = 1,mt
     mop(i)=abs(mt-i+1)
     if (mod(midx(i),2).eq.0) then
       nev=nev+1
       mevs(nev)=i
     else
       nod=nod+1
       mods(nod)=i
     endif
   enddo
#ifdef ALIASED
   allocate(nmxa(mt))
   do i = 1,mt
     m=abs(midx(i))
     mx=max(mta,jla-(m*jla)/(3*mta))+1
     if (mx.gt.nle) then
       nmxa(i)=-1
     else
       nmxa(i)=mx
     endif
   enddo
#endif
   return
   end subroutine get_dfs_dim
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine dfsvarini(jcap,ilevs,ntr,ncld,isiz,jsiz)
!-------------------------------------------------------------------------------
#ifdef MP
   use commpi, only : mype
#endif
   use constant, only : pi=>pi_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer          :: jcap,ilevs,ntr,ncld,isiz,jsiz
   integer          :: lapend
!
   call get_dfs_dim(jcap,ilevs,ntr,ncld,isiz,jsiz)
!
   VERTCORD='LORENZ'
!  VERTCORD='CHARNEY_PHILLIPS'
   allocate( aMSQUAR (mt)                                                      )
   allocate( JCOL2JS (0:jlha,2)                                                )
   allocate( SINLAT  (jgs:jge)                                                ,&
             TANLAT  (jgs:jge)                                                ,&
             COSLAT  (jgs:jge,3)                                               )
!
   allocate( AMATm   (mt,3,0:jlha,2)                                          ,&
             AMATm_  (mt,3,0:jlha,2)                                          ,&
             DMATm   (mt,3,0:jlha,2)                                          ,&
             DMATm_  (mt,3,0:jlha,2)                                          ,&
             DMATmi_ (mt,3,0:jlha,2,levsp)    )   ! should be global
!
   LAP_DIM=4
   lapend=LAP_DIM
#ifndef BARO_TEST
   lapend=LAP_DIM*3
#endif
   allocate( ZDMATmi_a_RE_v(mt,3,0:jlha,2,lapend)                             ,&
             ZDMATmi_a_RE_d(mt,3,0:jlha,2,lapend)                             ,&
             ZDMATmi_a_RE_t(mt,3,0:jlha,2,lapend)                             ,&
             ZDMATmi_a_RE_q(mt,3,0:jlha,2,lapend)                             ,&
             ZDMATmi_a_IM_v(mt,3,0:jlha,2,lapend)                             ,&
             ZDMATmi_a_IM_d(mt,3,0:jlha,2,lapend)                             ,&
             ZDMATmi_a_IM_t(mt,3,0:jlha,2,lapend)                             ,&
             ZDMATmi_a_IM_q(mt,3,0:jlha,2,lapend)                              )
!-------------------------------------------------------------------------------
! for 3d-model 
!-------------------------------------------------------------------------------
   allocate( sigma      (levs)                                                ,&
             sigmafull  (levs+1)                                              ,&
             tsigma     (levs)                                                ,&
             delsig     (levs)                                                ,&
             valpha     (levs)                                                ,&
             beta       (levs)                                                ,&
             gama1      (levs)                                                ,&
             gama2      (levs)                                                ,&
             TAVEXY     (levs)                                                ,&
             FAVEXY     (levs)                                                 )
   allocate( GV         (levs)                                                ,&
             PHS        (levs)                                                ,&
             TRS        (levs,levs)                                           ,&
             GRS        (levs,levs)                                           ,&
             GRSINV     (levs,levs)                                            )
   allocate( PHSDM      (mt,jlg)                                              ,&
             TRSDM      (mt,jlg,levsp)                                         )
   allocate( sigdot     (ib,jbw,levs+1)                                       ,&
             APSK       (ib,jbw,levs)                                         ,&
             DIVN       (ib,jbw)                                              ,&
             wR         (ib,jbw,levs)                                         ,&
             dwR        (ib,jbw,levs)                                         ,&
             ADT1       (ib,jbw,levs+1)                                        )
#ifdef HDSZ
   allocate( Z_VIS             (levs)                                         ,&
             T_VIS         (jbw,levs)                                         ,&
             T_REF      (ib,jbw,levs)                                          )
#endif
   allocate( EIGVAL     (levs)                                                ,&
             VERMAT     (levs,levs)                                           ,&
             EIGVEC     (levs,levs)                                           ,&
             AEIGVEC    (levs,levs)                                            )
!------------------------------------------------------------------------------
!
! spectral variables - tendency
!
   allocate(VJ(1:mt,1:jlg,   1:tlevsp+1) )
   DJ  =>   VJ(1:mt,1:jlg,tlvp:tlevsp+1)
   TJ  =>   VJ(1:mt,1:jlg,tltp:tlevsp+1)
#ifndef NISLQ
   QJ  =>   VJ(1:mt,1:jlg,tlqp:tlevsp+1)
#endif
   PSJ =>   VJ(1:mt,1:jlg,tlevsp+1)
#ifdef PERT
!
!  pertubation
!
   allocate(V0(1:mt,1:jlg,   1:tlevsp+1) )
   D0  =>   V0(1:mt,1:jlg,tlvp:tlevsp+1)
   T0  =>   V0(1:mt,1:jlg,tltp:tlevsp+1)
#ifndef NISLQ
   Q0  =>   V0(1:mt,1:jlg,tlqp:tlevsp+1)
#endif
   PSL0=>   V1(1:mt,1:jlg,tlevsp+1)
#endif
!
! n-1 time step
!
   allocate(V1(1:mt,1:jlg,   1:tlevsp+2) )
   D1  =>   V1(1:mt,1:jlg,tlvp:tlevsp+1)
   T1  =>   V1(1:mt,1:jlg,tltp:tlevsp+1)
#ifndef NISLQ
   Q1  =>   V1(1:mt,1:jlg,tlqp:tlevsp+1)
#endif
   PSL1=>   V1(1:mt,1:jlg,tlevsp+1)
   SF1 =>   V1(1:mt,1:jlg,tlevsp+2)
#ifdef MP
   allocate(SF1l (1:mt,1:jl)             )
#endif
!
! n time step
!
   allocate(V2(1:mt,1:jlg,   1:tlevsp+3) )
   D2  =>   V2(1:mt,1:jlg,tlvp:tlevsp+3)
   T2  =>   V2(1:mt,1:jlg,tltp:tlevsp+3)
#ifndef NISLQ
   Q2  =>   V2(1:mt,1:jlg,tlqp:tlevsp+3)
#endif
   PSL2=>   V2(1:mt,1:jlg,tlevsp+1)
   XPS2=>   V2(1:mt,1:jlg,tlevsp+2)
   LAP2=>   V2(1:mt,1:jlg,tlevsp+3)
!
! n+1 time step
!
   allocate(V3(1:mt,1:jlg,   1:tlevsp+3) )
   D3  =>   V3(1:mt,1:jlg,tlvp:tlevsp+3)
   T3  =>   V3(1:mt,1:jlg,tltp:tlevsp+3)
#ifndef NISLQ
   Q3  =>   V3(1:mt,1:jlg,tlqp:tlevsp+3)
#endif
   PSL3=>   V3(1:mt,1:jlg,tlevsp+1)
   XPS3=>   V3(1:mt,1:jlg,tlevsp+2)
   LAP3=>   V3(1:mt,1:jlg,tlevsp+3)
!
! lon,lat
!
   allocate(xlon(1:ib,1:jbw))
   allocate(xlat(1:ib,1:jbw))
!
! grid variables
!
   allocate(VOR(1:ib,1:jbw,  1:tlevs+3)   )
   DIO =>   VOR(1:ib,1:jbw,tlv:tlevs+3)
   TAI =>   VOR(1:ib,1:jbw,tlt:tlevs+3)
#ifndef NISLQ
   QAI =>   VOR(1:ib,1:jbw,tlq:tlevs+3)
#endif
   prs =>   VOR(1:ib,1:jbw,tlevs+1)
   xpsl=>   VOR(1:ib,1:jbw,tlevs+2)
   plap=>   VOR(1:ib,1:jbw,tlevs+3)
!
! uv wind
!
   allocate(VXC(1:ib,1:jbw,     1:2*levs+1) )
   VYC =>   VXC(1:ib,1:jbw,levs+1:2*levs+1)
   ypsl=>   VXC(1:ib,1:jbw,2*levs+1)
   allocate(spdmax(levs))
#ifdef CLM_CWF
   allocate(cgs(ib,jbw,levs))
#endif
#ifdef BARO_TEST
   allocate(SFAI(ib,jbw))
#endif
   allocate(APSN(ib,jbw))
#ifdef NISLQ
   allocate(QJ(1:mt,1:jlg,1:levhp))
#ifdef PERT
   allocate(Q0(1:mt,1:jlg,1:levhp))
#endif
   allocate(Q1(1:mt,1:jlg,1:levhp))
   allocate(Q2(1:mt,1:jlg,1:levhp))
   allocate(Q3(1:mt,1:jlg,1:levhp))
   allocate(QAI(1:ib,1:jbw,1:levh))
#endif /* NISLQ end */
!==========================================================!
!     coefficient of NEWTONIAN COOLING and DISSIPATION     !
!==========================================================!
   PI2     = pi*2
   PI180   = PI/180.d00
   DISSI   = 0./ ( (PI2* 3.0 )*(MTA*(MTA+1.)) )
   GENER   = 0./ (  PI2*15.0 )
   cooling = 0./ (  PI2*15.0 )
!
   press   = 0.0
   damp1   = 0.050d0
   damp2   = 0.050d0
!
   return
   end subroutine dfsvarini
!-------------------------------------------------------------------------------
   end module dfsvar
