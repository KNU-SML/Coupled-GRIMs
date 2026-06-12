#include <define.h>
   subroutine dyn_sph_driver
!-------------------------------------------------------------------------------
!
! subroutine: dyn_sph_driver         
!
! abstract: 
! - computes dynamic non-linear tendency terms of temp. div. ln(ps)
! - omputes predicted values of vorticity and moisture
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! references : 
!   hong et al. (2013, apjas): global/regional integrated model system (grims)
!   park et al. (2013, mwr): dfs dynamical core
!   byun and hong (2007, j. climate): single-column model (SMP)
!   kanamitsu et al. (2002, bams): ncep dynamical seasonal forecast system 2000
!
!-------------------------------------------------------------------------------
   use paramodel, only : JCAP1S,LONF2S,LATG2S,LEVSS,LEVHS,LCAPS,LCAP22S       ,&
                         LNT2S,levs_,levh_,lnt22_,ncpus_                      ,&
#ifndef DYNAMIC_ALLOC
                         ncpus=>ncpus_                                        ,&
#endif
#ifdef MP
                         lonf22p_,levhp_,lnt2p_                               ,&
                         lnt22p_,lln22p_,levsp_,lcap22p_,latg2p_              ,&
#else
                         lnt2_                                                ,&
#endif
                         jcap1_,lonf22_,lonf2_,latg2_,ntotal_,ngases_
   use comkindex
#ifdef NISLQ
   use nislq    , only : slq_q2,slq_psfc2,slq_u2,slq_v2,slq_w2
   use paramodel, only : LONF22S
#endif
   use constant, only  : cp_,g_,rd_,rerth_,hfus_,hvap_,qmin_
   use comfibm
   use comgrad
   use comgda
   use comfcst, only   : spdlat
#ifndef SMP
#ifdef PSPLIT
   use compsplit
#endif
#ifdef MP
   use commpi
   use compspec
#endif
#ifdef REDUCE_GRID
   use comreduce
#else
   use paramodel, only : lonf_
#endif
#endif  /* SMP */
#ifdef DG3
   use diag_3d_module, only : diag_3d_get, diag_3d_arrange
   use comgda
#endif
!
#ifdef SMP
   use comfcst, only : dtbdy,curtime,vvel,hour1,hour2,hours,houre,ftim1,ftim2 ,&
                       uo1,vo1,to1,qo1,wo1,ps1,uo2,vo2,to2,qo2,wo2,ps2        ,&
                       at1,aq1,dv1,at2,aq2,dv2
#ifdef CLM_CWF
   use comfcst, only : wdiv,hadq
#endif
#ifdef REV_FRC
   use comfcst, only : ac1, ac2
#endif
#ifdef RLX_FRC
   use comfcst, only : tscl1, tscl2
#endif
#endif /* SMP */
!
#ifdef CLM_CWF
   use comfcst, only : cgs
#endif
#ifdef RAS
   use phys_ras_module, only  : ras_setup
#endif
#ifdef RASV2
   use phys_ras2_module, only : ras2_setup
   use comfcst, only          : sig, sgb, prj, rasal, krmin, krmax, kfmax,     &
                                kctop, ncrnd, dsfc
#endif
#ifdef STOCH
   use stochastic
#endif
!
#include "abort.h"
!
! syn(1, 0*levs_+1, lan)  ze
! syn(1, 1*levs_+1, lan)  di
! syn(1, 2*levs_+1, lan)  te
! syn(1, 3*levs_+1, lan)  rq
! syn(1, 4*levs_+1, lan)  uln
! syn(1, 5*levs_+1, lan)  vln
! syn(1, 6*levs_+1, lan)  dpdphi
! syn(1, 6*levs_+2, lan)  dpdlam
! syn(1, 6*levs_+3, lan)  q
!
! dyn(1, 0*levs_+1, lan)  d(t)/d(phi)
! dyn(1, 1*levs_+1, lan)  d(rq)/d(phi)
! dyn(1, 2*levs_+1, lan)  d(t)/d(lam)
! dyn(1, 3*levs_+1, lan)  d(rq)/d(lam)
! dyn(1, 4*levs_+1, lan)  d(u)/d(lam)
! dyn(1, 5*levs_+1, lan)  d(v)/d(lam)
! dyn(1, 6*levs_+1, lan)  d(u)/d(phi)
! dyn(1, 7*levs_+1, lan)  d(v)/d(phi)
!
! anl(1, 0*levs_+1, lan)  z     dqdt
! anl(1, 0*levs_+2, lan)  uu    dudt
! anl(1, 1*levs_+2, lan)  vv    dvdt
! anl(1, 2*levs_+2, lan)  y     dtdt
! anl(1, 3*levs_+2, lan)  rt    drdt
!
   integer    ::  lotd   ,lotds                                               ,&
                  kdtphi ,kdtphis                                             ,&
                  kdrphi ,kdrphis                                             ,&
                  kdtlam ,kdtlams                                             ,&
                  kdrlam ,kdrlams                                             ,&
                  kdulam ,kdulams                                             ,&
                  kdvlam ,kdvlams                                             ,&
                  kduphi ,kduphis                                             ,&
                  kdvphi ,kdvphis 
!
   integer    ::  lota   ,lotas                                               ,&
                  lotat  ,lotats                                              ,&
                  kap    ,kaps                                                ,&
                  kau    ,kaus                                                ,&
                  kav    ,kavs                                                ,&
                  kat    ,kats                                                ,&
                  kar    ,kars
   integer    ::  klev
!
! local array
!
#ifdef SMP
#ifdef DBG
#ifndef CLM_CWF
   real                 ::  hadq(lnt22_,levs_)
#endif
   real                 ::  hadt(lnt22_,levs_)
   real                 ::  vadt(lnt22_,levs_), vadq(lnt22_,levh_)
   real                 ::  adbc(lnt22_,levs_)
   real                 ::  rlxt(lnt22_,levs_), rlxq(lnt22_,levh_)
#endif /* DBG end */
#define NCPUSS ncpus_
#else	/* ~SMP */
#ifdef MP
   real, allocatable, dimension(:,:,:)  ::  syf,grs,dyf,dgr,anf,gra
#endif
   real, allocatable, dimension(:)      ::  scs2
   real, allocatable, dimension(:,:,:)  ::  syn,syntop,dyn,anl,anltop
   real, allocatable, dimension(:,:,:,:)::  flp,flm
#ifdef HYBRID
   integer              ::   nvcn
   real                 ::   xvcn
#endif
!
#ifdef PSPLIT
#ifdef DG
#ifdef MP
#define NCPUSS latg2p_
#else
#define NCPUSS ncpus_
#endif
   real                 ::   tgmxl(NCPUSS),tgmnl(NCPUSS)
   integer              ::   igmxl(NCPUSS),kgmxl(NCPUSS)
   integer              ::   igmnl(NCPUSS),kgmnl(NCPUSS)
#endif /* DG end */
#endif /* PSPLIT end */
#endif /* SMP end */
!
#ifdef DG
#ifdef MP
#define NCPUSS latg2p_
#else
#define NCPUSS ncpus_
#endif
   real                 ::  cldt(LONF2S,levs_,NCPUSS)
   real                 ::  clcv(LONF2S,levs_,NCPUSS)
#endif
#ifdef DG3
   real                 ::  gda(nwgda,kdgda,NCPUSS)
#endif
#ifdef NISLQ
#ifdef HYBRID
   real, dimension(LONF22S,levs_+1)           ::  pdot
#else
   real, dimension(LONF22S,levs_+1)           ::  dot
#endif
   real, dimension(LONF22S,levh_  ,LATG2S)    ::  tmp1,tmp2,tmp3,slq
#endif /* NISLQ end */
!
#ifdef PSPLIT 
#ifdef DG3
   real                 ::  gda(nwgda,kdgda,NCPUSS)
#endif
#ifdef RAS
   logical,parameter    ::  ras=.true.
   integer,parameter    ::  nsphys=1
   real,parameter       ::  cp=cp_, alhl=hvap_, grav=g_, rgas=rd_
   real                 ::  sig(levs_+1), prj(levs_+1), prh(levs_)
   real                 ::  fpk(levs_),   hpk(levs_),   sgb(levs_)
   real                 ::  ods(levs_), rasal(levs_),  prns(levs_/2)
   real                 ::  rannum(200)
!
   call ras_setup(levs_, si, sl, del, cp, rgas, deltim                         &
                ,nsphys, thour                                                 &
                ,sig, sgb, prh, prj, hpk, fpk, ods, prns                       &
                ,rasal, lm, krmin, krmax, nstrp                                &
                ,ncrnd, rannum, afac, ufac)
#endif
#ifdef RASV2
   logical,parameter    ::  ras=.true.
   integer,parameter    ::  nsphys=1
   real,parameter       ::  cp=cp_, alhl=hvap_, grav=g_, rgas=rd_, hfus=hfus_
!
!  pdd is the lowest pressure above which downdraft is allowed
!
   real, parameter      ::  pdd=600.0
   real                 ::  rannum(200,5)
!
!-------------------------------------------------------------------------------
!
   frain = .5
   if(inistp.eq.1) frain = 1.
   dtp= deltim / frain
!
!  check fhour and thour
!
   call ras2_setup(levs_,  si,  sl, cp_, rd_, dtp, nsphys, fhour               &
                ,sig,   sgb, prj                                               &
                ,rasal, krmin, krmax, kfmax, ncrnd, rannum                     &
                ,mct,   kctop, deltim, dsfc, LONF2S*LATG2S)
#endif
!
   if(ngases_.ge.1) call rad_ozone_setup(idate,fhour,no3p,no3l)
!
#endif /* PSPLIT end */
#ifdef STOCH
#ifdef NLT_U
   call stc_forcing(stc_u)
#endif
#ifdef NLT_V
   call stc_forcing(stc_v)
#endif
#ifdef NLT_T
   call stc_forcing(stc_t)
#endif
#ifdef NLT_Q
   call stc_forcing(stc_q)
#endif
#endif /* STOCH end */
!
! K-loop index
!
#ifndef NISLQ
   lotd    =6*levs_+2*levh_
#else
   lotd    =6*levs_
#endif
   kdtphi  =1
#ifndef NISLQ
   kdrphi  =kdtphi+levs_
   kdtlam  =kdrphi+levh_
   kdrlam  =kdtlam+levs_
   kdulam  =kdrlam+levh_
#else
   kdtlam  =kdtphi+levs_
   kdulam  =kdtlam+levs_
#endif
   kdvlam  =kdulam+levs_
   kduphi  =kdvlam+levs_
   kdvphi  =kduphi+levs_
!
#ifndef NISLQ
   lota    =3*levs_+levh_+1
#else
   lota    =3*levs_+1
#endif
   lotat   =2*levs_
   kap     =1
   kau     =kap+1
   kav     =kau+levs_
   kat     =kav+levs_
#ifndef NISLQ
   kar     =kat+levs_
#endif
!
#ifndef NISLQ
   lotds   =6*LEVSS+2*LEVHS
#else
   lotds   =6*LEVSS
#endif
   kdtphis =1
#ifndef NISLQ
   kdrphis =kdtphis+LEVSS
   kdtlams =kdrphis+LEVHS
   kdrlams =kdtlams+LEVSS
   kdulams =kdrlams+LEVHS
#else
   kdtlams =kdtphis+LEVSS
   kdulams =kdtlams+LEVSS
#endif
   kdvlams =kdulams+LEVSS
   kduphis =kdvlams+LEVSS
   kdvphis =kduphis+LEVSS
!
#ifndef NISLQ
   lotas   =3*LEVSS+LEVHS+1
#else
   lotas   =3*LEVSS+1
#endif
   lotats  =2*LEVSS
   kaps    =1
   kaus    =kaps+1
   kavs    =kaus+LEVSS
   kats    =kavs+LEVSS
#ifndef NISLQ
   kars    =kats+LEVSS
#endif
!
#ifndef SMP
#ifdef MP
   allocate( syf(lonf22_,lotss,latg2p_),grs(lonf22p_,lots,latg2p_),            &
             dyf(lonf22_,lotds,latg2p_),dgr(lonf22p_,lotd,latg2p_),            &
             anf(lonf22_,lotas,latg2p_),gra(lonf22p_,lota,latg2p_) )
#endif
!
#ifdef MP
#define NCPUSS latg2_
#else
#define NCPUSS ncpus_
#endif
   allocate(  syn(LCAP22S,lotss,NCPUSS),syntop(2,JCAP1S,lotsts),               &
              dyn(LCAP22S,lotds,NCPUSS),scs2(JCAP1S),                          &
              anl(LCAP22S,lotas,NCPUSS),anltop(2,JCAP1S,lotats),               &
              flp(2,JCAP1S,lotas,NCPUSS),flm(2,JCAP1S,lotas,NCPUSS) )
#endif
!
! initialize local variables
!
   y=0.
#ifdef CLM_CWF
   hadq=0.
#endif
#ifdef SMP
#ifdef DBG
#ifndef CLM_CWF
   hadq=0.
#endif
   hadt=0.
   vadt=0.
   vadq=0.
   adbc=0.
   rlxt=0.
   rlxq=0.
#endif /* DBG end */
#else /* ~SMP */
#ifdef MP
   syf=0. ; grs=0.    ; dyf=0. ; dgr=0.  ; anf=0. ; gra=0. 
#endif
   syn=0. ; syntop=0. ; dyn=0. ; scs2=0. ; anl=0. ; anltop=0. ; flp=0. ; flm=0.
#ifdef PSPLIT
#ifdef DG
   tgmxl=0. ; tgmnl=0.
   igmxl=0  ; kgmxl=0  ; igmnl=0  ; kgmnl=0 
#endif /* DG end */
#endif /* PSPLIT end */
#endif /* SMP end */
!
#ifdef MP
   llstr=lwvstr(mype)
   llens=lwvlen(mype)
   jstr=latstr(mype)
   jend=latstr(mype)+latlen(mype)-1
   lons2=lonlen(mype)*2
   lats2=latlen(mype)
#else
   llstr=0
   llens=jcap1_
   lons2=lonf2_
   lats2=latg2_
#endif
#ifdef SMP
!=========================================================================
! SMP : read dynamic tendency terms (advection)
!=========================================================================
!
#ifdef DG3
#define LATX lan
#else
#define LATX lat
#endif
   lonsd2 = lons2
   lat = 1
   latdon = 0
   lan = lat - latdon
!
   write(6,*) 'dyn_sph_driver check bdy'
   do k = 1,levs_
     if (k.lt.3) then
       write(6,'(a4,i2)') ' k =',k
       write(6,220) ' uo1 vo1 to1 qo1 ',uo1(1,k),vo1(1,k),to1(1,k),qo1(1,k)
       write(6,220) ' at1 aq1 wo1 dv1 ',at1(1,k),aq1(1,k),wo1(1,k),dv1(1,k)
       write(6,230) ' ps1             ',ps1(1)
       write(6,230) ' tscl1           ',tscl1(1,k)
     endif
   enddo
220  format(a17,4e15.7)
230  format(a17,2e15.7)
!
! ... set wind & pressure ...
!
   do l = 1,lnt2_
     do k = 1,levs_
       x(l,k) = ftim2*uo1(l,k) + ftim1*uo2(l,k)
       w(l,k) = ftim2*vo1(l,k) + ftim1*vo2(l,k)
       vvel(l,k) = ftim2*wo1(l,k) + ftim1*wo2(l,k)
#ifdef CLM_CWF
       wdiv(l,k) = ftim2*dv1(l,k) + ftim1*dv2(l,k)
#endif
     enddo
     z(l) = ftim2*ps1(l) + ftim1*ps2(l)
   enddo
!
#ifndef REV_FRC
#ifdef DBG
   call sph_nonlinear_tend(lnt2_,vvel,te,rq,q,del,sl,rdel2,y,rt,vadt,adbc)
#else
   call sph_nonlinear_tend(lnt2_,vvel,te,rq,q,del,sl,rdel2,y,rt)
#endif
#endif
!
   do l = 1,lnt2_
     do k = 1,levs_
#ifdef RLX_FRC
#ifndef TOGA
       advtscl = ftim2*tscl1(l,k) + ftim1*tscl2(l,k) ! relaxation time scale (s)
#else
       advtscl = 3600.*24.  ! from Iacobellis & Somerville (2000)
#endif
#endif
       advq    = -1.0*(ftim2*aq1(l,k) + ftim1*aq2(l,k))
#ifdef RLX_FRC
#ifdef CLM_CWF
       hadq(l,k) = advq
#endif
       qobs = ftim2*qo1(l,k) + ftim1*qo2(l,k)
       advq = rt(l,k) + advq + (qobs-rm(l,k))/advtscl
#else
       advq = rt(l,k) + advq
#endif
#ifdef DBG
       vadq(l,k) = rt(l,k)
#ifdef RLX_FRC
       rlxq(l,k) = (qobs-rm(l,k))/advtscl
#endif
#endif
       rt(l,k) = rm(l,k) + 2.0*deltim*advq
#ifdef SMP_NUDGING_RAD
       rt(l,k) = qobs
#else
       rt(l,k) = max(rt(l,k),qmin_)
#endif
     enddo
!
     if (ntotal_.gt.1) then
       lev_qrs = levs_ + 1
       do k = levs_+1,levh_
!ld3v0   rt(l,k) = 0.0      ! no time change -> error
!ld3v1   rt(l,k) = rm(l,k)  ! zero advection (horizontal & vertical)
         advq = rt(l,k)
         if (k.ge.lev_qrs) advq = 0.0
         rt(l,k) = rm(l,k) + 2.0*deltim*advq ! only vertical advection
         rt(l,k) = max(rt(l,k),qmin_)
       enddo
     endif
!
     do k = 1,levs_
#ifdef RLX_FRC
#ifndef TOGA
       advtscl = ftim2*tscl1(l,k) + ftim1*tscl2(l,k) ! relaxation time scale (s)
#else
       advtscl = 3600.*24.  ! from Iacobellis & Somerville (2000)
#endif
#endif
       adtv = -1.0*(ftim2*at1(l,k) + ftim1*at2(l,k))  ! advection of Tv
#ifdef REV_FRC
       adia = ftim2*ac1(l,k) + ftim1*ac2(l,k)
       advt = y(l,k) + adtv + adia
#else
#ifdef RLX_FRC
       tobs = ftim2*to1(l,k) + ftim1*to2(l,k)
       advt = y(l,k) + adtv + (tobs-tem(l,k))/advtscl
#else
       advt = y(l,k) + adtv
#endif
#endif
#ifdef DBG
       hadt(l,k) = adtv
#ifdef RLX_FRC
       rlxt(l,k) = (tobs-tem(l,k))/advtscl
#endif
#ifdef REV_FRC
       adbc(l,k) = adia
#endif
#endif
#ifdef SMP_NUDGING_RAD
       y(l,k) = tobs
#else
       y(l,k) = tem(l,k) + 2.0*deltim*advt
#endif
     enddo
   enddo
!
#ifdef DBG
   write(6,'(A,i5,6f8.3)') 'BDY',kdt,thour,hours,houre,ftim1,ftim2,dtbdy
#ifdef REV_FRC
   write(6,'(7e13.5)')uo1(1,1),vo1(1,1),at1(1,1),ac1(1,1),aq1(1,1),            &
             wo1(1,1),exp(ps1(1))
   write(6,'(7e13.5)')uo2(1,1),vo2(1,1),at2(1,1),ac2(1,1),aq2(1,1),            &
             wo2(1,1),exp(ps2(1))
   write(6,'(6e13.5)')x(1,1),w(1,1),y(1,1),rt(1,1),vvel(1,1),exp(z(1))
#else
   write(6,'(6e13.5)')uo1(1,1),vo1(1,1),at1(1,1),aq1(1,1),wo1(1,1),exp(ps1(1))
   write(6,'(6e13.5)')uo2(1,1),vo2(1,1),at2(1,1),aq2(1,1),wo2(1,1),exp(ps2(1))
   write(6,'(6e13.5)')x(1,1),w(1,1),y(1,1),rt(1,1),vvel(1,1),exp(z(1))
#endif
!
   do k = 1,levs_
     write(95,'(f8.3,2i5,14e13.5)')                                            &
         curtime,kdt,k,x(1,k),w(1,k),y(1,k),rt(1,k),vvel(1,k),                 &
         z(1),hadq(1,k),vadq(1,k),rlxq(1,k),                                   &
         hadt(1,k),vadt(1,k),adbc(1,k),rlxt(1,k)
   enddo
#endif
!
   first=.false.
!
#else
!=========================================================================
! GMP : calculate dynamic tendency terms
!=========================================================================
!
#ifdef MP
!   
#ifndef NISLQ
   klev=3+ntotal_  ! ze,di,te,rq
#else
   klev=3          ! ze,di,te
#endif
   call mpnn2nk(ze ,lnt22p_,levs_, zea,lln22p_,levsp_,klev)
   call mpnn2n (q  ,lnt22p_,       qa ,lln22p_,       1)
#define QS qa
#define DPDPHIS dpdphia
#define DPDLAMS dpdlama
#define QLAPS qlapa
#define DIS dia
#define ZES zea
#define ULNS ulna
#define VLNS vlna
#else
#define QS q
#define DPDPHIS dpdphi
#define DPDLAMS dpdlam
#define QLAPS qlap
#define DIS di
#define ZES ze
#define ULNS uln
#define VLNS vln
#endif
   call sph_del_log_ps(QS,DPDPHIS,syntop(1,1,2*LEVSS+1),DPDLAMS,               &
               llstr,llens,lwvdef)
   call sph_solve_laplacian(QS,QLAPS,llstr,llens,lwvdef)
   call sph_divor2wind(DIS,ZES,ULNS,VLNS,syntop(1,1,1),                        &
               syntop(1,1,LEVSS+1),llstr,llens,lwvdef)
#undef QS
#undef DPDPHIS
#undef DPDLAMS
#undef QLAPS
#undef DIS
#undef ZES
#undef ULNS
#undef VLNS
!
   do n = 1,2
     do j = 1,JCAP1S
       do l = 1,lotats
         anltop(n,j,l)=0.0
       enddo
     enddo
   enddo
#ifdef PSPLIT
#ifdef DG
   tgmx=-1.e20
   tgmn= 1.e20
#endif
#endif
!
#ifdef MP
#define LNT2X lln22p_
#define UUS uua
#define VVS vva
#define YS ya
#define RTS rta
#define ZS za
#else
#define LNT2X lnt2_
#define UUS uu
#define VVS vv
#define YS y
#define RTS rt
#define ZS z
#endif
!
#ifdef ORIGIN_THREAD
!$doacross share(uu,uua,vv,vva,y,ya,rt,rta),local(j,k)
#endif
#ifdef OPENMP
!$omp parallel do private(j,k)
#endif
   do k = 1,LEVSS
     do j = 1,LNT2X
       UUS(j,k)=0.0
       VVS(j,k)=0.0
       YS(j,k)=0.0
     enddo
   enddo
#ifndef NISLQ
   do k = 1,LEVHS
     do j = 1,LNT2X
       RTS(j,k)=0.0
     enddo
   enddo
#endif
   do j = 1,LNT2X
     ZS(j,1)=0.e0
   enddo
#undef LNT2X
#undef UUS
#undef VVS
#undef YS
#undef RTS
#undef ZS
!
#ifndef MP
! compute latitude band limits in single mode
!
   last=mod(latg2_,ncpus)
   nggs=(latg2_-last)/ncpus
   if(last.ne.0)nggs=nggs+1
   inclat=ncpus
   lat1=1-ncpus
   lat2=0
   latdon=0
   !---------------------------------------------
   do 10000 ngg = 1,nggs
   !---------------------------------------------
   if((ngg.eq.nggs).and.(last.ne.0)) inclat=last
   lat1=lat1+ncpus
   lat2=lat2+inclat
#endif
!
#ifdef PSPLIT
#ifdef DG3
#ifdef MP
   lat1=jstr
   lat2=jend
   latdon=jstr-1
#define LATX lan
#else
#define LATX lat
#endif /* MP end */
   do lat = lat1,lat2
     lan=lat-latdon
     call diag_3d_get(LATX,nwgda*kdgda,gda(1,1,lan))
   enddo
#undef LATX
#endif /* DG3 end */
!
#ifdef DG
#ifdef MP
#define NCPUSS latg2p_
#else
#define NCPUSS ncpus
#endif /* MP end */
   do k = 1,NCPUSS
     tgmxl(k)=tgmx
     tgmnl(k)=tgmn
   enddo
#undef NCPUSS
#endif /* DG end */
#endif /* PSPLIT end */
!
#ifdef MP
   lat1=1
   lat2=latg2_
   latdon=0
#define ZES zea
#else
#define ZES ze
#endif
!
! first lat loop
!
#ifdef ORIGIN_THREAD
!$doacross share(syntop,syn,qtt,qvv,lat1,lat2,latdon,
!$&        colrad,ze,zea,llstr,llens,lwvdef,lcapdp,mype,lcapd)
!$&        local(lat,lan,llensd)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(syntop,syn,qtt,qvv,lat1,lat2,latdon)
!mic$1 shared(colrad,ze,zea,llstr,llens,lwvdef,lcapdp,mype,lcapd)
!mic$1 private(lat,lan,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,llensd)
#endif
   do lat = lat1,lat2
     lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
     llensd=lcapdp(lat,mype)
#else
     llensd=lcapd(lat)
#endif
#else
     llensd=llens
#endif
     call sph_sum_coeff(ZES,syn(1,kszs,lan),qtt(1,lat),llstr,llensd,           &
                          lwvdef,lotss)
     call sph_sum_coeff_top(syn(1,ksus,lan),syntop(1,1,1),                     &
                              qvv(1,lat),llstr,llensd,lwvdef,2*LEVSS+1)
   enddo
!
! compute merid. derivs. of temp. and moisture using qdd.
!
#ifdef ORIGIN_THREAD
!$doacross share(dyn,qdd,lat1,lat2,latdon,tea,te,
!$&        llstr,llens,lwvdef,lcapd,lcapdp,mype),
!$&        local(lat,lan,i,k,llensd)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(dyn,qdd)
!mic$1 shared(lat1,lat2,latdon,llens)
!mic$1 shared(tea,te,llstr,llens,lwvdef,lcapd,lcapdp,mype)
!mic$1 private(lat,lan,i,k,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,i,k,llensd)
#endif
   do lat = lat1,lat2
     lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
     llensd=lcapdp(lat,mype)
#else
     llensd=lcapd(lat)
#endif
#else
     llensd=llens
#endif
#ifndef NISLQ
     klev=LEVSS+LEVHS  ! te,rq
#else
     klev=LEVSS        ! te
#endif
     call sph_sum_coeff(ZES(1,ksts),dyn(1,kdtphis,lan),qdd(1,lat),             &
                 llstr,llensd,lwvdef,klev)
!
!   d(t)/d(phi)  d(rq)/d(phi) in s. hemi.
!
#ifndef NISLQ
     klev=LEVSS+LEVHS
#else
     klev=LEVSS
#endif
     do k = 1,klev
       do i = 1,llensd*2
           dyn(i+LCAPS,kdtphis-1+k,lan)=-dyn(i+LCAPS,kdtphis-1+k,lan)
       enddo
     enddo
   enddo
!
#ifdef ORIGIN_THREAD
!$doacross share(dyn,rcs2,syn,lat1,lat2,latdon,
!$&              llstr,llens,lwvdef,lcapd,lcapdp,mype),
!$&              local(lat,lan,j,k,l,kdd,kss,scs2,llensd)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(dyn,rcs2,llstr,llens,lwvdef)
!mic$1 shared(syn,lat1,lat2,latdon,lcapd,lcapdp,mype)
!mic$1 private(lat,lan,j,k,l,kdd,kss,scs2,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,j,k,l,kdd,kss,scs2,llensd)
#endif
!
   do lat = lat1,lat2
     lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
     llensd=lcapdp(lat,mype)
#else
     llensd=lcapd(lat)
#endif
#else
     llensd=llens
#endif
!
!    calculate t rq u v zonal derivs. by multiplication with i*l
!
     do l = 1,llensd
       j=lwvdef(llstr+l)
       scs2(l)=float(j)/rerth_ *rcs2(lat)
     enddo
!
!    calculate t rq u v zonal derivs. by multiplication with i*l
!    note scs2=rcs2*l/rerth_
!     
!    d(t)/d(lam)  d(rq)/d(lam)  d(u)/d(lam)  d(v)/d(lam) .......
!     
#ifndef NISLQ
     klev=3*LEVSS+LEVHS
#else
     klev=3*LEVSS
#endif
     do k = 1,klev
       kdd=kdtlams-1+k
       kss=ksts-1+k
       do j=1,llensd
         dyn(      2*j-1,kdd,lan)=-syn(      2*j  ,kss,lan)*scs2(j)
         dyn(      2*j  ,kdd,lan)= syn(      2*j-1,kss,lan)*scs2(j)
         dyn(LCAPS+2*j-1,kdd,lan)=-syn(LCAPS+2*j  ,kss,lan)*scs2(j)
         dyn(LCAPS+2*j  ,kdd,lan)= syn(LCAPS+2*j-1,kss,lan)*scs2(j)
       enddo
     enddo
   enddo
!
#undef ZES
!
#ifdef MP
   call mpnl2ny(syn,lcap22p_,latg2_,syf,lonf22_,latg2p_,lotss,1,lotss)
#ifndef NISLQ
   klev=4*levsp_+2*levhp_
#else
   klev=4*levsp_
#endif
   call mpnl2ny(dyn,lcap22p_,latg2_,dyf,lonf22_,latg2p_,lotds,1,klev)
!
   lat1=jstr
   lat2=jend
   latdon=jstr-1
#define SYNS syf
#define DYNS dyf
#else
#define SYNS syn
#define DYNS dyn
#endif
#ifdef ORIGIN_THREAD
!$doacross share(syf,syn,dyf,dyn,lat1,lat2,latdon,latdef,lcapd,lonfd),
!$&        local(lat,lan,k,j,lcapf,lonff)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$& shared(syf,syn,dyf,dyn)
!mic$& shared(lat1,lat2,latdon,latdef,lcapd,lonfd)
!mic$& private(lat,lan,k,j,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,k,j,lcapf,lonff)
#endif
   do lat = lat1,lat2
#ifdef REDUCE_GRID
     lcapf=lcapd(latdef(lat))
     lonff=lonfd(latdef(lat))
#else
     lcapf=jcap1_
     lonff=lonf_
#endif
     lan=lat-latdon
     call sph_wave2grid(SYNS(1,1,lan),SYNS(1,1,lan),lotss*2,                   &
                   lcapf,lonff,latdef(lat),1)
!
#ifndef NISLQ
     klev=4*LEVSS+2*LEVHS
#else
     klev=4*LEVSS
#endif
     call sph_wave2grid(DYNS(1,kdtphis,lan),DYNS(1,kdtphis,lan),2*klev,        &
                   lcapf,lonff,latdef(lat),1)
     do k = 1,LEVSS
       do j = 1,lonff*2
         DYNS(j,kduphis-1+k,lan)=                                              &
               DYNS(j,kdvlams-1+k,lan)-SYNS(j,kszs-1+k,lan)
         DYNS(j,kdvphis-1+k,lan)=                                              &
               -DYNS(j,kdulams-1+k,lan)+SYNS(j,ksds-1+k,lan)
       enddo
     enddo
   enddo
!
#undef SYNS
#undef DYNS
!
#ifdef MP
!
#ifndef NISLQ
   klev=5+ntotal_
#else
   klev=5
#endif
   call mpnk2nx(syf,lonf22_,lotss,                                             &
                grs,lonf22p_,lots,latg2p_,levsp_,levs_,kszs,ksz,klev)
   call mpx2nx (syf,lonf22_,lotss,                                             &
#ifdef PSPLIT
                grs,lonf22p_,lots,latg2p_,kspphis,kspphi,4)
#else
                grs,lonf22p_,lots,latg2p_,kspphis,kspphi,3)
#endif
#ifndef NISLQ
   klev=6+2*ntotal_
#else
   klev=6
#endif
   call mpnk2nx(dyf,lonf22_,lotds,                                             &
                dgr,lonf22p_,lotd,latg2p_,levsp_,levs_,1,1,klev)
!
   lat1=jstr
   lat2=jend
   latdon=jstr-1
#define SYNS grs
#define DYNS dgr
#define ANLS gra
#define LATX lan
#else
#define SYNS syn
#define DYNS dyn
#define ANLS anl
#define LATX lat
#endif
#ifdef ORIGIN_THREAD
!$doacross share(dgr,dyn,rcs2,grs,syn,gra,anl,lat1,lat2,latdon,spdlat,
!$&              rbs2,del,rdel2,ci,p1,p2,h1,h2,tov,lons2,lonfdp,mype),
!$&              local(lat,lan,j,k,lonsd2)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(dyn,dgr,rcs2)
!mic$1 shared(syn,grs,anl,gra,lat1,lat2,latdon,spdlat)
!mic$1 shared(rbs2,del,rdel2,ci,p1,p2,h1,h2,tov,lons2,lonfdp,mype)
!mic$1 private(lat,lan,j,k,lonsd2)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,j,k,lonsd2)
#endif
   do lat = lat1,lat2
     latrue=latdef(lat)
     lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
     lonsd2=lonfdp(lan,mype)*2
#else
     lonsd2=lonfd(latdef(lat))*2
#endif
#else
     lonsd2=lons2
#endif
!
#ifdef MP
     if( lonsd2.gt.0 ) then
#endif
       do k = 1,levs_
         do j = 1,lonsd2
#ifdef PSPLIT
           ANLS(j,kau-1+k,lan)=dudtm(j,k,LATX)
           ANLS(j,kav-1+k,lan)=dvdtm(j,k,LATX)
           ANLS(j,kat-1+k,lan)=dtdtm(j,k,LATX)
#else
           ANLS(j,kau-1+k,lan)=0.
           ANLS(j,kav-1+k,lan)=0.
           ANLS(j,kat-1+k,lan)=0.
#endif
         enddo
       enddo
!
       do k = 1,levh_
         do j = 1,lonsd2
#ifdef NISLQ
           slq(j,k,lan)=slq_q2(j,k,LATX)
#else
#ifdef PSPLIT
           ANLS(j,kar-1+k,lan)=drdtm(j,k,LATX)
#else
           ANLS(j,kar-1+k,lan)=0.0
#endif
#endif
         enddo
       enddo
       do j = 1,lonsd2
         ANLS(j,kap,lan)=0.0
       enddo
#ifdef STOCH 
!
! update grid-point values for u,v,t,q
!
   do k=1,levs_
   do i=1,lonsd2
     grid_u1(i,k,lan)=grid_u2(i,k,lan)
     grid_u2(i,k,lan)=SYNS(i,ksu+k-1,lan)
!
     grid_v1(i,k,lan)=grid_v2(i,k,lan)
     grid_v2(i,k,lan)=SYNS(i,ksv+k-1,lan)
!
     grid_t1(i,k,lan)=grid_t2(i,k,lan)
     grid_t2(i,k,lan)=SYNS(i,kst+k-1,lan)
!
     grid_q1(i,k,lan)=grid_q2(i,k,lan)
#ifdef NISLQ
     grid_q2(i,k,lan)=slq(i,k,lan)
#else
     grid_q2(i,k,lan)=SYNS(i,ksr+k-1,lan)
#endif
   enddo
   enddo
#endif /* STOCH end */
!
#ifdef HYBRID
       call sph_nonlinear_tend_hybrid(lonsd2,lan,                              &
               ak5,bk5,                                                        &
               SYNS(1,ksd,lan),SYNS(1,kst,lan),                                &
               SYNS(1,ksz,lan),SYNS(1,ksu,lan),SYNS(1,ksv,lan),                &
#ifndef NISLQ
               SYNS(1,ksr,lan),                                                &
#else
               slq(1,1,lan),                                                   &
#endif
               SYNS(1,kspphi,lan),SYNS(1,ksplam,lan),                          &
               SYNS(1,ksp,lan),                                                &
               rbs2(LATX),                                                     &
               spdlat(1,LATX),deltim,nvcn,xvcn,                                &
               DYNS(1,kdtphi,lan),DYNS(1,kdtlam,lan),                          &
#ifndef NISLQ
               DYNS(1,kdrphi,lan),DYNS(1,kdrlam,lan),                          &
#else
               tmp1(1,1,lan)     ,tmp2(1,1,lan)     ,                          &
#endif
               DYNS(1,kdulam,lan),DYNS(1,kdvlam,lan),                          &
               DYNS(1,kduphi,lan),DYNS(1,kdvphi,lan),                          &
               ANLS(1,kap,lan),ANLS(1,kat,lan),                                &
#ifndef NISLQ
               ANLS(1,kar,lan),ANLS(1,kau,lan),ANLS(1,kav,lan))
#else
               tmp3(1,1,lan)  ,ANLS(1,kau,lan),ANLS(1,kav,lan),pdot)
#endif
#else /* ~HYBRID */
       do k = 1,levs_
         do j=1,lonsd2
           SYNS(j,kst-1+k,lan)=SYNS(j,kst-1+k,lan)-tov(k)
         enddo
       enddo
!
       call sph_nonlinear_tend(lonsd2,lan,                                     &
               SYNS(1,ksd,lan),SYNS(1,kst,lan),                                &
               SYNS(1,ksz,lan),SYNS(1,ksu,lan),SYNS(1,ksv,lan),                &
#ifndef NISLQ
               SYNS(1,ksr,lan),                                                &
#else
               slq(1,1,lan),                                                   &
#endif
               SYNS(1,kspphi,lan),SYNS(1,ksplam,lan),                          &
               rbs2(LATX),del,rdel2,ci,p1,p2,h1,h2,tov,                        &
               spdlat(1,LATX),                                                 &
#ifdef CLM_CWF
               LATX, cgs,                                                      &
#endif

               DYNS(1,kdtphi,lan),DYNS(1,kdtlam,lan),                          &
#ifndef NISLQ
               DYNS(1,kdrphi,lan),DYNS(1,kdrlam,lan),                          &
#else
               tmp1(1,1,lan)     ,tmp2(1,1,lan)     ,                          &
#endif
               DYNS(1,kdulam,lan),DYNS(1,kdvlam,lan),                          &
               DYNS(1,kduphi,lan),DYNS(1,kdvphi,lan),                          &
               ANLS(1,kap,lan),ANLS(1,kat,lan),                                &
#ifndef NISLQ
               ANLS(1,kar,lan),ANLS(1,kau,lan),ANLS(1,kav,lan))
#else
               tmp3(1,1,lan)  ,ANLS(1,kau,lan),ANLS(1,kav,lan),dot)
#endif
#endif /* HYBRID end */
!
#ifdef STOCH
#if defined(NLT_U) || defined(NLT_V) || defined(NLT_T) || defined(NLT_Q)
! stochastic forcing in nonlinear tendency
!
        do k = 1,levs_
          do i = 1,lonsd2
#ifdef NLT_U
            ANLS(i,kau+k-1,lan)=ANLS(i,kau+k-1,lan)*stc_u(i,k,lan)
#endif
#ifdef NLT_V
            ANLS(i,kav+k-1,lan)=ANLS(i,kav+k-1,lan)*stc_v(i,k,lan)
#endif
#ifdef NLT_T
            ANLS(i,kat+k-1,lan)=ANLS(i,kat+k-1,lan)*stc_t(i,k,lan)
#endif
#ifdef NLT_Q
#ifndef NISLQ
            ANLS(i,kar+k-1,lan)=ANLS(i,kar+k-1,lan)*stc_q(i,k,lan)
#endif
#endif
          enddo
        enddo
#endif /* NLT end */
#endif /* STOCH end */
!
#ifdef NISLQ
!      dyn_grid for nislq
!
!      Psfc
!
       do i = 1,lonsd2
         slq_psfc2(i,LATX)=exp(SYNS(i,ksp,lan))
       enddo
!
!      U, V
!
       do k = 1,levs_
         do i=1,lonsd2
           slq_u2(i,levs_+1-k,LATX)=SYNS(i,ksu+k-1,lan)
           slq_v2(i,levs_+1-k,LATX)=SYNS(i,ksv+k-1,lan)
         enddo
       enddo
!
!      Pdot
!
       do k = 1,levs_+1
         do i = 1,lonsd2
#ifdef HYBRID
!
!          pdot:top to bottom
!
           slq_w2(i,k,LATX)=pdot(i,k)
#else
!
!          dot:bottom to top
!
           slq_w2(i,levs_+2-k,LATX)=-dot(i,k)*exp(SYNS(i,ksp,lan))
#endif
         enddo
       enddo
#endif /* NISLQ end */
!
#ifdef PSPLIT
#ifndef HYBRID
       do k = 1,levs_
         do j = 1,lonsd2
           SYNS(j,kst-1+k,lan)=SYNS(j,kst-1+k,lan)+tov(k)
         enddo
       enddo
#endif
       if(.not.stepone) then
         call phys_main_solver(lonsd2,                                         &
                 SYNS(1,ksplam,lan),SYNS(1,kspphi,lan),                        &
                 SYNS(1,ksu,lan),SYNS(1,ksv,lan),                              &
                 SYNS(1,ksp,lan),SYNS(1,kst,lan),                              &
                 SYNS(1,ksr,lan),SYNS(1,ksd,lan),                              &
                 SYNS(1,ksplap,lan),                                           &
                 dtdtm(1,1,LATX),drdtm(1,1,LATX),                              &
                 dudtm(1,1,LATX),dvdtm(1,1,LATX),                              &
#ifdef DG
                 tgmxl(lan),igmxl(lan),kgmxl(lan),                             &
                 tgmnl(lan),igmnl(lan),kgmnl(lan),                             &
#endif
#ifdef DG3
                gda(1,1,lan),                                                  &
#endif
#ifdef RAS
                 ras,levs_,cp,alhl,grav,rgas,                                  &
                 sig, sgb, prh, prj, hpk, fpk, ods, prns,                      &
                 rasal, lm, krmin, krmax, nstrp,                               &
                 ncrnd, rannum, afac, ufac,                                    &
#endif
#ifdef RASV2
                 ras,rgas, cp, grav, alhl,                                     &
                 sig,prj,sgb,rasal, rannum, dsfc(1,LATX),                      &
                 pdd,krmin, krmax, kfmax, ncrnd, mct,kctop,                    &
#endif
#ifndef VIC
                 LATX,1.0)
#else
                 LATX,1.0,idate)
#endif
       endif
!
#ifdef DG
       if(tgmxl(lan).gt.tgmx) then
         tgmx=tgmxl(lan)
         igmx=igmxl(lan)
         kgmx=kgmxl(lan)
         jgmx=lat
       else if(tgmnl(lan).lt.tgmn) then
         tgmn=tgmnl(lan)
         igmn=igmnl(lan)
         kgmn=kgmnl(lan)
         jgmn=lat
       endif
#endif
#endif
#ifdef MP
     endif
#endif
!
   enddo
!
#ifdef PSPLIT
#ifdef DG3
#ifdef MP
#define LATX lan
#else
#define LATX lat
#endif
   do lat = lat1,lat2
     lan=lat-latdon
     call diag_3d_arrange(LATX,nwgda*kdgda,gda(1,1,lan))
   enddo
#endif
#endif
#undef SYNS
#undef DYNS
#undef ANLS
#undef LATX
!
#ifdef MP
#ifndef NISLQ
   klev=3+ntotal_
#else
   klev=3
#endif
   call mpnx2nk(gra,lonf22p_,lota,                                             &
                anf,lonf22_ ,lotas,latg2p_,levs_,levsp_,kau,kaus,klev)
   call mpnx2x (gra,lonf22p_,lota,                                             &
                anf,lonf22_ ,lotas,latg2p_,kap,kaps,1)
   if( levlen(mype).lt.levsp_ ) then
     do lan = 1,lats2
       do i = 1,lonf22_
         anf(i,kaus+levsp_-1,lan)=0.0
         anf(i,kavs+levsp_-1,lan)=0.0
         anf(i,kats+levsp_-1,lan)=0.0
       enddo
     enddo
#ifndef NISLQ
     do lan = 1,lats2
       do i = 1,lonf22_
         do kc = 1,ntotal_
           anf(i,kars+levsp_*kc-1,lan)=0.0
         enddo
       enddo
     enddo
#endif
   endif
!
   lat1=jstr
   lat2=jend
   latdon=jstr-1
#define ANLS anf
#else
#define ANLS anl
#endif
#ifdef ORIGIN_THREAD
!$doacross share(latdon,lat1,lat2,anl,anf,latdef,lcapd,lonfd),
!$&        local(lat,lan,lcapf,lonff)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(latdon,lat1,lat2,anl,anf,latdef,lcapd,lonfd)
!mic$1 private(lat,lan,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,lcapf,lonff)
#endif
   do lat = lat1,lat2
     lan=lat-latdon
#ifdef REDUCE_GRID
     lcapf=lcapd(latdef(lat))
     lonff=lonfd(latdef(lat))
#else
     lcapf=jcap1_
     lonff=lonf_
#endif
     call sph_wave2grid(ANLS(1,1,lan),ANLS(1,1,lan),2*lotas,                   &
                  lcapf,lonff,latdef(lat),-1)
   enddo
#undef ANLS
!
#ifdef MP
   call mpny2nl(anf,lonf22_ ,latg2p_,anl,lcap22p_,latg2_ ,lotas,kaps,lotas)
!
   lat1=1
   lat2=latg2_
   latdon=0
#define ZS za
#else
#define ZS z
#endif
#ifdef ORIGIN_THREAD
!$doacross share(latdon,lat1,lat2,anl,flp,flm,llens,lcapdp,mype),
!$&        local(lat,lan,llensd)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(latdon,lat1,lat2,anl,flp,flm,llens,lcapdp,mype)
!mic$1 private(lat,lan,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,llensd)
#endif
   do lat = lat1,lat2
     lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
     llensd=lcapdp(lat,mype)
#else
     llensd=lcapd(lat)
#endif
#else
     llensd=llens
#endif
!
     call sph_grid2wave(flp(1,1,1,lan),flm(1,1,1,lan),anl(1,1,lan),            &
                          llensd,lotas)
   enddo
! 
!  no multi-threads should be given in following loop
!  otherwise results will be non-produceable
!
   do lat = lat1,lat2
     lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
     llensd=lcapdp(lat,mype)
#else
     llensd=lcapd(lat)
#endif
#else
     llensd=llens
#endif
!
#define DEFAULT
#ifdef FL2I
#undef DEFAULT
     call sph_comp_coeff1(flp(1,1,kaps,lan),flm(1,1,kaps,lan),ZS,qww(1,lat),   &
               llstr,llensd,lwvdef,lotas)
#endif
#ifdef DEFAULT
     call sph_comp_coeff2(flp(1,1,kaps,lan),flm(1,1,kaps,lan),ZS,qww(1,lat),   &
               llstr,llensd,lwvdef,lotas)
#endif
!
     call sph_sum_wind_coeff(flp(1,1,kaus,lan),flm(1,1,kaus,lan),              &
                 flp(1,1,kavs,lan),flm(1,1,kavs,lan),                          &
                 anltop(1,1,1),anltop(1,1,LEVSS+1),                            &
                 qvv(1,lat),wgt(lat),                                          &
                 llstr,llensd,lwvdef,LEVSS)
   enddo
#undef ZS
!
#ifndef MP
   latdon=latdon+(lat2-lat1+1)
!---------------------------------------------------
10000 continue
!---------------------------------------------------
#endif
#ifdef ORIGIN_THREAD
!$doacross share(spdmax,spdlat,lats2),local(k,lat)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(spdmax,spdlat,lats2)
!mic$1 private(k,lat)
#endif
#ifdef OPENMP
!$omp parallel do private(k,lat)
#endif
   do k = 1,levs_
     spdmax(k) = 0.0
     do lat = 1,lats2
       spdmax(k)=max(spdmax(k),spdlat(k,lat))
     enddo
     spdmax(k)=sqrt(spdmax(k))
   enddo
#ifdef MP
   call mpgetspd(spdmax)
#endif
#ifndef NOPRINT
   if( iope ) then
     write(6,100)(spdmax(k),k=1,levs_)
100  format(' global checked speed maxima for all layers ',                    &
        :/' spdmx(001:010)=',10f5.0,:/' spdmx(011:020)=',10f5.0,               &
        :/' spdmx(021:030)=',10f5.0,:/' spdmx(031:040)=',10f5.0,               &
        :/' spdmx(041:050)=',10f5.0,:/' spdmx(051:060)=',10f5.0,               &
        :/' spdmx(061:070)=',10f5.0,:/' spdmx(071:080)=',10f5.0,               &
        :/' spdmx(081:090)=',10f5.0,:/' spdmx(091:100)=',10f5.0,               &
        :/' spdmx(101:110)=',10f5.0,:/' spdmx(111:120)=',10f5.0,               &
        :/' spdmx(121:130)=',10f5.0,:/' spdmx(131:140)=',10f5.0,               &
        :/' spdmx(141:150)=',10f5.0,:/' spdmx(151:160)=',10f5.0,               &
        :/' spdmx(161:170)=',10f5.0,:/' spdmx(171:180)=',10f5.0,               &
        :/' spdmx(181:190)=',10f5.0,:/' spdmx(191:200)=',10f5.0)
   endif
#endif
!
#ifndef DCMIP
   do k =1,levs_
     if(spdmax(k).eq.0.) then
       if (iope) write(6,*)'run failure.  spdmax=0.'
       call MPABORT
     endif
   enddo
#endif
!
!  input : w=d(u)/d(t) x=d(v)/d(t)
!  output: uln=d(di)/d(t) vln=d(ze)/d(t)
!
#ifdef MP
#define UUS uua
#define VVS vva
#define XS xa
#define WS wa
#else
#define UUS uu
#define VVS vv
#define XS x
#define WS w
#endif
   call sph_wind2divor(UUS,VVS,XS,WS,anltop(1,1,1),                            &
               anltop(1,1,LEVSS+1),llstr,llens,lwvdef)
#undef UUS
#undef VVS
#undef XS
#undef WS
!
#ifdef MP
#ifndef NISLQ
   klev=3+ntotal_  ! y,rt,x,w
#else
   klev=3          ! y,x,w
#endif
   call mpn2nn (za ,lln22p_,        z ,lnt22p_,      1)
   call mpnk2nn(za(1,kys),lln22p_,levsp_,z(1,ky),lnt22p_,levs_,klev)
#endif
!
!  subtract off linear dependence on divergence
!
   do k = 1,levs_
     do j = 1,levs_
       do i = 1,LNT2S
#ifndef HYBRID
         y(i,k)=y(i,k)-bm(k,j)*di(i,j)
#else
         y(i,k)=y(i,k)+bm(k,j)*di(i,j)
#endif
       enddo
     enddo
   enddo
#ifdef HYBRID
   do j = 1,levs_
     do i = 1,LNT2S
       z(i,1)=z(i,1)+sv(j)*di(i,j)
     enddo
   enddo
#endif
!
! add topog. contriba into x.
! integrate vorticity amd moisture in time
!
   do k = 1,levs_
     do i = 1,LNT2S
       x(i,k)=x(i,k)+gz(i)
       w(i,k)=zem(i,k)+2.*deltim*w(i,k)
     enddo
   enddo
#ifndef NISLQ
   do k = 1,levh_
     do i = 1,LNT2S
       rt(i,k)= rm(i,k)+2.*deltim* rt(i,k)
     enddo
   enddo
#endif /* NISLQ end */
!
   if( iope ) then
     do k = 1,levs_
       w(1,k)=0.e0
       w(2,k)=0.e0
     enddo
   endif
#ifdef PSPLIT
#ifdef DG
   if(iope) then
     write(6,'("phys_main_driver t range ",                                    &
                  2(4x,f6.1," @i,k,lat ",3i4))')                               &
           tgmx,igmx,kgmx,jgmx,tgmn,igmn,kgmn,jgmn
   endif
#endif	/* DG */
#endif /* PSPLIT */
!=========================================================================
! END - SMP or GMP
!=========================================================================
#endif 		/* SMP */
!
#ifndef SMP
#ifdef MP
   deallocate( syf,grs,dyf,dgr,anf,gra )
#endif
   deallocate( syn,syntop,dyn,scs2,anl,anltop,flp,flm )
#endif
   return
   end subroutine dyn_sph_driver
!
