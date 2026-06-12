#include <define.h>
!-------------------------------------------------------------------------------
   subroutine gmp_start_setup(n1,                                              &
#ifdef RMP
              nrsmi1,nrsmi2,nrflip,                                            &
              nrsmo1,nrsmo2,nrflop,nrsfli,nrsflx,nrinit,nrpken,                &
#endif
#ifdef LFM
              nlfmsgi,nlfmsfi,nlfmsgo,nlfmsfo,klenp,weight,fwindow,            &
#endif
              kpfix,ksfcx,komlx,ksig,ksfc,kpost,krestart,klfm,krsm)
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [gmp_start_setup]
!      |
!      |--- A-1) [dyn_sigma_setup] *
!      |    A-2) [sph_matrix_init] - [sph_matrix_init_sub] *
!      |--- B-1) [dyn_hybrid_setup] *
!      |    B-2) [sph_matrix_init_hybrid] *
!      |
!      |--- [solve_sph] *
!      |--- [dyn_sincos_lat] *
!      |--- [lsm_init_module]
!               |--- [lsm_init_dfkt]
!               |--- [lsm_init_df] 
!               |--- [lsm_init_kt]
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only : LONF2S,LATG2S,lnt2_,lonf_,latg_,lonf2_,latg2_,        &
                         levm1_,lnut2_,levs_,jcap_,jcap1_,                     &
                         nwater_,igen_,mtnvar_
#ifdef SMP
   use paramodel, only : lnt22_
   use comfcst, only    :  dtbdy, curtime, vvel  ,                             &
                           hour1,hour2,hours,houre,ftim1,ftim2
#ifdef CLM_CWF
   use comfcst, only    :  wdiv, hadq
#endif
#endif
   use constant, only  : pi_,rerth_,                                           &
                         rhoair0_, rhoh2o_, rhosnow_, cliq_, cvap_
   use comsfc          ! cv,cvb,cvt
#ifdef DFS
   use dfsvar, only  : tov=>TAVEXY,coslat,                                     &
                       VERMAT,EIGVAL,EIGVEC,AEIGVEC,                           &
                       sl=>sigma,si=>sigmafull,del=>delsig,                    &
                       AMATm,DMATm,aMSQUAR,DMATm_,AMATm_,jcol2js,              &
                       jlha,jls,mt,jgs,iope,ib,jb,latdef,                      &
                       jbw,iba,jbwa,xlond=>xlon,xlatd=>xlat
#endif
#ifdef MP
   use commpi          ! mype
   use paramodel, only : lnt2p_,lonf2p_,lln2p_,latg2p_,levp1_
#endif
   use comfibm
   use comcon
   use comgpd
#ifdef REDUCE_GRID
   use comreduce
#endif
#ifdef WSM3
   use module_mp_wsm3
#endif
#ifdef WSM5
   use module_mp_wsm5
#endif
#ifdef WSM6
   use module_mp_wsm6
#endif
#ifdef WDM5
   use module_mp_wdm5
#endif
#ifdef WDM6
   use module_mp_wdm6
#endif
#ifdef REDUCE_GRID
   use module_sph_reduced_grid, only : sph_reduced_grid_init
#endif
#ifdef DGP
   use diag_point_module, only       : diag_point_prepare
#endif
   use lsm_init_module, only         : dfkt_init, lsm_init_df, lsm_init_kt
   use module_trans, only            : dyn_trans2model_grid
   use module_sph_legendre
#ifdef NISLQ
   use nislq, only                   : nislq_init
#endif
#ifdef STOCH
   use stochastic, only              : stc_init
#endif
#ifndef SWRMDC
#ifndef ICECLOUD
   use rdparm, only   : imbx, lp1
#else
   use rdparm8, only   : imbx, lp1
#endif
#else
   use rdparm99, only   : imbx, lp1
#endif
#ifdef RRTMGLW
   use module_ra_rrtmg_lw, only : rrtmg_lwinit
#endif
#ifdef RRTMGSW
   use module_ra_rrtmg_sw, only : rrtmg_swinit
#endif
   use comconsts,   only        : kbmaxin,kbmin,kmaxin,levshcin,kpblmaxin
   use comconsts,   only        : ictin,icbin
#ifdef DCMIP
   use dcmip_grims, only : dcmip_grims_init
#endif /* DCMIP end */
!-------------------------------------------------------------------------------
#ifndef HYBRID
   implicit none
!
#endif
#include "abort.h"
#ifdef MP
   integer              ::  ilist(100)
   real                 ::  rlist(100)
#ifdef DFS
   real                 ::  grid1d(lonf_,latg_)
   real                 ::  grid2d(lonf_,latg_)
#else
   real                 ::  rdex(lnt2_)
#endif
   real                 ::  grid1(lonf2_,latg2_)
   real                 ::  grid2(lonf2_,latg2_)
   real                 ::  gridv(lonf2_,latg2_,mtnvar_)
#endif
   integer              ::  ndigit
   data ndigit/1/
#ifndef DFS
#ifdef MP
   real                 ::  tmpq(lnt2p_)
#endif
   real                 ::  tmpqtt(lnt2_),tmpqdd(lnt2_),tmpqvv(lnut2_)
#endif
!
#ifdef SMP
#ifndef CLM_CWF
   integer              ::  lsmask
   real                 ::  rad
#endif
#endif
!
#ifdef LFM
   real weight(*)
   real weix(1000)
#endif
!
   real endhour,dt80,dlat,cvmint
   integer i,labl,ldebug,ncpus,jcap,levs,n1,nmtnv,lev,j,jj,k                   &
          ,intpfix,intsfcx,intomlx,intsig,intsfc,intflx,intrestart             &
          ,intlfm,intrsm,kpfix,ksig,ksfc,kpost,ncpus1,klfm,krsm
       integer krestart,ksfcx,komlx,ncldb1
       integer ind,l,ll,n,maxi,lat,llstr,llens,len
       real    fact, qmaxall,qttcut, wcsa
#ifdef RMP
   integer              ::  nrsmi1,nrsmi2,nrflip,                              &
                            nrsmo1,nrsmo2,nrflop,nrsfli,nrsflx,nrinit,nrpken
#endif
   namelist/namsmf/ con,num,labl,endhour,ldebug,filta,icen,igen,icen2          &
                   ,ienst,iensi,runid,usrid,ncpus                              &
                   ,intpfix,intsfcx,intomlx,intsig,intsfc,intflx,intrestart    &
#ifdef LFM
                   ,intlfm,intrsm,ndigit,critfs,filtwin
   real                 ::  critfs,filtwin
#else
                   ,intlfm,intrsm,ndigit
#endif
   integer              ::  ios
#if defined (RRTMGSW) || defined (RRTMGLW) || defined (GSFCSW)
   integer              ::  ids,ide, jds,jde, kds,kde,                         &
                            ims,ime, jms,jme, kms,kme,                         &
                            its,ite, jts,jte, kts,kte
#endif
   logical, parameter   ::  allowed_to_read = .true.
!-------------------------------------------------------------------------------
   limlow=1
   jcap=jcap_
   levs=levs_
   filta= 0.92
   dt80=939.14e0
   percut=27502.e0
   icen=7
   igen=igen_
   icen2=0
   ienst=0
   iensi=0
   runid=0
   usrid=0
   call dyn_cpu_number(ncpus)
!
!  define unit numbers
!
!  input
   n1    = 11
   nmtnv = 24
#ifdef LFM
   nlfmsgi=30
   nlfmsfi=31
#endif
#ifdef RMP
   nrsmi1= 30
   nrsmi2= 31
   nrflip= 32
#endif
#ifdef LFM
   nlfmsgo= 70
   nlfmsfo= 71
#endif
#ifdef RMP
   nrsmo1= 70
   nrsmo2= 71
   nrflop= 72
   nrsflx= 73
   nrsfli= 74
   nrpken= 77
   nrinit= 78
#endif
!
!  output interval control
!
   intpfix=24
   intsfcx=24
   intomlx=24
   intsig=24
   intsfc=24
   intflx=24
   intrestart=24
   intlfm=24
   intrsm=24
#if defined (RRTMGSW) || defined (RRTMGLW) || defined (GSFCSW)
!
! define the dimension constants
!  
   ids = 1
   ide = imbx
   jds = 1
   jde = 1
   kds = 1
   kde = levs
!  
   ims = 1
   ime = imbx
   jms = 1
   jme = 1
   kms = 1
   kme = lp1
!
   its = 1
   ite = imbx
   jts = 1
   jte = 1
   kts = 1
   kte = levs
#endif
!
!.... cmean,clstp control time averaging of convective clds in kuo
!
   clstp=99.
!
!...  averaging interval for conv cld approx 3 hrs (num of timesteps)
!
#ifndef NOPRINT
   if(iope) write(6,100) jcap, levs
100  format (1x,'gmp_start_setup jcap levs ',i4,i4,' created Jan 2011')
#endif
   filtb =(1.e0-filta) * 0.5e0
#ifdef DCMIP
!
! dcmip
!
#ifdef HYBRID
   call dcmip_grims_init(ak5,bk5)
#else
   call dcmip_grims_init(si,sl)
#endif
#endif /* DCMIP end */
!
   if( iope ) then
#ifdef SMP
     open(unit=n1,file='basic.info',form='unformatted',status='old')
     read(n1) dtbdy, xlon, xlat, lsmask, ci, si, del, sl, cl, rpi
#ifdef DBG
     write(6,*) ' gmp_start_setup basic info'
     write(6,*) dtbdy, xlon, xlat, lsmask, ci, si, del, sl, cl, rpi
#endif
     close(n1)
     rad = pi_/180.
     do i = 1, lonf2_
       do j = 1, latg2_
         xlon(i,j) = xlon(i,j)*rad
         xlat(i,j) = xlat(i,j)*rad
         write(6,'(2F7.3)') xlon(i,j), xlat(i,j)
       enddo
     enddo
#else /* ~SMP */
#ifdef DFS
#ifdef BARO_TEST
     read(n1,*)si(:)
     do k = 1,levs_
       sl(k)=(si(k)+si(k+1))*0.5
       del(k)=si(k)-si(k+1)
     enddo
     read(n1,*)TOV(:)
#else
#ifndef HYBRID
     call dyn_sigma_setup(ci,si,del,sl,cl,rpi,tov,n1)
#else
     call dyn_hybrid_setup(ak5,bk5,ci,si,del,sl,cl,rpi,tov,n1)
#endif
#endif /* BARO_TEST end */
#else /* SPH */
#ifndef HYBRID
     call dyn_sigma_setup(ci,si,del,sl,cl,rpi,n1)
#else
     call dyn_hybrid_setup(ak5,bk5,ci,si,del,sl,cl,rpi,n1)
#endif
#endif /* DFS end */
#endif /* SMP */
   endif
#ifdef MP
   call mpbcastr(ci,levp1_)
   call mpbcastr(si,levp1_)
   call mpbcastr(del,levs_)
   call mpbcastr(sl,levs_)
   call mpbcastr(cl,levs_)
   call mpbcastr(rpi,levm1_)
   call mpbcastr(tov,levs_)
#ifdef HYBRID
   call mpbcastr(ak5,levp1_)
   call mpbcastr(bk5,levp1_)
#endif
#endif
   sl1=sl(1)
#ifndef DFS
   do lev=1,levs_
     tov(lev)=300.0
   enddo
#endif
#ifdef SMP
   rcs2(1) = 1.0
   wgt(1) = 1.0
#else
#ifdef DFS
!
! calculate matrix and dissipation coef.
!
   call dfs_constant
!
! calculate matrix for simi-implicit
!
#ifdef HYBRID
   call dfs_vertical_descr_hybrid(VERMAT, EIGVAL, EIGVEC, AEIGVEC,levs_,ak5,bk5)
#else
   call dfs_vertical_descr(VERMAT, EIGVAL, EIGVEC, AEIGVEC,levs_)
#endif
#ifdef DBG
   if (iope) then
     write(6,*)'EIGVAL=',EIGVAL
     !write(6,*)'TOV=',TOV(:)
   endif
#endif
!
! calculate coslat, sinlat
!
   call dfs_sincos_lat
!
! aditionally colatitudes
!
   dlat=pi_/latg_
   do j = 1,latg2_
     colrad(j)=(j-0.5)*dlat
   enddo
   do j = 1,latg2_
     rcs2(j)=coslat(j,3)
   enddo
   do j = 1,jb
     jj=latdef(jls+j-jgs)
     colrab(j)=colrad(jj)
     rbs2  (j)=rcs2  (jj)
   enddo
!
! Matrices for Laplacian and inversion
!
   call dfs_semi_matrix(AMATm,DMATm,aMSQUAR,DMATm_,AMATm_,jcol2js,             &
                     mt,jlha)
!CALL test_fft
#else
#ifndef HYBRID
   call sph_matrix_init(levs_,si,sl,tov,am,bm,sv,gv,cm)
#else
   call sph_matrix_init_hybrid(levs_,ak5,bk5,am,bm,sv,gv,cm)
#endif 
   call sph_gaussian_lat(latg2_, colrad, wgt, wgtcs, rcs2)
#endif
#endif			/* SMP */
#ifndef DFS
#ifdef MP
   do j = 1,latlen(mype)
     jj=latdef(latstr(mype)+j-1)
     colrab(j)=colrad(jj)
     wgb   (j)=wgt   (jj)
     wgbcs (j)=wgtcs (jj)
     rbs2  (j)=rcs2  (jj)
   enddo
#else
   do j = 1,latg2_
     jj=latdef(j)
#ifndef SMP
     colrab(j)=colrad(jj)
     wgb   (j)=wgt   (jj)
     wgbcs (j)=wgtcs (jj)
#endif
     rbs2  (j)=rcs2  (jj)
   enddo
#endif /* MP end */
#ifndef SMP
   call sph_poly_epsilon1(eps, jcap_)
#endif
#endif /* ~DFS end */
#ifdef NISLQ
!
! initialize nislq
!
#ifdef MP
   call nislq_init(nrow,myrow,colrad,rbs2)
#else
   call nislq_init(1,0,colrad,rbs2)
#endif
#endif /* NISLQ end */
#ifdef STOCH
   call stc_init(idate)
#endif
!
!     rpi(k) = (sl(k+1)/sl(k))**rk  from dyn_sigma_setup  k=1...levm1_
!
   do k = 1,levm1_
     rpirec(k) = 1.e0/rpi(k)
   enddo
   do k = 1,levs_
     rdel2(k)=0.5e0/del(k)
   enddo
#ifndef DFS
   ind=0
   do ll=1,jcap1_
     n=ll-2
     maxi=jcap1_+1-ll
     do i = 1,maxi
       ind=ind+1
       n=n+1
       ndex(ind*2-1) = n
       ndex(ind*2  ) = n
       fact=float(n*(n+1))
       snnp1(ind*2-1) = fact
       snnp1(ind*2  ) = fact
     enddo
   enddo
#ifdef MP
   call spcshfli(snnp1,lnt2_,1,jcap_,lwvdef)
   do n = 1,lnt2_
     rdex(n)=ndex(n)
   enddo
   call spcshfli(rdex,lnt2_,1,jcap_,lwvdef)
   do n = 1,lnt2_
     ndex(n)=rdex(n)
   enddo
#endif
#endif                  /* not DFS */
!
!  initialize cv, cvt and cvb
!
   do j = 1,LATG2S
     do i = 1,LONF2S
       cv (i,j) = 0.e0
       cvt(i,j) = 0.e0
       cvb(i,j) = 100.e0
     enddo
   enddo
!
!  temporarily set some cons and nums 
!
   num(1:28)=0
   num( 1) = 11
   num( 2) = 11
   num( 3) = 51
   num( 4) = 52
   num( 5) =  0
   num( 6) =  1
   num( 7) =  0
   num( 8) =  1
   num( 9) =  8
   num(10) = 15
   num(11) =  1
   num(12) = 23
   num(13) =  0
   num(14) = 55
   num(15) =  0
   num(16) = 11
   num(17) = 51
   num(18) =  4
   num(19) =  2
   num(20) =  6
   num(21) = 15
   num(22) = 10
   num(23) =  1
   num(24) =  0
   num(25) =  0
   num(26) =  0
   num(27) =  0
   num(28) =  0
   con(1)=0.
   con(3)=0.        ! gmp dyn_digital_filter initialization in hour or not (0.)
   con(4)=1.        ! dtswav in hour for gmp
   con(5)=3.        ! dtlwav in hour for gmp
   con(6)=0.
   con(7)=12.
   con(17)=120.     ! forecast ending hour
#ifdef RMP
   con(11)=400.     ! rmp deltim
   con(12)=21600.   ! rmp nesting period in second
   con(13)=6.0      ! rmp initialization step in hour or not (0.)
   con(14)=1.       ! rmp dtswav in hour
   con(15)=1.       ! rmp dtlwav in hour
   con(16)=0.       ! rmp start forecast period
   con(17)=36.      ! rmp ending forecast period
   con(18)=0.       ! rmp local diffusion (1) or not (0)
   con(19)=0.       ! rmp lateral boundary relaxation 1 blending 0
#endif
   num(31)=1
   num(32)=0
   num(30)=0
   endhour=0.
!
!  replace num con from namelist input
!
   if( iope ) then
     rewind(95)
     read(95,namsmf)
!
#ifdef MP
     ilist(1)=labl
     ilist(2)=ldebug
     ilist(3)=icen
     ilist(4)=igen
     ilist(5)=icen2
     ilist(6)=ienst
     ilist(7)=iensi
     ilist(8)=ncpus
     ilist(9)=intpfix
     ilist(10)=intsfcx 
     ilist(11)=intomlx 
     ilist(12)=intsig
     ilist(13)=intsfc
     ilist(14)=intflx
     ilist(15)=intrestart
     ilist(16)=intlfm
     ilist(17)=intrsm
     ilist(18)=ndigit
     rlist(1)=endhour
     rlist(2)=filta
     rlist(3)=runid
     rlist(4)=usrid
#ifdef LFM
     rlist(5)=critfs
     rlist(6)=filtwin
#endif
#endif
   endif
#ifdef MP
!
   call mpbcastr(con,1700)
   call mpbcasti(num,1700)
   call mpbcasti(ilist,18)
#ifdef LFM
   call mpbcastr(rlist,6)
#else
   call mpbcastr(rlist,4)
#endif
   labl   =ilist(1)
   ldebug =ilist(2)
   icen   =ilist(3)
   igen   =ilist(4)
   icen2  =ilist(5)
   ienst  =ilist(6)
   iensi  =ilist(7)
   ncpus  =ilist(8)
   intpfix=ilist(9)
   intsfcx=ilist(10)
   intomlx=ilist(11)
   intsig =ilist(12)
   intsfc =ilist(13)
   intflx =ilist(14)
   intrestart =ilist(15)
   intlfm =ilist(16)
   intrsm =ilist(17)
   ndigit =ilist(18)
   endhour=rlist(1)
   filta  =rlist(2)
   runid  =rlist(3)
   usrid  =rlist(4)
#ifdef LFM
   critfs =rlist(5)
   filtwin=rlist(6)
#endif
!
#endif
!
!  adjust constants with fhour from sig file
!
#if !defined(BARO_TEST) && !defined(DCMIP)
   if( iope ) then
     close(n1)
     open (unit=n1,file='sigit',form='unformatted',err=999)
     go to 998
 999 continue
     write(6,*)'error opening sigit in gmp_start_setup'
     call MPABORT
 998 continue
     rewind n1
     read(n1)
     read(n1) fhour
   endif
#else
   fhour=0.0
#endif

#ifdef MP
   call mpbcastr(fhour,1)
#endif
   if(fhour.eq.0.) then
     if(num(5).eq.-1) num(5)=2
     if(num(5).eq.-2) num(5)=1
   else
     if(num(5).eq.-1) num(5)=0
     if(num(5).eq.-2) num(5)=0
     if(con(17)-fhour.lt.con(7)) then
       con(7)=con(17)-fhour
     endif
   endif
   if(con(1).le.0.) con(1)=dt80 *80./jcap
   if(num(7).le.0.and.con(7).ne.0.) then
     num(7)=3600.*con(7)/con(1)+0.99
     con(1)=nint(3600.*con(7)/num(7))
   else
     con(7)=num(7)*con(1)/3600.
   endif
   if(num(32).eq.0) num(32)=num(7)
   if(num(1).gt.0) con(6)=num(1)
   dtswav=con(4)
   dtlwav=con(5)
#ifndef DCMIP
   if(dtswav.gt.float(intflx)) then
     write(6,*)'dtswav.gt.intflx. dtswav=',dtswav,' intflx=',intflx
     call MPABORT
   endif
   if(dtlwav.gt.float(intflx)) then
     if (iope)                                                                 &
       write(6,*)'dtlwav.gt.intflx. dtlwav=',dtlwav,' intflx=',intflx
     call MPABORT
   endif
!
!  check dtswav and dtlwav if they are reasonable
!
   if(mod(int(dtswav*3600.),int(con(1))).ne.0) then
     if(iope) then
       write(6,*)'dtswav must be a multiple of timestep'
       write(6,*)'dtswav,con(1)=',dtswav*3600,con(1)
     endif
     call MPABORT
   endif
   if(mod(24*60,int(dtswav*60.)).ne.0) then
     if(iope)                                                                  &
       write(6,*)'24*60 must be must be a multiple of dtswav*60'
     call MPABORT
   endif
   if(mod(int(dtlwav*3600.),int(con(1))).ne.0) then
     if(iope)                                                                  &
       write(6,*)'dtlwav must be a multiple of timestep'
     call MPABORT
   endif
   if(mod(24*60,int(dtlwav*60.)).ne.0) then
     if(iope)                                                                  &
       write(6,*)'24*60 must be must be a multiple of dtlwav*60'
     call MPABORT
   endif
!
   cowave=0.
   dtwave=0.
!
!  cvmint - maximum conv. cld accumulation time interval in hours
!
   cvmint = con(4)
   dtcvav =  min (cvmint,  max (dtswav,dtlwav))
#endif /* ~DCMIP */
!
#ifdef DBG
   if(iope) write(6,201)(num(i),i=1,28)
201   format(1h0,'num=',28(1x,i2))
   if(iope) write(6,*)'con'
   if(iope) write(6,*)(con(i),i=1,10)
#endif
!
!  compute the output and cycling frequency
!
   kpfix=3600*intpfix/con(1)+0.5
   ksfcx=3600*intsfcx/con(1)+0.5
   komlx=3600*intomlx/con(1)+0.5
   ksig=3600*intsig/con(1)+0.5
   ksfc=3600*intsfc/con(1)+0.5
   kpost=3600*intflx/con(1)+0.5
   krestart=3600*intrestart/con(1)+0.5
   krsm=3600*intrsm/con(1)+0.5
   klfm=3600*intlfm/con(1)+0.5
!
   dk=num(9)
   dk=dk*(10.e0)**num(10)
   tk=num(20)
   tk=tk*(10.e0)**num(21)
   if(num(20).eq.0)tk=dk
#ifndef NOPRINT
   if(iope) write(6,105)con(1),filta,dk,tk
105   format(2x,'deltim = ',f5.0,1x,'filta = ',f4.2,1x,'dk = ',e8.2,1x,    &
            'tk = ',e8.2)
#endif
#ifdef DGP
   npoint=num(1300)
   if(npoint.lt.0.or.npoint.gt.nptken) then
     write(6,*)'DGP points disabled - grid points exceed ',nptken
     npoint=0
   endif
   isave=0
   itnum=0
   if(npoint.ne.0) then
     isave=1
     itnum=1
     isshrt=num(1301)
     ilshrt=num(1302)
     ikfreq=num(1303)
     call diag_point_prepare(con,colrad,lonf_,latg2_,n1)
   endif
#endif
   ncpus1=ncpus+1
   ncldb1=ncpus*lonf2_/lonf2_+1
! 
!     call  sph_comp_index  to set common/comind/ for subs. transi,sph_matrix_trans.
#if !defined(SMP) && !defined(DFS)
   call  sph_comp_index
#ifdef DBG
   if(iope) write(6,*)' done sph_comp_index '
#endif
#endif
! 
   call funct_svp_init
#ifdef DBG
   if(iope) write(6,*)' done funct_svp_init '
#endif
   call funct_dew_point_temp_init
#ifdef DBG
   if(iope) write(6,*)' done funct_dew_point_temp_init '
#endif
   call funct_pot_temp_init
#ifdef DBG
   if(iope) write(6,*)' done funct_pot_temp_init '
#endif
   call funct_moist_adiabat_init
#ifdef DBG
   if(iope) write(6,*)' done funct_moist_adiabat_init '
#endif
#if !defined(DFS) && !defined(SMP)
   call sph_poly_funct_init
#ifdef DBG
   if(iope) write(6,*)' done sph_poly_funct_init '
#endif
   call sph_poly_epsilon2(epsi,jcap_)
#ifdef DBG
   if(iope) write(6,*)' done sph_poly_epsilon2 '
#endif
   call sph_derivative_init(epsi)
#ifdef DBG
   if(iope) write(6,*)' done sph_derivative_init '
#endif
#endif			          /* not SMP  .and. not DFS */
#ifdef REDUCE_GRID
   qmaxall=0.0
#endif
!
!--for the rad initialization--------------------
!
#ifdef RRTMGLW
   call rrtmg_lwinit(                                                          &
                      allowed_to_read ,                                        &
                      ids, ide, jds, jde, kds, kde,                            &
                      ims, ime, jms, jme, kms, kme,                            &
                      its, ite, jts, jte, kts, kte                             )
#endif
#ifdef RRTMGSW
   call rrtmg_swinit(                                                          &
                      allowed_to_read ,                                        &
                      ids, ide, jds, jde, kds, kde,                            &
                      ims, ime, jms, jme, kms, kme,                            &
                      its, ite, jts, jte, kts, kte                             )
#endif
!
!--for the mps initialization--------------------
!
   if (allowed_to_read) then
#ifdef WSM3
     if (nwater_.eq.3) call wsm3init(rhoair0_, rhoh2o_, rhosnow_, cliq_, cvap_,&
                           allowed_to_read)
#endif
#ifdef WSM5
     if (nwater_.eq.5) call wsm5init(rhoair0_, rhoh2o_, rhosnow_, cliq_, cvap_,&
                           allowed_to_read)
#endif
#ifdef WSM6
     if (nwater_.eq.6) call wsm6init(rhoair0_, rhoh2o_, rhosnow_, cliq_, cvap_,&
                           allowed_to_read)
#endif
#ifdef WDM5
     if (nwater_.eq.8) call wdm5init(rhoair0_, rhoh2o_, rhosnow_, cliq_, cvap_,&
                           allowed_to_read)
#endif
#ifdef WDM6
     if (nwater_.eq.9) call wdm6init(rhoair0_, rhoh2o_, rhosnow_, cliq_, cvap_,&
                           allowed_to_read)
#endif
   endif
!
#ifndef DFS
#ifndef SMP
   do lat = 1,latg2_
!
     call sph_poly_funct(tmpqtt,tmpqvv,colrad,lat)
#ifdef REDUCE_GRID
     do n = 1,lnt2_
       qmaxall=max(qmaxall,abs(tmpqtt(n)))
     enddo
#endif
     call solve_sph(tmpqtt,tmpqvv,tmpqdd,epsi)
     wcsa=rcs2(lat)/rerth_
#ifdef MP
     call mpsf2p (tmpqtt,lnt2_,tmpq,lnt2p_,1)
     call mpnn2n (tmpq,lnt2p_,qtt(1,lat),lln2p_,1)
     call mpsf2p (tmpqdd,lnt2_,tmpq,lnt2p_,1)
     call mpnn2n (tmpq,lnt2p_,qdd(1,lat),lln2p_,1)
     do n = 1,lln2p_
       qww(n,lat)=qtt(n,lat)*wgt(lat)
       qdd(n,lat)=qdd(n,lat)*wcsa
     enddo
     llstr=lwvstr(mype)
     llens=lwvlen(mype)
#else
     do n = 1,lnt2_
       qtt(n,lat)=tmpqtt(n)
       qww(n,lat)=tmpqtt(n)*wgt(lat)
       qdd(n,lat)=tmpqdd(n)*wcsa
     enddo
     llstr=0
     llens=jcap1_
#endif
     len=2*jcap1_
     j=len
     do l = 0,jcap_
       tmpqvv(2*l+1)=tmpqvv(j+1)
       tmpqvv(2*l+2)=tmpqvv(j+2)
       j=j+len
       len=len-2
     enddo
     do ll = 1,llens
       j=lwvdef(llstr+ll)
       l=ll-1
       qvv(2*l+1,lat)=tmpqvv(2*j+1)
       qvv(2*l+2,lat)=tmpqvv(2*j+2)
     enddo
!
   enddo
#ifdef DBG
   if(iope) write(6,*)' done qtt qdd qww and qvv '
#endif
#endif
#endif              /* not DFS */
#ifdef REDUCE_GRID
   qttcut=qmaxall/(10.**ndigit)
   if( ndigit.eq.0 ) qttcut=0.0
   if( iope ) then
   write(6,*)' reduce grid is on with ',ndigit,' digit accuracy.'
   endif
!
   do lat = 1,latg2_
     call sph_poly_funct(tmpqtt,tmpqvv,colrad,lat)
!
!  dyn_trans2model_grid, dyn_trans2output_grid depend on sph_reduced_grid
!
     call sph_reduced_grid_init(tmpqtt,lnt2_,jcap_,qttcut,                     &
                   lcapd(lat),lonfd(lat))
#ifdef DBG
     if( iope ) then
     write(6,*)' --- lat=',lat,' needs local jcap and lonf as ',               &
                   lcapd(lat)-1,lonfd(lat)
     endif
#endif
   enddo
#ifdef MP
!
! --- mpgf2p and mpgp2f depend on preduceg
!
   call preduceg
   if( iope ) write(6,*)' done preduceg '
#endif
!
#endif
!
#ifndef DCMIP
! ------- finish setting of reduce grid --------------
!
#ifdef MP
#define HPRIMES gridv
#else
#define HPRIMES hprime
#endif
   if( iope ) then
     read(nmtnv) HPRIMES
#ifdef DBG
#ifndef NOPRINT
     call print_maxmin_six(HPRIMES,lonf2_,latg2_,1,latg2_,"mtnvar")
#endif
#endif
     call dyn_trans2model_grid(HPRIMES,mtnvar_)
   endif
#ifdef MP
   call mpgf2p(HPRIMES,lonf2_,latg2_,hprime,lonf2p_,latg2p_,mtnvar_)
#endif
#ifdef DBG
   if( iope ) write(6,*)' done hprime '
#endif
#undef HPRIMES
#endif /* ~DCMIP end */
!
#ifdef MP
#define XLONS grid1
#define XLATS grid2
#else
#define XLONS xlon
#define XLATS xlat
#endif
#ifndef SMP
   if( iope ) then
     call dyn_lon_lat(XLONS,XLATS,colrad,lonf_,latg_)
     call dyn_trans2model_grid(XLONS,1)
     call dyn_trans2model_grid(XLATS,1)
   endif
#ifdef MP
   call mpgf2p(grid1,lonf2_,latg2_,xlon,lonf2p_,latg2p_,1)
   call mpgf2p(grid2,lonf2_,latg2_,xlat,lonf2p_,latg2p_,1)
#endif
#ifdef DFS
   do j = 1,latg2_
     do i = 1,lonf_
       grid1d(i,        j)=grid1(i      ,j)
       grid1d(i,latg_+1-j)=grid1(lonf_+i,j)
       grid2d(i,        j)=grid2(i      ,j)
       grid2d(i,latg_+1-j)=grid2(lonf_+i,j)
     enddo
   enddo
#ifdef MP
   call mpgf2pd(grid1d,iba,jbwa,xlond,ib,jbw,1)
   call mpgf2pd(grid2d,iba,jbwa,xlatd,ib,jbw,1)
#else
   xlond=grid1d
   xlatd=grid2d
#endif
#endif /* DFS end */
#ifdef DBG
   if(iope) write(6,*)' done xlon xlat '
#endif
#undef XLONS
#undef XLATS
   do j = 1,latg2_
     sinlat(j) = cos(colrad(j))
   enddo
#endif
#ifdef MP
#define SINLABS grid1
#define COSLABS grid2
#else
#define SINLABS sinlab
#define COSLABS coslab
#endif
   if( iope ) then
#ifdef SMP
     call dyn_sincos_lat(SINLABS,COSLABS,xlon,xlat,lonf_,latg_)
#else
     call dyn_sincos_lat(SINLABS,COSLABS,colrad,lonf_,latg_)
#endif
     call dyn_trans2model_grid(SINLABS,1)
     call dyn_trans2model_grid(COSLABS,1)
   endif
#ifdef MP
   call mpgf2p(grid1,lonf2_,latg2_,sinlab,lonf2p_,latg2p_,1)
   call mpgf2p(grid2,lonf2_,latg2_,coslab,lonf2p_,latg2p_,1)
#endif
#ifdef DBG
   if(iope) write(6,*)' done sinlab coslab '
#endif
#undef SINLABS
#undef COSLABS
!
   call dfkt_init
   if(iope) write(6,*)' done dfkt_init'
#ifndef DFS
#ifndef SMP
   call sph_fft_trans
#ifdef DBG
   if(iope) write(6,*)' done sph_fft_trans '
#endif
#endif
   call lsm_init_df
#ifdef DBG
   if(iope) write(6,*)' done lsm_init_df '
#endif
#endif			/* DFS */
   call lsm_init_kt
#ifdef DBG
   if(iope) write(6,*)' done lsm_init_kt '
#endif
#ifndef DCMIP
   call rad_read_init
#endif
#ifdef DBG
   if(iope) write(6,*)' done rad_read_init '
#endif
!
! Define top layer for search of the downdraft originating layer
! and the maximum thetae for updraft
!
   kbmaxin = levs_
   kbmin   = levs_
   kmaxin  = levs_
   kpblmaxin  = levs_
   levshcin  = levs_
   do k = 1,levs_
     if(sl(k).gt.0.45) kbmaxin = k + 1
     if(sl(k).gt.0.70) kbmin   = k + 1
#ifdef NCEP2010
     if(sl(k).gt.0.04) kmaxin  = k + 1
#else
     if(sl(k).gt.0.05) kmaxin  = k + 1
#endif
#ifdef GRIMSCV
     if(sl(k).gt.0.60) levshcin  = k + 1
#else
     if(sl(k).gt.0.70) levshcin  = k + 1
#endif
     if(sl(k).gt.0.60) kpblmaxin  = k + 1
   enddo
   if(iope) write(6,*)' kbmax,kbm,kmax in gmp_start_setup ', kbmaxin,kbmin,kmaxin
   if(iope) write(6,*)' k-max for gwdo reference top in gmp_start_setup ', kpblmaxin
   if(iope) write(6,*)' k-max for shallow clouds top in gmp_start_setup ',levshcin
!
! -- calculate pressure (originated from rad_main_solver)
! -- find cloud tops (originated from rad_lw_nasa_driver)
!   level index separating high and middle clouds (ict)
!   level index separating middle and low clouds   (icb)
!  ict       =15 for 28 layer   !about 440mb for mls standard atmos.
!  icb       =19 for 28 layer   !about 702mb for mls standard atmos.
!
   ictin = 1
   icbin = 1
   do k = 1,levs_
     if(si(k).gt.0.45) ictin = levs_-k+1
     if(si(k).gt.0.70) icbin = levs_-k+1
   enddo
   if(iope) write(6,*)' ictin,icbin in gmp_start_setup ', ictin,icbin
#ifndef DCMIP
   call rad_initialize
#ifdef DBG
   if(iope) write(6,*)' done rad_initialize '
#endif
#endif
!
!  filtwin in the unit of hour
!
#ifdef LFM
   sechr=60.*60.
   klenp=nint(filtwin*sechr/con(1))+1
   nlenp=(klenp-1)/2+1
   critfl=9.0e10
   clancz=1.0
   tinc=con(1)
   critsc=critfs*sechr
#ifndef NOPRINT
   if( iope )                                                                  &
      write(6,*) 'klenp,nlenp,tinc,critfs,critsc,critfl,clancz=',              &
                  klenp,nlenp,tinc,critfs,critsc,critfl,clancz
#endif
   call filtcof(nlenp,tinc,critsc,critfl,clancz,weix)
   fwindow=filtwin
#endif
!
!  expand weights
!
#ifdef LFM
   weight(nlenp)=weix(1)
   do k=2,nlenp
      weight(nlenp+k-1)=weix(k)
      weight(nlenp-k+1)=weix(k)
   enddo
#ifndef NOPRINT
   do k = 1,klenp
     if(iope) write(6,*)k,' weight=',weight(k)
     call flush(6)
   enddo
#endif
#endif
!
   return
   end subroutine gmp_start_setup
!
!-------------------------------------------------------------------------------
#ifdef DFS
   subroutine dyn_sigma_setup(ci, si, del, sl, cl, rpi, tavexy,n1)
!-------------------------------------------------------------------------------
   use dfsvar, only    : mtg,jlg
#else
   subroutine dyn_sigma_setup(ci, si, del, sl, cl, rpi, n1)                          
#endif
   use constant, only  : rd_,cp_
   use paramodel, only : levs_,levp1_,levm1_
#include <abort.h>
!-------------------------------------------------------------------------------
!
! subprogram:    setsig      sets up model sigma structure.                     
!                                                                               
! abstract: sets up model sigma structure based on vertical                     
!   sigma spacing defined in the subroutine.                                    
!                                                                               
! program history log:                                                          
!   88-04-05  joseph sela                                                       
!                                                                               
! usage:    call setsig (ci, si, del, sl, cl, rpi)                              
!                                                                               
!   output argument list:                                                       
!     ci       - array of 1.0-si at each level.                                 
!     si       - array of sigma value at each level.                            
!     del      - array of sigma spacing at each layer.                          
!     sl       - array of sigma at midpoint of sigma layers.                    
!     cl       - array of 1.0-sl at each layer midpoint.                        
!     rpi      - array of pi ratios needed in thermodynamic equation.           
!                                                                               
!   output files:                                                               
!     output   - printout file.                                                 
!                                                                               
!-------------------------------------------------------------------------------
   real rk,rk1,rkr                                                           
#ifdef DFS
   real ci(levp1_), si(levp1_),del(levs_), sl(levs_), cl(levs_), rpi(levm1_),  &
        tavexy(levs_),dum(mtg*jlg)
#else
   real ci(levp1_), si(levp1_),del(levs_), sl(levs_), cl(levs_), rpi(levm1_)
#endif
   integer idate(4)                                                          
!-------------------------------------------------------------------------------
#ifndef DCMIP
#ifndef NOPRINT
   print 98, n1                                                             
   98  format (1x,'begin setsig - getting sigs from unit',i4)                  
#endif
   close(n1)
   open (unit=n1,file='sigit',form='unformatted',err=999)
   go to 998
   999 continue
   print *,'error opening sigit in gmp_start_setup'
   call MPABORT
  998 continue
   rewind n1
   read(n1)                                                                 
   read(n1) fhour,idate,si,sl                                               
#ifdef DFS
   read(n1)                  ! topo
   read(n1)dum(1:mtg*jlg),tavexy(1:levs_)  ! ps & ave_t
#endif
   rewind n1                                                                
#else /* DCMIP */
#ifdef DFS
   tavexy(1:levs_)=300.
#endif
#endif /* ~DCMIP end */
!
   do li = 1,levp1_                                                          
     ci(li) = 1.e0 - si(li)                                                   
   enddo
   do le = 1,levs_                                                           
     cl(le) = 1.e0 - sl(le)                                                   
     del(le) = si(le) - si(le+1)                                               
   enddo
!
! compute pi ratios for temp. matrix.                                       
!
   rk = rd_/cp_                                                              
   do le = 1,levm1_                                                          
     base = sl(le+1)/sl(le)                                                    
     rpi(le) = base**rk                                                        
   enddo
#ifdef DBG
   do le = 1,levp1_                                                          
     print 100, le, ci(le), si(le)                                             
100  format (1h , 'level=', i2, 2x, 'ci=', f6.3, 2x, 'si=', f6.3)              
   enddo
   print 97                                                                  
97 format (1h0)                                                              
   do le = 1,levs_                                                           
      print 101, le, cl(le), sl(le), del(le)                                    
101   format (1h , 'layer=', i2, 2x, 'cl=', f6.3, 2x, 'sl=', f6.3, 2x,         &
             'del=', f6.3)                                                            
   enddo
   print 102, (rpi(le), le=1,levm1_)                                         
102   format (1h0, 'rpi=', (18(1x,f6.3)) )                                      
#endif
   return                                                                    
   end subroutine dyn_sigma_setup
!
!-------------------------------------------------------------------------------
#ifdef DFS
   subroutine dyn_hybrid_setup(ak5,bk5,ci, si, del, sl, cl, rpi, tavexy,n1)
!-------------------------------------------------------------------------------
   use dfsvar, only    : mtg,jlg
#else
   subroutine dyn_hybrid_setup(ak5,bk5,ci, si, del, sl, cl, rpi, n1)
!-------------------------------------------------------------------------------
#endif
   use paramodel, only : levs_,levp1_,levm1_
   use constant,  only : cp_,rd_,akapa_
!-------------------------------------------------------------------------------
!
! subprogram:    dyn_hybrid_setup      sets up model hybrid structure.    
!                                                                               
! abstract: sets up model sigma structure based on vertical                     
!   sigma spacing defined in the subroutine.                                    
!                                                                               
! program history log:                                                          
!   88-04-05  joseph sela                                                       
!                                                                               
! usage:    call hybrid_setup (ci, si, del, sl, cl, rpi)                              
!                                                                               
!   output argument list:                                                       
!     ci       - array of 1.0-si at each level.                                 
!     si       - array of sigma value at each level.                            
!     del      - array of sigma spacing at each layer.                          
!     sl       - array of sigma at midpoint of sigma layers.                    
!     cl       - array of 1.0-sl at each layer midpoint.                        
!     rpi      - array of pi ratios needed in thermodynamic equation.           
!                                                                               
!   output files:                                                               
!     output   - printout file.                                                 
!                                                                               
!-------------------------------------------------------------------------------
   real     ::  rk1,rkr                                                           
   real     ::  ci(levp1_), si(levp1_),                                        &
                del(levs_), sl(levs_), cl(levs_), rpi(levm1_)                            
   real     ::  ak5(levp1_),bk5(levp1_)
   real     ::  ak5x(levp1_),bk5x(levp1_)
#ifdef DFS
   real     ::  tavexy(levs_), dum(mtg*jlg)
#endif
   integer  ::  idate(4)                                                          
!-------------------------------------------------------------------------------
#ifndef DCMIP
#ifndef NOPRINT
   print 98, n1                                                             
   98  format (1x,'begin hybrid_setup - getting sigs from unit',i4)                  
#endif
   close(n1)
   open (unit=n1,file='sigit ',form='unformatted',err=999)
   go to 998
   999 continue
   print *,'error opening sigit in gmp_start_setup'
#ifdef MP
#ifdef RMP
   call rmpabort
#else
   call mpabort
#endif
#else
   call abort
#endif
   998 continue
   rewind n1
   read(n1)                                                                 
   read(n1) fhour,idate,ak5x,bk5x
#ifdef DFS
   read(n1)                             ! topo
   read(n1)dum(1:mtg*jlg),tavexy(1:levs_)             ! ps & ave_t
#endif
!
! though hybrid layers from file (ak5 and bk5) are bottom to top 
! (same as sigma file), the gfs codes assume they are top to bot.
! therefore ak5 and bk5 are reversed here.
!
   rewind n1                                                                
!
! [t2b] ak5,bk5
! [b2t] ak5x,bk5x
!
   do l = 1,levp1_
     ak5(l)=ak5x(levp1_+1-l)/1000. ! Pa -> cb
     bk5(l)=bk5x(levp1_+1-l)
   enddo
#else /* DCMIP */
#ifdef DFS
   tavexy(1:levs_)=300.
#endif
   do l = 1,levp1_
     ak5(l)=ak5(l)/1000.   ! Pa -> cb
   enddo
#endif /* ~DCMIP end */
!
! [b2t] sl
!
   do l = 1,levp1_     
     si(l)=ak5(levp1_+1-l)/100.+bk5(levp1_+1-l)
!
! assume ps=100(cb) for reference. might be a source of error
! in radiative heating calculation
!
   enddo
!
   rk1 = akapa_ + 1.
   rkinv=1./akapa_
!
   do l = 1,levs_
     dif = si(l)**rk1 - si(l+1)**rk1
     dif = dif / (rk1*(si(l)-si(l+1)))
     sl(l) = dif**rkinv
   enddo
!
   do li = 1,levp1_                                                          
     ci(li) = 1.e0 - si(li)                                                   
   enddo
!
   do le = 1,levs_                                                           
     cl(le) = 1.e0 - sl(le)                                                   
     del(le) = si(le) - si(le+1)                                               
   enddo
!
! compute pi ratios for temp. matrix.                                       
!
   do le = 1,levm1_                                                          
     base = sl(le+1)/sl(le)                                                    
     rpi(le) = base**akapa_                                                        
   enddo
#ifdef DBG
   do le = 1,levp1_                                                          
     print 100, le, ci(le), si(le)                                             
     100 format (1h , 'level=', i2, 2x, 'ci=', f6.3, 2x, 'si=', f6.3)              
   enddo
   print 97                                                                  
   97  format (1h0)                                                              
   do le = 1,levs_                                                           
     print 101, le, cl(le), sl(le), del(le)                                    
     101 format (1h , 'layer=', i2, 2x, 'cl=', f6.3, 2x, 'sl=', f6.3, 2x,   &
    'del=', f6.3)                                                            
   enddo
   print 102, (rpi(le), le=1,levm1_)                                         
   102 format (1h0, 'rpi=', (18(1x,f6.3)) )                                      
#endif
!
   return                                                                    
   end subroutine dyn_hybrid_setup                                                                      
!
#ifndef DFS
!-------------------------------------------------------------------------------
   subroutine sph_matrix_init(kmx,si,sl,tov,am,bm,sv,gv,cm)                           
!-------------------------------------------------------------------------------
!
! subroutine:    sph_matrix_init  
!                                                                               
! abstract: 
!   computes the 3 matrices am and bm and sv.                           
!   the matrices represent the linearized gravity wave terms                    
!   of the equations involved in the semi-implicit time integration.            
!   the matrices are dependent only on the vertical structure.                  
!   am is the divergence equations linear dependence on temperature.           
!   bm is the temperature equations linear dependence on divergence.           
!   sv is the continuity equations linear dependence on divergence.            
!                                                                               
! program history log:                                                          
!   1988-04-06  joseph sela                                                       
!   1993-02-23  mark iredell           compact vertical formulation                       
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call sph_matrix_init(km,sl,tov,am,bm,sv)                                     
!   input argument list:                                                        
!     km       - integer number of vertical levels.                             
!     sl       - real (km) sigma level values.                                  
!     tov      - real (km) reference temperatures.                              
!                                                                               
!   output argument list:                                                       
!     am       - real (km,km) such that dd(k)/dt = ... + am(k,j)*t(j)           
!     bm       - real (km,km) such that dt(k)/dt = ... + bm(k,j)*d(j)           
!     sv       - real (km) such that dq/dt = ... + sv(j)*d(j)                   
!                                                                               
!-------------------------------------------------------------------------------
   use constant,  only : cp_,rd_,rerth_
   use paramodel, only : levs_
!-------------------------------------------------------------------------------
   real             ::  si(kmx+1),sl(kmx),tov(kmx) 
   real             ::  am(kmx,kmx),bm(kmx,kmx),sv(kmx),gv(kmx),cm(kmx,kmx)
   real, parameter  ::  rd=rd_,cp=cp_,rerth=rerth_
   real, parameter  ::  rocp=rd/cp,raa=rd/(rerth**2)
!
! local array 
!
   real             ::  cd(levs_,levs_+1),ci(levs_+1,levs_),                   &
                        cq(levs_+1,levs_),cql(levs_,levs_)    
   real             ::  rnu(levs_),rmu(levs_),ti(2:levs_),dt(levs_+1,levs_)
!
   km=kmx                                                                    
   call sph_matrix_init_sub(km,si,sl,cd,ci,cq,cql)                                          
   do j = 1,km                                                              
     sv(j)=cq(km+1,j)                                                        
     gv(j)=raa*tov(j)                                                        
     dt(1,j)=0.                                                              
     dt(km+1,j)=0.                                                           
   enddo
!
   do k = 1,km-1                                                            
     ti(k+1)=0.                                                              
     del=si(k+1)-si(k)                                                       
     rsl=log(si(k+1)/si(k))/del                                              
     rnu(k)=(1.-rsl*si(k))/del                                               
     rmu(k)=(rsl*si(k+1)-1.)/del                                             
   enddo
!
   rnu(km)=0.                                                                
   rmu(km)=1./si(km)                                                         
   do k = 2,km                                                              
     do j = 1,km                                                            
       ti(k)=ti(k)+ci(k,j)*tov(j)                                            
       dt(k,j)=(1-si(k))*cq(km+1,j)-cq(k,j)                                  
     enddo
   enddo
!
   do j = 1,km                                                              
     do k = 1,km                                                            
       am(k,j)=-raa*cql(k,j)                                                 
       bm(k,j)=(rocp-1)*tov(k)*cq(km+1,j)                                      &
                +rocp*tov(k)*(rnu(k)*dt(k+1,j)+rmu(k)*dt(k,j))                 
     enddo
     bm(j,j)=bm(j,j)-tov(j)                                                  
     do k = 1,km                                                            
       do i=2,km                                                            
         bm(k,j)=bm(k,j)-cd(k,i)*ti(i)*dt(i,j)                                 
       enddo
     enddo
   enddo
#ifdef ORIGIN_THREAD
!$doacross share(km,gv,sv,am,bm,cm),local(j,k,i)                                
#endif
#ifdef OPENMP
!$omp parallel do private(j,k,i)
#endif
   do j = 1,km                                                             
     do k = 1,km                                                            
         cm(k,j)=gv(k)*sv(j)                                                   
     enddo
     do k = 1,km                                                            
       do i = 1,km                                                            
         cm(k,j)=cm(k,j)+am(k,i)*bm(i,j)                                       
       enddo
     enddo
   enddo
!
   return                                                                    
   end subroutine sph_matrix_init
!
!-------------------------------------------------------------------------------
   subroutine sph_matrix_init_sub(kmx,si,sl,cd,ci,cq,cql)                                   
!-------------------------------------------------------------------------------
   use paramodel, only  :   levs_
   use constant, only   :   rd_,cp_
!-------------------------------------------------------------------------------
   implicit none
!
   real, parameter      ::  rd=rd_,cp=cp_,rocp=rd/cp
   integer, intent(in)  ::  kmx
   real, intent(in)     ::  si(kmx+1),sl(kmx)                                               
   real, intent(out)    ::  cd(kmx,kmx+1),ci(kmx+1,kmx)
   real, intent(out)    ::  cq(kmx+1,kmx),cql(kmx,kmx)
!
! local array                                                                   
!
   real                 ::  del(levs_),slk(levs_),alfa(2:levs_),beta(levs_-1)
   integer              ::  km,k,kd,ki,kt,kz
!-------------------------------------------------------------------------------
   km=kmx                                                                    
   do k = 1,km                                                                 
     del(k)=si(k+1)-si(k)                                                    
     slk(k)=sl(k)**rocp                                                      
   enddo                                                                     
!
   do k = 2,km                                                                 
     alfa(k)=0.5*(1.-slk(k-1)/slk(k))/rocp                                   
   enddo                                                                     
!
   do k = 1,km-1                                                               
     beta(k)=0.5*(slk(k+1)/slk(k)-1.)/rocp                                   
   enddo                                                                     
!                                                                             
   do kd = 1,km                                                                
     do ki=1,km+1                                                            
       cd(kd,ki)=0                                                           
     enddo                                                                   
   enddo                                                                     
!
   do kd = 1,km                                                                
     cd(kd,kd)=-1/del(kd)                                                    
     cd(kd,kd+1)=1/del(kd)                                                   
   enddo                                                                     
!                                                                             
   do ki = 1,km+1                                                              
     do kd=1,km                                                              
       ci(ki,kd)=0                                                           
     enddo                                                                   
   enddo                                                                     
!
   do ki = 2,km                                                                
     ci(ki,ki-1)=0.5                                                         
     ci(ki,ki)=0.5                                                           
   enddo                                                                     
!                                                                             
   do ki = 1,km+1                                                              
     do kd = ki,km                                                             
       cq(ki,kd)=0                                                           
     enddo                                                                   
     do kd = 1,ki-1                                                            
       cq(ki,kd)=del(kd)                                                     
     enddo                                                                   
   enddo                                                                     
!                                                                             
   cql(1,1)=del(1)-beta(1)*si(2)                                             
   do kt = 2,km-1                                                              
      cql(1,kt)=del(kt)-alfa(kt)*si(kt)-beta(kt)*si(kt+1)                     
   enddo                                                                     
!
   cql(1,km)=del(km)-alfa(km)*si(km)                                         
   do kz = 2,km                                                                
     cql(kz,1)=cql(1,1)+beta(1)                                              
     do kt=2,kz-1                                                            
       cql(kz,kt)=cql(1,kt)+alfa(kt)+beta(kt)                                
     enddo                                                                   
     cql(kz,kz)=cql(1,kz)+alfa(kz)                                           
     do kt=kz+1,km                                                           
       cql(kz,kt)=cql(1,kt)                                                  
     enddo                                                                   
   enddo                                                                     
!                                                                             
   return                                                                    
   end subroutine sph_matrix_init_sub
!
!-------------------------------------------------------------------------------
   subroutine sph_matrix_init_hybrid(kmx,ak5,bk5,am,bm,sv,gv,cm)
!-------------------------------------------------------------------------------
   use paramodel, only : levs_
   use constant, only  : cp_,rd_,rerth_,akapa_
#ifdef DBG
   use comio   ,only   : iope
#endif
!-------------------------------------------------------------------------------
   real     ::  pk5ref(levs_+1),beta,dpkref(levs_),                            &
                tref(levs_),psref,factor,                                      &
                alfaref(levs_),                                                &
                vecm(levs_), yecm(levs_,levs_),tecm(levs_,levs_)
   real     ::  am(levs_,levs_),bm(levs_,levs_),sv(levs_),                     &
                gv(levs_),cm(levs_,levs_)
   real     ::  ak5(levs_+1),bk5(levs_+1)
   integer  ::  k,kk,kmx,j,irow,icol,icolbeg,icolend
!-------------------------------------------------------------------------------
   do k = 1,levs_
     tref(k)=300.
   enddo
   psref=80.
   beta=1.
! 
   do k = 1,levs_+1
     pk5ref(k)=ak5(k)+bk5(k)*psref
#ifdef DBG
     if(iope) print 100,k,ak5(k),bk5(k),pk5ref(k)
#endif
   enddo
100   format('k=',i2,2x,'ak5=',e10.3,2x,'bk5=',e10.3,2x,'pk5ref=',e10.3)
! 
   do k = 1,levs_
     dpkref(k)=pk5ref(k+1)-pk5ref(k)
     gv(k)=beta*rd_*tref(k)/(rerth_*rerth_)
#ifdef DBG
     if(iope) print 110,k,dpkref(k),gv(k)
#endif
   enddo
110   format('k=',i2,2x,' in am_bm dpkref=',e11.4,2x,'gv=',e11.4)
#ifdef DBG
   if(iope) print*,'-  calculate alfaref  watch alfaref(1)     '
#endif
  
   alfaref(1)=log(2.) ! could also be=1.  but watch for layer values
   do k = 2,levs_
     alfaref(k)=1.-(pk5ref(k)/dpkref(k))*log(pk5ref(k+1)/pk5ref(k))
#ifdef DBG
     if(iope) print 210,k,k,k,k,k
210  format('alfa(',i2,')=1.-(pk5(',i2,')/dpk(',i2,'))*log(pk5(',i2,          &
       '+1)/pk5(',i2,'))')
#endif
   enddo
! 
#ifdef DBG
   if(iope) then
     print 125,alfaref(1)
125   format('worry --- alfaref(1)=',e10.3)
     do k = 1,levs_
       print 130,k,alfaref(k)
     enddo
130   format('k=',i2,2x,'alfaref',e16.8)
   endif
#endif
!
!  print*,'---- begin matrices computation -----------'
! 
!  print 144
144   format(1x,'begin yecm computation')
   yecm=0.
   do irow = 1,levs_
     yecm(irow,irow)=alfaref(irow)*rd_
     icolbeg=irow+1
     if(icolbeg.le.levs_)then
       do icol=icolbeg,levs_
         yecm(irow,icol)=rd_*log( pk5ref(icol+1)/pk5ref(icol) )
       enddo
     endif
   enddo
150    format('yecm(',i2,',',i2,')=rd_*log( pk5ref(',i2,                       &
              '+1)/pk5ref(',i2,'))')
!  print*,'-----------------1234567------------------'
160    format('yecm=',4(1x,e10.3))
! 
   tecm=0.
! 
   do irow = 1,levs_
!    print*,' doing row ...............................',irow
     tecm(irow,irow)=akapa_*tref(irow)*alfaref(irow)
     icolend=irow-1
 
     do icol = 1,icolend
       factor=(akapa_*tref(irow)/                                              &
                         dpkref(irow))*log(pk5ref(irow+1)/pk5ref(irow))
       tecm(irow,icol)=factor*dpkref(icol)
     enddo
   enddo
165    format('irow=',i2,2x,'factor=',e16.8,2x,'icolend=',i2)
166    format('factor=(akapa_*tref/dpkref(',i2,'))*log(pk5ref(',i2,            &
       '+1)/pk5ref(',i2,'))')
167    format('innerlup irow=',i2,2x,'icol=',i2,2x,'tecm(ir,ic)=',e12.4)
!
!  print*,'4444444  print yecm      44444444444444444'
! 
!  do irow=1,levs_
!       print*,'yecm row=',irow,'levs_=',levs_
!       print 1700,(yecm(irow,j),j=1,levs_/2)
!       print 1701,(yecm(irow,j),j=levs_/2+1,levs_)
!  enddo
!1700   format('  a  ',10(1x,e10.3))
!1701   format('  b  ',10(1x,e10.3))
! 
!  print*,'5555555  print tecm      55555555555555555'
! 
!  do irow=1,levs_
!    print*,'tecm row=',irow,'levs_=',levs_
!    print 1700,(tecm(irow,j),j=1,levs_/2)
!    print 1701,(tecm(irow,j),j=levs_/2+1,levs_)
!  enddo
! 
!  print*,'666666666666666666666666666666666666666666'
!  print 171
!171   format(1x,'begin vvec dcomputation')
!
   do icol = 1,levs_
     vecm(icol)=dpkref(icol)/psref
!    print 175,icol,vecm(icol)
   enddo
!175    format('icol=',i2,2x,'vecm=',e16.8)
 
   do j = 1,levs_
     sv(j)=vecm(levs_+1-j)
     do k = 1,levs_
       am(k,j)=yecm(levs_+1-k,levs_+1-j)
       bm(k,j)=tecm(levs_+1-k,levs_+1-j)
     enddo
   enddo
! 
   do j = 1,levs_
     do k = 1,levs_
       am(k,j)=am(k,j)*beta/(rerth_*rerth_)
     enddo
   enddo
!
   do j = 1,levs_
     do k = 1,levs_
          cm(k,j)=gv(k)*sv(j)
     enddo
     do k = 1,levs_
       do i = 1,levs_
            cm(k,j)=cm(k,j)+am(k,i)*bm(i,j)
       enddo
     enddo
   enddo
!
   return
   end subroutine sph_matrix_init_hybrid
!
!-------------------------------------------------------------------------------
   subroutine solve_sph(qlnt,qlnv,qdert,epsi)
!-------------------------------------------------------------------------------
!                                            
!   part between guards made into sr ggozri.     
!   7 dec 1990      m. rozwodoski               
!                                                    
!   compute pln derivatives in ibm order.           
!                                                  
!-------------------------------------------------------------------------------
   use paramodel, only : jcap1_,jcap2_,lnt2_,lnut2_,twoj1_
   use comfcst, only   : dxa, dxb
!-------------------------------------------------------------------------------
   real                 ::  qlnt(lnt2_)   
   real                 ::  qlnv(lnut2_) 
   real                 ::  qdert(lnt2_)
   real                 ::  epsi(jcap2_,jcap1_)      
!-------------------------------------------------------------------------------
   lp0 = 0                                                                   
   lp1 = 2                                                                   
   len = twoj1_                                                              
   do i = 1,jcap1_                                                        
     do ll = 1,len                                                           
       qdert(ll+lp0) = qlnv(ll+lp1) * dxb(ll+lp0)                         
     enddo
     lp1 = lp1 + len + 2                                                       
     lp0 = lp0 + len                                                           
     len = len - 2                                                             
   enddo
   !
   lend = lnt2_ - 4                                              
   do ll=1,lend                                                          
     qdert(ll+2) = qdert(ll+2) + qlnt(ll) * dxa(ll+2)                   
   enddo
   !                                                   
   return                                                                    
   end subroutine solve_sph
!
#endif /* ~DFS end */
!-------------------------------------------------------------------------------
#ifdef SMP
   subroutine dyn_sincos_lat(sinlat,coslat,xlon,xlat,lon,lat)                             
#else
   subroutine dyn_sincos_lat(sinlat,coslat,colrad,lon,lat)                             
#endif
!-------------------------------------------------------------------------------
   integer, intent(in)  ::  lon,lat
   real, intent(out)    ::  sinlat(lon,lat),coslat(lon,lat)
#ifdef SMP
   real, intent(in)     ::  xlat(lon,lat),xlon(lon,lat)
!-------------------------------------------------------------------------------
   do j = 1, lat
     do i = 1, lon
       sinlat(i,j) = sin(xlat(i,j))
       coslat(i,j) = sqrt(1.e0 - sinlat(i,j)*sinlat(i,j))
     enddo
   enddo
#else
   real,intent(in)     ::   colrad(lat/2)
!-------------------------------------------------------------------------------
!
!  get normal sinlat and coslat 
!
   do j = 1,lat/2
     sinlaj = cos(colrad(j))
     coslaj = sqrt(1.e0 - sinlaj*sinlaj)
     do i = 1,lon
       sinlat(i,j) = sinlaj
       coslat(i,j) = coslaj
     enddo
   enddo
   do j = lat/2+1,lat
     jj=lat+1-j
     do i = 1,lon
       sinlat(i,j) = -sinlat(i,jj)
       coslat(i,j) =  coslat(i,jj)
     enddo
   enddo
#endif
!
   return                                                                    
   end subroutine dyn_sincos_lat 
!
