#include "define.h"
   subroutine gmp_start
!-------------------------------------------------------------------------------
!
! subroutine: gmp_start         
!
! abstract: make global initilization
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
#ifdef MP
   use paramodel, only : jcap_,lonf_,latg_,lonf2_,latg2_,lnt2_,                &
#else
   use paramodel, only : jcap1_,lonf_,latg_,lonf2_,latg2_,lnt2_,               &
#endif
                         LONF2S,LATG2S,levh_,levs_
   use comsfc        ! sfcfcs,sfcftyp,albedo,tsea,slmsk
#ifdef COS_BELL
   use cosbell, only   : init_cosbell
#endif
#ifdef DFS
   use dfsvar, only    : psl1,t1,d1,v1,q1,sf1,                                 &
#ifdef MP
                                          sf1l,                                &
#endif
#ifdef PERT
                         d0,                                                   &
#endif
                         psl2,t2,d2,v2,q2,                                     &
                         ntotal,iope,lnt2,nls,nle,                             &
                         sl=>sigma,si=>sigmafull,del=>delsig,                  &
                         coslat,gamma_v,gamma_d,gamma_t,gamma_q,               &
                         mt,levsp,levhp,jbw,ib,jl,jlg,                         &
                         levs,levh,                                            &
#if defined(BARO_TEST) || defined(COS_BELL)
                         mta,vor,dio,tai,qai,prs,sfai,tavexy,                  &
#endif
#if defined(BIN_DBG)
#ifndef BARO_TEST
                         mta,vor,dio,tai,qai,prs,tavexy,                       &
#endif
                         mta,xpsl,ypsl,vxc,vyc,                                &
                         AMATm,DMATm,AMATm_,DMATm_,aMSQUAR,jcol2js
#endif
                         mta,jb,jgs,jba
#else
   use module_sph_semi_implicit, only : sph_semi_gwave
#endif /* DFS end */
#ifdef MP
   use commpi       ! mype,mycol,myrow,master,lntlen,lonlen,latlen,lntstr
#endif
   use comfibm
   use comcon
   use comgpd
   use comfgsm
#ifdef EXPLICIT_CLOUDINESS
#ifdef CPS_QINI
   use comfcst, only                 : qcicps,qrscps,taucld,cldwp,cldip
#endif
#endif
#ifdef LFM
   use comlfm
#endif
#ifdef DG3
   use diag_3d_module, only          : diag_3d_index, diag_3d_zero_out
#endif
   use module_file_write, only       : file_write_bin
#ifdef NISLQ
   use paramodel, only : ntotal_
   use nislq    , only : slq_q1,slq_q2
#endif
#ifdef HORA
   use comdyn, only : hora_n2,workt,workd,workz,workr
#ifndef DFS
   use paramodel, only : LNT22S
#endif
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#include <abort.h>
!
#ifndef DFS
#ifdef COS_BELL
   real,dimension(LONF2S,LATG2S,levs_)  ::  vor,tai,qai
   real,dimension(LONF2S,LATG2S)        ::  sfai,prs
   common /comcosbell/vor,tai,qai,prs,sfai
#endif
#endif /* ~DFS end */
#if defined(EXPLICIT_CLOUDINESS) && defined(CPS_QINI)
   integer                              :: ioerr
#endif
!
   character(len=128)                   ::  fno,fname
   integer                              ::  itread,nvar,lotg,k,j,m
   integer                              ::  nfstep,iday,l,i
!-------------------------------------------------------------------------------
#ifndef OMP
   itread=0
#else
   itread=omp_get_num_threads()
#endif
#ifdef MPI_DBG
#ifdef MP
   write(fname,'(A,3(I3.3,A))')'fcst_std-p',mype,'-r',myrow,'-c',mycol,'.out'
#else
   write(fname,'(A)')'fcst_std.out'
#endif
   close(96)
   open(96,file=trim(fname),status='unknown',form='formatted')
#endif
!
#ifdef DFS
   lnts2=lnt2
   lons2=ib*2
   lats2=jb
#else
#ifdef MP
   iope=mype.eq.master
   call mpdimset(jcap_,levs_,lonf_,latg_)
   lnts2=lntlen(mype)*2
   lons2=lonlen(mype)*2
   lats2=latlen(mype)
   lnoffset=lntstr(mype)*2
#else
   iope=.true.
!
   do l = 1,jcap1_
     lwvdef(l)=l-1
   enddo
!
   do l = 1,latg2_
     latdef(l)=l
   enddo
!
   lnts2=lnt2_
   lons2=lonf2_
   lats2=latg2_
   lnoffset=0
#endif
#endif /* MP end */
!
   call gmp_start_setup(n1,                                                    &
#ifdef RMP
    nrsmi1,nrsmi2,nrflip,                                                      &
    nrsmo1,nrsmo2,nrflop,nrsfli,nrsflx,nrinit,nrpken,                          &
#endif
#ifdef LFM
    nlfmsgi,nlfmsfi,nlfmsgo,nlfmsfo,klenp,weight,filtwin,                      &
#endif
    kpfix,ksfcx,komlx,ksig,ksfc,kpost,krestart,klfm,krsm)
!
#ifdef LFM
#ifndef NOPRINT
   if(iope)                                                                    &
     write(6,*) 'nlfmsgi,nlfmsfi,nlfmsgo,nlfmsfo,klenp,wght,fltwin=',          &
            nlfmsgi,nlfmsfi,nlfmsgo,nlfmsfo,klenp,(weight(i),i=1,10),          &
            filtwin
#endif
#endif /* LFM end */
!
!  fhour is forecast hour of the start of the segment
!  thour is forecast hour rounded to nearest full hour
!  shour is forecast time in seconds measured from the start of the segment
!
!  shour is initialized to zero in step1
!
!  thour=fhour+int(shour/3600.+0.5e0)
!
   thour=fhour
#ifdef DFI
!
!  dyn_digital_filter:  do digital filter initialization setup
!         nummax and numsum saved and passed in comver
!
   if( con(3).ne.0.0 ) then
     nummax=nint(con(3)*3600./con(1)/2.)
     numsum=-nummax-1
#ifndef NOPRINT
     if(iope) write(6,*)' do gmp digital filter initialization '
#endif
     print*, 'call dyn_digital_filter'
     call dyn_digital_filter(0,con(3),shour/3600.,solsec,lugi)
   else
     nummax=0
     numsum=-1
   endif
#endif
!
!  zero-out arrays
!
   call dyn_flux_zero_out(0)
   call dyn_tmax_zero_out
!
#ifdef DG3
   call diag_3d_index
   call diag_3d_zero_out
#endif
!
!  read the t-dt or t=0 sigma file
!
#ifndef DCMIP
   if(iope) then
     open (unit=n1,file='sigit',form='unformatted',err=999)
     go to 998
999  continue
     write(6,*) 'error opening sigit in rmp_start'
     call MPABORT
998  continue
   endif
#endif
!
#ifndef DFS
#ifdef COS_BELL
!
!  u,v,z,ps,zs
!
   call init_cosbell(VOR,TAI,QAI,PRS,SFAI,rm,ib,jbw,colrad,0.0)
   call file_write_bin(333,qai,ib,jbw,1    ,0)
#else
#ifdef DCMIP
   fhour=0.
   idate(4)=2000    ! year
   idate(2)=1       ! month
   idate(3)=1       ! day
   idate(1)=0       ! hour
#ifdef NISLQ
   call dcmip_update(fhour,gz,qm,dim,zem,tem,rm,z00,slq_q1)
#else
   call dcmip_update(fhour,gz,qm,dim,zem,tem,rm,z00)
#endif
   call dcmip_write
#else /* ~DCMIP */
#ifndef HYBRID
   call file_read_sigma(n1,fhour,idate,gz,qm,tem,dim,zem,rm,sl,si,z00)
#else
   call file_read_sigma(n1,fhour,idate,gz,qm,tem,dim,zem,rm,ak5,bk5,z00)
#endif
#endif /* DCMIP end */
#endif /* COS_BELL end */
#else /* DFS */
#ifdef BARO_TEST
   call dfs_baro_test_init(VOR, DIO, TAI, PRS, SFAI)
!
   sfai(:,:)=0.0
   fhour=0.0
   idate(1)=0
   idate(2)=1
   idate(3)=1
   idate(4)=2000
!
   nvar=3
   lotg=nvar*levs_
   call dfs_fft_driver(1,dio,ib,jbw,lotg,D1,mt,jl,lotg,jlg,                    &
                         levsp,levs_,nvar,coslat,1)
   call dfs_fft_driver(1,prs,ib,jbw,1,psl1,mt,jl,1,jlg,                        &
                         1,1,0,coslat,1)
   call file_write_bin(333,PRS,ib,jbw,1    ,0)
   call file_write_bin(333,dio,ib,jbw,3*levs_,0)
#else /* ~BARO_TEST */
#ifdef COS_BELL
!
!  u,v,z
!
   call init_cosbell(VOR,TAI,QAI,PRS,SFAI,Q1,ib,jbw,colrad,0.0)
   call file_write_bin(333,qai,ib,jbw,1    ,0)
#else
#ifdef DCMIP
   fhour=0.
   idate(4)=2000    ! year
   idate(2)=1       ! month
   idate(3)=1       ! day
   idate(1)=0       ! hour
#ifdef NISLQ
   call dcmip_update(fhour,sf1,psl1,d1,v1,t1,q1,z00,slq_q1)
#else
   call dcmip_update(fhour,sf1,psl1,d1,v1,t1,q1,z00)
#endif /* NISLQ end */
   call dcmip_write
#else /* ~DCMIP */
#ifndef HYBRID
   call file_read_sigma(n1,fhour,idate,sf1,psl1,t1,d1,v1,q1,sl,si,z00)
#else
   call file_read_sigma(n1,fhour,idate,sf1,psl1,t1,d1,v1,q1,ak5,bk5,z00)
#endif
#endif /* DCMIP end */
#endif /* COS_BELL end */
#endif /* BARO_TEST end */
#endif /* ~DFS end */
#ifdef HORA
!
! U(-1)=U(0) for high-order time filter
!
#ifdef DFS
   call hora_n2(t1 ,d1 ,v1 ,q1,workt,workd,workz,workr,mt*jlg,levsp,levhp)
#else
   call hora_n2(tem,dim,zem,rm,workt,workd,workz,workr,LNT22S,levs_,levh_)
#endif
#endif /* HORA end */
#ifdef NISLQ
#ifndef DCMIP
!
! nislq moisture at n-1 time (wave to grid)
!
#ifdef DFS
   call dfs_fft_driver(-1,slq_q1,ib,jbw,levh,q1,mt,jlg,levhp,jlg,              &
                          levsp,levs,ntotal_,coslat,1)
#else
   call sph_fft_driver(-1,slq_q1,rm,ntotal_)
#endif
!
!   remove negative values (qmin=0.)
!   do j=1,LATG2S
!     do k=1,levh_
!     do i=1,LONF2S
   slq_q1=max(slq_q1,0.)
!     enddo
!     enddo
!   enddo
!
#endif /* ~DCMIP end */
#endif /* NISLQ end */
!
   if( iope ) then
     rewind n1
#ifndef NOPRINT
     write(6,9877)n1,itread,fhour
9877 format(1h ,'n1,itread,fhour after tread',2(i4,1x),f10.1)
     write(6,*)' input t=t0 full values'
#endif
   endif
#ifndef NOPRINT
!
#ifdef DBG
#ifdef DFS
   call dyn_comp_rms(psl1,d1,t1,v1,del,q1)
#else
   call dyn_comp_rms(qm,dim,tem,zem,del,rm)
#endif /* DFS end */
#endif /* DBG end */
#endif /* ~NOPRINT end */
!
   if(fhour.eq.0.) then
#ifdef DFS
     do j = 1,jlg
       do m = 1,mt
         psl2(m,j)=psl1(m,j)
       enddo
     enddo
!
     do k = 1,levsp
       do j = 1,jlg
         do m = 1,mt
           t2(m,j,k)=t1(m,j,k)
           d2(m,j,k)=d1(m,j,k)
           v2(m,j,k)=v1(m,j,k)
         enddo
       enddo
     enddo
#ifndef NISLQ
     do k = 1,levhp
       do j = 1,jlg
         do m = 1,mt
           q2(m,j,k)=q1(m,j,k)
         enddo
       enddo
     enddo
#endif
#else
     do i = 1,lnts2
       q(i)=qm(i)
     enddo
!
     do k = 1,levs_
       do i = 1,lnts2
         te(i,k)=tem(i,k)
         di(i,k)=dim(i,k)
         ze(i,k)=zem(i,k)
       enddo
     enddo
#ifndef NISLQ
     do k = 1,levh_
       do i = 1,lnts2
         rq(i,k)=rm(i,k)
       enddo
     enddo
#endif
#endif
!
#ifdef NISLQ
!
! nislq moisture at n time (n=n-1)
!
     slq_q2(:,:,:)=slq_q1(:,:,:)
#endif
!
#if defined(BARO_TEST) || defined(COS_BELL)
     deltim=con(1)
     limlow=1
     inistp=0
     isave=1
     stepone=.false.
#else
     nfstep=2
     deltim=con(1)/2.e0**nfstep
     limlow=2-nfstep
     inistp=1
#ifdef DGP
     isave=1
#endif
     stepone=.true.
#endif
   else
     if(iope) then
       close(n1)
       open (unit=n1,file='sigitdt',form='unformatted',err=889)
       go to 888
889    continue
       write(6,*) 'error opening sigitdt in gmp_start'
       call MPABORT
888    continue
     endif
#ifdef DFS
#ifndef HYBRID
     call file_read_sigma(n1,fhour,idate,sf1,psl2,t2,d2,v2,q2,sl,si,z00)
#else
     call file_read_sigma(n1,fhour,idate,sf1,psl2,t2,d2,v2,q2,ak5,bk5,z00)
#endif
#else /* SPH */
#ifndef HYBRID
     call file_read_sigma(n1,fhour,idate,gz,q,te,di,ze,rq,sl,si,z00)
#else
     call file_read_sigma(n1,fhour,idate,gz,q,te,di,ze,rq,ak5,bk5,z00)
#endif
#endif /* DFS end */
!
#ifdef NISLQ
!
! nislq moisture at n time (wave to grid)
!
#ifdef DFS
     call dfs_fft_driver(-1,slq_q2,ib,jbw,levh,q2,mt,jlg,levhp,jlg,            &
                          levsp,levs,ntotal_,coslat,1)
#else
     call sph_fft_driver(-1,slq_q2,rq,ntotal_)
#endif
!
! remove negative values (qmin=0.)
!      do j=1,LATG2S
!        do k=1,levh_
!        do i=1,LONF2S
     slq_q2=max(slq_q2,0.)
!        enddo
!        enddo
!      enddo
!
#endif
!
     if( iope ) then
       rewind n1
#ifndef NOPRINT
       write(6,9878)n1,itread,fhour
9878   format(1h ,'n1,itread,fhour after tread',2(i4,1x),f10.1)
       write(6,*) ' input t=t0+dt full values'
#endif
     endif
!
     deltim=con(1)
     limlow=1
     inistp=0
#ifdef DGP
     isave=0
#endif
     stepone=.false.
   endif
!
#ifdef DBG
#ifdef DFS
   call dyn_comp_rms(psl2,d2,t2,v2,del,q2)
#else /* SPH */
   if(iope) print*, 'q in gmp_start', (q(i),i=1,20)
   call dyn_comp_rms(q,di,te,ze,del,rq)
#endif
#endif
!
   shour=0.0
   dtpost=0.0
   maxstp=num(7)
   dthr = con(1)/3600.e0
   hdthr = 0.5 * dthr
#ifndef NO_PHYSICS
!
!  set initial solhr
!
   solhr=fhour+idate(1)
   iday=solhr/24.e0
   solhr=solhr-iday*24.e0
   solsec=solhr*3600.
!
   if( iope ) then
     write(6,*) ' initial solhr = ',solhr
   endif
!
!   read fixed fields from fixfld prog
!
   fno='sfci'
   call sfc_read_file_driver(n1,fno,sfcftyp,                                   &
              labs,idate(4),idate(2),idate(3),idate(1),fhour,                  &
              sfcfcs,LONF2S,LATG2S,0)
#ifdef DBG
   if ( iope ) then
     write(6,*) 'just after sfc_read_file_driver'
     call print_maxmin_six(tsea,lonf2_,latg2_,1,1,'tsea')
     call print_maxmin_six(albedo,lonf2_,latg2_,1,1,'albedo')
     call print_maxmin_six(slmsk,lonf2_,latg2_,1,1,'slmsk')
     call print_maxmin_six(hml0,lonf2_,latg2_,1,1,'hml0')
   endif
#endif
#ifdef EXPLICIT_CLOUDINESS
#ifdef CPS_QINI
!
!   read initial radiation properties
!
   qcicps = 0.
   qricps = 0.
   taucld = 0.
   cldwp = 0.
   cldip = 0.
!
   if(iope) then
     close(44)
     open (unit=44,file='cpscldi',form='unformatted',status='old',iostat=ioerr)
   endif
   call mpbcasti(ioerr,1)
!
   if(ioerr.eq.0) call file_read_cpscldi(qcicps,qrscps,taucld,cldwp,cldip)
!
#ifdef DBG
   if ( iope ) then
     call print_maxmin_six(qcicps,lonf2_,latg2_,1,1,'qcicps')
   endif
#endif
!
#endif
#endif
!
#endif /* ~NO_PHYSICS end */
!
#ifdef DFI
   limlow=limlow-nummax
#endif
#ifdef LFM
   call dyn_lfm_initialize(ipstep,fhour)
   icstep=nint(fhour*3600./deltim)+1
#ifndef NOPRINT
   if(iope) then
     write(6,*) 'ipstep from dyn_lfm_initialize =',ipstep
     write(6,*) 'icstep computed from fhour=',icstep
   endif
#endif
#endif
!
   if( .not. stepone ) then
     avprs0=0.0
#ifdef DFS
     call dfs_mass_adjustment(psl2,avprs0)
#else
     call sph_mass_adjustment(q,avprs0)
#endif
   endif
#if defined(DFS) && defined(MP)
!
! update local for sf1
!
   call mpm2mn(SF1 ,jlg,SF1l ,jl,mt,1)
#endif
#ifndef SMP
#ifdef DFS
#ifdef PERT
!
! for perturbation diffusion
!
   do k = 1,(3+ntotal)*levs_
     d0(1:mt, 1:jl, k) = d1(1:mt, 1:jl, k)
   end do
#endif
   call dfs_diffusion_coef(deltim,gamma_v,mta,1)
   call dfs_diffusion_coef(deltim,gamma_d,mta,2)
   call dfs_diffusion_coef(deltim,gamma_t,mta,3)
   call dfs_diffusion_coef(deltim,gamma_q,mta,4)
#else
   call sph_semi_gwave(deltim,am,bm,gv,sv,cm)
#endif /* DFS end */
#ifdef DBG
   if(iope) write(6,*) ' gmp_start:  sph_semi_gwave for deltim= ',deltim
#endif
#endif /* ~SMP end */
!
#ifdef DGP
   if(npoint.gt.0) then
     isave = 1
     itnum = 1
     if(fhour.eq.0.) then
       if(istep.eq.1) then
         if (ikfreq.gt.1) isave = 0
         if (ikfreq.eq.1) itnum = 2
       end if
     end if
   end if
#endif
!
#ifdef CHEM
   call start_chem
#endif
!
   return
   end subroutine gmp_start
