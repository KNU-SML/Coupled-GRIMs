#include "define.h"
   subroutine mtn_solver(im,jm,nm,nr,nw)
!-------------------------------------------------------------------------------
#ifdef SMP_RA2SFC
   use paramodel, only : jcap_=>jcapscm_,lonf_=>lonfscm_,latg_=>latgscm_      ,&
#else
   use paramodel, only : jcap_          ,lonf_          ,latg_                ,&
#endif
                         imn=>imn_      ,jmn=>jmn_      ,mtnres_              ,&
                         mtnvar=>mtnvar_                                      ,&
                         ngases_,nwater_,ntotal_,levs_
   use constant, only  : pi_
#ifdef DFS
   use dfsvar,   only  : mt,jl,jlg,get_dfs_dim,ib,jbw,jbwa,jge,jle,            &
                         sinlat,coslat
#endif
   use module_file_write, only : file_write_bin
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  im,jm,nm,nr,nw
!
   integer              ::  zavg(imn,jmn),zvar(imn,jmn)
   integer              ::  zmax(imn,jmn),zslm(imn,jmn)
!
   real                 ::  cosclt(latg_),wgtclt(latg_),rclt(latg_),xlat(latg_)
   real                 ::  slm(lonf_,latg_),oro(lonf_,latg_)!,var(lonf_,latg_)
#ifdef DFS
   real                 ::  xlatd(latg_)
   real, allocatable    ::  orsdfs(:,:),WNM2D1(:,:)
#endif
   real                 ::  ors((jcap_+1)*(jcap_+2))
   real                 ::  orsmth((jcap_+1)*(jcap_+2))
   real                 ::  orf(lonf_,latg_)
!
   real,pointer         ::  var(:,:)
   real,pointer         ::  var4(:,:),oa(:,:,:)
   real,pointer         ::  ol(:,:,:)
#ifdef FLOW_BLOCKING
   real,pointer         ::  omax(:,:)
#endif
!  real                 ::  var4(lonf_,latg_),oa(lonf_,latg_,4)
!  real                 ::  ol(lonf_,latg_,4)
   real                 ::  work(lonf_,latg_),work1(lonf_,latg_),              &
                            work2(lonf_,latg_),work3(lonf_,latg_),             &
                            work4(lonf_,latg_),                                &
                            work5(lonf_,latg_),work6(lonf_,latg_),             &
                            glat(jmn)
#ifdef FLOW_BLOCKING
   real,target          ::  hprime(lonf_,latg_,11)
#else
   real,target          ::  hprime(lonf_,latg_,10)
#endif
!  real,allocatable     ::  hprime(:,:,:)
   integer              ::  ist(lonf_,latg_),ien(lonf_,latg_)
   integer              ::  jst(latg_),jen(latg_)
   integer              ::  iwork(lonf_,latg_,4)
!
!   equivalence (var (1,1),hprime(1,1,1))
!   equivalence (var4(1,1),hprime(1,1,2))
!   equivalence (oa(1,1,1),hprime(1,1,3))
!   equivalence (oa(1,1,2),hprime(1,1,4))
!   equivalence (oa(1,1,3),hprime(1,1,5))
!   equivalence (oa(1,1,4),hprime(1,1,6))
!   equivalence (ol(1,1,1),hprime(1,1,7))
!   equivalence (ol(1,1,2),hprime(1,1,8))
!   equivalence (ol(1,1,3),hprime(1,1,9))
!   equivalence (ol(1,1,4),hprime(1,1,10))
#ifdef FLOW_BLOCKING
!   equivalence (omax(1.1),hprime(1,1,11))
#endif
!
   logical lsmth,lenh
   real     ::  rk,enhc
!
   data lsmth/.true./,rk/0.7e-6/,lenh/.false./,enhc/1.2/
!
   integer  ::  i,j,ie,iw,jn,js,l,m,nn,n,nmax,ns,nmin,ivar
   real     ::  phi,dlat,degrad,rlat,slma,oroa,vara,width,fac
!
!  restricted data 
!
   integer  ::  numi(latg_),irowsep
!
!  data  numi/latg_*lonf_/
!
!  data numi  /
!    .   30,  30,  30,  40,  48,  56,  60,  72,  72,  80,  90,  90,
!    .   96, 110, 110, 120, 120, 128, 144, 144, 144, 144, 154, 160,
!    .  160, 168, 168, 180, 180, 180, 180, 180, 180, 192, 192, 192,
!    .  192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192 /
!
   numi=lonf_
!
!  do i = latg_/2+1,latg_
!    numi(i) = numi(latg_+1-i)
!  enddo
!
!  pointing to memory addresses of hprime
!
   var=>hprime(:,:,1)
   var4=>hprime(:,:,2)
   oa=>hprime(:,:,3:6)
   ol=>hprime(:,:,7:10)
#ifdef FLOW_BLOCKING
   omax=>hprime(:,:,11)
#endif
!
#ifdef DFS
   call get_dfs_dim(jcap_,levs_,ngases_,nwater_,lonf_,latg_)
   allocate(orsdfs(mt,jl))
   irowsep=0
#else
   irowsep=0
#endif
!
!  set constants and zero fields
!
   degrad = 180./pi_
!
!- open(unit=11,form='formatted',err=900) ! average 
!- open(unit=12,form='formatted',err=900) ! variance 
!- open(unit=13,form='formatted',err=900) ! maximum 
!- open(unit=14,form='formatted',err=900) ! sea-land-lake-mask 
!
   read(11,11) zavg
   read(12,11) zvar
   read(13,11) zmax
   read(14,12) zslm
11    format(20i4)
12    format(80i1)
!
   call mtn_gaussian_lat(jm,cosclt,wgtclt)
!
   do j = 1,jm/2
     rclt(j)=acos(cosclt(j))
   enddo
!
   do j = 1,jm/2
     phi = rclt(j) * degrad
     xlat(j) = 90. - phi
     xlat(jm-j+1) =  phi - 90.
   enddo
#ifdef DFS
! 
! calculate coslat
!
   call dfs_sincos_lat
   print*,'sinlat=',sinlat
   dlat=180./jm
   rlat=pi_/jm
!
   do j = 1,jm/2
     xlatd(j) = 90.-(j-0.5)*dlat
     xlatd(jm-j+1) = -xlatd(j)
     !print*,'xlat=',xlatd(j),xlatd(jm-j+1)
   enddo
#endif
!
!     compute mountain data : oro slm var oc
!
   call mtn_make_orography(zavg,zvar,zslm,oro,slm,var,var4,glat,               &
                           ist,ien,jst,jen,im,jm,imn,jmn,xlat,numi)
!
!     compute mountain data : oa ol
!
   call mtn_make_statistics(zavg,zmax,var,glat,oa,ol,iwork,                    &
#ifdef FLOW_BLOCKING
               omax,                                                           &
#endif
               work1,work2,work3,work4,                                        &
               work5,work6,                                                    &
               ist,ien,jst,jen,im,jm,imn,jmn,xlat,numi)
#ifdef AQUA_PLANET
!
   do j = 1,jm
     do i = 1,im
        slm(i,j)=0.0
        oro(i,j)=0.0
        var(i,j)=0.0
     enddo
   enddo
#endif
!
!  enhancement
!
   if(lenh) then
     do j = 1,jm
       do i = 1,im
         oro(i,j)=oro(i,j)+var(i,j)*enhc
       enddo
     enddo
   endif
!
!  remove isolated points
!
   do j = 2,jm-1
     jn=j-1
     js=j+1
     do i = 1,im
       iw=mod(i+im-2,im)+1
       ie=mod(i,im)+1
       slma=slm(iw,jn)+slm(i,jn)+slm(ie,jn)+                                   &
              slm(iw,j )          +slm(ie,j )+                                 &
              slm(iw,js)+slm(i,js)+slm(ie,js)
       if(slm(i,j).eq.0..and.slma.eq.8.) then
         slm(i,j)=1.
         oroa=(oro(iw,jn)+oro(i,jn)+oro(ie,jn)+                                &
                  oro(iw,j )          +oro(ie,j )+                             &
                  oro(iw,js)+oro(i,js)+oro(ie,js))/8.
         vara=(var(iw,jn)+var(i,jn)+var(ie,jn)+                                &
                  var(iw,j )          +var(ie,j )+                             &
                  var(iw,js)+var(i,js)+var(ie,js))/8.
         print '("sea ",2f8.0," modified to land",2f8.0," at ",2i8)',          &
                     oro(i,j),var(i,j),oroa,vara,i,j
         oro(i,j)=oroa
         var(i,j)=vara
       elseif(slm(i,j).eq.1..and.slma.eq.0.) then
         slm(i,j)=0.
         oroa=(oro(iw,jn)+oro(i,jn)+oro(ie,jn)+                                &
                  oro(iw,j )          +oro(ie,j )+                             &
                  oro(iw,js)+oro(i,js)+oro(ie,js))/8.
         vara=(var(iw,jn)+var(i,jn)+var(ie,jn)+                                &
                  var(iw,j )          +var(ie,j )+                             &
                  var(iw,js)+var(i,js)+var(ie,js))/8.
         print '("land",2f8.0," modified to sea ",2f8.0," at ",2i8)',          &
                    oro(i,j),var(i,j),oroa,vara,i,j
         oro(i,j)=oroa
         var(i,j)=vara
       endif
     enddo
   enddo
!
!     zero over ocean
!
   do i = 1,im
     do j = 1,jm
       if(slm(i,j).eq.0.) then
         var4(i,j) = 0.
         oa(i,j,1) = 0.
         oa(i,j,2) = 0.
         oa(i,j,3) = 0.
         oa(i,j,4) = 0.
         ol(i,j,1) = 0.
         ol(i,j,2) = 0.
         ol(i,j,3) = 0.
         ol(i,j,4) = 0.
#ifdef FLOW_BLOCKING
         omax(i,j) = 0.
#endif
       endif
     enddo
   enddo
!
!  spectrally truncate orography
!
   call mtn_trans_wave_grid(1,oro,ors,0,0.,im,jm,nm,nr)
!
!  spectral smoothing
!
   if(lsmth) then
     l=1
     do m = 1,jcap_+1
       do nn = m,jcap_+1
         n=nn
!
!  hoskins(1990)
!
!        nnx=(n-1)*(n-1)*n*n
!        fac=exp(-rk*float(nnx))
!
!  navarra(1994)
!
!        fac=1.-float(n-1)/float(jcap_+1)
!
!  empirical (kanamitsu, 1998)
!
!        beta=1.
!        alfa=0.6
!           fac=1.-alfa*(float(n-1)/float(jcap_+1))**beta
!
!  empirical (ebisuzaki, 1998)
!
         width=0.6
         nmax=sqrt(float((jcap_)*(jcap_+1)))
         nmin=(1.-width)*nmax
         ns=sqrt(float((n-1)*n))
         if(ns.lt.nmin) then
           fac=1.
         else
           fac=(float(nmax)-float(ns))/float(nmax-nmin)
         endif
!
!  empirical (navarra, 1994)
!
!        alfa=32.
!        beta=8.
!        fac=exp(-alfa*(n/jcap_)**(2*beta))
!
         orsmth(l  )=fac*ors(l  )
         orsmth(l+1)=fac*ors(l+1)
         l=l+2
       enddo
     enddo
   endif
!
!  output fields
!
   write(51) slm
   call file_write_bin(63,slm,im,jm,1,irowsep)
   call mtn_quick_print(slm,im,jm,1)
!
   if(mtnvar.eq.1) then
     write(53) var
   else
      write(53) hprime
   endif
!
   call file_write_bin(63,hprime,im,jm,mtnvar,irowsep)
!
   do ivar = 1,mtnvar
     call mtn_quick_print(hprime(1,1,ivar),im,jm,1.)
   enddo
!
#ifdef DFS
   call mtn_trans_wave_grid(-101,orf,ors,0,0.,im,jm,nm,nr)
#else
   call mtn_trans_wave_grid(-1,orf,ors,0,0.,im,jm,nm,nr)
#endif
   call file_write_bin(63,oro,im,jm,1,irowsep)
   call file_write_bin(63,orf,im,jm,1,irowsep)
   write(55) oro
   call mtn_quick_print(oro,im,jm,0.001)
!
   write(52) orf
   call mtn_quick_print(orf,im,jm,0.001)
#ifdef DFS
   call dfs_fft_driver(1,orf,lonf_,latg_,1,orsdfs,mt,jl,1,jl,                  &
                       1,1,1,coslat,1)
#ifdef ALIASED
   call dfs_cut_alias(orsdfs,mt,0,jlg-1,1)
#endif
   write(54) orsdfs
#else
   write(54) ors
   write(61) orsmth
#endif
!
#ifdef DFS
   call mtn_trans_wave_grid(-101,orf,orsmth,0,0.,im,jm,nm,nr)
   call dfs_fft_driver(1,orf,lonf_,latg_,1,orsdfs,mt,jl,1,jl,                  &
                       1,1,1,coslat,1)
   write(61) orsdfs
   call dfs_fft_driver(-1,oro,lonf_,latg_,1,orsdfs,mt,jl,1,jl,                 &
                        1,1,1,coslat,1)
#else
   call mtn_trans_wave_grid(-1,orf,orsmth,0,0.,im,jm,nm,nr)
#endif
!
#ifdef DFS
   write(62) oro
   call file_write_bin(63,oro,im,jm,1,irowsep)
#else
   write(62) orf
   call file_write_bin(63,orf,im,jm,1,irowsep)
#endif
!
   return
   end subroutine mtn_solver
!-------------------------------------------------------------------------------
