#include "define.h"
   subroutine mtn_solver_scm(im,jm,nm,nr,nw)
!-------------------------------------------------------------------------------
   use paramodel, only : jcap_,lonf_,latg_,                                    &
                         imn=>imn_,jmn=>jmn_,mtnres_,                          &
                         mtnvar=>mtnvar_,nnw=>lnt2_,                           &
                         ngases_,nwater_,ntotal_,levs_
   use constant, only : pi_
#ifdef DFS
   use dfsvar,   only : mt,jl,jlg,get_dfs_dim,ib,jbw,coslat
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
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
#ifdef SMP
   real                 ::  xlon, ylat,                                        &
                            delxn,xnsum,xland,xwatr,xl1,xs1,xw1,xw2,xv2,xw4,   &
                            height,slm1p,oro1p,var1p
   integer              ::  lsmask,ipos,jpos,istr,iend,jstr,jend
#endif
   real                 ::  ors(nnw)
   real                 ::  orsmth(nnw)
   real                 ::  orf(lonf_,latg_)
!
   real, pointer        ::  var(:,:)
   real, pointer        ::  var4(:,:),oa(:,:,:)
   real, pointer        ::  ol(:,:,:)
!   real                 ::  var4(lonf_,latg_),oa(lonf_,latg_,4)
!   real                 ::  ol(lonf_,latg_,4)
   real                 ::  work(lonf_,latg_),work1(lonf_,latg_),              &
                            work2(lonf_,latg_),work3(lonf_,latg_),             &
                            work4(lonf_,latg_),                                &
                            work5(lonf_,latg_),work6(lonf_,latg_),             &
                            glat(jmn)
   real, target         ::  hprime(lonf_,latg_,10)
!   real                 ::  hprime(lonf_,latg_,10)
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
!  add msk
!
   var=>hprime(:,:,1)
   var4=>hprime(:,:,2)
   oa=>hprime(:,:,3:6)
   ol=>hprime(:,:,7:10)
#ifdef DFS
   call get_dfs_dim(jcap_,levs_,ngases_,nwater_,lonf_,latg_)
   print*,'hoon:mt,jl=',mt,jl
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
#ifdef SMP
!
#define GXLON _gxlon_
#define GYLAT _gylat_
#define GLSMSK _lsmsk_
!                           ... observation position
   xlon = GXLON
   ylat = GYLAT
   lsmask = GLSMSK
!
   delxn = 360./imn      ! mountain data resolution
   ipos = int(xlon/delxn) + 1
   jpos = int((ylat+90.0)/delxn) + 1
   print *, 'POS ',imn,jmn,ipos,jpos,delxn
   istr = ipos-4
   iend = ipos+4
   jstr = jpos-4
   jend = jpos+4
!
   xnsum = 0.0
   xland = 0.0
   xwatr = 0.0
   xl1 = 0.0
   xs1 = 0.0
   xw1 = 0.0
   xw2 = 0.0
   xv2 = 0.0
   xw4 = 0.0
!
   do i = istr,iend
     do j = jstr,jend
       xland = xland + float(zslm(i,j))
       xwatr = xwatr + float(1-zslm(i,j))
       xnsum = xnsum + 1.
       height = float(zavg(i,j))
       if(height.lt.-990.) height = 0.0
       xl1 = xl1 + height * float(zslm(i,j))
       xs1 = xs1 + height * float(1-zslm(i,j))
       xw1 = xw1 + height
       xw2 = xw2 + height ** 2
       xv2 = xv2 + float(zvar(i,j)) ** 2
       print 1967, i,j,xland,xwatr,xnsum,height,xl1,xs1,xw1,xw2,xv2
1967   format(2i5,9e13.5)
     enddo
   enddo
!
   if(xnsum.gt.1.) then
     slm1p = float(nint(xland/xnsum))
     if (float(lsmask).eq.slm1p) then
       if(slm1p.ne.0.) then
         oro1p = xl1 / xland
       else
         oro1p = xs1 / xwatr
       endif
     else
       print *, 'Land-Sea mask is not consistent with "define.h" ! '
       call abort
     endif
     var1p=sqrt(max((xv2+xw2)/xnsum-(xw1/xnsum)**2,0.))
     print 1966, xnsum,slm1p,oro1p,var1p
1966  format(4e13.5)
   endif
!
   do j = 1,latg_
     do i = 1,lonf_
       oro(i,j) = oro1p
       var(i,j) = var1p
     enddo
   enddo
!
   if(mtnvar.eq.1) then
     write(53) var
   else
     write(53) hprime
   endif
!
   write(62) oro
#endif
!
   return
   end subroutine mtn_solver_scm
!-------------------------------------------------------------------------------
