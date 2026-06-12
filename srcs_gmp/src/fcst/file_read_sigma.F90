#include "define.h"
#ifndef HYBRID
   subroutine file_read_sigma(n,fhour,idate,gz,q,te,di,ze,rq,sl,si,z00)
#else
   subroutine file_read_sigma(n,fhour,idate,gz,q,te,di,ze,rq,ak5,bk5,z00)
#endif
!-------------------------------------------------------------------------------
   use paramodel, only : igen_,jcap_,latg_,levh_,levp1_,levs_,lnt22_,          &
                         lnt2_,lonf_,LNT22S,                                   &
                         nwater_,ngases_,ntotal_,                              &
                         nwmass_,icloud_,kcloud_,igases_,kgases_
#ifdef DFS
   use dfsvar, only : mt,jl,lnt2a,lnt2,tavexy,mtg,jlg,levsp,levhp,             &
#ifdef MP
                      levs,levh,ntotal,                                        &
#endif
                      ib,jbw,dio,coslat,iope
#endif
#ifdef MP
   use commpi
   use paramodel, only : lnt22p_
#endif
   use comio
   use constant, only  : rerth_,g_
   use comfcst, only   : kdum, kdum2, kens
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! subprogram:    file_read_sigma       reads sigma level spectral coefficients.
!
! abstract: reads a complete set of sigma level spectral coefficients
!   at a single time to be used to start the model forecast.
!   the subroutine compares si and sl (the models vertical
!   structure) computed in dyn_sigma_setup with the si and sl of the
!   input coefficients in order to make sure the coefficients
!   were generated under the same vertical structure.
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call file_read_sigma (n,fhour,idate,gz,q,te,di,ze,rq,sl,si,z00)
!   input argument list:
!     n        - fortran unit number for file containing the
!                complete set of sigma level coefficients.
!     sl       - sigma layers     computed in dyn_sigma_setup.
!     si       - sigma interfaces computed in dyn_sigma_setup.
!
!   output argument list:
!     fhour    - forecast hour of the set of coefficients read
!                from unit n.
!     idate    - idate(1)=initial hour (gmt) of forecast from
!                         which coefficients were made.
!                idate(2)=month (1-12).
!                idate(3)=day of the month.
!                idate(4)=year of the century.
!     gz       - laplacian of topography.
!                gz is multiplied by the constant snnp1 array.
!                gz is then multiplied by the gravity constant and
!                divided by the square of the radius of the earth.
!     q        - ln(psfc)          coefficients.
!     te       - temperature       coefficients.
!     di       - divergence        coefficients.
!     ze       - vorticity         coefficients.
!     rq       - specific humidity coefficients.
!     z00      - mean topography.
!                z00 is set equal to gz(1) after gz is read.
!
!   input files:
!     unit n   - complete set of sigma level spectral coefficients.
!
!   output files:
!     output   - print file.
!
!-------------------------------------------------------------------------------
#include "abort.h"
#ifdef MP
   real, allocatable  ::  spec(:),coef1(:,:),coef2(:,:)
#endif
!
#ifdef DFS
#define lnt2_ lnt2a
#ifdef MP
#define lnt22p_ lnt2
#define lnt22_ lnt2a
#else
#define lnt22_ lnt2a
#endif
#endif
#ifdef MP
#define LNT22S lnt22p_
#else
#define LNT22S lnt22_
#endif
!
   integer  ::  n,idate(4)
   real     ::  fhour,z00
#ifdef DFS
   real     ::  gz(mt*jlg),q(mt*jlg)
   real     ::  te(mt*jlg,levsp),di(mt*jlg,levsp)
   real     ::  ze(mt*jlg,levsp),rq(mt*jlg,levhp)
#ifdef MP
   real   , target , dimension(mt*jl ,2)              ::  WORKL2
   real   , target , dimension(mt*jl ,3*levs +levh )  ::  WORKL3
   real   , target , dimension(mt*jlg,2)              ::  WORKG2
   real   , target , dimension(mt*jlg,3*levsp+levhp)  ::  WORKG3
   real   , pointer, dimension(:)                     ::  gz1,q1
   real   , pointer, dimension(:,:)                   ::  te1,di1,ze1,rq1
#endif
#else /* SPH */
   real     ::  gz(LNT22S),q(LNT22S),te(LNT22S,levs_),di(LNT22S,levs_)
   real     ::  ze(LNT22S,levs_),rq(LNT22S,levh_)
#endif /* DFS end */
   real     ::  si(levp1_),sl(levs_)
   real     ::  xi(levp1_),xl(levs_)
#ifdef HYBRID
   real     ::  ak5(levp1_),bk5(levp1_)
   real     ::  xak(levp1_),xbk(levp1_)
#endif
!
   integer  ::  k,i,j,nc,kc,nt,kt,nw,kw
   integer  ::  iwater,igases
   real     ::  dummy(kdum), dummy2(kdum2),ensemble(kens)
!
   real     ::  work(LNT22S,levs_)
   real     ::  waves,xlayers,trun,order,realform,gencode,rlond,rlatd,rlonp
   real     ::  rlonr,rlatr,gases,subcen,ppid,slid,vcid,vmid,vtid
   real     ::  pdryini,water,ga2,rlatp
!
   integer, save  ::   ifp
   data ifp/0/
!-------------------------------------------------------------------------------
#if defined(DFS) && defined (MP)
   gz1=>WORKL2(1:mt*jl,1)
   q1 =>WORKL2(1:mt*jl,2)
   te1=>WORKL3(1:mt*jl,       1:  levs)
   di1=>WORKL3(1:mt*jl,1*levs+1:2*levs)
   ze1=>WORKL3(1:mt*jl,2*levs+1:3*levs)
   rq1=>WORKL3(1:mt*jl,3*levs+1:3*levs+levh)
#else
   if(ifp.eq.0) then
     call sph_comp_index
     ifp=1
   endif
#endif
!
!  end addition
!
#ifdef MP
   allocate (spec(lnt22_))
   allocate (coef1(lnt22_,levs_))
   allocate (coef2(lnt22_,levs_))
#endif
!
!     spectral data file format
!     lab
!     hour,idate(4),si(levp1_),sl(levs_)
!     zln q te di ze
!
#ifdef MP
   if (iope) then
#endif
#ifdef ASSIGN
     call assign('assign -R')
#endif
     rewind n
     read(n)lab
#ifdef PRINT
     print 3000,lab,n
3000 format(1x,'file_read_sigma lab  ',4a8,' n=',i3)
#endif
#ifdef MP
   endif
#endif
!
#ifdef MP
   if( iope ) then
!
#endif
#ifndef HYBRID
   read(n,err=201)fhour,idate,(xi(k),k=1,levp1_),(xl(k),k=1,levs_)             &
#else
   read(n,err=201)fhour,idate,(xak(k),k=1,levp1_),(xbk(k),k=1,levp1_)          &
#endif
          ,dummy,waves,xlayers,trun,order,realform,gencode                     &
          ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,water                           &
          ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid                &
          ,pdryini,dummy2,gases
   igases=nint(gases)
   iwater=nint(water)
#ifndef NOPRINT
   print *,'file_read_sigma unit,fhour,idate=',n,fhour,idate
   print *,' number of water input = ',iwater
   print *,' number of gases input = ',igases
   print *,' number of wave  = ',LNT22S
#endif
   goto 202
201 continue
   rewind n
   read(n) lab
#ifdef MP
!
#endif
#ifndef HYBRID
   read(n,err=201)fhour,idate,(xi(k),k=1,levp1_),(xl(k),k=1,levs_)
#else
   read(n,err=201)fhour,idate,(xak(k),k=1,levp1_),(xbk(k),k=1,levp1_)
#endif
   do i = 1,kdum
     dummy(i)=0.
   enddo
   waves=jcap_
   xlayers=levs_
   trun=1.
   order=2.
   realform=1.
   gencode=igen_
   rlond=lonf_
   rlatd=latg_
   rlonp=lonf_
   rlatp=latg_
   rlonr=lonf_
   rlatr=latg_
   gases=0.
   water=1.
   pdryini=0.
   subcen=0.
!
   do i = 1,kens
     ensemble(i)=0.
   enddo
!
   ppid=0.
   slid=0.
   vcid=0.
   vmid=0.
   vtid=0.
!
   do k = 1,kdum2
     dummy2(k)=0.
   enddo
!
   igases=nint(gases)
   iwater=nint(water)
#ifndef NOPRINT
   print *,'file_read_sigma old format unit,fhour,idate=',n,fhour,idate
   print *,' number of water input = ',iwater
   print *,' number of gases input = ',igases
#endif
202 continue
!
#ifndef DRY_MODEL
   if(iwater.le.0) then
     print *,'moisture not predicted'
     call MPABORT
   endif
#else
   iwater=0
#endif
!
#ifdef MP
   endif
!
   call mpbcastr(fhour,1)
   call mpbcasti(idate,4)
#ifndef HYBRID
   call mpbcastr(xi,levp1_)
   call mpbcastr(xl,levs_)
#else
   call mpbcastr(xak,levp1_)
   call mpbcastr(xbk,levp1_)
#endif
   call mpbcasti(igases,1)
   call mpbcasti(iwater,1)
#endif /* MP end */
!
   ga2=g_/(rerth_*rerth_)
#ifdef MP
   if( iope ) then
     read(n)(spec(i),i=1,lnt2_)
#ifdef DFS
     z00=spec(jcap_+1)
#ifndef NOPRINT
     write(6,*)' file_read_sigma gz z00=',z00
#endif
   endif
   call mpbcastr(z00,1)
   call mpsf2p(spec,mtg,jlg,gz1,mt,jl,1)
#else /* SPH */
     z00=spec(1)
     call spcshfli(spec,lnt22_,1,jcap_,lwvdef)
     do j = 1,lnt2_
       spec(j)=spec(j)*snnp1(j)*ga2
     enddo
     call spcshflo(spec,lnt22_,1,jcap_,lwvdef)
   endif
   call mpbcastr(z00,1)
   call mpsf2p(spec,lnt22_,gz,lnt22p_,1)
#endif /* DFS end */
#else /* ~MP */
   read(n)(gz(i),i=1,lnt2_)
#ifdef DFS
   z00=gz(jcap_+1)
#else
   z00=gz(1)
   do j = 1,lnt2_
     gz(j)=gz(j)*snnp1(j)*ga2
   enddo
#endif /* DFS end */
#endif /* MP end */
!
#ifdef MP
   if( iope ) then
#ifdef DFS
     read(n)(spec(i),i=1,lnt2_),tavexy(:)
     write(6, *)' file_read_sigma q '
     call print_maxmin_seven(spec,lnt2_,lnt22_,1,1,1,'q file_read_sigma')
   endif
   call mpbcastr(tavexy,levs_)
   call mpsf2p(spec,mtg,jlg,q1,mt,jl,1)
#else /* SPH */
     read(n)(spec(i),i=1,lnt2_)
     print *,' file_read_sigma q '
   endif
   call mpsf2p(spec,lnt22_,q,lnt22p_,1)
#endif /* DFS end */
#else /* ~MP */
#ifdef DFS
   read(n)(q(i),i=1,lnt2_),tavexy(:)
#else
   read(n)(q(i),i=1,lnt2_)
#endif
#endif /* MP end */
!
#ifdef MP
   if( iope ) then
     do k = 1,levs_
       read(n)(coef1(i,k),i=1,lnt2_)
     enddo
#ifndef NOPRINT
     call print_maxmin_seven(coef1,lnt2_,lnt22_,levs_,1,levs_,                 &
                                                          'te file_read_sigma')
#endif
   endif
#ifdef DFS
   call mpsf2p(coef1,mtg,jlg,te1,mt,jl,levs_)
#else
   call mpsf2p(coef1,lnt22_,te,lnt22p_,levs_)
#endif
#else				/* not MP */
   do k = 1,levs_
     read(n)(te(i,k),i=1,lnt2_)
   enddo
#endif
!
#ifdef MP
   if( iope ) then
     do k = 1,levs_
       read(n)(coef1(i,k),i=1,lnt2_)
       read(n)(coef2(i,k),i=1,lnt2_)
     enddo
#ifndef NOPRINT
     call print_maxmin_seven(coef1,lnt2_,lnt22_,levs_,1,levs_,                 &
                                                          'di file_read_sigma')
     call print_maxmin_seven(coef2,lnt2_,lnt22_,levs_,1,levs_,                 &
                                                          'vo file_read_sigma')
#endif
   endif
!
#ifdef DFS
   call mpsf2p(coef1,mtg,jlg,di1,mt,jl,levs_)
   call mpsf2p(coef2,mtg,jlg,ze1,mt,jl,levs_)
#else /* SPH */
   call mpsf2p(coef1,lnt22_,di,lnt22p_,levs_)
   call mpsf2p(coef2,lnt22_,ze,lnt22p_,levs_)
#endif /* DFS end */
#else /* ~MP */
   do k = 1,levs_
     read(n)(di(i,k),i=1,lnt2_)
     read(n)(ze(i,k),i=1,lnt2_)
   enddo
#endif /* MP end */
!
#ifdef MP
   if( iope ) then
#endif
   do k = 1,levs_
#ifdef MP
     read(n)(coef1(i,k),i=1,lnt2_)
   enddo
#ifndef NOPRINT
   call print_maxmin_seven(coef1,lnt2_,lnt22_,levs_,1,levs_,                   &
                                                          'rq file_read_sigma')
#endif
   endif
#ifdef DFS
   call mpsf2p(coef1,mtg,jlg,rq1,mt,jl,levs_)
#else
   call mpsf2p(coef1,lnt22_,rq,lnt22p_,levs_)
#endif
#else /* ~MP */
     read(n)(rq(i,k),i=1,lnt2_)
   enddo
!
   call print_maxmin_seven(rq,lnt2_,lnt22_,levs_,1,levs_,                      &
                                              'rq aft read in file_read_sigma')
#endif /* MP end */
   if(iwater.ge.2) then
     do nw = 2,iwater
       kw = (nw-1)*levs_+1
#ifdef MP
       if( iope ) then
         do k = 1,levs_
           read(n)(coef1(i,k),i=1,lnt2_)
         enddo
       endif
       if(nwater_.ge.nw) then
#ifdef DFS
         call mpsf2p(coef1,mtg,jlg,rq1(1,kw),mt,jl,levs_)
#else
         call mpsf2p(coef1,lnt22_,rq(1,kw),lnt22p_,levs_)
#endif
       endif
#else /* ~MP */
       do k = 1,levs_
         read(n)(work(i,k),i=1,lnt2_)
       enddo
       if(nwater_.ge.nw) then
         do k = 1,levs_
           do i = 1,lnt2_
             rq(i,k+kw-1) = work(i,k)
           enddo
         enddo
       endif
#endif /* MP */
     enddo
   endif
!
   if(nwater_.gt.iwater) then
     do nw = iwater+1,nwater_
       kw = (nw-1)*levs_+1
       do k = kw,kw+levs_-1
#if defined(DFS) && defined(MP)
         rq1(1:mt*jl,k)= 0.0
#else
         do i = 1,LNT22S
           rq(i,k) = 0.0
         enddo
#endif
       enddo
     enddo
   endif
!
   if(igases.ge.1) then
     do nt = 1,igases
       kt = (nt-1)*levs_+kgases_
#ifdef MP
       if( iope ) then
         do k = 1,levs_
           read(n)(coef1(i,k),i=1,lnt2_)
         enddo
       endif
       if(ngases_.ge.nt) then
#ifdef DFS
         call mpsf2p(coef1,mtg,jlg,rq1(1,kt),mt,jl,levs_)
#else
         call mpsf2p(coef1,lnt22_,rq(1,kt),lnt22p_,levs_)
#endif
       endif
#else /* ~MP */
       do k = 1,levs_
         read(n)(work(i,k),i=1,lnt2_)
       enddo
       if(ngases_.ge.nt) then
         do k = 1,levs_
           do i = 1,lnt2_
             rq(i,k+kt-1) = work(i,k)
           enddo
         enddo
       endif
#endif /* MP end */
     enddo
   endif
!
   if(ngases_.gt.igases) then
     do nt = igases+1,ngases_
       kt = (nt-1)*levs_ + kgases_
       do k = kt,kt+levs_-1
#if defined(DFS) && defined(MP)
         rq1(:,k)= 0.0
#else
         do i = 1,LNT22S
           rq(i,k) = 0.0
         enddo
#endif
       enddo
     enddo
   endif
!
#if defined(DFS) && defined(MP)
   call mpmn2m (WORKL2,jl,     WORKG2,jlg,mt      ,2)
   call mpmn2mz(WORKL3,jl,levs,WORKG3,jlg,mt,levsp,3+ntotal)
   do i = 1,mt*jlg
     gz(i)=WORKG2(i,1)
     q (i)=WORKG2(i,2)
   enddo
!
   do k = 1,levsp
     do i = 1,mt*jlg
       te(i,k)=WORKG3(i,k)
       di(i,k)=WORKG3(i,levsp+k)
       ze(i,k)=WORKG3(i,2*levsp+k)
     enddo
   enddo
!
   do k = 1,levhp
     do i = 1,mt*jlg
       rq(i,k)=WORKG3(i,3*levsp+k)
     enddo
   enddo
#endif
#ifndef NOPRINT
#ifdef MP
   if( iope ) then
#ifndef DFS
     call print_maxmin_six(rq,LNT22S,levh_,1,levh_,                            &
                                           'rq global coef in file_read_sigma')
#endif
#endif /* MP end */
#ifndef HYBRID
   do k = 1,levs_
     xl(k)=xl(k)-sl(k)
   enddo
!
   do k = 1,levp1_
     xi(k)=xi(k)-si(k)
   enddo
#else
   do k = 1,levp1_
     xak(k)=xak(levp1_+1-k)/1000.
     xbk(k)=xbk(levp1_+1-k)
   enddo
   do k = 1,levp1_
     xak(k)=xak(k)-ak5(k)
     xbk(k)=xbk(k)-bk5(k)
   enddo
#endif
#ifndef HYBRID
   print 100,(xi(k),k=1,levp1_)
   print 100,(xl(k),k=1,levs_)
#else
   print 100,(xak(k),k=1,levp1_)
   print 100,(xbk(k),k=1,levp1_)
#endif
100 format(1x, 12(e9.3))
   print 101,n,fhour,idate
101 format (1x, 'if above two rows not zero,inconsistency in sig.def'          &
              ,' on n=',i2,2x,f6.1,2x,4(i6))
#ifdef MP
   endif
#endif
#endif /* ~NOPRINT end */
!
#ifdef MP
   deallocate (spec)
   deallocate (coef1)
   deallocate (coef2)
#endif
!
   return
   end subroutine file_read_sigma
