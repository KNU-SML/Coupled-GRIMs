#include "define.h"
#ifndef HYBRID
   subroutine file_write_sigma(n,fhour,idate,q,te,di,ze,rq,sl,si,gz,           &
#else
   subroutine file_write_sigma(n,fhour,idate,q,te,di,ze,rq,ak5,bk5,gz,         &
#endif
#ifdef DFS
                     tavexy,                                                   &
#else
                     z00,                                                      &
#endif
#ifdef SMP
                     vvel,                                                     &
#endif
                     itpdt)
!-------------------------------------------------------------------------------
!
! subprogram:    file_write_sigma      writes sigma level spectral coefficients.
!
! abstract: writes a complete set of forecast sigma level
!   spectral coefficients for all model variables.
!
! program history log:
!   1988-04-29  joseph sela
!   1988-11-02  mark rozwodoski  changed second record to selalabel.
!   2000-01-01  hann-ming henry juang  mpi                                     
!   2000-01-01  song-you hong          physcis options                         
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2006-04-01  hoon park              double-fourier spectral (dfs)           
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                 
! usage:    call file_write_sigma (n,fhour,idate,z,q,te,di,ze,rq,sl,si,gz,z00)
!   input argument list:
!     n        - fortran unit number for file to be written to.
!     fhour    - forecast hour of the set of coefficients.
!     idate    - idate(1)=initial hour (gmt) of forecast.
!                idate(2)=month (1-12).
!                idate(3)=day of the month.
!                idate(4)=year of the century.
!     q        - ln(psfc)          coefficients.
!     te       - temperature       coefficients.
!     di       - divergence        coefficients.
!     ze       - vorticity         coefficients.
!     rq       - specific humidity coefficients.
!     sl       - sigma layers     computed in dyn_sigma_setup.
!     si       - sigma interfaces computed in dyn_sigma_setup.
!     gz       - laplacian of topography.
!     z00      - mean topography.
!     itpdt    - flag for time level t, t+dt or initialized (1, 2 or 3)
!                4 for sigit, 5 for sigitdt file names
!
!   output argument list:
!     z        - topography spectral coefficients.
!
!   output files:
!     unit n   - complete set of forecast sigma level spectral
!                coefficients.
!     output   - print file.
!
!-------------------------------------------------------------------------------
#ifdef DFS
   use dfsvar, only : lnt2,mtg,jlg,lnt2a,iope,mt,jl,levsp,levhp,levs,levh,ntotal
#else
   use paramodel, only : LNT22S
#ifdef MP
   use paramodel, only : lnt22p_
   use commpi
#endif
#endif
   use paramodel, only : lnt2_,lnt22_,levs_,ntotal_,jcap_,lonf_,latg_,         &
                         ngases_,nwater_,levp1_,levh_
   use comio
   use comfcst, only   : kdum, kdum2, kens
#ifdef NISLQ
   use nislq, only     : slq_q1,slq_q2
#ifdef DFS
   use dfsvar, only    : ib,jbw,coslat
#endif
#endif
#ifdef DCMIP
#ifndef DFS
   use dcmip_grims, only : gphf
#endif
#endif /* DCMIP end */
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#include "abort.h"
!
#ifdef MP
   real, allocatable  ::  coef1(:,:),coef2(:,:),coef3(:,:)
#endif
#ifdef DFS
#define LNT22S lnt2
#define lnt22_ lnt2a
#define lnt2_  lnt2a
#endif /* DFS end */
!
#ifdef SMP
   real     ::  vvel(lnt22_,levs_)
#endif
   real     ::  gz(LNT22S),q(LNT22S),spec(lnt22_),                             &
                te(LNT22S,levs_),di(LNT22S,levs_),ze(LNT22S,levs_),            &
                rq(LNT22S,levh_),si(levp1_),sl(levs_)
   integer  ::  idate(4),itpdt,n
#ifdef HYBRID
   real     ::  ak5(levp1_),bk5(levp1_),ak5x(levp1_),bk5x(levp1_)
#endif
   real     ::  fhour
#ifdef DFS
   real     ::  tavexy(levs_)
#else
   real     ::  z00
#endif
   real     ::  orog(lnt22_)
   real     ::  dummy(kdum),dummy2(kdum2),ensemble(kens)
!
   character*128 fni(5)
   integer nchi(5),ncho,k,i,nt,kt,kt1,kt2
   data fni/'sig','sigp','sigi','sigit','sigitdt'/
   data nchi/3,4,4,5,7/
   character*128 fno
   real waves,xlayers,trun,order,realform
   real vtid,vmid,vcid,slid,ppid,subcen,gencode
   real rlond,rlatd,rlonp,rlatp,rlonr,rlatr,gases,water,pdryini
!-------------------------------------------------------------------------------
#ifdef NISLQ
!
! nislq grid transform
!
   if(itpdt.eq.1 .or. itpdt.eq.4) then
#ifdef DFS
     call dfs_fft_driver(1,slq_q1,ib,jbw,levh,rq,mt,jlg,levhp,jlg,             &
                            levsp,levs,ntotal_,coslat,1)
#else
     call sph_fft_driver(1,slq_q1,rq,ntotal_)
#endif
   else if (itpdt.eq.2 .or. itpdt.eq.5) then
#ifdef DFS
     call dfs_fft_driver(1,slq_q2,ib,jbw,levh,rq,mt,jlg,levhp,jlg,             &
                            levsp,levs,ntotal_,coslat,1)
#else
     call sph_fft_driver(1,slq_q2,rq,ntotal_)
#endif
   endif
#endif
#ifdef NISLQ_GRIB
!
! nislq_write
!
   call file_write_slq(idate,itpdt,fhour)
#endif /* NISLQ_GRIB end */
#ifndef DCMIP
#ifndef DFS
!
!  get topography from history file
!
   if( iope ) then
     close(n)
     open(unit=n,file='./sigit',form='unformatted',err=700)
     go to 701
700 continue
     write(6,*) ' error in opening file sigit at file_write_sigma'
     call MPABORT
701 continue
     rewind n
     read(n)
#ifndef NOPRINT
     write(6,*) ' read lab '
#endif
     read(n)
#ifndef NOPRINT
     write(6,*) ' read fhour idate  '
#endif
     read(n)( orog(i),i=1,lnt2_)
#ifndef NOPRINT
     write(6,*) ' read gz '
#endif
     close(n)
   endif
#endif          /* not DFS */
#endif /* ~DCMIP end */
!
   if( iope ) then
     if(itpdt.lt.4) then
       call file_name(fni(itpdt),nchi(itpdt),fhour,fno,ncho)
     else
       fno=fni(itpdt)
       ncho=nchi(itpdt)
     endif
!
#ifdef ASSIGN
     call assign('assign -R')
#endif
     close(n)
     open(unit=n,file=fno(1:ncho),form='unformatted',err=900)
     go to 901
900  continue
     write(6,*) ' error in opening file at file_write_sigma ',fno(1:ncho)
     call MPABORT
901  continue
#ifndef NOPRINT
     write(6,*) ' file ',fno(1:ncho),' opened. unit=',n
#endif
   endif
!
#ifdef MP
   allocate (coef1(lnt22_,levs_))
   allocate (coef2(lnt22_,levs_))
   allocate (coef3(lnt22_,levh_))
#endif
!
   if( iope ) then
     rewind(n)
     write(n) lab
#ifndef NOPRINT
     write(6,'(1x,A,a32,A,i3)')' file_write_sigma lab ',lab,' n= ',n
#endif
     do k = 1,kdum
       dummy(k)=0.
     enddo
     waves=jcap_
     xlayers=levs_
     trun=1.
     order=2.
     realform=1.
     gencode=igen
     rlond=lonf_
     rlatd=latg_
     rlonp=lonf_
     rlatp=latg_
     rlonr=lonf_
     rlatr=latg_
     water=nwater_
     gases=ngases_
     pdryini=1.
     subcen=icen2
     ensemble(1)=ienst
     ensemble(2)=iensi
     ppid=0.
     slid=0.
     vcid=0.
     vmid=0.
     vtid=0.
     do k = 1,kdum2
       dummy2(k)=0.
     enddo
#ifdef HYBRID
     do k = 1,levp1_
       ak5x(k)=ak5(levp1_+1-k)*1000. ! cb -> Pa
       bk5x(k)=bk5(levp1_+1-k)
     enddo
#endif
     if(ngases_.eq.0.and.nwater_.eq.1) then
#ifndef HYBRID
       write(n)fhour,idate,si,sl                                               &
#else
       write(n)fhour,idate,ak5x,bk5x                                           &
#endif
                 ,dummy,waves,xlayers,trun,order,realform,gencode              &
                 ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,water                    &
                 ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid         &
                 ,dummy2
     else
#ifndef HYBRID
       write(n)fhour,idate,si,sl                                               &
#else
       write(n)fhour,idate,ak5x,bk5x                                           &
#endif
          ,dummy,waves,xlayers,trun,order,realform,gencode                     &
          ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,water                           &
          ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid                &
          ,pdryini,dummy2,gases
     endif
#ifndef NOPRINT
     write(6,*) ' file_write_sigma fhour idate ... '
#endif
   endif
!
! gz
!
#ifdef DFS
#ifdef MP
   call mpm2f(gz,mt,spec,mtg,jlg)
!
#define OROG spec
#else
#define OROG gz
#endif
#endif /* DFS end */
   if (iope) then
#ifdef DCMIP
#ifndef DFS
     OROG=gphf
#endif
#endif /* DCMIP end */
     write(n)( OROG(i),i=1,lnt2_)
#ifndef NOPRINT
     write(6,*) ' file_write_sigma OROG '
#endif
   endif
!
!  q
!
#ifdef DFS
#ifdef MP
   call mpm2f(q,mt,spec,mtg,jlg)
!
   if (iope) then
     write(n)(spec(i),i=1,lnt2_),tavexy(1:levs_)
#ifndef NOPRINT
     write(6,*) ' file_write_sigma q '
#endif
   endif
!
#else
   write(n)( q(i),i=1,lnt2_),tavexy(1:levs_)
#endif
#else /* SPH */
#ifdef MP
   call mpsp2f(q,lnt22p_,spec,lnt22_,1)
!
   if( iope ) then
     write(n)(spec(i),i=1,lnt2_)
#ifndef NOPRINT
     write(6,*) ' file_write_sigma q '
#endif
   endif
#else
   write(n)( q(i),i=1,lnt2_)
#endif
#endif /* DFS end */
!
!  te
!
#ifdef DFS
#ifdef MP
   call mpmz2f(te,mt,jlg,levsp,coef1,mtg,jlg,levs,levs)
!
   if( iope ) then
     do k = 1,levs_
       write(n)(coef1(i,k),i=1,lnt2_)
     enddo
#ifndef NOPRINT
     write(6,*) ' file_write_sigma te '
#endif
   endif
#else
   do k = 1,levs_
     write(n)(te(i,k),i=1,lnt2_)
   enddo
#endif
#else /* SPH */
#ifdef MP
   call mpsp2f(te,lnt22p_,coef1,lnt22_,levs_)
!
   if( iope ) then
     do k = 1,levs_
       write(n)(coef1(i,k),i=1,lnt2_)
     enddo
#ifndef NOPRINT
     write(6,*) ' file_write_sigma te '
#endif
   endif
#else
   do k = 1,levs_
     write(n)(te(i,k),i=1,lnt2_)
   enddo
#endif
#endif /* DFS end */
!
!  di ze
!
#ifdef DFS
#ifdef MP
   call mpmz2f(di,mt,jlg,levsp,coef1,mtg,jlg,levs,levs)
   call mpmz2f(ze,mt,jlg,levsp,coef2,mtg,jlg,levs,levs)
!
   if( iope ) then
     do k = 1,levs_
       write(n)(coef1(i,k),i=1,lnt2_)
       write(n)(coef2(i,k),i=1,lnt2_)
     enddo
#ifndef NOPRINT
     write(6,*) ' file_write_sigma di ze '
#endif
   endif
#else
   do k = 1,levs_
     write(n)(di(i,k),i=1,lnt2_)
     write(n)(ze(i,k),i=1,lnt2_)
   enddo
#endif
#else /* SPH */
#ifdef MP
   call mpsp2f(di,lnt22p_,coef1,lnt22_,levs_)
   call mpsp2f(ze,lnt22p_,coef2,lnt22_,levs_)
!
   if( iope ) then
     do k = 1,levs_
       write(n)(coef1(i,k),i=1,lnt2_)
       write(n)(coef2(i,k),i=1,lnt2_)
     enddo
#ifndef NOPRINT
     write(6,*) ' file_write_sigma di ze '
#endif
   endif
#else
   do k = 1,levs_
     write(n)(di(i,k),i=1,lnt2_)
     write(n)(ze(i,k),i=1,lnt2_)
   enddo
#endif /* MP end */
#endif /* DFS end */
!
!  rq
!
#ifdef MP
#ifdef DFS
!  call mpmz2f(rq,mt,jlg,levhp,coef3,mtg,jlg,levs,levh)
!
   do nt = 1,ntotal_
     kt1 = (nt-1)*levsp+1
     kt2 = (nt-1)*levs+1
     call mpmz2f(rq(1,kt1),mt,jlg,levsp,coef3(1,kt2),mtg,jlg,levs,levs)
   enddo
#else
!  call mpsp2f(rq,lnt22p_,coef3,lnt22_,levh_)
!
   do nt = 1,ntotal_
     kt = (nt-1)*levs_+1
     call mpsp2f(rq(1,kt),lnt22p_,coef3(1,kt),lnt22_,levs_)
   enddo
#endif
#endif /* MP end */
!
   do nt = 1,ntotal_
     kt = (nt-1)*levs_+1
#ifdef DFS
#ifdef MP
     if( iope ) then
       do k = 1,levs_
         do i = 1,lnt2_
           coef1(i,k)=coef3(i,kt+k-1)
         enddo
       enddo
       do k = 1,levs_
         write(n)(coef1(i,k),i=1,lnt2_)
       enddo
#ifndef NOPRINT
       write(6,*) ' file_write_sigma rq '
#endif
     endif
#else
     do k = 1,levs_
       write(n)(rq(i,kt+k-1),i=1,lnt2_)
     enddo
#endif /* MP end */
#else /* SPH */
#ifdef MP
     if( iope ) then
       do k = 1,levs_
         do i = 1,lnt2_
           coef1(i,k)=coef3(i,kt+k-1)
         enddo
       enddo
       do k = 1,levs_
         write(n)(coef1(i,k),i=1,lnt2_)
       enddo
#ifndef NOPRINT
       write(6,*) ' file_write_sigma rq '
#endif
     endif
#else
     do k = 1,levs_
       write(n)(rq(i,kt+k-1),i=1,lnt2_)
     enddo
#endif /* MP end */
#endif /* DFS end */
   enddo
#ifdef SMP
!
! vvel
!
   do k = 1,levs_
     write(n)(vvel(i,k),i=1,lnt2_)
   enddo
#endif
!
#ifndef NOPRINT
#ifdef MP
   if( iope ) then
#endif
     close(n)
     write(6,3001)fhour,idate,n
     write(6,*)'number of water = ',water
     write(6,*)'number of gases = ',gases
3001  format(1x,'file_write_sigma fhour=',f10.1,2x,4i6,2x,'n=',i2)
#ifdef MP
   endif
#endif
#endif
!
#ifdef MP
   deallocate (coef1)
   deallocate (coef2)
   deallocate (coef3)
#endif
!
   return
   end subroutine file_write_sigma
!-------------------------------------------------------------------------------
