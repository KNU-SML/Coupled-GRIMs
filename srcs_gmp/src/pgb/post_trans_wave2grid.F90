#include <define.h>
   subroutine post_trans_wave2grid
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [post_trans_wave2grid]
!           |
!           |-- [post_trans_wave_grid] *
!           |        |-- [post_gaussian_lat] *
!           |-- [post_wave2grid] *
!
!-------------------------------------------------------------------------------
   end subroutine post_trans_wave2grid
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_trans_wave_grid(idir,grid,wave,mlt,fac,                     &
                                   imax,jmax,maxwv,iromb)
!-------------------------------------------------------------------------------
!
! subprogram:  post_trans_wave_grid     spherical transform
!
! abstract: transforms a field between grid and spectral domains.
!   this versatile routine will:
!     ...transform grid to spectral or transform spectral to grid;
!     ...pass a gaussian grid or an equally-space grid;
!     ...pass a triangular truncation or a rhomboidal truncation;
!     ...optionally transform with derivatives of legendre functions;
!     ...optionally transform from grid dividing by coslat**2;
!     ...optionally multiply spectral field by complex factors.
!
! usage: call post_trans_wave_grid(idir,grid,wave,mlt,fac,imax,jmax,maxwv,iromb)
!
!   input argument list:
!     idir     - must be one of 1,2,3,4,101,102,103,104,-1,-2,-101,-102,
!                idir> 0 to transform grid to wave,
!                idir < 0 to transform wave to grid,
!                abs(idir) < 100 for gaussian grid,
!                abs(idir)> 100 for equally-spaced grid,
!                idir is odd for normal transform,
!                idir is even for legendre derivative transform,
!                idir last digit> 2 to transform with cosine latitude.
!     grid     - if idir> 0, real (imax,jmax) field to transform.
!                grid starts at north pole and greenwich meridian.
!     wave     - if idir < 0, complex (kmax) field to transform,
!                where kmax=(maxwv+1)*(iromb+1)*maxwv+2)/2.
!                wave starts at the global mean component and then
!                contains the zonally symmetric components.
!     mlt      - multiplication option on wave values.
!                mlt = 0 for no multiplication,
!                mlt = 1 to multiply wave by fac
!                mlt = -1 to multiply wave by sqrt(-1)*fac
!     fac      - if mlt.ne.0, real (kmax) field to multiply wave.
!     imax     - longitudinal dimension of the grid
!     jmax     - latitudinal dimension of the grid
!     maxwv    - spectral truncation of the wave
!     iromb    - iromb = 0 for triangular truncation
!                iromb = 1 for rhomboidal truncation
!
!   output argument list:
!     grid     - if idir < 0, real (imax,jmax) field output.
!                grid starts at north pole and greenwich meridian.
!     wave     - if idir> 0, complex (kmax) field output,
!                where kmax=(maxwv+1)*(iromb+1)*maxwv+2)/2.
!                wave starts at the global mean component and then
!                contains the zonally symmetric components.
!
!   subprograms called:
!     unique:
!     post_gaussian_lat   - compute gaussian latitudes
!     post_equall_colat   - compute equally-spaced latitudes
!     (post_legendre)     - compute legendre polynomials
!     fftfax     - fft (library call can be substituted)
!     rfftmlt    - fft (library call can be substituted)
!
!  ::: structure ::: This file contains ... 
!
!     [post_trans_wave_grid] * ----- [post_gaussian_lat] *
!
!-------------------------------------------------------------------------------
   use paramodel, only : jcap_,latg_,lonf_, nepx=>jcap1_, ntrigs=>lonf2_, &
                         npnm=>lnt_
!-------------------------------------------------------------------------------
#undef MP
   save
   real                 ::  grid(imax,jmax)
   complex              ::  wave((maxwv+1)*((iromb+1)*maxwv+2)/2)
   real                 ::  fac((maxwv+1)*((iromb+1)*maxwv+2)/2)
!
   integer              ::  ngg
   real                 ::  gg(2*(1+(lonf_+1)/2),latg_)
   complex              ::  ww(npnm),w2(-1:1),ws
   real                 ::  trigs(ntrigs)
   integer              ::  ifax(100),indir
!
   real                 ::  cosclt(latg_),wgtclt(latg_)
   real                 ::  wrkfft(2*lonf_*latg_)
   integer              ::  is(npnm)
   real                 ::  pnm(npnm)
   real                 ::  ex(0:nepx),px(-1:nepx)
!-------------------------------------------------------------------------------
!
   if (idir.gt.0) then
     indir=1
   else
     indir=-1
   endif
!
   kmax=(maxwv+1)*((iromb+1)*maxwv+2)/2
   ipd=1-mod(abs(indir),2)
   isd=1-2*ipd
   icd=(mod(abs(indir),10)-1)/2
   jump=2*(1+(imax+1)/2)
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
   call fftfax (lonf_,ifax,trigs)
#endif
#ifdef DEFAULT
   call fax (ifax, lonf_,3)
   call fftrig (trigs,lonf_,3)
#endif
!
   if(abs(idir).lt.100) then
     call post_gaussian_lat(jmax,cosclt,wgtclt)
   else
#ifdef DFS
     call post_dfs_lat(jmax,cosclt,wgtclt)
#else
     call post_equall_colat(jmax,cosclt,wgtclt)
#endif
   endif
!
   k=0
   do m = 0,maxwv
     do n = m,iromb*m+maxwv
       k=k+1
       is(k)=isd*(1-2*mod(n-m,2))
     enddo
   enddo
!
   if(indir.gt.0) then
!
     do j = 1,jmax
       do i = 1,imax
         gg(i,j)=grid(i,j)
       enddo
     enddo
!
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
     call rfftmlt(gg,wrkfft,trigs,ifax,1,jump,imax,jmax,-1)
#endif
#ifdef DEFAULT
     call fft99m (gg,wrkfft,trigs,ifax,1,jump,imax,jmax,-1)
#endif
!
     do k = 1,kmax
       wave(k)=0.
     enddo
!
     do j = 1,(jmax+1)/2
       jr=jmax+1-j
       call post_legendre(ipd,cosclt(j),maxwv,iromb,ex,px,pnm)
       wj=wgtclt(j)
       if(icd.ne.0.and.cosclt(j).lt.1.) wj=wj/(1.-cosclt(j)**2)
       k=0
       do m = 0,maxwv
         w2(1)=wj*cmplx(gg(2*m+1,j)+gg(2*m+1,jr),                              &
                        gg(2*m+2,j)+gg(2*m+2,jr))
         w2(-1)=wj*cmplx(gg(2*m+1,j)-gg(2*m+1,jr),                             &
                         gg(2*m+2,j)-gg(2*m+2,jr))
         do n = m,iromb*m+maxwv
           k=k+1
           wave(k)=wave(k)+w2(is(k))*pnm(k)
         enddo
       enddo
     enddo
!
     if(mlt.ne.0) then
       ws=cmplx(1.,0.)
       if(mlt.lt.0) ws=cmplx(0.,1.)
       do k = 1,kmax
         wave(k)=wave(k)*ws*fac(k)
       enddo
     endif
!
   else
!
     if(mlt.ne.0) then
       ws=cmplx(1.,0.)
       if(mlt.lt.0) ws=cmplx(0.,1.)
       do k = 1,kmax
         ww(k)=wave(k)*ws*fac(k)
       enddo
     else
       do k = 1,kmax
         ww(k)=wave(k)
       enddo
     endif
!
     do j = 1,(jmax+1)/2
       call post_legendre(ipd,cosclt(j),maxwv,iromb,ex,px,pnm)
       k=0
       do m = 0,maxwv
         w2(1)=0.
         w2(-1)=0.
         do n = m,iromb*m+maxwv
           k=k+1
           w2(1)=w2(1)+ww(k)*pnm(k)
           w2(-1)=w2(-1)+ww(k)*pnm(k)*is(k)
         enddo
         gg(2*m+1,jmax+1-j)=real(w2(-1))
         gg(2*m+2,jmax+1-j)=aimag(w2(-1))
         gg(2*m+1,j)=real(w2(1))
         gg(2*m+2,j)=aimag(w2(1))
       enddo
     enddo
!
     do i = 2*maxwv+3,jump
       do j = 1,jmax
         gg(i,j)=0.
       enddo
     enddo
!
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
     call rfftmlt(gg,wrkfft,trigs,ifax,1,jump,imax,jmax,1)
#endif
#ifdef DEFAULT
     call fft99m (gg,wrkfft,trigs,ifax,1,jump,imax,jmax,1)
#endif
!
     do j = 1,jmax
       do i = 1,imax
         grid(i,j)=gg(i,j)
       enddo
     enddo
   endif
!
   return
   end subroutine post_trans_wave_grid
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_gaussian_lat(k,a,w)
!-------------------------------------------------------------------------------
   use paramodel, only : npk=>latg2_
!-------------------------------------------------------------------------------
!
! subprogram:  post_gaussian_lat   compute gaussian latitudes
!
! abstract: computes cosines of colatitude and gaussian weights
!   on the gaussian latitudes.  the k gaussian latitudes are at
!   the zeroes of the legendre polynomial of order k.
!
! program history log:
!   1992-04-16  iredell
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!
! usage:    call post_gaussian_lat(k,a,w)
!
!   input argument list:
!     k        - number of latitudes.
!
!   output argument list:
!     a        - real (k) cosines of colatitude.
!     w        - real (k) gaussian weights.
!
!   remarks: fortran 9x extensions are used.
!
!-------------------------------------------------------------------------------
#undef MP
   save
   real                 ::  a(k),w(k)
   real                 ::  pk(npk),pkm1(npk),pkm2(npk)
   integer,parameter    ::  kz=50
   real                 ::  bz(kz)
   data bz        / 2.4048255577,  5.5200781103,                               &
     8.6537279129, 11.7915344391, 14.9309177086, 18.0710639679,                &
    21.2116366299, 24.3524715308, 27.4934791320, 30.6346064684,                &
    33.7758202136, 36.9170983537, 40.0584257646, 43.1997917132,                &
    46.3411883717, 49.4826098974, 52.6240518411, 55.7655107550,                &
    58.9069839261, 62.0484691902, 65.1899648002, 68.3314693299,                &
    71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711,                &
    84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819,                &
    96.6052679510, 99.7468198587, 102.888374254, 106.029930916,                &
    109.171489649, 112.313050280, 115.454612653, 118.596176630,                &
    121.737742088, 124.879308913, 128.020877005, 131.162446275,                &
    134.304016638, 137.445588020, 140.587160352, 143.728733573,                &
    146.870307625, 150.011882457, 153.153458019, 156.295034268 /
!-------------------------------------------------------------------------------
   pi=4.*atan(1.)
   c=(1.-(2./pi)**2)*0.25
   eps=1.e-14
   kh=k/2
   r=1./sqrt((k+0.5)**2+c)
!
   do j = 1,min(kh,kz)
     a(j)=cos(bz(j)*r)
   enddo
   do j = kz+1,kh
     a(j)=cos((bz(kz)+(j-kz)*pi)*r)
   enddo
!
   spmax=1.
   kount=0
!
   do while(spmax.gt.eps)
     kount=kount+1
     spmax=0.
     do j = 1,kh
       pkm1(j)=1.
       pk(j)=a(j)
     enddo
     do n = 2,k
       do j = 1,kh
         pkm2(j)=pkm1(j)
         pkm1(j)=pk(j)
         pk(j)=((2*n-1)*a(j)*pkm1(j)-(n-1)*pkm2(j))/n
       enddo
     enddo
     do j = 1,kh
       sp=pk(j)*(1.-a(j)**2)/(k*(pkm1(j)-a(j)*pk(j)))
       a(j)=a(j)-sp
       spmax=max(spmax,abs(sp))
     enddo
     if(kount.gt.20) then
       print *,'chgr_gaussian_lat not converging'
       call abort
     endif
   enddo
!
   do j = 1,kh
     w(j)=(2.*(1.-a(j)**2))/(k*pkm1(j))**2
     a(k+1-j)=a(j)
     w(k+1-j)=w(j)
   enddo
!
   if(k.ne.kh*2) then
     j=kh+1
     a(j)=0.
     w(j)=1./k**2
     do n = 2,k,2
       w(j)=w(j)*n**2/(n-1)**2
     enddo
   endif
!
   return
   end subroutine post_gaussian_lat
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
#ifdef SMP
   subroutine post_wave2grid(nflds,igases,iwater,ss,f                         ,&
#else
   subroutine post_wave2grid(nflds,igases,iwater,trig,ifax                    ,&
                             eps,epstop,ss,sstop,coslat,sinlat,f              ,&
#endif
                             io2,io22,johf)
!-------------------------------------------------------------------------------
!
! subprogram:    post_wave2grid        transform spectral to grid
!
! abstract: transforms spectral to gridded data on a latitude pair
!           and computes dry temperature and surface pressure
!           and winds and gradients without a cosine latitude factor.
!           subprogram post_read_coeff should be called already
!           to read spectral data and initialize utility fields.
!           this subprogram can be called from a multiprocessed segment.
!
! program history log:
!   91-10-31  mark iredell
!
! usage:    call post_wave2grid(jcap,nc,nctop,levs,trig,ifax,eps,epstop,ss,sstop,
!    &                io2,io22,coslat,sinlat,f)
!
!   input argument list:
!     jcap     - integer spectral truncation
!     nc       - integer number of spectral coefficients
!     nctop    - integer number of spectral coefficients over top
!     levs     - integer number of levels
!     trig     - real (io2) trigonometric quantities for the fft
!     ifax     - integer (100) factors for the fft
!     eps      - real ((jcap+1)*(jcap+2)/2) sqrt((n**2-l**2)/(4*n**2-1))
!     epstop   - real (jcap+1) sqrt((n**2-l**2)/(4*n**2-1)) over top
!     ss       - real (nc,6*levs+6) spectral coefs
!     sstop    - real (nctop,6*levs+6) spectral coefs over top
!     io2    - integer number of valid data points per latitude pair
!     io22   - integer longitude dimension of data (>=io2+4)
!     coslat   - real cosine of latitude of the latitude pair
!     sinlat   - real sine of latitude of the northern latitude
!
!   output argument list:
!     f        - real (io22,6*levs+6) gridded data
!                (:,1:levs)             vorticity
!                (:,levs+1:2*levs)      divergence
!                (:,2*levs+1:3*levs)    temperature
!                (:,3*levs+1:4*levs)    specific humidity
!                (:,4*levs+1)           d(lnps)/dx
!                (:,4*levs+2)           d(lnps)/dy
!                (:,4*levs+3:5*levs+2)  zonal wind
!                (:,5*levs+3:6*levs+2)  meridional wind
!                (:,6*levs+3)           surface pressure
!                (:,6*levs+4)           orography
!                (:,6*levs+5)           d(orog)/dx
!                (:,6*levs+6)           d(orog)/dy
!
! subprograms called:
!   post_sph_poly          compute associated legendre polynomials
!   post_synth_coeff1      synthesize fourier from spectral coefficients
!   rfftmlt      fast fourier transform
!
!-------------------------------------------------------------------------------
   use paramodel, only  :  jcap => jcap_ , levs  => levs_ , &
                           nc   => lnt22_, nctop => twoj1_
   use constant, only   :  rd_,rv_
!-------------------------------------------------------------------------------
   real, parameter  ::  fv= rv_ / rd_ -1.
   integer          ::  io2,io22,johf
#ifdef SMP
   real             ::  ss(nc,nflds),f(io22,nflds)
#else
   real             ::  trig(io2)
   integer          ::  ifax(100)
   real             ::  eps((jcap+1)*(jcap+2)/2),epstop(jcap+1)
   real             ::  ss(nc,nflds),sstop(nctop,nflds),f(io22,nflds)
   real             ::  pln((jcap+1)*(jcap+2)/2),plntop(jcap+1)
   real             ::  wfft(io22,2*(5*levs+igases*levs+6+iwater*levs))
   integer          ::  mp(5*levs+igases*levs+6+iwater*levs)
#endif
   logical          ::  first
   data first/.true./
!-------------------------------------------------------------------------------
#ifdef SMP
   print*,' post_wave2grid ',nflds,igases,iwater
!
   do n = 1,nflds
     do i = 1,io22
       f(i,n)=0.0    ! initialization
     enddo
   enddo
!
   do n = 1,nflds
     do i = 1,io2
       f(i,n)=ss(i,n)
     enddo
   enddo
#else
!
!  transform spectral coefficients to fourier coefficients
!
   if (first) then
     first=.false.
     print*,' post_wave2grid : nflds,iwater,igases,trig(1),ifax(1),coslat,sinlat '
     print*,nflds,iwater,igases,trig(1),ifax(1),coslat,sinlat
   endif
!
   call post_sph_poly(jcap,sinlat,coslat,eps,epstop,pln,plntop)
   call post_synth_coeff1(jcap,io22/2,nc,nctop,nflds,pln,plntop,ss,sstop,f)
!
!  transform fourier coefficients to gridded data
!
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
   call rfftmlt(f,wfft,trig,ifax,1,io22/2,io2/2,2*nflds,1)
#endif
#ifdef DEFAULT
   call fft99m (f,wfft,trig,ifax,1,io22/2,io2/2,2*nflds,1)
#endif
!
!  move southern hemisphere latitude after northern hemisphere latitude
!
   do k = 1,nflds
     do i = 1,io2/2
       f(io2/2+i,k)=f(io22/2+i,k)
     enddo
   enddo
#endif
! 
!  compute dry temperature from virtual temperature
!  and surface pressure from log surface pressure
!  and divide gradients and winds by cosine of latitude.
!
   do k = 1,levs
     do i = 1,io2
       f(i,4*levs+k+6)=f(i,4*levs+k+6)/(1.+fv*f(i,5*levs+k+6))
#ifndef SMP
       if(coslat.ne.0.) then
         f(i,2*levs+6+k)=f(i,2*levs+6+k)/coslat
         f(i,3*levs+6+k)=f(i,3*levs+6+k)/coslat
       else
         f(i,2*levs+6+k)=0.
         f(i,3*levs+6+k)=0.
       endif
#endif
     enddo
   enddo
!
   do i = 1,io2
     f(i,4)=exp(f(i,4))
   enddo
!
#ifndef SMP
   if(coslat.ne.0.) then
     do i = 1,io2
       f(i,5)=f(i,5)/coslat
       f(i,6)=f(i,6)/coslat
       f(i,2)=f(i,2)/coslat
       f(i,3)=f(i,3)/coslat
     enddo
   else
     do i = 1,io2
       f(i,5)=0.
       f(i,6)=0.
       f(i,2)=0.
       f(i,3)=0.
     enddo
   endif
#endif
#ifdef DBG
   call print_maxmin_seven(f(1,1),io2,io22,nflds,1,nflds,'after post_wave2grid')
#endif
!
   return
   end subroutine post_wave2grid
!-------------------------------------------------------------------------------
