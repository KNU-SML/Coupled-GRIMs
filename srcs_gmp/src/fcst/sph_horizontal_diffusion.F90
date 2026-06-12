#include "define.h"
   subroutine sph_horizontal_diffusion(rt,w,deltim,qm,sl,x,y,lnt2s,lnoffset)
#ifdef NCEP_GFS_DIFFUSION
#define GFS_FACTOR 1.0
#else
#define GFS_FACTOR 0.1
#endif
!-------------------------------------------------------------------------------
!
! subroutine:    sph_horizontal_diffusion   
!
! abstract: horizontal diffusion of temperature, moisture,
!           vorticity and divergence. the implicit linear equation
!           is solved using the lapse rates of globally averaged
!           temperature and moisture to transform the laplacian
!           from constant pressure to constant sigma surfaces.
!           for the t126 operational model, sph_horizontal_diffusion invokes
!           second order leith diffusion only above wavenumber 69
!           with a time scale of 13080 seconds at wavenumber 126.
!           for higher resolutions, full diffusion for all wave number is applied.
!
! program history log:
!   1991-03-15  mark iredell           development
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!   2011-05-01  song-you hong          modifications for high resolutions as in gfs
!
! usage:    call sph_horizontal_diffusion (rt,w,deltim,qm,sl,x,y)
!   input argument list:
!     rt       - specific humidity
!     w        - vorticity
!     deltim   - timestep
!     qm       - ln(psfc)
!     sl       - sigma layer values
!     x        - divergence
!     y        - temperature
!
!   output argument list:
!     rt       - specific humidity
!     w        - vorticity
!     x        - divergence
!     y        - temperature
!
! remarks: local variables that can be modified to change diffusion are:
!          rtnp- reciprocal of time scale of diffusion at wavenumber np
!          np  - wavenumber at which rtnp diffusion is applied
!          n0  - maximum wavenumber for zero diffusion
!          jdel- order of diffusion
!
!-------------------------------------------------------------------------------
   use paramodel, only : ngases_,levs_,levh_,latg_,jcap_,LNT22S,lnt2_,ntotal_
   use constant, only  : rerth_
#ifdef DFS
   use dfsvar, only : iope
#endif
   use comio
#ifdef MP
   use commpi
#endif
!-------------------------------------------------------------------------------
   real                 ::  sl(levs_),qm(LNT22S)
   real                 ::  w(LNT22S,levs_),x(LNT22S,levs_)
   real                 ::  y(LNT22S,levs_),rt(LNT22S,levh_)
!
! local
!
   real,save,allocatable  ::  dn(:),rtrd(:),rthk(:)
   integer,save         ::  ifirst
   data ifirst/0/
!
   real,parameter       ::  difcof=3.e15
   integer,parameter    ::  lefres=80
   integer,save         ::  n0,jdel
   real,save            ::  fshk,slrd0,rtrd1
   real                 ::  fact,rfact,rfactrd
   integer              ::  kk,n
!
   if(.not.allocated(dn)) allocate(dn(lnt2_))
#ifdef SPH_GFS_DIFFUSION
   if(.not.allocated(rtrd)) allocate(rtrd(levs_))
   if(.not.allocated(rthk)) allocate(rthk(levs_))
#endif
   if(ifirst.eq.0) then
     ifirst=1
#ifdef SPH_GFS_DIFFUSION
     if(jcap_.gt.170) then
!
!  maximum wavenumber for zero diffusion
!
       n0=0
!
!  reciprocal of time scale of diffusion at reference wavenumber np
!
       rtnp=(jcap_/170.)**4*1.1/3600.*GFS_FACTOR
       np=jcap_
!
!  order of diffusion (even power to raise del)
!
       jdel=8
!
!  extra height-dependent diffusion factor per scale height
!
       fshk=2.2
     else
#endif
       n0=0.55*jcap_
       rtnp=difcof/(rerth_**4)*float(lefres*(lefres+1))**2
       np=jcap_
       jdel=2
       fshk=1.0
#ifdef SPH_GFS_DIFFUSION
     endif
#endif
!
#ifndef NOPRINT
     if( iope ) write(6,6)rtnp,np,n0,jdel
6      format(' horizontal diffusion parameters'/                              &
     '   effective ',6pf10.3,' microhertz at wavenumber ',i4/                  &
     '   maximum wavenumber for zero diffusion ',i4/                           &
     '   order of diffusion ',i2)
#endif
#ifdef SPH_GFS_DIFFUSION
!
!  sigma level at which to begin rayleigh damping
!
     slrd0=0.002
!
!  reciprocal of time scale per scale height
!  above beginning sigma level for rayleigh damping
!
     rtrd1=1./(5*86400)
!
     do k = 1,levs_
       if(sl(k).lt.slrd0) then
         rtrd(k)=rtrd1*log(slrd0/sl(k))
       else
         rtrd(k)=0
       endif
       rthk(k)=(sl(k))**log(1/fshk)
     enddo
#endif
!
     jdelh=jdel/2
     npd=max(np-n0,0)
     dn1=2.*rtnp/float(npd*(npd+1))**jdelh
     i=0
     do nm = 0,jcap_
       do m = 0,jcap_-nm
         nd=max(nm+m-n0,0)
         dn(i+1)=dn1*float(nd*(nd+1))**jdelh
         dn(i+2)=dn(i+1)
         i=i+2
       enddo
     enddo
#ifdef MP
     call spcshfli(dn,lnt2_,1,jcap_,lwvdef)
#endif
#ifdef DBG
     if(iope) write(6,*) 'dn computation completed'
#endif
   endif
#ifdef SPH_GFS_DIFFUSION
!
   istr=1
   if(lnoffset.eq.0) istr=3
   do k = 1,levs_
     do i = istr,lnt2s
       ii=i+lnoffset
       w(i,k)= w(i,k)/(1.+deltim*dn(ii)*rthk(k)+deltim*rtrd(k))
       x(i,k)= x(i,k)/(1.+deltim*dn(ii)*rthk(k)+deltim*rtrd(k))
       y(i,k)= y(i,k)/(1.+deltim*dn(ii)*rthk(k))
     enddo
   enddo
#ifndef NISLQ
   do n = 1,ntotal_
     do k = 1,levs_
       kk=(n-1)*levs_+k 
       do i = istr,lnt2s
         ii=i+lnoffset
         rt(i,kk)=rt(i,kk)/(1.+deltim*dn(ii)*rthk(k))
       enddo
     enddo
   enddo
#endif
#else /* ~SPH_GFS_DIFFUSION */
!
   istr=1
   if(lnoffset.eq.0) istr=3
   do k = 1,levs_
     do i = istr,lnt2s
       ii=i+lnoffset
       w(i,k)= w(i,k)/(1.+deltim*dn(ii))
       x(i,k)= x(i,k)/(1.+deltim*dn(ii))
       y(i,k)= y(i,k)/(1.+deltim*dn(ii))
      enddo
   enddo
#ifndef NISLQ
   do k = 1,levh_
     do i = istr,lnt2s
       ii=i+lnoffset
       rt(i,k)=rt(i,k)/(1.+deltim*dn(ii))
     enddo
   enddo
#endif
#endif /* SPH_GFS_DIFFUSION end */
#ifdef DBG
   if (iope) write(6,*) 'diffusion computation completed'
#endif
!
   return
   end subroutine sph_horizontal_diffusion
