#include <define.h>
   subroutine post_pgb_3d
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [post_pgb_3d]
!           |
!           |-- [post_get_rh] *
!           |-- [post_get_omega] *
!           |-- [post_get_height] *
!
!-------------------------------------------------------------------------------
   end subroutine post_pgb_3d
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine post_get_rh(im,ix,km,sl,ps,q,t,qs,r)
!-------------------------------------------------------------------------------
!
! abstract: 
!   calculates relative humidity as a function of pressure,
!   specific humidity and temperature on the sigma layers.
!   saturation specific humidity is calculated from saturation vapor
!   pressure which is returned from a lookup table routine fpvs.
!
! program history log:
!   1992-10-31  iredell                development
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call post_get_rh(im,ix,km,sl,ps,q,t,qs,r)
!
!   input argument list:
!     im       - integer number of points
!     ix       - integer first dimension of upper air data
!     km       - integer number of levels
!     sl       - real (km) sigma values
!     ps       - real (im) surface pressure in kpa
!     q        - real (ix,km) specific humidity in kg/kg
!     t        - real (ix,km) temperature in k
!
!   output argument list:
!     qs       - real (im,km) saturated specific humidity in kg/kg
!     r        - real (im,km) relative humidity in percent
!
! subprograms called:
!   (fpvs)   - function to compute saturation vapor pressure
!
!-------------------------------------------------------------------------------
   use constant, only : rd_,rv_
!-------------------------------------------------------------------------------
   real            ::  sl(im,km),ps(im)
   real            ::  qs(im,km),r(im,km),q(ix,km),t(ix,km)
   real,parameter  ::  eps=rd_/rv_,epsm1=rd_/rv_-1.
!-------------------------------------------------------------------------------
   do k = 1,km
     do i = 1,im
#ifdef ICE
       es=fpvs(t(i,k))
#else
       es=fpvs0(t(i,k))
#endif
       qs(i,k)=eps*es/(sl(i,k)*ps(i)+epsm1*es)
       r(i,k)=min(max(q(i,k)/qs(i,k),0.),1.)*100.
     enddo
   enddo
!
   return
   end subroutine post_get_rh
!
!-------------------------------------------------------------------------------
   subroutine post_get_omega(im,ix,km,si,sl,ps,psx,psy,d,u,v,o,os)
!-------------------------------------------------------------------------------
!
! abstract: 
!   calculates pressure vertical velocity omega as a function
!   of surface pressure, surface pressure gradients, and divergence
!   and wind components on the sigma surfaces.  the formula for omega
!   is derived from the continuity equation
!     o=(sig*v.grad(lnps)-sum((d+v.grad(lnps))*dsig))*ps*1.e3
!   where the sum is taken from the top of the atmosphere.
!
! program history log:
!   1992-10-31  iredell
!
! usage:    call omega(im,ix,km,si,sl,ps,psx,psy,d,u,v,o,os)
!
!   input argument list:
!     im       - integer number of points
!     ix       - integer first dimension of upper air data
!     km       - integer number of levels
!     si       - real (km+1) sigma interface values
!     sl       - real (km) sigma values
!     ps       - real (im) surface pressure in kpa
!     psx      - real (im) zonal gradient of log pressure in 1/m
!     psy      - real (im) merid gradient of log pressure in 1/m
!     d        - real (ix,km) divergence in 1/s
!     u        - real (ix,km) zonal wind in m/s
!     v        - real (ix,km) merid wind in m/s
!
!   output argument list:
!     o        - real (im,km) pressure vertical velocity in pa/s
!     os       - real (im) surface pressure tendency in pa/s
!
!-------------------------------------------------------------------------------
   real     ::  si(im,km+1),sl(im,km)
   real     ::  ps(im),psx(im),psy(im)
   real     ::  d(ix,km),u(ix,km),v(ix,km)
   real     ::  o(im,km),os(im)
   real     ::  sum(im)
! 
   do i = 1,im
     sum(i)=0.
   enddo
!
   do k = km,1,-1
     do i = 1,im
       vgradp=u(i,k)*psx(i)+v(i,k)*psy(i)
       gradpv=vgradp+d(i,k)
       sum(i)=sum(i)+gradpv*(sl(i,k)-si(i,k+1))
       o(i,k)=(vgradp*sl(i,k)-sum(i))*ps(i)*1.e3
       sum(i)=sum(i)+gradpv*(si(i,k)-sl(i,k))
     enddo
   enddo
!
   do i = 1,im
     os(i)=-sum(i)*ps(i)*1.e3
   enddo
!
   return
   end subroutine post_get_omega
!
!-------------------------------------------------------------------------------
   subroutine post_get_height(im,ijm,km,si,sl,zs,t,q,z,zi)
!-------------------------------------------------------------------------------
!
! abstract: 
!   calculates geopotential heights on both the sigma interfaces
!   and the sigma full levels as a function of orography, temperature
!   and moisture.  virtual temperature is calculated from temperature
!   and moisture and the hydrostatic equation is integrated
!     dz=rd/g*tv*dlnp
!
! program history log:
!   1992-10-31  iredell
!
! usage:    call post_get_height(im,ijm,km,si,sl,zs,t,q,z,zi)
!
!   input argument list:
!     im       - integer number of points
!     ijm       - integer first dimension of upper air data
!     km       - integer number of levels
!     si       - real (km+1) sigma interface values
!     sl       - real (km) sigma values
!     zs       - real (im) orography is m
!     t        - real (ijm,km) temperature in k
!     q        - real (ijm,km) specific humidity in kg/kg
!
!   output argument list:
!     z        - real (ijm,km) heights on the full levels in m
!     zi       - real (ijm,km) heights on the interfaces in m
!
!-------------------------------------------------------------------------------
   use constant, only : g_,rd_,rv_
!-------------------------------------------------------------------------------
   integer              ::  km,ijm,im
   real                 ::  si(ijm,km+1),sl(ijm,km),zs(im),t(ijm,km),q(ijm,km)
   real                 ::  z(ijm,km),zi(ijm,km)
   real,parameter       ::  rog=rd_/g_,fvirt=rv_/rd_-1.
!-------------------------------------------------------------------------------
   do i = 1,im
     zi(i,1)=zs(i)
   enddo
!
#ifdef DBG
   print*, 'im,km,ijm',im,km,ijm
#endif
!
   do k = 1,km-1
     do i = 1,im
       ca=rog*log(si(i,k)/sl(i,k))
       cb=rog*log(sl(i,k)/si(i,k+1))
       tv=t(i,k)*(1.+fvirt*q(i,k))
       z(i,k)=zi(i,k)+ca*tv
       zi(i,k+1)=z(i,k)+cb*tv
     enddo
   enddo
!
   do i = 1,im
     ca=rog*log(si(i,km)/sl(i,km))
     tv=t(i,km)*(1.+fvirt*q(i,km))
     z(i,km)=zi(i,km)+ca*tv
   enddo
!
   return
   end subroutine post_get_height
!-------------------------------------------------------------------------------
