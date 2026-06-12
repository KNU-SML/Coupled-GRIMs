#include <define.h>
   subroutine phys_moist_adiabat(ilev,klev,k1,k2,                              &
                                 prsl,prsik,prslk,tenv,qenv,                   &
                                 klcl,kbot,ktop,tcld,qcld,rd,rv,               &
                                      ids,ide, jds,jde, kds,kde,               &
                                      ims,ime, jms,jme, kms,kme,               &
                                      its,ite, jts,jte, kts,kte)
!-------------------------------------------------------------------------------
!
! subprogram: phys_moist_adiabat       
!
! abstract: 
! - compute moist adiabatic cloud soundings
! - atmospheric columns of temperature and specific humidity
!   are examined by this routine for conditional instability.
!   the test parcel is chosen from the layer between layers k1 and k2
!   that has the warmest potential wet-bulb temperature.
!   excess cloud temperatures and specific humidities are returned
!   where the lifted parcel is found to be buoyant.
!   fast inlinable functions are invoked to compute
!   dewpoint and lifting condensation level temperatures,
!   equivalent potential temperature at the lcl, and
!   temperature and specific humidity of the ascending parcel.
!
! program history log:
!   1983-11-01  phillips
!   1991-05-07  iredell                arguments changed, code tidied
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:  call phys_moist_adiabat(ilev,klev,k1,k2,                             &
!                                 prsl,prslk,prsik,tenv,qenv,                  &
!                                 klcl,kbot,ktop,tcld,qcld,rd,rv,              &
!                                      ids,ide, jds,jde, kds,kde,              &
!                                      ims,ime, jms,jme, kms,kme,              &
!                                      its,ite, jts,jte, kts,kte)
!
!   input argument list:
!     ilev         - integer number of atmospheric columns
!     klev         - integer number of sigma levels in a column
!     k1           - integer lowest level from which a parcel can originate
!     k2           - integer highest level from which a parcel can originate
!     prsl         - real (ilev,klev) pressure values
!     prslk,prsik  - real (ilev,klev) pressure values to the kappa
!     tenv         - real (ilev,klev) environment temperatures
!     qenv         - real (ilev,klev) environment specific humidities
!
!   output argument list:
!     klcl     - integer (ilev) level just above lcl (klev+1 if no lcl)
!     kbot     - integer (ilev) level just above cloud bottom
!     ktop     - integer (ilev) level just below cloud top
!              - note that kbot(i) gt ktop(i) if no cloud.
!     tcld     - real (ilev,klev) of excess cloud temperatures.
!                (parcel t minus environ t, or 0. where no cloud)
!     qcld     - real (ilev,klev) of excess cloud specific humidities.
!                (parcel q minus environ q, or 0. where no cloud)
!
! subprograms called:
!     ftdp     - function to compute dewpoint temperature
!     ftlcl    - function to compute lcl temperature
!     fthe     - function to compute equivalent potential temperature
!     ftma     - function to compute parcel temperature and humidity
!
! remarks: all functions are inlined by fpp.
!          nonstandard automatic arrays are used.
!
!--------------------------------------------------------------------------------
   integer              ::  ilev,klev,k1,k2,kl
   real                 ::  prsl(ilev,klev),prslk(ilev,klev),prsik(ilev,klev)
   real                 ::  tenv(ilev,klev),qenv(ilev,klev)
   integer              ::  klcl(ilev),kbot(ilev),ktop(ilev)
   real                 ::  tcld(ilev,klev),qcld(ilev,klev)
!
!  local arrays
!
   real                 ::  slkma(ilev)
   real                 ::  thema(ilev)
!
   integer,parameter    ::  nx=151,ny=121
   real                 ::  tbtma(nx,ny),tbqma(nx,ny)
   real                 ::  eps,epsm1,ftv
   common/comma/ c1xma,c2xma,c1yma,c2yma,tbtma,tbqma
!
!  compute parameters
!   
   eps=rd/rv
   epsm1=rd/rv-1.
   ftv=rv/rd-1.
!
!  determine warmest potential wet-bulb temperature between k1 and k2.
!  compute its lifting condensation level.
!
   do i = its,ilev
     slkma(i)=0.
     thema(i)=0.
   enddo
!
   do k = k1,k2
     do i = its,ilev
       pv=prsl(i,k)*qenv(i,k)/(eps-epsm1*qenv(i,k))
       tdpd=tenv(i,k)-ftdp(pv)
       if(tdpd.gt.0.) then
         tlcl=ftlcl(tenv(i,k),tdpd)
         slklcl=prslk(i,k)/prsik(i,1)*tlcl/tenv(i,k)
       else
         tlcl=tenv(i,k)
         slklcl=prslk(i,k)/prsik(i,1)
       endif
       thelcl=fthe(tlcl,slklcl*prsik(i,1))
       if(thelcl.gt.thema(i)) then
         slkma(i)=slklcl
         thema(i)=thelcl
       endif
     enddo
   enddo
!
!  set cloud temperatures and humidities wherever the parcel lifted up
!  the moist adiabat is buoyant with respect to the environment.
!
   do i = its,ilev
     klcl(i)=klev+1
     kbot(i)=klev+1
     ktop(i)=0
   enddo
!
   do k = kts,klev
     do i = its,ilev
       tcld(i,k)=0.
       qcld(i,k)=0.
     enddo
   enddo
!
   do k = k1,klev
     do i = its,ilev
       if(prslk(i,k)/prsik(i,1).le.slkma(i)) then
         klcl(i)=min(klcl(i),k)
!
! insert ftma
!        tma=ftma(thema(i),prslk(i,k),qma)
!
         xj=min(max(c1xma+c2xma*thema(i),1.),float(nx))
         yj=min(max(c1yma+c2yma*prslk(i,k),1.),float(ny))
         jx=min(xj,nx-1.)
         jy=min(yj,ny-1.)
         ftx1=tbtma(jx,jy)+(xj-jx)*                                            &
                     (tbtma(jx+1,jy)-tbtma(jx,jy))
         ftx2=tbtma(jx,jy+1)+(xj-jx)*                                          &
                     (tbtma(jx+1,jy+1)-tbtma(jx,jy+1))
         ftma1=ftx1+(yj-jy)*(ftx2-ftx1)
         qx1=tbqma(jx,jy)+(xj-jx)*(tbqma(jx+1,jy)-tbqma(jx,jy))
         qx2=tbqma(jx,jy+1)+(xj-jx)*(tbqma(jx+1,jy+1)-tbqma(jx,jy+1))
         qma=qx1+(yj-jy)*(qx2-qx1)
         tma=ftma1
!
         tvcld=tma*(1.+ftv*qma)
         tvenv=tenv(i,k)*(1.+ftv*qenv(i,k))
         if(tvcld.gt.tvenv) then
           kbot(i)=min(kbot(i),k)
           ktop(i)=max(ktop(i),k)
           tcld(i,k)=tma-tenv(i,k)
           qcld(i,k)=qma-qenv(i,k)
         endif
       endif
     enddo
   enddo
!
   return
   end subroutine phys_moist_adiabat
