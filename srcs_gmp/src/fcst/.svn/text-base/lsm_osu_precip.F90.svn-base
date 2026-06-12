#include <define.h>
   subroutine lsm_osu_precip(ims2,imx2,kmx,                                    &
          rhscnpy,rhsmc,ai,bi,ci,smc,slimsk,                                   &
          canopy,precip,runoff,snowmt,                                         &
#ifndef HYDRO
          zsoil,soiltyp,sigmaf,delt,lat)
#else
          zsoil,soiltyp,sigmaf,delt,lat,hydrow)
#endif
!-------------------------------------------------------------------------------
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only : LONF2S
   use varsfc, only : lsoil_
#ifdef DFS
   use dfsvar, only : iope
#endif
   use comio
!-------------------------------------------------------------------------------
#include "abort.h"
#ifndef HYDRO
   real,parameter       ::  scanop=2.0,rhoh2o=1000.
#else
   real,parameter       ::  scanop=.5,rhoh2o=1000.
#endif
   real,parameter       ::  ctfil1=.5,ctfil2=1.-ctfil1
   real,parameter       ::  rffact=.15
!
   integer              ::  ims2,imx2,kmx,lat
   real                 ::  rhscnpy(imx2),rhsmc(imx2,kmx)
   real                 ::  ai(imx2,kmx),bi(imx2,kmx),ci(imx2,kmx)
   real                 ::  smc(imx2,kmx),canopy(imx2),precip(imx2)
   real                 ::  sigmaf(imx2),slimsk(imx2),runoff(imx2)
   real                 ::  zsoil(imx2,kmx),snowmt(imx2)
   integer              ::  soiltyp(imx2)
!
! local array
!
   logical              ::  flag(LONF2S)
   real                 ::  osu_funct_thsat
   real                 ::  prcp(LONF2S),inf(LONF2S)
   real                 ::  infmax(LONF2S)
   real                 ::  tsat(LONF2S),dsat(LONF2S)
   real                 ::  ksat(LONF2S)
   real                 ::  smsoil(LONF2S,lsoil_)
   real                 ::  cnpy(LONF2S)
   real                 ::  dew(LONF2S)
#ifdef HYDRO
   real                 ::  hydrow(LONF2S,lsoil_)
#endif
!-------------------------------------------------------------------------------
   im = ims2
   km = kmx
   latd = 44
   lond = 353
   delt2 = delt * 2.
!
!  precipitation rate is needed in unit of kg m-2 s-1
!
   do i = 1,im
     prcp(i) = rhoh2o * (precip(i)/delt+snowmt(i))
     runoff(i) = 0.
     cnpy(i) = canopy(i)
   enddo
!
!  update canopy water content
!
   do i = 1,im
     if(slimsk(i).eq.1.) then
       rhscnpy(i) = rhscnpy(i) + sigmaf(i) * prcp(i)
       canopy(i) = canopy(i) + delt * rhscnpy(i)
       canopy(i) = max(canopy(i),0.)
       prcp(i) = prcp(i) * (1. - sigmaf(i))
       if(canopy(i).gt.scanop) then
         drip = canopy(i) - scanop
         canopy(i) = scanop
         prcp(i) = prcp(i) + drip / delt
       endif
!
!  calculate infiltration rate
!
       inf(i) = prcp(i)
       tsat(i) = osu_funct_thsat(soiltyp(i))
#ifdef HYDRO
       hydrow(i,1)=hydrow(i,1)-prcp(i)
#else
       infmax(i) = (-zsoil(i,1)) *                                           &
                 ((tsat(i) - smc(i,1)) / delt - rhsmc(i,1)) * rhoh2o
       infmax(i) = max(rffact*infmax(i),0.)
       if(inf(i).gt.infmax(i)) then
         runoff(i) = inf(i) - infmax(i)
         inf(i) = infmax(i)
       endif
       inf(i) = inf(i) / rhoh2o
       rhsmc(i,1) = rhsmc(i,1) - inf(i) / zsoil(i,1)
#endif
     endif
   enddo
!
!  we currently ignore the effect of rain on sea ice
!
   do i = 1,im
     flag(i) = slimsk(i).eq.1.
   enddo
!
!  solve the tri-diagonal matrix
!
#ifndef HYDRO
   do k = 1,km
     do i = 1,im
       if(flag(i))  then
         rhsmc(i,k) = rhsmc(i,k) * delt
         ai(i,k) = ai(i,k) * delt
         bi(i,k) = 1. + bi(i,k) * delt
         ci(i,k) = ci(i,k) * delt
       endif
     enddo
   enddo
#else
   do i = 1,im
     if(flag(i))  then
       rhsmc(i,1) = rhsmc(i,1) * delt
       rhsmc(i,2) = rhsmc(i,2) * delt
       ai(i,1) = 1. + ai(i,1) * delt /                                       &
                   (-zsoil(i,1)*rhoh2o)
       bi(i,1) = bi(i,1) * delt /                                            &
                   (-zsoil(i,1)*rhoh2o)
       ci(i,1) = smc(i,1) + (ci(i,1) - hydrow(i,1))*delt /                   &
                   (-zsoil(i,1)*rhoh2o)
       ai(i,2) = ai(i,2) * delt /                                            &
                   (-(zsoil(i,2)-zsoil(i,1))*rhoh2o)
       bi(i,2) = 1. + bi(i,2) * delt /                                       &
                   (-(zsoil(i,2)-zsoil(i,1))*rhoh2o)
       ci(i,2) = smc(i,2) + (ci(i,2) - hydrow(i,2))*delt /                   &
                   (-(zsoil(i,2)-zsoil(i,1))*rhoh2o)
     endif
   enddo
#endif
#ifndef HYDRO
!
!  implicit solver
!
!  forward elimination
!
   do i = 1,im
     if(flag(i)) then
       ci(i,1) = -ci(i,1) / bi(i,1)
       rhsmc(i,1) = rhsmc(i,1) / bi(i,1)
     endif
   enddo
   do k = 2,km
     do i = 1,im
       if(flag(i)) then
         cc = 1. / (bi(i,k) + ai(i,k) * ci(i,k-1))
         ci(i,k) = -ci(i,k) * cc
         rhsmc(i,k) = (rhsmc(i,k) - ai(i,k) * rhsmc(i,k-1)) * cc
       endif
     enddo
   enddo
!  backward substituttion
   do i = 1,im
     if(flag(i)) then
       ci(i,km) = rhsmc(i,km)
     endif
   enddo
   do k = km-1,1
     do i = 1,im
       if(flag(i)) then
         ci(i,k) = ci(i,k) * ci(i,k+1) + rhsmc(i,k)
       endif
     enddo
   enddo
!
!  update soil moisture
!
   do k = 1,km
     do i = 1,im
       if(flag(i)) then
         smsoil(i,k) = smc(i,k) + ci(i,k)
         smsoil(i,k) = max(smsoil(i,k),0.)
         tdif = max(smsoil(i,k) - tsat(i),0.)
         runoff(i) = runoff(i) - rhoh2o * tdif * zsoil(i,k) / delt
         smsoil(i,k) = smsoil(i,k) - tdif
       endif
     enddo
   enddo
#else
!
! explicit solver
!
   do i=1,im
     if(flag(i)) then
       if (km.ne.2) then
         if (iope) write(6,*)'This code is only for two soil levels'
         call MPABORT
       endif
!
! acr
! solving the following system of two equations and two unknowns 
! with x1=top layer soil moisture content and x2=bottom layer 
! soil moisture content:
! ai(i,1)*x1 + bi(i,1)*x2 = ci(i,1)
! ai(i,2)*x1 + bi(i,2)*x2 = ci(i,1)
!
       smsoil(i,1) = ci(i,1)/ai(i,1)-bi(i,1)/ai(i,1) *                       &
                  (ci(i,2) - ai(i,2)*ci(i,1)/ai(i,1)) /                      &
                  (bi(i,2) - bi(i,1)*ai(i,2)/ai(i,1))
       smsoil(i,2) = (ci(i,2)-ai(i,2)*ci(i,1)/ai(i,1)) /                     &
                  (bi(i,2)-bi(i,1)*ai(i,2)/ai(i,1))
     endif
   enddo
!
!  update soil moisture
!
   do k = 1,km
     do i = 1,im
       if(flag(i)) then
         if(k.gt.1) then
           dz=zsoil(i,k)-zsoil(i,k-1)
         else
           dz=zsoil(i,k)
         endif
!
!  handling negative soil moisture
!
         if(smsoil(i,k).lt.0.) then
           runoff(i)=runoff(i)-smsoil(i,k)*rhoh2o*dz/delt
           rnof=smsoil(i,k)
           smsoil(i,k)=0.
         else
           rnof=0.
         endif
!
!  handling over-saturation of soil moisture
!
         tdif = max(smsoil(i,k) - tsat(i),0.)
         runoff(i) = runoff(i)-rhoh2o*tdif*dz/delt
         smsoil(i,k) = smsoil(i,k) - tdif
       endif
     enddo
   enddo
#endif
   do k = 1,km
     do i = 1,im
       if(flag(i)) then
         smc(i,k) = smsoil(i,k)
       endif
     enddo
   enddo
   do i = 1,im
     if(flag(i)) then
       canopy(i) = canopy(i)
     endif
   enddo
!
   return
   end
