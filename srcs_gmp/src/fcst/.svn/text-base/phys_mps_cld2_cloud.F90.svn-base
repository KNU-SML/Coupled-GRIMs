#include <define.h>
   subroutine phys_mps_cld2_cloud (im ,ix  ,km  ,deltim ,prsi ,prsl ,          &
                         q ,cwm ,t   ,lat    ,                                 &
                        tp ,qp  ,psp ,tp1    ,qp1  ,psp1 ,rcs)
!-------------------------------------------------------------------------------
!
! subroutine: phys_mps_cld2_cloud
! 
! abstract:
!   subroutine for grid-scale condensation & evaporation        
!   for the mrf model at ncep.                                 
!
! program history log:
!   1995-01   q. zhao      created
!   1996-01   s.-y  hong   modified
!   1998-09   h.-l. pan    modified
!   1999-08   s. moorthi   modified 
!   2000-04   s.-y. hong   modified
!
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS,levs_
   use constant, only  : cp=>cp_,g=>g_,eliv=>hsub_,elwv=>hvap_,                &
                         rd_,rv_,psat_,ttp_
!-------------------------------------------------------------------------------
   real,parameter       ::  h1=1.e0,h2=2.e0,h1000=1000.0
   real,parameter       ::  d00=0.e0,d125=.125e0,d5=0.5e0
   real,parameter       ::  a1=psat_,eps=rd_/rv_,epsm1=eps-1.
   real,parameter       ::  epsq=2.e-12,tm10=ttp_-10.,r=rd_
   real,parameter       ::  cpr=cp*r,rcpr=h1/(cpr)
   real,parameter       ::  rcp=h1/cp
   integer              ::  im, km
!
   real                 ::  prsl(ix,km), prsi(ix,km+1)
   real                 ::  q(ix,km), t(ix,km)
   real                 ::  tp(ix,km),  qp(ix,km),  psp(ix), cwm(ix,km)
   real                 ::  tp1(ix,km), qp1(ix,km), psp1(ix)
!
!  local variables and arrays
!
   real                 ::  qi(ILOTS),qint(ILOTS)
   real                 ::  u(levs_)
   integer              ::  iw(ILOTS,levs_)
   logical              ::  lprnt
!
!-----------------prepare constants for later uses------------------------------
!
   dt      = 2.0*deltim
   rdt     = h1/dt
   c0      = 1.5e-4
   c1      = 300.0
   c2      = 0.5
   us      = h1
   cclimit = 1.0e-3
   climit  = 1.0e-20
   u00b    = 0.85
   u00t    = 0.85
   rcsi = 1.0 / rcs
!-------------------------------------------------------------------------------
   do k = 1,km
     u(k) = u00b + (u00t-u00b)*float(k-1)/float(km-1)
     u(k)  = min(0.99, u(k)+(1.0-u(k))*(1.0-rcsi))
   enddo
!------------------set ice id iw(i,k)=0 at top of the domain--------------------
   do  i = 1,im
     iw(i,km) = d00
   enddo
!-------------------------------------------------------------------------------
   do k = km,1,-1
!------------------qw, qi and qint----------------------------------------------
     do i = 1,im                                    
       tmt0  = t(i,k)-273.16                                                
       tmt15 = min(tmt0,-15.)                                            
       qik   = max(q(i,k),epsq)
       cwmik = max(cwm(i,k),climit)
       ai    = 0.008855
       bi    = 1.0
       if (tmt0 .lt. -20.0) then
         ai = 0.007225
         bi = 0.9674
       end if
!
!  the global qsat computation is done in cb
!
       pres    = prsl(i,k)
!         qw      = fpvs(t(i,k))
       qw      = fpvs0(t(i,k))
       qw      = eps * qw / (pres + epsm1 * qw)
       qw      = max(qw,epsq)
       qi(i)   = qw *(bi+ai*min(tmt0,0.))
       qint(i) = qw *(1.-0.00032*tmt15*(tmt15+15.))
       if (tmt0 .le. -40.) qint(i) = qi(i)
!-------------------ice-water id number iw--------------------------------------
       if(tmt0.lt.-15.0) then
         u00ik = u(k)
         fi    = qik - u00ik*qi(i)    
         if(fi.gt.d00.or.cwmik.gt.climit) then                    
           iw(i,k) = 1                                                   
         else                                                           
           iw(i,k) = 0                                                   
         end if                                                         
       end if
!
       if(tmt0.ge.0.0) then
         iw(i,k) = 0
       end if
!
       if (tmt0 .lt. 0.0 .and. tmt0 .ge. -15.0) then
         iw(i,k) = 0
         if (k .lt. km) then
           if (iw(i,k+1) .eq. 1 .and. cwmik .gt. climit) iw(i,k) = 1
         endif
       end if
     enddo
!
!--------------condensation and evaporation of cloud----------------------------
!
     do i = 1,im
!
!------------------------at, aq and dp/dt---------------------------------------
!
       qik   = max(q(i,k),epsq)
       cwmik = max(cwm(i,k),climit)
       iwik  = iw(i,k)
       u00ik = u(k)
       tik   = t(i,k)
       pres  = prsl(i,k) * h1000
       pp0   = psp(i)*prsl(i,k)/prsi(i,1) * h1000
       at    = (tik-tp(i,k)) * rdt
       aq    = (qik-qp(i,k)) * rdt
       ap    = (pres-pp0)    * rdt
!
!----------------the satuation specific humidity--------------------------------
!
       fiw   = float(iwik)
       elv   = (h1-fiw)*elwv    + fiw*eliv
       qc    = (h1-fiw)*qint(i) + fiw*qi(i)
!
!----------------the relative humidity------------------------------------------
!
       if(qc.le.1.0e-10) then
         rqik=d00 
       else
         rqik = qik/qc
       endif
!
!----------------cloud cover ratio ccrik----------------------------------------
!
       if (rqik .lt. u00ik) then
          ccrik = d00
       elseif(rqik.ge.us) then
          ccrik = us
       else
          rqikk  = min(us,rqik)
          ccrik = h1-sqrt((us-rqikk)/(us-u00ik))
       endif
!
!   if no cloud exists then evaporate any existing cloud condensate
!----------------evaporation of cloud water-------------------------------------
!
       e0 = d00
       if (ccrik.le.cclimit.and.cwmik.gt.climit)  then 
         e0 = max(qc*(u00ik-rqik)*rdt, 0.0)
         e0 = min(cwmik*rdt,   e0)
         e0 = max(0.0,e0)
       end if
!
!   if cloud cover > 0.2 condense water vapor in to cloud condensate
!-----------the eqs. for cond. has been reorganized to reduce cpu---------------
!
       cond = d00
       if (ccrik .gt. 0.20 .and. qc .gt. epsq) then
         us00   = us  - u00ik 
         ccrik1 = 1.0 - ccrik
         aa     = eps*elv*pres*qik
         ab     = ccrik*ccrik1*qc*us00
         ac     = ab + 0.5*cwmik
         ad     = ab * ccrik1
         ae     = cpr*tik*tik
         af     = ae * pres
         ag     = aa * elv
         ai     = cp * aa
         cond   = (ac-ad)*(af*aq-ai*at+ae*qik*ap)/(ac*(af+ag))
         condi  = (qik   -u00ik   *qc*1.0)*rdt
         cond   = min(cond, condi)
         cond = max(cond, d00)
       end if
       cone0    = (cond-e0) * dt
       cwm(i,k) = cwm(i,k) + cone0
       t(i,k)   = t(i,k)   + elv*rcp*cone0
       q(i,k)   = q(i,k)   - cone0
     enddo                                  ! end of i-loop!
   enddo                                    ! end of k-loop!
!
!----------------store t, q, ps for next time step------------------------------
!
   do k = 1,km
     do i = 1,im
       tp(i,k)  = tp1(i,k)
       qp(i,k)  = qp1(i,k)
       tp1(i,k) = t(i,k)
       qp1(i,k) = max(q(i,k),epsq)
     enddo
   enddo
   do i = 1,im
     psp(i)  = psp1(i)
     psp1(i) = prsi(i,1)
   enddo
!-------------------------------------------------------------------------------
   return
   end subroutine phys_mps_cld2_cloud
