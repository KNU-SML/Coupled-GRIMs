#include <define.h>
   subroutine phys_cps_kuo(delt,delprsi,prsl,prsi,prslk,prsik,                 &
                           qn,q,t,rn,                                          &
                           jcap,kbot,ktop,icps,rd,rv,cp,t0c,g,hvap,            &
                             ids,ide, jds,jde, kds,kde,                        &
                             ims,ime, jms,jme, kms,kme,                        &
                             its,ite, jts,jte, kts,kte)
!-------------------------------------------------------------------------------
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only : ilot=>ILOTS,klot=>levs_
!-------------------------------------------------------------------------------
   real                ::  delprsi(ims:ime,kms:kme),                           &
                           prsl(ims:ime,kms:kme),prslk(ims:ime,kms:kme),       &
                           prsi(ims:ime,kms:kme+1),prsik(ims:ime,kms:kme+1)
   real                ::  q(ims:ime,kms:kme),qn(ims:ime,kms:kme)
   real                ::  t(ims:ime,kms:kme),rn(ims:ime)
   integer             ::  kbot(ims:ime),ktop(ims:ime),icps(ims:ime)
!
!  bounds of parcel origin
!
   integer,parameter    ::  kliftl=2,kliftu=3
!
!  sigma below which to compute moisture convergence
!
   real,parameter       ::  sigdhq=0.65
!
!  local variables and arrays
!
!  this is checkout demonstration
!
   integer              ::  index2(ilot),klcl2(ilot),kbot2(ilot),ktop2(ilot)
   real                 ::  adq2(ilot),aqc2(ilot),atc2(ilot),rn2(ilot)
   real                 ::  delprsi2(ilot,klot),prsl2(ilot,klot),prsi2(ilot,klot)
   real                 ::  prslk2(ilot,klot),prsik2(ilot,klot)
   real                 ::  sdq(ilot),sd2(ilot),sdq2(ilot),sdqc2(ilot)
   real                 ::  sq2(ilot),sqs2(ilot),sqc2(ilot),stc2(ilot)
   real                 ::  dq2(ilot*klot),q2(ilot*klot),qc2(ilot*klot)
   real                 ::  t2(ilot*klot),tc2(ilot*klot)
!
#ifdef KUM
   real                 ::  bdq2(ilot), q3(ilot*klot)
#endif
   real                ::  cpoel,elocp,el2orc,eps,epsm1
!-------------------------------------------------------------------------------
!
!  physical parameters
!
   cpoel=cp/hvap
   elocp=hvap/cp
   el2orc=hvap*hvap/(rv*cp)
   eps=rd/rv
   epsm1=rd/rv-1.
!
! critical relative humidity for computing beta
!
   rhcrit = 1.
#ifdef KUM
   rhcrit = 1.
#endif
!
! condensation factor to compute sdqc2
!
   gamma = 1.
#ifdef KUM
   gamma = 1.
#endif
!
! constant beta
!
#ifdef KUM
   beta0 = 0.
   beta0 = 0.1
!
#endif
!-------------------------------------------------------------------------------
!  initialize arrays and compute moisture convergence.
!  compress fields to points with minimum moisture convergence
!  and acceptable boundary layer temperatures.
!
   do i = ims,ime
     rn(i)=0.
     kbot(i)=kte+1
     ktop(i)=0
     icps(i)=0
     sdq(i)=0.
   enddo
!
   do k = kts,kte
     do i = ims,ime
       if(prsl(i,k).ge.prsi(i,1)*sigdhq) then
         sdq(i)=sdq(i)+(q(i,k)-qn(i,k))*(delprsi(i,k)/prsi(i,1))
       endif
     enddo
   enddo
!
!  crconv is in units of g times meters of rain per second.
!
   crconv=2.e-6
   crconv=crconv*(jcap/80.)**2
   if(jcap.eq.126) crconv=6.e-6
   gmrn=delt*crconv
!
   im2=0
   do i = ims,ime
     if(sdq(i)*prsi(i,1).gt.gmrn.and.t(i,2).gt.max(t0c+5.,t(i,3))) then
       im2=im2+1
       index2(im2)=i
     endif
   enddo
   if(im2.eq.0) return
!
   do k = kts,kte
#ifdef CRAY_THREAD
!fpp$ select(vector)
#endif
     do i = 1,im2
       ik=(k-1)*im2+i
       prsl2(i,k)=prsl(index2(i),k)
       prsi2(i,k)=prsi(index2(i),k)
       prslk2(i,k)=prslk(index2(i),k)
       prsik2(i,k)=prsik(index2(i),k)
       delprsi2(i,k)=delprsi(index2(i),k)
       dq2(ik)=q(index2(i),k)-qn(index2(i),k)
       q2(ik)=qn(index2(i),k)
       t2(ik)=t(index2(i),k)
     enddo
   enddo
!
!-------------------------------------------------------------------------------
!  compute moist adiabat and determine cloud boundaries.
!  sum fields within cloud.  restore humidity.
!
   call phys_moist_adiabat(im2,kte,kliftl,kliftu,                              &
                           prsl2,prsik2,prslk2,t2,q2,                          &
                           klcl2,kbot2,ktop2,tc2,qc2,rd,rv,                    &
                                   ids,ide, jds,jde, kds,kde,                  &
                                   ims,ime, jms,jme, kms,kme,                  &
                                   its,ite, jts,jte, kts,kte)
 
   kbm2=kte+1
   ktm2=0
   kbx2=0
!
   do i = 1,im2
     kbm2=min(kbm2,kbot2(i))
     ktm2=max(ktm2,ktop2(i))
!
!  uncomment next line to turn on evaporation of falling rain
!       if(kbot2(i).le.ktop2(i)) kbx2=max(kbx2,kbot2(i)-1)
!
     sd2(i)=0.
     sdq2(i)=0.
     sq2(i)=0.
     sqs2(i)=0.
     sqc2(i)=0.
     stc2(i)=0.
#ifdef KUM
     adq2(i) = 0.
     bdq2(i) = 0.
#endif
   enddo
!
   if(ktm2.lt.kbm2) return
!
   do k = kbm2,ktm2
     do i = 1,im2
       ik=(k-1)*im2+i
       if(k.ge.kbot2(i).and.k.le.ktop2(i)) then
#ifdef ICE
         pvs2 = fpvs(t2(ik))
#else
         pvs2 = fpvs0(t2(ik))
#endif
         qs2=eps*pvs2/(prsl2(i,k)+epsm1*pvs2)
         sd2(i)=sd2(i)+(delprsi2(i,k)/prsi(i,1))
         sdq2(i)=sdq2(i)+dq2(ik)*(delprsi2(i,k)/prsi(i,1))
         sq2(i)=sq2(i)+q2(ik)*(delprsi2(i,k)/prsi(i,1))
         sqs2(i)=sqs2(i)+qs2*(delprsi2(i,k)/prsi(i,1))
         sqc2(i)=sqc2(i)+qc2(ik)*(delprsi2(i,k)/prsi(i,1))
         stc2(i)=stc2(i)+tc2(ik)*(delprsi2(i,k)/prsi(i,1))
       endif
     enddo
   enddo
!
   do ik = 1,im2*kte
#ifdef KUM
     q3(ik) = q2(ik)
#endif
     q2(ik)=q2(ik)+dq2(ik)
   enddo
!
!-------------------------------------------------------------------------------
!  cloud must extend over 0.3 of ps and have moisture convergence.
!  compute partitioning of heating and moistening and rainfall.
!  evaporate rain below cloud base.
!  expand fields back again.
!
   do i = 1,im2
     rn2(i)=0.
     sdqmax=sqc2(i)+cpoel*stc2(i)
     sdqc2(i)=min( gamma*sdq2(i), sdqmax )
!
     if (sd2(i).gt.0.3.and.sdqc2(i).gt.0..and.sqs2(i).gt.0.) then
       adq2(i)=-sdqc2(i)/sdq2(i)
#ifdef KUM
       if ( adq2(i) .lt. -1. ) then
         adq2(i) = -1.
         bdq2(i) = -(sdqc2(i) - sdq2(i))/sq2(i)
       end if
#endif
       beta = rhcrit - sq2(i)/sqs2(i)
#ifdef KUM
       beta = beta0
#endif
       beta=min(beta,sqc2(i)/sdqc2(i))
       beta=max(beta,1.-cpoel*stc2(i)/sdqc2(i))
       if(beta.le.0.) then
         aqc2(i)=0.
         betanow = 1.
#ifdef KUM
         betanow = 1. - beta
#endif
         atc2(i) = betanow*sdqc2(i)/(cpoel*stc2(i))
#ifdef KUM
         bdq2(i) = bdq2(i) + beta*sdqc2(i)/sq2(i)
#endif
       elseif(beta.ge.1.) then
         aqc2(i)=sdqc2(i)/sqc2(i)
         atc2(i)=0.
       else
         aqc2(i)=beta*sdqc2(i)/sqc2(i)
         atc2(i)=(1.-beta)*sdqc2(i)/(cpoel*stc2(i))
       endif
     else
       kbot2(i)=kte+1
       ktop2(i)=0
     endif
!
   enddo
!
   do k = kbm2,ktm2
     do i = 1,im2
       ik=(k-1)*im2+i
       if(k.ge.kbot2(i).and.k.le.ktop2(i)) then
         dpovg=delprsi2(i,k)/g
#ifndef KUM
         qchg2=aqc2(i)*qc2(ik)+adq2(i)*dq2(ik)
#else
         qchg2=aqc2(i)*qc2(ik)+adq2(i)*dq2(ik) + bdq2(i)*q3(ik)
#endif
         tchg2=atc2(i)*tc2(ik)
         q2(ik)=q2(ik)+qchg2
         t2(ik)=t2(ik)+tchg2
         rn2(i)=rn2(i)+dpovg*cpoel*tchg2
       endif
     enddo
   enddo
!
   do k = kbx2,1,-1
     do i = 1,im2
       if(rn2(i).gt.0..and.k.lt.kbot2(i)) then
         ik=(k-1)*im2+i
#ifdef ICE
         pvs2=fpvs(t2(ik))
#else
         pvs2=fpvs0(t2(ik))
#endif
         qs2=eps*pvs2/(prsl2(i,k)+epsm1*pvs2)
         qchg2=qs2-q2(ik)
         if(qchg2.gt.0.) then
           dpovg=delprsi2(i,k)/g
           qchg2=qchg2/(1.+el2orc*qs2/t2(ik)**2)
           qchg2=qchg2*(1.-exp(-0.32*sqrt(2.*delt*rn2(i))))
           rnchg2=min(dpovg*qchg2,rn2(i))
           qchg2=rnchg2/dpovg
           q2(ik)=q2(ik)+qchg2
           t2(ik)=t2(ik)-elocp*qchg2
           rn2(i)=rn2(i)-rnchg2
         endif
       endif
     enddo
   enddo
!
   do i = 1,im2
     rn(index2(i))=rn2(i)
     kbot(index2(i))=kbot2(i)
     ktop(index2(i))=ktop2(i)
     if(rn2(i).gt.0.) icps(index2(i))=1
   enddo
!
   do k = 1,ktm2
#ifdef CRAY_THREAD
!fpp$ select(vector)
#endif
     do i = 1,im2
       ik=(k-1)*im2+i
       q(index2(i),k)=q2(ik)
       t(index2(i),k)=t2(ik)
     enddo
   enddo
!
   return
   end subroutine phys_cps_kuo
