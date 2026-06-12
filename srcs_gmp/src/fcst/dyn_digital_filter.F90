#include <define.h>
!-------------------------------------------------------------------------------
   subroutine dyn_digital_filter(icall,hrini,chour,solsec,n1)
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [dyn_digital_filter]
!      |
!      |--- [dyn_digital_filter_flux] *
!      !
!      |--- [dyn_3d_digital_filter] <-- [dyn_3d_module] 
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only     :  LONF2S,LONF2S,LATG2S,LNT2S,levh_,levs_
   use constant, only      :  pi_
   use varsfc, only        :  lalbd_,nsoil_,nsoil_
#ifdef DFS
   use dfsvar, only        :  iope
#endif
   use dyn_3d_module, only :  dyn_3d_digital_filter
   use comsfc              !  sfcftyp,sfcfcs
   use comfibm
   use comio
!-------------------------------------------------------------------------------
   real                   ::  dthour,dshour,dchour,dsolsec,totsum
   real                   ::  qs(LNT2S)                                       ,&
                              tes(LNT2S,levs_),rqs(LNT2S,levh_)               ,&
                              dis(LNT2S,levs_),zes(LNT2S,levs_)
   character(len=128)     ::  fno
   integer                ::  ncho
!-------------------------------------------------------------------------------
   if(numsum.ge.nummax) return
!
   if( icall.eq.0 ) then
#ifndef NOPRINT
     if (iope) then
       print *,' initial dyn_digital_filter '
       print *,' ini time is ',hrini,' hour.'
     endif
#endif
     do k = 1,levs_
       do i = 1,LNT2S
         dis(i,k) = 0.0
         zes(i,k) = 0.0
         tes(i,k) = 0.0
       enddo
     enddo
!
     do k = 1,levh_
       do i = 1,LNT2S
         rqs(i,k) = 0.0
       enddo
     enddo
!
     do i = 1,LNT2S
       qs(i) = 0.0
     enddo 
!
     totsum=0.0
   endif
!
   numsum=numsum+1
#ifndef NOPRINT
   if (iope) then
     print *,' ---- in dyn_digital_filter ---- numsum nummax ',numsum,nummax
   endif
#endif
   if( numsum.ne.0 ) then
     sc = pi_ / nummax
     sx= numsum*sc
     tx= numsum*pi_
     wx= tx/ ( nummax+1 )
     digfil= sin(wx)/wx * sin(sx)/tx
   else
     digfil = 1.0/nummax
   endif
   totsum = totsum + digfil
!
!------------------------do summation with window---
!
! first lat loop
#ifdef ORIGIN_THREAD
!$doacross share(di ,ze ,te ,rq , q,
!$&              dis,zes,tes,rqs, qs,
!$&              digfil),
!$&        local(j,k)
#endif
#ifdef OPENMP
!$omp parallel do private(j,k)
#endif
!      autoscope
!
! .......obtain full field values
!
   do k = 1,levs_
     do j = 1,LNT2S
       dis(j,k) = dis(j,k) + digfil*di(j,k)
       zes(j,k) = zes(j,k) + digfil*ze(j,k)
       tes(j,k) = tes(j,k) + digfil*te(j,k)
     enddo
   enddo
   !
   do k = 1,levh_
     do j = 1,LNT2S
       rqs(j,k) = rqs(j,k) + digfil*rq(j,k)
     enddo
   enddo
   !
   do j = 1,LNT2S
     qs(j) = qs(j) + digfil*q(j)
   enddo
!
! save
!
   if( numsum.eq.0 ) then
     dthour=thour
     dshour=shour
     dchour=chour
     dsolsec=solsec
#ifndef NOPRINT
     if (iope) then
       print *,' numsum=0, save thour= ',dthour
     endif
#endif
     call file_name('sfc',3,thour,fno,ncho)
     call sfc_read_file_driver(n1,fno,sfcftyp,                                 &
                 labs,idate(4),idate(2),idate(3),idate(1),                     &
                 thour,sfcfcs,LONF2S,LATG2S,1)
   endif
!
! restore
!
   if( numsum.eq.nummax ) then
#ifndef NOPRINT
     if (iope) then
        print *,' numsum=nummax reassign perturbation '
        print *,' with normalized factor=',totsum,' at hour=',dthour
     endif
#endif
     hrini=0
     thour=dthour
     shour=dshour
     chour=dchour
     solsec=dsolsec
     fno='sfci'
     call sfc_read_file_driver(n1,fno,sfcftyp,                                 &
                 labs,idate(4),idate(2),idate(3),idate(1),                     &
                 thour,sfcfcs,LONF2S,LATG2S,0)
!
#ifdef ORIGIN_THREAD
!$doacross share(di ,ze ,te ,rq , q,
!$&              dis,zes,tes,rqs, qs,
!$&              dim,zem,tem, rm, qm,
!$&              totsum),
!$&        local(j,k)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(di ,ze ,te ,rq , q )
!mic$1 shared(dis,zes,tes,rqs, qs)
!mic$1 shared(dim,zem,tem, rm, qm)
!mic$1 shared(totsum)
!mic$1 private(j,k) 
#endif
#ifdef OPENMP
!$omp parallel do private(j,k)
#endif
!      autoscope
!
     do k = 1,levs_
       do j = 1,LNT2S
         di (j,k) = dis(j,k) / totsum
         ze (j,k) = zes(j,k) / totsum
         te (j,k) = tes(j,k) / totsum
         dim(j,k) = dis(j,k) / totsum
         zem(j,k) = zes(j,k) / totsum
         tem(j,k) = tes(j,k) / totsum
       enddo
     enddo
!
     do k = 1,levh_
       do j = 1,LNT2S
         rq (j,k) = rqs(j,k) / totsum
         rm (j,k) = rqs(j,k) / totsum
       enddo
     enddo
!
     do j = 1,LNT2S
       qm(j) = qs(j) / totsum
       q (j) = qs(j) / totsum
     enddo
!
     do l = 1,LATG2S
       do j = 1,LONF2S
         raintot(j,l)=0.5*raintot(j,l)
       enddo
     enddo
     call dyn_digital_filter_flux(0.5,dusfc,dvsfc,dtsfc,dqsfc,dlwsfc,ulwsfc,   &
               raincps,gflux,dugwd,dvgwd,psmean)
#ifdef DG3
     call diag_3d_digital_filter(0.5)
#endif
     do i = 1,LONF2S
       do j = 1,LATG2S
         do k = 1,25
           fluxr(i,j,k)=0.5*fluxr(i,j,k)
         enddo
       enddo
     enddo
     do i = 1,LONF2S
       do j = 1,LATG2S
         cvavg(i,j)=0.5*cvavg(i,j)
       enddo
     enddo
   endif
!
   return
   end subroutine dyn_digital_filter
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine dyn_digital_filter_flux(fac,dusfc,dvsfc,dtsfc,dqsfc,             &
                     dlwsfc,ulwsfc,raincps,gflux,                               &
                     dugwd,dvgwd,psmean,dtfulx)
!-------------------------------------------------------------------------------
   use paramodel, only : LONF2S,LATG2S
!-------------------------------------------------------------------------------
   real  ::  dusfc(LONF2S,LATG2S),dvsfc(LONF2S,LATG2S)
   real  ::  dtsfc(LONF2S,LATG2S),dqsfc(LONF2S,LATG2S)
   real  ::  dlwsfc(LONF2S,LATG2S),ulwsfc(LONF2S,LATG2S)
   real  ::  raincps(LONF2S,LATG2S),gflux(LONF2S,LATG2S)
   real  ::  dugwd(LONF2S,LATG2S),dvgwd(LONF2S,LATG2S)
   real  ::  psmean(LONF2S,LATG2S)
!-------------------------------------------------------------------------------
   do l = 1,LATG2S
     do j = 1,LONF2S
       dusfc(j,l)=fac*dusfc(j,l)
       dvsfc(j,l)=fac*dvsfc(j,l)
       dtsfc(j,l)=fac*dtsfc(j,l)
       dqsfc(j,l)=fac*dqsfc(j,l)
       dlwsfc(j,l)=fac*dlwsfc(j,l)
       ulwsfc(j,l)=fac*ulwsfc(j,l)
       raincps(j,l)=fac*raincps(j,l)
       gflux(j,l)=fac*gflux(j,l)
       dugwd(j,l)=fac*dugwd(j,l)
       dvgwd(j,l)=fac*dvgwd(j,l)
       psmean(j,l)=fac*psmean(j,l)
     enddo
   enddo
!
   return
   end subroutine dyn_digital_filter_flux
