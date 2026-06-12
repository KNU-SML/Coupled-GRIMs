#include "define.h"
MODULE module_ocean_slab_wrf

   contains

   SUBROUTINE ocean_mixed_layer_wrf(I,J,TML,T0ML,H,H0,HUML,                    &
           HVML,TSK,HFX,LH,GSW,GLW,UGLW,TMOML,                                 &
           UAIR,VAIR,UST,F,EMISS,STBOLT,G,DT,OML_GAMMA,                        & 
           ids,ide, jds,jde, kds,kde,                                          &
           ims,ime, jms,jme, kms,kme,                                          &
           its,ite, jts,jte, kts,kte                            )

!-------------------------------------------------------------------------------
!
! program history log:
!   2009-01-01  song-you hong          implementation from WRF
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! reference :
!   kim and hong (2010 JGR)
!
!------------------------------------------------------------------------------
   IMPLICIT NONE
!
   INTEGER, INTENT(IN   )    ::      I, J
   INTEGER, INTENT(IN   )    ::      ids,ide, jds,jde, kds,kde, &
                                     ims,ime, jms,jme, kms,kme, &
                                     its,ite, jts,jte, kts,kte

   REAL,    INTENT(INOUT)    ::   TML, H, H0, HUML, HVML, TSK

   REAL,    INTENT(IN      )    :: T0ML, HFX, LH, GSW, GLW, UGLW,  &
                                   UAIR, VAIR, UST, F, EMISS, TMOML
   REAL,    INTENT(IN) :: STBOLT, G, DT, OML_GAMMA
!
! Local
!
   REAL :: rhoair, rhowater, Gam, alp, BV2, A1, A2, B2, u, v, wspd, &
           hu1, hv1, hu2, hv2, taux, tauy, tauxair, tauyair, q, hold, &
           hsqrd, thp, cwater, ust2
   CHARACTER(LEN=120) :: time_series
!----------------------------------------------------------------
!
   hu1=huml
   hv1=hvml
   rhoair=1.
   rhowater=1000.
   cwater=4200.
!
! Deep ocean lapse rate (K/m) - from Rich
!
   Gam=oml_gamma
!
!  if(i.eq.1 .and. j.eq.1 .or. i.eq.105.and.j.eq.105) print *, 'gamma = ', gam
!  Gam=0.14
!  Gam=5.6/40.
!  if(i.eq.1 .and. j.eq.1 ) print *, 'gamma = ', gam
!  Gam=5./100.
!
! Thermal expansion coeff (/K)
!
!  alp=.0002
!  temp dependence (/K)
!
   alp=max((tml-273.15)*1.e-5, 1.e-6)
   BV2=alp*g*Gam
   thp=t0ml-Gam*(h-h0)

   A1=(tml-thp)*h - 0.5*Gam*h*h
   if(h.ne.0.)then
     u=hu1/h
     v=hv1/h
   else
     u=0.
     v=0.
   endif
!
!  time step
!
   q=(-hfx-lh+gsw+glw-stbolt*emiss*tml*tml*tml*tml)/(rhowater*cwater)
!  q=(-hfx-lh+gsw+glw-uglw)/(rhowater*cwater)
   wspd=max(sqrt(uair*uair+vair*vair),0.1)
!  wspd=sqrt(uair*uair+vair*vair)
!
   if (wspd .lt. 1.e-10 ) then
#ifdef DBG
     print *, 'OML::i,j,wspd are ', i,j,wspd
#endif
     wspd = 1.e-10
   endif
!
! limit ust to 1.6 to give a value of ust for water of 0.05
!       ust2=min(ust, 1.6)
! new limit for ust: reduce atmospheric ust by half for ocean
!
   ust2=0.5*ust
   tauxair=ust2*ust2*uair/wspd
   taux=rhoair/rhowater*tauxair
   tauyair=ust2*ust2*vair/wspd
   tauy=rhoair/rhowater*tauyair
!
! note: forward-backward coriolis force for effective time-centering
!
   hu2=hu1+dt*( f*hv1 + taux)
   hv2=hv1+dt*(-f*hu2 + tauy)
!
! consider the flux effect
!
   A2=A1
!  A2=A1+q*dt
!
   huml=hu2
   hvml=hv2

   hold=h
   B2=hu2*hu2+hv2*hv2
   hsqrd=-A2/Gam + sqrt(A2*A2/(Gam*Gam) + 2.*B2/BV2)
   h=sqrt(max(hsqrd,0.0))
!
! limit to positive h change
!
   if(h.lt.hold) h=hold
!
!  if(h.ne.0.)then
! no change unless tml is warmer than layer mean temp tmol or tsk-5 (see omlinit)
!
   if(tml.ge.tmoml .and. h.ne.0.)then
     tml=max(t0ml - Gam*(h-h0) + 0.5*Gam*h + A2/h, tmoml)
     u=hu2/h
     v=hv2/h
   else
     tml=t0ml
     u=0.
     v=0.
   endif
!
   tsk=tml
!
! ww: output point data
!     if( (i.eq.190 .and. j.eq.115) .or. (i.eq.170 .and. j.eq.125) ) then
!        write(jtime,fmt='("TS ",f10.0)') float(itimestep)
!        CALL wrf_message ( TRIM(jtime) )
!        write(time_series,fmt='("OML",2I4,2F9.5,2F8.2,2E15.5,F8.3)') &
!              i,j,u,v,tml,h,taux,tauy,a2
!        CALL wrf_message ( TRIM(time_series) )
!     end if
!
   END SUBROUTINE ocean_mixed_layer_wrf
!
!-------------------------------------------------------------------------------
    SUBROUTINE ocean_ml_wrf_init(                                              &
                           tsk,tml,t0ml,hml,h0ml,hml0,huml,hvml,               &
                           tmoml,its,itf,jts,jtf)
!-------------------------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------------------------
!  LOGICAL , INTENT(IN)      ::      allowed_to_read
!  LOGICAL , INTENT(IN)      ::      start_of_simulation
!  INTEGER, INTENT(IN   )    ::      ids,ide, jds,jde, kds,kde, &
!                                    ims,ime, jms,jme, kms,kme, &
!                                    its,ite, jts,jte, kts,kte
!
   INTEGER, INTENT(IN   )    ::      its,itf,jts,jtf
!
   REAL,    DIMENSION( its:itf, jts:jtf ),  INTENT(INOUT)    ::   TSK
!
   REAL,    DIMENSION( its:itf, jts:jtf )                                    , &
            INTENT(INOUT)    :: T0ML,TML,HML,H0ML,hml0,HUML,HVML,TMOML
!
!  LOCAR VAR
!
!  INTEGER                   ::      L,J,I,itf,jtf
   INTEGER                   ::      L,J,I
   CHARACTER*1024 message
!-------------------------------------------------------------------------------

!   IF(start_of_simulation .AND. &
!      PRESENT(tml) .AND. PRESENT(t0ml) ) THEN
!
   DO J = jts,jtf
     DO I = its,itf
       TML(I,J)=TSK(I,J)
       T0ML(I,J)=TSK(I,J)
!
! MAY HAVE INPUT OF HML BUT FOR NOW SET HERE
!
       H0ML(I,J)=hml0(i,j)
       HML(I,J)=hml0(i,j)
       HUML(I,J)=0.
       HVML(I,J)=0.
       TMOML(I,J)=TSK(I,J)-5.

!      IF(TMOML(I,J).GT.200. .and. TMOML(I,J).LE.201.) TMOML(I,J)=TSK(I,J)
     ENDDO
   ENDDO
!   ENDIF
!  ENDIF
!
   END SUBROUTINE ocean_ml_wrf_init
!
!-------------------------------------------------------------------------------
   SUBROUTINE sst_skin_update(glw,uglw,gsw,hfx,qfx,tsk,ust,emiss,              &
                dtw1,sstsk,dt,stbolt,                                          &
                ids, ide, jds, jde, kds, kde,                                  &
                ims, ime, jms, jme, kms, kme,                                  &
                its, ite, jts, jte, kts, kte                       )
#ifdef DBG
#ifdef DFS
   use dfsvar, only : iope
#else
   use comio, only : iope
#endif
#endif
!-------------------------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------------------------
!
!  modify the 2d - module to 1d module....hong
!
!-------------------------------------------------------------------------------
   INTEGER , INTENT(IN)           :: ids, ide, jds, jde, kds, kde,             &
                                     ims, ime, jms, jme, kms, kme,             &
                                     its, ite, jts, jte, kts, kte
!
!  REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: glw, gsw, uglw
!  REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(IN) :: hfx, qfx
!  REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(IN) :: ust, emiss
!  REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(INOUT  ) :: dtw1         ! warm temp difference (C)
!  REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(INOUT  ) :: sstsk,TSK        ! skin sst (K)
   REAL,     INTENT(IN   ) :: glw, gsw, uglw
   REAL,     INTENT(IN) :: hfx, qfx
   REAL,     INTENT(IN) :: ust, emiss
   REAL,     INTENT(INOUT  ) :: dtw1         ! warm temp difference (C)
   REAL,     INTENT(INOUT  ) :: sstsk,TSK    ! skin sst (K)
   REAL,     INTENT(IN )   ::   DT           ! model time step
   REAL,     INTENT(IN )   ::   STBOLT       ! Stefan-Boltzmann constant (W/m^2/K^4)
!
! Local
!
   REAL :: lw, sw, q, qn, zeta, dep, dtw3, skinmax, skinmin
   REAL :: fs, con1, con2, con3, con4, con5, zlan, q2, ts, phi, qn1
   REAL :: usw, qo, swo, us, tb, dtc, dtw, alw, dtwo, delt, f1
   INTEGER :: i, j, k
!-------------------------------------------------------------------------------
   INTEGER , PARAMETER :: n=1152
   REAL , PARAMETER :: z1=3.,an=.3,zk=.4,rho=1.2,rhow=1025.,cw=4190.
   REAL , PARAMETER :: g=9.8,znuw=1.e-6,zkw=1.4e-7,sdate=1201.6667
!-------------------------------------------------------------------------------
!  parameter(g=9.8,delt=900.,znuw=1.e-6,zkw=1.4e-7)
!
!  Input arguments
!  (all fluxes are positive downwards)
!     real qo      ! LH + SH + LW (W/m^2), + down
!     real swo      ! Net shortwave flux (W/m^2), + down
!     real u       ! Wind speed (m/s)
!     real us      ! Atmospheric friction velocity (m/s)
!     real tb      ! Bulk temperature (deg C)
!     real dtwo    ! Warm layer temp. diff. from previous time (deg C)
!  Local variables
!     real lw
!     real sw
!     real q       ! LH + SH + LW
!     real qn      ! Q + R_s - R(-d)
!     real zeta    ! -z / L
!     real dep     ! Skin layer depth (m)
!     real dtw3
!  Output variables
!     real dtw     ! Warm layer temp. diff. (deg C)
!     real dtc     ! Cool skin temp. diff. (deg C)
!     real ts      ! Skin temperature (deg C)
!      q=lh+sh+lwo
!
   skinmax=-9999.
   skinmin=9999.
!
!  do i = its,ite
!  do j = jts,jte
!
!  if(xland(i,j).ge.1.5) then
!  if(xland.eq.0) then
!  qo=glw(i,j)-emiss(i,j)*stbolt*(sstsk(i,j)**4)-qfx(i,j)-hfx(i,j)
!
   qo=glw-uglw-qfx-hfx
   swo=gsw
   us=max(ust,0.01)
   tb=tsk-273.15
   dtwo=dtw1
   delt=dt
!
   q=qo/(rhow*cw)
   sw=swo/(rhow*cw)
!
! TEMPORARY KLUDGE
!
!  f1=1.-0.28*exp(-71.5*z1)-0.27*exp(-2.8*z1)-0.45*exp(-0.07*z1)
!
   f1=1.                   -0.27*exp(-2.8*z1)-0.45*exp(-0.07*z1)
!
! cool skin
!
   dtc=0.0
!
! tb in C
!
   alw=1.e-5*max(tb,1.)
   con4=16.*g*alw*znuw**3/zkw**2
   usw=sqrt(rho/rhow)*us
   con5=con4/usw**4
!
! otherwise, iterations would be needed for the computation of fs
! iteration impact is less than 0.03C
!
   q2=max(1./(rhow*cw),-q)
   zlan=6./(1.+(con5*q2)**0.75)**0.333
   dep=zlan*znuw/usw                    ! skin layer depth (m)
   fs=0.065+11.*dep-(6.6e-5/dep)*(1.-exp(-dep/8.e-4))
   fs=max(fs,0.01)          ! fract. of solar rad. absorbed in sublayer
   dtc=dep*(q+sw*fs)/zkw            ! cool skin temp. diff (deg C)
   dtc=min(dtc,0.)
!
! warm layer (X. Zeng)
!
   dtw=0.0
!
! tb in C
!
   alw=1.e-5*max(tb,1.)
   con1=sqrt(5.*z1*g*alw/an)
   con2=zk*g*alw
   qn=q+sw*f1
   usw=sqrt(rho/rhow)*us
!
!  does not change when qn is positive
!
   if(dtwo.gt.0..and.qn.lt.0.) then
     qn1=sqrt(dtwo)*usw**2/con1
     qn=max(qn,qn1)
   endif
!
   zeta=z1*con2*qn/usw**3
!
   if(zeta.gt.0.) then
     phi=1.+5.*zeta
   else
     phi=1./sqrt(1.-16.*zeta)
   endif
!
   con3=zk*usw/(z1*phi)
!
! use all SW flux
!
   dtw=(dtwo+(an+1.)/an*(q+sw*f1)*                                             &
                          delt/z1)/(1.+(an+1.)*con3*delt)
   dtw=max(0.,dtw)
   dtwo=dtw
   ts = tb + dtw + dtc
!
   skinmax=amax1(skinmax,ts-tb)
   skinmin=amin1(skinmin,ts-tb)
   sstsk=ts+273.15       ! convert ts (in C) to sstsk (in K)
   dtw1=dtw              ! dtw always in C
!
   tsk = sstsk
!
!  end do
!  end do
#ifdef DBG
#ifdef OPENMP
!$omp master
#endif
   if(iope) print *, 'check skin sst skinmax = ', skinmax, '  skinmin = ', skinmin
#ifdef OPENMP
!$omp end master
#endif
#endif
!
   return

   END SUBROUTINE sst_skin_update
!-------------------------------------------------------------------------------
END MODULE module_ocean_slab_wrf


