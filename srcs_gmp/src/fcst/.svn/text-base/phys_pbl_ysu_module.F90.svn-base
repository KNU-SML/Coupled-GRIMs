!
#undef WRF
!
#define GAMMA_PBPS
#define ENT_PBPS
#define GAMMA_PT
#define ENT_PT
#define GAMMA_EL
#define KTKE_YSU
#define MLNOH
#define TKE_OUTPUT
#define VDIFQ_OUTPUT
#ifndef WRF
#undef TKE_OUTPUT
#undef VDIFQ_OUTPUT 
#include <define.h>
#endif
!
MODULE module_ysutke
!
#ifdef WRF
   USE MODULE_MODEL_CONSTANTS
#else
   USE CONSTANT, ONLY : G=>G_,R_D=>RD_,R_V=>RV_,CP=>CP_                 &
                       ,XLS=>HSUB_,XLV=>HVAP_,P608=>RVORDM1_
#endif
!
!-----------------------------------------------------------------------
!
      INTEGER :: ITRMX=5 ! ITERATION COUNT FOR MIXING LENGTH COMPUTATION
      REAL,PARAMETER :: PI=3.1415926,VKARMAN=0.4
!
!-----------------------------------------------------------------------
!***  QNSE MODEL CONSTANTS
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: EPSQ2L=0.01
      REAL,PARAMETER :: C0=0.55,CEPS=C0**3,BLCKDR=0.0063,CN=0.75        &
     &                 ,AM1=8.0,AM2=2.3,AM3=35.0,AH1=1.4,AH2=-0.01      &
     &                 ,AH3=1.29,AH4=2.44,AH5=19.8                      &
     &                 ,ARIMIN=0.127,BM1=2.88,BM2=16.0,BH1=3.6,BH2=16.0 &
     &                 ,BH3=720.0,EPSKM=1.E-3
      REAL,PARAMETER :: CAPA=R_D/CP
      REAL,PARAMETER :: RLIVWV=XLS/XLV,ELOCP=2.72E6/CP
      REAL,PARAMETER :: EPS1=1.E-12,EPS2=0.
      REAL,PARAMETER :: EPSL=0.32,EPSRU=1.E-7,EPSRS=1.E-7               &
     &                 ,EPSTRB=1.E-24
      REAL,PARAMETER :: EPSA=1.E-8,EPSIT=1.E-4,EPSU2=1.E-4,EPSUST=0.07
      REAL,PARAMETER :: ALPH=0.30,BETA=1./273.,EL0MAX=1000.,EL0MIN=1.   &
     &                 ,ELFC=0.23*0.5,GAM1=0.2222222222222222222        &
     &                 ,PRT=1.
      REAL,PARAMETER :: A1=0.659888514560862645                         &
     &                 ,A2X=0.6574209922667784586                       &
     &                 ,B1=11.87799326209552761                         &
     &                 ,B2=7.226971804046074028                         &
     &                 ,C1=0.000830955950095854396
      REAL,PARAMETER :: A2S=17.2693882,A3S=273.16,A4S=35.86
      REAL,PARAMETER :: ELZ0=0.,ESQ=5.0,EXCM=0.001                      &
     &                 ,FHNEU=0.8,GLKBR=10.,GLKBS=30.                   &
     &                 ,QVISC=2.1E-5,RFC=0.191,RIC=0.505,SMALL=0.35     &
     &                 ,SQPR=0.84,SQSC=0.84,SQVISC=258.2,TVISC=2.1E-5   &
     &                 ,USTC=0.7,USTR=0.225,VISC=1.5E-5                 &
     &                 ,WOLD=0.15,WWST=1.2,ZTMAX=1.,ZTFC=1.,ZTMIN=-5.
!
#ifdef WRF
      REAL,PARAMETER :: SEAFC=0.98,PQ0SEA=PQ0*SEAFC
#else
      REAL,PARAMETER :: SEAFC=0.98,PQ0=379.90516,PQ0SEA=PQ0*SEAFC
#endif
!
      REAL,PARAMETER :: BTG=BETA*G,CZIV=SMALL*GLKBS                     &
#ifdef WRF
     &                 ,ESQHF=0.5*5.0,GRRS=GLKBR/GLKBS                  &
#else
     &                 ,EP_1=R_V/R_D-1.,ESQHF=0.5*5.0,GRRS=GLKBR/GLKBS  &
#endif
     &                 ,RB1=1./B1,RTVISC=1./TVISC,RVISC=1./VISC         &
     &                 ,ZQRZT=SQSC/SQPR
!
      REAL,PARAMETER :: ADNH= 9.*A1*A2X*A2X*(12.*A1+3.*B2)*BTG*BTG      &
     &                 ,ADNM=18.*A1*A1*A2X*(B2-3.*A2X)*BTG              &
     &                 ,ANMH=-9.*A1*A2X*A2X*BTG*BTG                     &
     &                 ,ANMM=-3.*A1*A2X*(3.*A2X+3.*B2*C1+18.*A1*C1-B2)  &
     &                                *BTG                              &
     &                 ,BDNH= 3.*A2X*(7.*A1+B2)*BTG                     &
     &                 ,BDNM= 6.*A1*A1                                  &
     &                 ,BEQH= A2X*B1*BTG+3.*A2X*(7.*A1+B2)*BTG          &
     &                 ,BEQM=-A1*B1*(1.-3.*C1)+6.*A1*A1                 &
     &                 ,BNMH=-A2X*BTG                                   &
     &                 ,BNMM=A1*(1.-3.*C1)                              &
     &                 ,BSHH=9.*A1*A2X*A2X*BTG                          &
     &                 ,BSHM=18.*A1*A1*A2X*C1                           &
     &                 ,BSMH=-3.*A1*A2X*(3.*A2X+3.*B2*C1+12.*A1*C1-B2)  &
     &                                *BTG                              &
     &                 ,CESH=A2X                                        &
     &                 ,CESM=A1*(1.-3.*C1)                              &
     &                 ,CNV=EP_1*G/BTG                                  &
     &                 ,ELFCS=VKARMAN*BTG                               &
     &                 ,FZQ1=RTVISC*QVISC*ZQRZT                         &
     &                 ,FZQ2=RTVISC*QVISC*ZQRZT                         &
     &                 ,FZT1=RVISC *TVISC*SQPR                          &
     &                 ,FZT2=CZIV*GRRS*TVISC*SQPR                       &
     &                 ,FZU1=CZIV*VISC                                  &
     &                 ,PIHF=0.5*PI                                     &
     &                 ,RFAC=RIC/(FHNEU*RFC*RFC)                        &
     &                 ,RQVISC=1./QVISC                                 &
     &                 ,RRIC=1./RIC                                     &
     &                 ,USTFC=0.018/G                                   &
     &                 ,WNEW=1.-WOLD                                    &
     &                 ,WWST2=WWST*WWST
!
!-----------------------------------------------------------------------
!***  FREE TERM IN THE EQUILIBRIUM EQUATION FOR (L/Q)**2
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: AEQH=9.*A1*A2X*A2X*B1*BTG*BTG                   &
     &                      +9.*A1*A2X*A2X*(12.*A1+3.*B2)*BTG*BTG       &
     &                 ,AEQM=3.*A1*A2X*B1*(3.*A2X+3.*B2*C1+18.*A1*C1-B2)&
     &                      *BTG+18.*A1*A1*A2X*(B2-3.*A2X)*BTG
!
!-----------------------------------------------------------------------
!***  FORBIDDEN TURBULENCE AREA
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: REQU=-AEQH/AEQM                                 &
     &                 ,EPSGH=1.E-9,EPSGM=REQU*EPSGH
!
!-----------------------------------------------------------------------
!***  NEAR ISOTROPY FOR SHEAR TURBULENCE, WW/Q2 LOWER LIMIT
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: UBRYL=(18.*REQU*A1*A1*A2X*B2*C1*BTG             &
     &                         +9.*A1*A2X*A2X*B2*BTG*BTG)               &
     &                        /(REQU*ADNM+ADNH)                         &
     &                 ,UBRY=(1.+EPSRS)*UBRYL,UBRY3=3.*UBRY
!
      REAL,PARAMETER :: AUBH=27.*A1*A2X*A2X*B2*BTG*BTG-ADNH*UBRY3       &
     &                 ,AUBM=54.*A1*A1*A2X*B2*C1*BTG -ADNM*UBRY3        &
     &                 ,BUBH=(9.*A1*A2X+3.*A2X*B2)*BTG-BDNH*UBRY3       &
     &                 ,BUBM=18.*A1*A1*C1           -BDNM*UBRY3         &
     &                 ,CUBR=1.                     -     UBRY3         &
     &                 ,RCUBR=1./CUBR
!
!-----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                          SUBROUTINE MIXLEN                            &
!----------------------------------------------------------------------
!   ******************************************************************
!   *                                                                *
!   *                   LEVEL 2.5 MIXING LENGTH                      *
!   *                                                                *
!   ******************************************************************
!
#ifdef MLNOH
     &(LMH,U,V,T,THE,Q,CWM,Q2,Z,USTAR,CORF,EPSHOL                      &
#else
     &(LMH,U,V,T,THE,Q,CWM,Q2,Z,USTAR,CORF                             &
#endif
     &,S2,GH,RI,EL,PBLH,LPBL,LMXL,CT                                   &
#ifdef GAMMA_EL
     &,HGAMU,HGAMV,HGAMT,PBLFLG                                        &
     &,ZFACENTK,UFXPBL,VFXPBL,HFXPBL                                   &
#endif
     &,IDS,IDE,JDS,JDE,KDS,KDE                                         &
     &,IMS,IME,JMS,JME,KMS,KME                                         &
     &,ITS,ITE,JTS,JTE,KTS,KTE)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: LMH
!
      INTEGER,INTENT(IN) :: LMXL,LPBL
!
      REAL,DIMENSION(KTS:KTE),INTENT(IN) :: CWM,Q,Q2,T,THE,U,V
!
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
!
      REAL,INTENT(IN) :: PBLH
!
      REAL,DIMENSION(KTS+1:KTE),INTENT(OUT) :: EL,RI,GH,S2
!
      REAL,INTENT(INOUT) :: CT
!
      REAL,INTENT(IN)    :: CORF,USTAR
#ifdef MLNOH
      REAL,INTENT(INOUT) :: EPSHOL
#endif
#ifdef GAMMA_EL
!
      REAL,DIMENSION(KTS+1:KTE),INTENT(IN) :: ZFACENTK
      REAL,INTENT(IN) :: HGAMU,HGAMV,HGAMT,UFXPBL,VFXPBL,HFXPBL
      REAL :: SUK,SVK
!
      LOGICAL,INTENT(IN) :: PBLFLG
#endif
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K,LPBLM
!
      REAL :: A,ADEN,B,BDEN,AUBR,BUBR,BLMX,EL0,ELOQ2X,GHL,S2L          &
     &       ,QOL2ST,QOL2UN,QDZL,RDZ,SQ,SREL,SZQ,TEM,THM,VKRMZ,RLAMBDA &
     &       ,RLB,RLN,F
!
      REAL,DIMENSION(KTS:KTE) :: Q1,EN2
!
      REAL,DIMENSION(KTS+1:KTE) :: DTH,ELM,REL
#ifdef MLNOH
      REAL,PARAMETER :: ELCBL=0.77
      REAL :: CKP
#endif
!----------------------------------------------------------------------
!**********************************************************************
!-----------------------------------------------------------------------
      DO K=KTS,KTE
        Q1(K)=0.
      ENDDO
!
      DO K=KTS+1,KTE
        DTH(K)=THE(K)-THE(K-1)
      ENDDO
!
      DO K=KTS+2,KTE
        IF(DTH(K)>0..AND.DTH(K-1)<=0.)THEN
          DTH(K)=DTH(K)+CT
          EXIT
        ENDIF
      ENDDO
!
      CT=0.
!----------------------------------------------------------------------
!***  COMPUTE LOCAL GRADIENT RICHARDSON NUMBER
!----------------------------------------------------------------------
      DO K=KTE,KTS+1,-1
        RDZ=2./(Z(K+1)-Z(K-1))
        S2L=((U(K)-U(K-1))**2+(V(K)-V(K-1))**2)*RDZ*RDZ   ! S**2
#ifdef GAMMA_EL
#ifdef GAMMA_PBPS
        IF(PBLFLG.AND.K.LE.LPBL)THEN
          SUK=(U(K)-U(K-1))*RDZ
          SVK=(V(K)-V(K-1))*RDZ
          S2L=(SUK-HGAMU/PBLH)*SUK+(SVK-HGAMV/PBLH)*SVK
        ENDIF
#endif
#endif
        S2L=MAX(S2L,EPSGM)
        S2(K)=S2L
!
        TEM=(T(K)+T(K-1))*0.5
        THM=(THE(K)+THE(K-1))*0.5
!
        A=THM*P608
        B=(ELOCP/TEM-1.-P608)*THM
!
        GHL=(DTH(K)*((Q(K)+Q(K-1)+CWM(K)+CWM(K-1))*(0.5*P608)+1.)      &
     &     +(Q(K)-Q(K-1)+CWM(K)-CWM(K-1))*A                            &
     &     +(CWM(K)-CWM(K-1))*B)*RDZ                       ! dTheta/dz
#ifdef GAMMA_EL
#ifdef GAMMA_PBPS
        IF(PBLFLG.AND.K.LE.LPBL)THEN
          GHL=GHL-HGAMT/PBLH
        ENDIF
#endif
#endif
!
        IF(ABS(GHL)<=EPSGH)GHL=EPSGH
!
        EN2(K)=GHL*G/THM                                   ! N**2
!
        GH(K)=GHL
        RI(K)=EN2(K)/S2L
      ENDDO
!
!----------------------------------------------------------------------
!***  FIND MAXIMUM MIXING LENGTHS AND THE LEVEL OF THE PBL TOP
!----------------------------------------------------------------------
!
      DO K=KTE,KTS+1,-1
        S2L=S2(K)
        GHL=GH(K)
!
        IF(GHL>=EPSGH)THEN
          IF(S2L/GHL<=REQU)THEN
            ELM(K)=EPSL
          ELSE
            AUBR=(AUBM*S2L+AUBH*GHL)*GHL
            BUBR= BUBM*S2L+BUBH*GHL
            QOL2ST=(-0.5*BUBR+SQRT(BUBR*BUBR*0.25-AUBR*CUBR))*RCUBR
            ELOQ2X=1./QOL2ST
            ELM(K)=MAX(SQRT(ELOQ2X*Q2(K)),EPSL)
          ENDIF
        ELSE
          ADEN=(ADNM*S2L+ADNH*GHL)*GHL
          BDEN= BDNM*S2L+BDNH*GHL
          QOL2UN=-0.5*BDEN+SQRT(BDEN*BDEN*0.25-ADEN)
          ELOQ2X=1./(QOL2UN+EPSRU)       ! repsr1/qol2un
          ELM(K)=MAX(SQRT(ELOQ2X*Q2(K)),EPSL)
        ENDIF
      ENDDO
!
!----------------------------------------------------------------------
      DO K=LPBL,LMH,-1
        Q1(K)=SQRT(Q2(K))
      ENDDO
!----------------------------------------------------------------------
      SZQ=0.
      SQ =0.
!
      DO K=KTE,KTS+1,-1
        QDZL=(Q1(K)+Q1(K-1))*(Z(K)-Z(K-1))
        SZQ=(Z(K)+Z(K-1)-Z(LMH)-Z(LMH))*QDZL+SZQ
        SQ=QDZL+SQ
      ENDDO
!
!----------------------------------------------------------------------
!***  COMPUTATION OF ASYMPTOTIC L IN BLACKADAR FORMULA
!----------------------------------------------------------------------
!
      EL0=MIN(ALPH*SZQ*0.5/SQ,EL0MAX)
      EL0=MAX(EL0            ,EL0MIN)
!
!----------------------------------------------------------------------
!***  ABOVE THE PBL TOP
!----------------------------------------------------------------------
!
      LPBLM=MIN(LPBL+1,KTE)
!
      DO K=KTE,LPBLM,-1
        EL(K)=(Z(K+1)-Z(K-1))*ELFC
        REL(K)=EL(K)/ELM(K)
      ENDDO
!
!----------------------------------------------------------------------
!***  INSIDE THE PBL
!----------------------------------------------------------------------
!
#ifdef MLNOH
      EPSHOL=MIN(EPSHOL,0.0)
      CKP=ELCBL*((1.0-8.0*EPSHOL)**(1./3.))
#endif
      IF(LPBL>LMH)THEN
        DO K=LPBL,LMH+1,-1
          VKRMZ=(Z(K)-Z(LMH))*VKARMAN
#ifdef MLNOH
          IF(PBLFLG) THEN
            VKRMZ=CKP*(Z(K)-Z(LMH))*VKARMAN
            EL(K)=VKRMZ/(VKRMZ/EL0+1.)
          ELSE
            EL(K)=VKRMZ/(VKRMZ/EL0+1.)
          ENDIF
#else
          EL(K)=VKRMZ/(VKRMZ/EL0+1.)
#endif
          REL(K)=EL(K)/ELM(K)
        ENDDO
      ENDIF
!
      DO K=LPBL-1,LMH+2,-1
        SREL=MIN(((REL(K-1)+REL(K+1))*0.5+REL(K))*0.5,REL(K))
        EL(K)=MAX(SREL*ELM(K),EPSL)
      ENDDO
!
!----------------------------------------------------------------------
!***  MIXING LENGTH FOR THE QNSE MODEL IN STABLE CASE
!----------------------------------------------------------------------
!
      F=MAX(CORF,EPS1)
      RLAMBDA=F/(BLCKDR*USTAR)
      DO K=KTE,KTS+1,-1
        IF(EN2(K)>=0.0)THEN                          ! Stable case
          VKRMZ=(Z(K)-Z(LMH))*VKARMAN
          RLB=RLAMBDA+1./VKRMZ
          RLN=SQRT(2.*EN2(K)/Q2(K))/CN
!         EL(K)=MIN(1./(RLB+RLN),ELM(K))
          EL(K)=1./(RLB+RLN)
        ENDIF
      ENDDO
!
!----------------------------------------------------------------------
      END SUBROUTINE MIXLEN
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                          SUBROUTINE PRODQ2                            &
!----------------------------------------------------------------------
!   ******************************************************************
!   *                                                                *
!   *            LEVEL 2.5 Q2 PRODUCTION/DISSIPATION                 *
!   *                                                                *
!   ******************************************************************
!
     &(LMH,DTTURBL,USTAR,S2,RI,Q2,EL,Z,AKM,AKH                         &
#ifdef GAMMA_PBPS
     &,UXK,VXK,THXK,THVXK                                              &
     &,HGAMU,HGAMV,HGAMT                                               &
     &,HPBL,PBLFLG,KPBL                                                &
#endif
#ifdef ENT_PBPS
     &,ZFACENTK,UFXPBL,VFXPBL,HFXPBL                                   &
#endif
#ifdef TKE_OUTPUT
     &,PS1D,PB1D,EPS1D                                                 &
#endif
     &,IDS,IDE,JDS,JDE,KDS,KDE                                         &
     &,IMS,IME,JMS,JME,KMS,KME                                         &
     &,ITS,ITE,JTS,JTE,KTS,KTE)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: LMH
!
      REAL,INTENT(IN) :: DTTURBL,USTAR
!
      REAL,DIMENSION(KTS+1:KTE),INTENT(IN) :: S2,RI,AKM,AKH,EL
!
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
!
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: Q2
#ifdef GAMMA_PBPS
!
      REAL,DIMENSION(KTS:KTE),INTENT(IN) :: UXK,VXK,THXK,THVXK
      REAL,INTENT(IN) :: HGAMU,HGAMV,HGAMT,HPBL
!
      INTEGER,INTENT(IN) :: KPBL
      LOGICAL,INTENT(IN) :: PBLFLG
!
#endif
#ifdef ENT_PBPS
      REAL,DIMENSION(KTS+1:KTE),INTENT(IN) :: ZFACENTK
      REAL,INTENT(IN) :: UFXPBL,VFXPBL,HFXPBL
!
#endif
#ifdef TKE_OUTPUT
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: PS1D,PB1D,EPS1D
!
#endif
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K
!
      REAL :: S2L,Q2L,DELTAZ,AKML,AKHL,EN2,PR,BPR,DIS,RC02
#ifdef GAMMA_PBPS
      REAL :: SUK,SVK,GTHVK,GOVRTHVK,PRU,PRV
#endif
#ifdef ENT_PBPS
      REAL :: ZFACENTL
#endif
#ifdef TKE_OUTPUT
      REAL :: Q2_SFC
#endif
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
      RC02=2.0/(C0*C0)
      main_integration: DO K=KTS+1,KTE
        DELTAZ=0.5*(Z(K+1)-Z(K-1))
        S2L=S2(K)
        Q2L=Q2(K)
#ifdef GAMMA_PBPS
        SUK=(UXK(K)-UXK(K-1))/DELTAZ
        SVK=(VXK(K)-VXK(K-1))/DELTAZ
        GTHVK=(THVXK(K)-THVXK(K-1))/DELTAZ
        GOVRTHVK=G/(0.5*(THVXK(K)+THVXK(K-1)))
#endif
        AKML=AKM(K)
        AKHL=AKH(K)
#ifdef ENT_PBPS
        ZFACENTL=ZFACENTK(K)
#endif
        EN2=RI(K)*S2L                                           !N**2
!
!***  TURBULENCE PRODUCTION TERM
!
#ifdef GAMMA_PBPS
        IF(PBLFLG.AND.K.LE.KPBL)THEN
          PRU=(AKML*(SUK-HGAMU/HPBL))*SUK
          PRV=(AKML*(SVK-HGAMV/HPBL))*SVK
#ifdef ENT_PBPS
          PRU=(AKML*(SUK-HGAMU/HPBL)-UFXPBL*ZFACENTL)*SUK
          PRV=(AKML*(SVK-HGAMV/HPBL)-VFXPBL*ZFACENTL)*SVK
#endif
        ELSE
          PRU=AKML*SUK*SUK
          PRV=AKML*SVK*SVK
        ENDIF
        PR=PRU+PRV
#else
        PR=AKML*S2L
#endif
#ifdef TKE_OUTPUT
        PS1D(K)=PR
#endif
!
!***  BUOYANCY PRODUCTION
!
#ifdef GAMMA_PBPS
        IF(PBLFLG.AND.K.LE.KPBL)THEN
          BPR=(AKHL*(GTHVK-HGAMT/HPBL))*GOVRTHVK
#ifdef ENT_PBPS
          BPR=(AKHL*(GTHVK-HGAMT/HPBL)-HFXPBL*ZFACENTL)*GOVRTHVK
#endif
        ELSE
          BPR=AKHL*GTHVK*GOVRTHVK
        ENDIF
#else
        BPR=AKHL*EN2
#endif
#ifdef TKE_OUTPUT
        PB1D(K)=-BPR
#endif
!
!***  DISSIPATION
!
        DIS=CEPS*(0.5*Q2L)**1.5/EL(K)
#ifdef TKE_OUTPUT
        EPS1D(K)=-DIS
#endif
!
        Q2L=Q2L+2.0*(PR-BPR-DIS)*DTTURBL
        Q2(K)=AMAX1(Q2L,EPSQ2L)
!----------------------------------------------------------------------
!***  END OF PRODUCTION/DISSIPATION LOOP
!----------------------------------------------------------------------
!
      ENDDO main_integration
!
!----------------------------------------------------------------------
!***  LOWER BOUNDARY CONDITION FOR Q2
!----------------------------------------------------------------------
!
      Q2(KTS)=AMAX1(RC02*USTAR*USTAR,EPSQ2L)
!----------------------------------------------------------------------
!
      END SUBROUTINE PRODQ2
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                           SUBROUTINE VDIFQ                            &
!   ******************************************************************
!   *                                                                *
!   *               VERTICAL DIFFUSION OF Q2 (TKE)                   *
!   *                                                                *
!   ******************************************************************
     &(LMH,DTDIF,Q2,EL,Z                                               &
#ifdef KTKE_YSU
     &,AKHK                                                            &
#endif
#ifdef TKE_OUTPUT
     &,PT1D                                                            &
#endif
#ifdef GAMMA_PT
     &,HGAME,HPBL,PBLFLG,KPBL                                          &
#endif
#ifdef ENT_PT
     &,EFXPBL                                                          &
#endif
#ifdef VDIFQ_OUTPUT
     &,XKZE1D,EFLX_L1D,EFLX_NL1D                                       &
#endif
     &,IDS,IDE,JDS,JDE,KDS,KDE                                         &
     &,IMS,IME,JMS,JME,KMS,KME                                         &
     &,ITS,ITE,JTS,JTE,KTS,KTE)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: LMH
!
      REAL,INTENT(IN) :: DTDIF
!
      REAL,DIMENSION(KTS+1:KTE),INTENT(IN) :: EL
#ifdef KTKE_YSU
      REAL,DIMENSION(KTS+1:KTE),INTENT(IN) :: AKHK
#endif
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
!
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: Q2
#ifdef TKE_OUTPUT
!
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: PT1D
#endif
#ifdef GAMMA_PT
!
      REAL,DIMENSION(KTS:KTE),INTENT(IN)   :: HGAME
      REAL,INTENT(IN) :: HPBL
      INTEGER,INTENT(IN) :: KPBL
      LOGICAL,INTENT(IN) :: PBLFLG
!
#endif
#ifdef ENT_PT
      REAL,INTENT(IN) :: EFXPBL
!
#endif
#ifdef VDIFQ_OUTPUT
!
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: XKZE1D,EFLX_L1D,EFLX_NL1D
#endif
!
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K
!
      REAL :: ADEN,AKQS,BDEN,BESH,BESM,CDEN,CF,DTOZS,ELL,ELOQ2,ELOQ4   &
     &       ,ELQDZ,ESH,ESM,ESQHF,GHL,GML,Q1L,RDEN,RDZ
#ifdef ENT_PT
      REAL :: ZAK
#endif
!
      REAL,DIMENSION(KTS+2:KTE) :: AKQ,CM,CR,DTOZ,RSQ2
#ifdef ENT_PT
      REAL,DIMENSION(KTS+1:KTE) :: ZFACENTK
#endif
#ifdef TKE_OUTPUT
!
      REAL,DIMENSION(KTS:KTE) :: Q2_SAVE
#endif
!
      REAL,PARAMETER :: C_K=1.0
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!***
!***  VERTICAL TURBULENT DIFFUSION
!***
!----------------------------------------------------------------------
      ESQHF=0.5*ESQ
#ifdef TKE_OUTPUT
!
      DO K=KTS,KTE
        Q2_SAVE(K)=Q2(K)
      ENDDO
#endif
#ifdef ENT_PT
      DO K=KTS+1,KTE
        ZAK=0.5*(Z(K)+Z(K-1)) !ZAK OF VDIFQ = ZA(K-1) OF YSU2D
        ZFACENTK(K)=(ZAK/HPBL)**3.0
      ENDDO
#endif
!
      DO K=KTE,KTS+2,-1
        DTOZ(K)=(DTDIF+DTDIF)/(Z(K+1)-Z(K-1))
        AKQ(K)=C_K*(AKHK(K)/(Z(K+1)-Z(K-1))+AKHK(K-1)/(Z(K)-Z(K-2)))
        CR(K)=-DTOZ(K)*AKQ(K)
      ENDDO
!
        AKQS=C_K*AKHK(KTS+1)/(Z(KTS+2)-Z(KTS))
#ifdef VDIFQ_OUTPUT
        XKZE1D(KTS)=AKQS*(Z(KTS+1)-Z(KTS))
      DO K=3,KTE
        XKZE1D(K-1)=AKQ(K)*(Z(K)-Z(K-1))
      ENDDO
      DO K=KTS,KTE-1
        EFLX_L1D(K)=-0.5*XKZE1D(K)/(Z(K+1)-Z(K))*(Q2(K+1)-Q2(K))
#ifdef GAMMA_PT
        IF(PBLFLG.AND.K.LT.KPBL) THEN
          EFLX_NL1D(K)=XKZE1D(K)*HGAME(K)/HPBL
        ENDIF
#endif
      ENDDO
!
#endif
      CM(KTE)=DTOZ(KTE)*AKQ(KTE)+1.
      RSQ2(KTE)=Q2(KTE)
!
      DO K=KTE-1,KTS+2,-1
        CF=-DTOZ(K)*AKQ(K+1)/CM(K+1)
        CM(K)=-CR(K+1)*CF+(AKQ(K+1)+AKQ(K))*DTOZ(K)+1.
        RSQ2(K)=-RSQ2(K+1)*CF+Q2(K)
#ifdef GAMMA_PT
      IF(PBLFLG.AND.K.LT.KPBL) THEN
        RSQ2(K)=RSQ2(K)-DTOZ(K)*(2.0*HGAME(K)/HPBL)*AKQ(K+1)*(Z(K+1)-Z(K)) &
                       +DTOZ(K)*(2.0*HGAME(K-1)/HPBL)*AKQ(K)*(Z(K)-Z(K-1))
#ifdef ENT_PT
        RSQ2(K)=RSQ2(K)-DTOZ(K)*2.0*EFXPBL*ZFACENTK(K+1)                   &
                       +DTOZ(K)*2.0*EFXPBL*ZFACENTK(K)
#endif
      ENDIF
#endif
      ENDDO
!
      DTOZS=(DTDIF+DTDIF)/(Z(KTS+2)-Z(KTS))
      CF=-DTOZS*AKQ(LMH+2)/CM(LMH+2)
!
#ifndef GAMMA_PT
      Q2(LMH+1)=(DTOZS*AKQS*Q2(LMH)-RSQ2(LMH+2)*CF+Q2(LMH+1))          &
     &        /((AKQ(LMH+2)+AKQS)*DTOZS-CR(LMH+2)*CF+1.)
#else
      IF(PBLFLG.AND.((LMH+1).LT.KPBL)) THEN
        Q2(LMH+1)=(DTOZS*AKQS*Q2(LMH)-RSQ2(LMH+2)*CF+Q2(LMH+1)         &
            -DTOZS*(2.0*HGAME(LMH+1)/HPBL)*AKQ(LMH+2)*(Z(LMH+2)-Z(LMH+1)) &
            +DTOZS*(2.0*HGAME(LMH)/HPBL)*AKQS*(Z(LMH+1)-Z(LMH)))
#ifdef ENT_PT
        Q2(LMH+1)=Q2(LMH+1)-DTOZS*2.0*EFXPBL*ZFACENTK(LMH+2)           &
                           +DTOZS*2.0*EFXPBL*ZFACENTK(LMH+1)
#endif
        Q2(LMH+1)=Q2(LMH+1)/((AKQ(LMH+2)+AKQS)*DTOZS-CR(LMH+2)*CF+1.)
      ELSE
        Q2(LMH+1)=(DTOZS*AKQS*Q2(LMH)-RSQ2(LMH+2)*CF+Q2(LMH+1))        &
       &        /((AKQ(LMH+2)+AKQS)*DTOZS-CR(LMH+2)*CF+1.)
      ENDIF
#endif
!
      DO K=LMH+2,KTE
        Q2(K)=(-CR(K)*Q2(K-1)+RSQ2(K))/CM(K)
      ENDDO
#ifdef TKE_OUTPUT
!
#ifndef GAMMA_PT
      DO K=KTS+1,KTE
        PT1D(K)=(Q2(K)-Q2_SAVE(K))/(2.0*DTDIF)
      ENDDO
#else
      DO K=KTS+1,KTE
        PT1D(K)=(Q2(K)-Q2_SAVE(K))/(2.0*DTDIF)
      ENDDO
#endif
#endif
!----------------------------------------------------------------------
!
      END SUBROUTINE VDIFQ
!
!----------------------------------------------------------------------
END MODULE module_ysutke
