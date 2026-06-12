#include "define.h"
   SUBROUTINE DFS_BARO_TEST_INIT( VOR, DIO, TAI, PRS, SFAI )
!-------------------------------------------------------------------------------
   use constant, only : pi_, rd_,rerth_,omega_,g_
   use dfsvar, only   : ibA,IB1A,jgs,levs,jge,jbwA,sigma,ib,jbw,iope
!-------------------------------------------------------------------------------
!                                                                    
! abstract :    Baroclinic-wave test case                                      
!                                                                    
!     Mon. Wea. Rev. Vol 132, 133-153 (2004)                         
!     by F. X. Giraldo and T. E. Rosmond                             
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!                                                                    
!-------------------------------------------------------------------------------
   implicit none
!
#ifdef MP
   real     ::  VOR(ib,jbw,levs), DIO(ib,jbw,levs)
   real     ::  TAI(ib,jbw,levs),SFAI(ib,jbw)
   real     ::  PRS(ib,jbw)
!
   real     ::  VORA(0:IB1A,jgs:jge,levs), DIOA(0:IB1A,jgs:jge,levs)
   real     ::  TAIA(0:IB1A,jgs:jge,levs),SFAIA(0:IB1A,jgs:jge)
   real     ::  PRSA(0:IB1A,jgs:jge)
#define VOR VORA
#define DIO DIOA
#define TAI TAIA
#define PRS PRSA
#define SFAI SFAIA
#else
   real     ::  VOR(0:IB1A,jgs:jge,levs), DIO(0:IB1A,jgs:jge,levs)
   real     ::  TAI(0:IB1A,jgs:jge,levs),SFAI(0:IB1A,jgs:jge)
   real     ::  PRS(0:IB1A,jgs:jge)
#endif
!
   INTEGER  ::  k,j,i
   real     ::  cosalong,a_in,rr,rra,drdlon,exppp,u,dudlon,dudlat,si6,Ac
   real     ::  Bc,firs,seco,thir,c23m,along,duzodlat,uzo,RCOSr,sinalati2
   real     ::  sinalong,alati,sigma_v,COSalati,SINalati,cosalati2,drdlat
   real     ::  SINalati_c,along_c,COSalong_c,SINalong_c,alati_c,COSalati_c
   real     ::  RGASd,u0,sigma_0,dlat
   real     ::  pi,pi2,radius,aOmega,Omega,aOmega2,gravi
!
!------ surface pressure ------------------------------
!
   do j = jgs,jge
     do i = 0,IB1A
       PRS(i,j)= 0.0d0        ! surface pressure=1000 hPa
     end do 
   end do 
#ifdef MP
#undef PRS
   call mpgf2pd(PRSA,iba,jbwa,PRS,ib,jbw,1)
#endif
!
!------ Vorticity and Divergence, Surface geopotential-
   if (iope) then
!
     PI    = pi_
     PI2   = PI*2
     radius = rerth_
     Omega = omega_
     gravi = g_
     aOmega = radius*Omega            ! Earth-Parameters
     aOmega2= aOmega**2               ! Earth-Parameters
!
     RGASd= rd_
     RR= radius/10.d0
     u0= 35.d0
     sigma_0= 0.252d0
     alati_c= PI2/9.d0
#ifdef IBMSP
     COSalati_c= COS(dble(alati_c))
     SINalati_c= SIN(dble(alati_c))
#else
     COSalati_c= DCOS(dble(alati_c))
     SINalati_c= DSIN(dble(alati_c))
#endif
     along_c= PI /9.d0
#ifdef IBMSP
     COSalong_c= COS(dble(along_c))
     SINalong_c= SIN(dble(along_c))
#else
     COSalong_c= DCOS(dble(along_c))
     SINalong_c= DSIN(dble(along_c))
#endif
     dlat=PI/JBWA
     do k = 1,levs
       sigma_v=(sigma(k)-sigma_0)*PI/2.d0
       do j = jgs,jge
         alati = 0.5*PI-dlat*(jgs-j+0.5d0)
#ifdef IBMSP
         COSalati = COS(dble(alati))
         SINalati = SIN(dble(alati))
         cosalati2= COS(dble(alati*2.d0))
         sinalati2= SIN(dble(alati*2.d0))
#else
         COSalati = DCOS(dble(alati))
         SINalati = DSIN(dble(alati))
         cosalati2= DCOS(dble(alati*2.d0))
         sinalati2= DSIN(dble(alati*2.d0))
#endif
         RCOSr= 1.d0/(COSalati*radius)
#ifdef IBMSP
         uzo= u0* COS(dble(sigma_v))**1.5d0 * sinalati2**2.d0
         duzodlat= u0* COS(dble(sigma_v))**1.5d0 * sinalati2*cosalati2 *4.d0
#else
         uzo= u0* DCOS(dble(sigma_v))**1.5d0 * sinalati2**2.d0
         duzodlat= u0* DCOS(dble(sigma_v))**1.5d0 * sinalati2*cosalati2 *4.d0
#endif
         do i = 0,IB1A
           along= (PI2/IBA)*(i      )
#ifdef IBMSP
           sinalong= SIN(dble(along-along_c))
           cosalong= COS(dble(along-along_c))
#else
           sinalong= DSIN(dble(along-along_c))
           cosalong= DCOS(dble(along-along_c))
#endif
!
           VOR(i,j,k)= -duzodlat/radius + RCOSr*sinalati*uzo  ! zonal-mean
           DIO(i,j,k)=  0.d0                                  ! zonal-mean
!
           a_in= SINalati_c*SINalati                                           &
                +COSalati_c*COSalati*cosalong
#ifdef IBMSP
           rr = radius* ACOS(dble(a_in))
#else
           rr = radius* DACOS(dble(a_in))
#endif
           rra= rr/radius
#ifdef IBMSP
           drdlon= ( COSalati_c*COSalati*sinalong )                            &
                 / ( (1.d0/radius)*SIN(dble(rra)) )
           drdlat= ( SINalati_c*COSalati                                       &
                  -COSalati_c*SINalati*cosalong )                              &
                 / (-(1.d0/radius)*SIN(dble(rra)) )
#else
           drdlon= ( COSalati_c*COSalati*sinalong )                            &
                 / ( (1.d0/radius)*DSIN(dble(rra)) )
           drdlat= ( SINalati_c*COSalati                                       &
                  -COSalati_c*SINalati*cosalong )                              &
                 / (-(1.d0/radius)*DSIN(dble(rra)) )
#endif
!
           u= DEXP(-(rr/RR)**2.d0)
           exppp= u*(-2)*(rr/RR)*(1.d0/RR)
           dudlon= exppp*drdlon
           dudlat= exppp*drdlat 
!
!------ Perturbed Vor, Diverg---
!
           VOR(i,j,k)= VOR(i,j,k)+(-dudlat/radius+RCOSr*SINalati*u)
           DIO(i,j,k)= DIO(i,j,k)+( dudlon*RCOSr  )
!
           si6= (SINalati**2)**3.d0
           Ac= -2.d0*si6*(COSalati**2+1/3.d0)+10.d0/63.d0
           Bc= aOmega*(1.6d0*COSalati**3.d0*(SINalati**2+2/3.d0)-PI/4)
!
!------ Perturbed temperature---
!
           firs= 0.75d0*sigma(k)*(PI*u0/RGASd)
#ifdef IBMSP
           seco= SIN(dble(sigma_v))*( COS(dble(sigma_v))**0.5d0 )
           thir= 2.d0*u0*( COS(dble(sigma_v))**1.5d0 )*Ac + Bc
#else
           seco= DSIN(dble(sigma_v))*( DCOS(dble(sigma_v))**0.5d0 )
           thir= 2.d0*u0*( DCOS(dble(sigma_v))**1.5d0 )*Ac + Bc
#endif
!
           TAI(i,j,k)= firs*seco*thir
!
!------ Surface geopotential ---
!
           if(k.eq.1) then
#ifdef IBMSP
             c23m= ( COS(dble( (1-sigma_0)*(PI/2.d0) )) )**1.5d0
#else
             c23m= ( DCOS(dble( (1-sigma_0)*(PI/2.d0) )) )**1.5d0
#endif
!
             SFAI(i,j)= u0*c23m*( u0*c23m* Ac + Bc ) 
           endif
!
         end do ! i 
       end do ! j
     end do ! k
!
!-----nondimensionalize--------------------
!
!    do k = 1,levs
!      do j = jgs,jge
!        do i = 0,IB1A
!          VOR(i,j,k)= VOR(i,j,k)/Omega 
!          DIO(i,j,k)= DIO(i,j,k)/Omega 
!          TAI(i,j,k)= TAI(i,j,k)/TSCALE
!          if(k.eq.1) SFAI(i,j) = SFAI(i,j) /FSCALE
!        end do ! i 
!      end do ! j
!    end do ! k
!
   endif         ! iope
#ifdef MP
#undef VOR
#undef DIO
#undef TAI
#undef SFAI
   call mpgf2pd(VORA,iba,jbwa,VOR,ib,jbw,levs)
   call mpgf2pd(DIOA,iba,jbwa,DIO,ib,jbw,levs)
   call mpgf2pd(TAIA,iba,jbwa,TAI,ib,jbw,levs)
   call mpgf2pd(SFAIA,iba,jbwa,SFAI,ib,jbw,1)
#endif
!
   RETURN
   END
