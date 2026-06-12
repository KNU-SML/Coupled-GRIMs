#include "define.h"
   SUBROUTINE dfs_ideal_temp ( PSL2 )
!-------------------------------------------------------------------------------
   use constant, only : rd_,g_,cp_,pi_
   use dfsvar, only : mt,jl,ib,ib1,jb,jb1,levsp,sigma,                         &
#ifdef HDSZ
                      coslat, T_REF, T_VIS, Z_VIS
#else
                      coslat
#endif
!-------------------------------------------------------------------------------
!                                                                    
! abstract :    Reference Temperature field ; Held-Suarez (1994, BAMS)         
!                                           ; Williamson et al.(1998, MWR)     
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
!-------------------------------------------------------------------------------
   real     ::  PSL2(MT,JL)
#ifdef HDSZ
!---local-variables--------------------------------------------------!
   real     ::  SUFPRS(IB,JBw)
   real     ::  T0,gammad,gammai,conA,cos2,theta,t_rad,sin2,p_i,gravit,Rgas
   real     ::  Rgasd,Rgasi,p_d,akap,delTh,phi0,delphi,peq2,prs,p,pp_d,pp_i
   real     ::  zss,sss,ass,cos4,amax,arigh,aleft,peq,ppl,delT,aaa
   integer  ::  i,j,k
!--------------------------------------------------------------------!
!
   CALL dfs_fft_driver( -1, SUFPRS, PSL2, JL, +1 )
!
!------------------------------temperature relaxation-------
!
!--beg--Held-Suarez-Williamson profile
!
   T0= 200.d0
   gammad= 2.000d0     ! [K/Km]
   gammai=-3.345d0
   gravit= g_          ! [m/s^2]
   Rgas  = rd_  /1000.d0        ! /Km -> /m
   Rgasd = Rgas*gammad/gravit
   Rgasi = Rgas*gammai/gravit
!
!-----
!         p_d=  0.00001d0/1000.d0
!         ppl=    2.d-13/1000.d0
!
   p_d= 100.d0/1000.d0
   ppl=   2.d0/1000.d0
   peq= p_d
   delphi= 15.d0
   phi0  = 60.d0
   conA  = 2.65d0/delphi
   delT  = 60.d0
   delTh = 10.d0
   akap= rd_ /cp_
!
!------
!
   do k = 1,levsp
     do j = -JB,JB-1
       theta= (90.d0/JB)*(j+0.5)
       t_rad= theta*pi_/180.d0
       cos2 = dcos(t_rad)**2
       sin2 = dsin(t_rad)**2
       p_i = peq-(peq-ppl)*(1.d0+dtanh(conA*(dabs(theta)-phi0)))/2.d0
       do i = 0,IB1
!          prs= (psuf(k)+psuf(k+1))/2.d0
         prs= sigma(k)*DEXP(SUFPRS(i,j)) ! normalized by 1000 hPa
         p  = prs                        ! should be vectorized
         IF(p.LE.p_d) THEN  !--stratosphere---
           pp_d=  p/p_d
           pp_i=  p/p_i
           aleft= ( min(1.d0,pp_d) )**Rgasd
           arigh= ( min(1.d0,pp_i) )**Rgasi
           T_REF(i,j,k)= T0*aleft + T0*(arigh-1)
           T_REF(i,j,k)= T_REF(i,j,k)   ! /TSCALE
         ELSE               !--troposphere---
           aaa= 315.d0-delT*sin2-delTh*dlog10(p)*cos2
           aaa= aaa*p**akap
           T_REF(i,j,k)= max(200.d0,aaa)
           T_REF(i,j,k)= T_REF(i,j,k) ! /TSCALE
         ENDIF              !--temp-end----
       end do ! i
     end do ! j
   end do ! levsp
!
!--end--Held-Suarez-Williamson profile
!
!------------------------------cooling coefficient----------
!
   ass= 1.d0/(40.d0*pi_*2)
   sss= 1.d0/( 4.d0*pi_*2)
   zss= 1.d0/( 1.d0*pi_*2)
!
   do k = 1,levsp
     amax= max( 0.d0,(sigma(k)-0.7d0)/0.3d0 )
     Z_VIS(  k)= zss* amax
     do j = -JB,JB1
       cos4= COSLAT(J,1)**4
       T_VIS(j,k)= ass+(sss-ass)* amax *cos4
     end do ! j
   end do ! k
#endif
!
   RETURN
   END SUBROUTINE dfs_ideal_temp
