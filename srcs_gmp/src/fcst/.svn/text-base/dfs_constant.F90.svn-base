!
   SUBROUTINE dfs_constant()
!-------------------------------------------------------------------------------
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
   use constant, only : pi=>pi_,rerth_,Omega_,akap=>akapa_
   use dfsvar, only   : pi2,PI180,sigmafull,delsig,tsigma,                     &
                        sigma,levs,gama1,beta,valpha,mt,                       &
                        cooling,damp1,damp2,DISSI,gama2,GENER
!                       ZDMATmi_a_RE,ZDMATmi_a_IM,jlha,lap_dim
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real     ::  ay,ayr
   integer  ::  k,nk,n1,m,lapend
!-------------------------------------------------------------------------------
!
   PI2   = PI*2
   PI180 = PI/180.d00
!
!  LAP_DIM= 4
!  lapend=LAP_DIM
!  if (LAP_DIM.eq.4) lapend=8
!  allocate(ZDMATmi_a_RE(mt,3,0:jlha,2,lapend))
!  allocate(ZDMATmi_a_IM(mt,3,0:jlha,2,lapend))
!
! DELTA SIGMA 
!
   tsigma(1)= delsig(1)
!
   do k = 2,levs
     tsigma(k)= tsigma(k-1)+delsig(k)
   end do
!
! to avoid confusing alpha in MAT200S
!    ALPHA -> vALPHA here
!    vALPHA(rlev) 
!      BETA(rlev) 
!
   do k = 1,levs-1
     ay  = ( sigma(k+1)/sigma(k  ) ) ** akap
     ayr = ( sigma(k  )/sigma(k+1) ) ** akap
     valpha(k)   = ( ay -1.)/(2*akap)
     beta(k+1) = ( 1.-ayr)/(2*akap)
   end do
!
! GAMA1(rlev)
! GAMA2(rlev) 
!
   do k = 1,levs
     nk=1
     n1=1
     if(k.eq.  1 ) n1=0
     if(k.eq.levs) nk=0
     ay  =   sigma(k)/ sigma(k-n1)
     ayr =   sigma(k)/ sigma(k+nk)
     gama1(k)  =   ay  ** akap
     gama2(k)  =   ayr ** akap
   end do
!
! coefficient of NEWTONIAN COOLING and DISSIPATION
!
   DISSI   = 0./ ( (PI2* 3.0 )*(MT*(MT+1.)) )
   GENER   = 0./ (  PI2*15.0 )
   cooling = 0./ (  PI2*15.0 )
!
!      press   = 0.0
!
   damp1   = 0.050d0
   damp2   = 0.050d0
!
!      do k=levs,levs
!      do M=-MT,MT
!         EKMAN(M,k)= 1.d0/ (  PI2*10.0 )
!      end do
!      end do
!
   return
   end
