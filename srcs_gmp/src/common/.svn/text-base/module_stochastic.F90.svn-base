#include <define.h>
   module stochastic
!-------------------------------------------------------------------------------
#ifdef STOCH
   implicit none
   public stc_forcing
!-------------------------------------------------------------------------------
!
! arrays for random number
!
#if defined(NLT_U) || defined(DYN_U) || defined(PHY_U)
   real, allocatable, dimension(:,:,:)  ::  stc_u
#endif
#if defined(NLT_V) || defined(DYN_V) || defined(PHY_V)
   real, allocatable, dimension(:,:,:)  ::  stc_v
#endif
#if defined(NLT_T) || defined(DYN_T) || defined(PHY_T)
   real, allocatable, dimension(:,:,:)  ::  stc_t
#endif
#if defined(NLT_Q) || defined(DYN_Q) || defined(PHY_Q)
   real, allocatable, dimension(:,:,:)  ::  stc_q
#endif
!
! grid values at n, n-1 time steps
!
#ifndef DFS
   real, allocatable, dimension(:,:,:)  ::  grid_u1,grid_v1,grid_t1,grid_q1,&
                                            grid_u2,grid_v2,grid_t2,grid_q2
#else
   real, allocatable, dimension(:,:,:)  ::  dfsg_u1,dfsg_v1,dfsg_t1,dfsg_q1,&
                                            dfsg_u2,dfsg_v2,dfsg_t2,dfsg_q2
#endif
   contains
!-------------------------------------------------------------------------------
   subroutine stc_init(idate)
!-------------------------------------------------------------------------------
   use paramodel, only : LONF2S,LATG2S,levs_
#ifdef DFS
   use dfsvar   , only : ib,jbw,levs
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variabels
!
   integer, intent(in) , dimension(4)  ::  idate
!  real   , intent(in)                 ::  fhour
!
! local variables
!
   integer                             ::  krsize, iseed
   real                                ::  seed0, me
   integer, allocatable, dimension(:)  ::  nrnd
   data krsize/20/
!
! allocation
!
#if defined(NLT_U) || defined(DYN_U) || defined(PHY_U)
   allocate(  stc_u(LONF2S,levs_,LATG2S)) ; stc_u=0.
#endif
#if defined(NLT_V) || defined(DYN_V) || defined(PHY_V)
   allocate(  stc_v(LONF2S,levs_,LATG2S)) ; stc_v=0.
#endif
#if defined(NLT_T) || defined(DYN_T) || defined(PHY_T)
   allocate(  stc_t(LONF2S,levs_,LATG2S)) ; stc_t=0.
#endif
#if defined(NLT_Q) || defined(DYN_Q) || defined(PHY_Q)
   allocate(  stc_q(LONF2S,levs_,LATG2S)) ; stc_q=0.
#endif
!
#ifndef DFS
   allocate(grid_u1(LONF2S,levs_,LATG2S)) ; grid_u1=0.
   allocate(grid_v1(LONF2S,levs_,LATG2S)) ; grid_v1=0.
   allocate(grid_t1(LONF2S,levs_,LATG2S)) ; grid_t1=0.
   allocate(grid_q1(LONF2S,levs_,LATG2S)) ; grid_q1=0.
   allocate(grid_u2(LONF2S,levs_,LATG2S)) ; grid_u2=0.
   allocate(grid_v2(LONF2S,levs_,LATG2S)) ; grid_v2=0.
   allocate(grid_t2(LONF2S,levs_,LATG2S)) ; grid_t2=0.
   allocate(grid_q2(LONF2S,levs_,LATG2S)) ; grid_q2=0.
#else
   allocate(dfsg_u1(ib,jbw,levs)        ) ; dfsg_u1=0.
   allocate(dfsg_v1(ib,jbw,levs)        ) ; dfsg_v1=0.
   allocate(dfsg_t1(ib,jbw,levs)        ) ; dfsg_t1=0.
   allocate(dfsg_q1(ib,jbw,levs)        ) ; dfsg_q1=0.
   allocate(dfsg_u2(ib,jbw,levs)        ) ; dfsg_u2=0.
   allocate(dfsg_v2(ib,jbw,levs)        ) ; dfsg_v2=0.
   allocate(dfsg_t2(ib,jbw,levs)        ) ; dfsg_t2=0.
   allocate(dfsg_q2(ib,jbw,levs)        ) ; dfsg_q2=0.
#endif
!
! random seed by date and fhour
!
   call random_seed(size=krsize)
   allocate(nrnd(krsize))
   seed0=idate(1)+idate(2)+idate(3)+idate(4)
!  me=0.0
!  iseed = mod(100.0*sqrt(fhour*3600+100.0*me),1.0E9) + 1 + seed0
   iseed = 1 + seed0
   nrnd(1:krsize) = iseed
!
! random seed
!
#ifdef IBM
   call random_seed(generator=2)
#else
   call random_seed(put=nrnd)
#endif
!
   deallocate(nrnd)
!
   return
   end subroutine stc_init
!-------------------------------------------------------------------------------
   subroutine stc_forcing(random_forcing)
!-------------------------------------------------------------------------------
   use paramodel, only : LONF2S,LATG2S,lonf2_,latg2_,levs_
#ifdef STOCH_SC
   use paramodel, only : lonf_
   use comfgrid , only : colrad
   use constant , only : pi_
#endif
   use comfver  , only : kdt,deltim,fhour,thour,shour
#ifdef DFS
   use dfsvar   , only : ib,jbw,iba,jbwa,iope,sl=>sigma
#else
   use comio    , only : iope
   use comfver  , only : sl
#ifdef MP
   use paramodel, only : lonf2p_,latg2p_
#endif
#endif /* DFS end */
#ifdef DFS
#define NX ib
#define NY jbw
#define NXG iba
#define NYG jbwa
#define MPGF2P mpgf2pd
#else
#define NX LONF2S
#define NY LATG2S
#define NXG lonf2_
#define NYG latg2_
#define MPGF2P mpgf2p
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! parameters
!
   real   , parameter   ::  Imax0=2.0  ! max interval
#ifdef STOCH_SC
   real   , parameter   ::  SCdeg=10.  ! degree of spatial correlation
#endif
#ifdef STOCH_TIME
   real   , parameter   ::  Tmax=10.   ! time (day) of max interval
#endif
#ifdef STOCH_TIME2
   real   , parameter   ::  T2_fh=120. ! starting time (hour) of stochasity
#endif
#ifdef STOCH_INT
   real   , parameter   ::  HourINT=3. ! update interval in time (hour)
#endif
!
! passing variables
!
   real  , dimension(LONF2S,levs_,LATG2S)  ::  random_forcing
!
! local variables
!
   integer                                 ::  i,j,k
   real                                    ::  fh,dummy
   real  , dimension(NXG,NYG,levs_)        ::  g_random_forcing
   real  , dimension(NX ,NY ,levs_)        ::  temp
   real  , dimension(levs_)                ::  interval,Imax
#ifdef STOCH_SC
   integer, parameter                      ::  SCnx=360/SCdeg,SCny=SCnx/2
   integer                                 ::  ilon,jlat
   real  , dimension(lonf2_)               ::  londeg
   real  , dimension(latg2_)               ::  latdeg
   real  , dimension(SCnx,SCny,levs_)      ::  SCref
#endif
!
   g_random_forcing=0.
   fh=deltim*kdt/3600.
#ifdef STOCH_SC
   SCref=0.
!
! longitude with SC
!
   londeg(1)      =0.
   londeg(lonf_+1)=londeg(1)
   do i = 2,lonf_
     londeg(      i)=londeg(i-1)+360./lonf_
     londeg(lonf_+i)=londeg(i)
   enddo
!
! latitude with SC
!
   do j = 1,latg2_
     latdeg(j)=colrad(j)*180./pi_
   enddo
#endif /* STOCH_SC end */
!
! no stochastic forcing for forward step
!
   random_forcing=1.
   if(kdt.le.2) return
#ifdef STOCH_DBG
   if(iope) then
     print *,'STOCH_DBG> kdt=',kdt
     print *,'STOCH_DBG> fhour,thour,shour,fh=',fhour,thour,shour,fh
   endif
#endif
!
#ifdef STOCH_INT
   if(mod(fh,hourINT).eq.0.) then  ! update stochastic forcing every hourINT
#endif
!
! maximum interval
!
   do k = 1,levs_
#ifdef STOCH_LEVEL
     Imax(k)=Imax0*exp(sl(k)-1.)
#else
     Imax(k)=Imax0
#endif
#ifdef STOCH_DBG
     if(iope) print *,'STOCH_DBG> k,sl,Imax=',k,sl(k),Imax(k)
#endif
   enddo
!
! interval
!
   do k = 1,levs_
#if defined(STOCH_TIME)
     if(fh/24.le.Tmax) then
       interval(k)=Imax(k)*exp((fh/24.-Tmax)/3.)
     else
       interval(k)=Imax(k)
     endif
#elif defined(STOCH_TIME2)
     if(fh.le.T2_fh) then
       interval(k)=0.0
     else
       interval(k)=Imax(k)
     endif
#else
     interval(k)=Imax(k)
#endif /* STOCH_TIME end */
#ifdef STOCH_DBG
     if(iope) print *,'STOCH_DBG> k,sl,interval=',k,sl(k),interval(k)
#endif
   enddo
!
   if(iope) then
     !
     ! spin-up for random number
     !
     do k = 1,1000
       call stc_interval_uniform(dummy,1.)
     enddo
#ifdef STOCH_SC
     do k = 1,levs_
       do j = 1,SCny
         do i = 1,SCnx
           call stc_interval_uniform(SCref(i,j,k),interval(k))
         enddo
       enddo
     enddo
     do k = 1,levs_
       do j = 1,NYG
         do i = 1,NXG
           ilon=londeg(i)/SCdeg+1
           jlat=latdeg(j)/SCdeg+1
           if(i.gt.NXG/2) then
             jlat=SCny-jlat+1
           endif
           g_random_forcing(i,j,k)=SCref(ilon,jlat,k)
         enddo
       enddo
     enddo
#else
     do k = 1,levs_
       do j = 1,NYG
         do i = 1,NXG
           call stc_interval_uniform(g_random_forcing(i,j,k),interval(k))
         enddo
       enddo
     enddo
#endif
#ifdef STOCH_DBG
   print *,'fhour,max,min=',fhour,maxval(g_random_forcing),minval(g_random_forcing)
#endif
   endif
!
#ifdef MP
! full to partial
!
   call MPGF2P(g_random_forcing,NXG,NYG,temp,NX,NY,levs_)
#else
   temp=g_random_forcing
#endif /* MP end */
!   random_forcing=0.
#ifndef DFS
   do j = 1,LATG2S
     do k = 1,levs_
       do i = 1,LONF2S
         random_forcing(i,k,j)=temp(i,j,k)
       enddo
     enddo
   enddo
#else
   do j = 1,jbw/2
     do k = 1,levs_
       do i = 1,ib
         random_forcing(i   ,k,j)=temp(i,j      ,k)
         random_forcing(ib+i,k,j)=temp(i,jbw-j+1,k)
       enddo
     enddo
   enddo
#endif /* ~DFS end */
!
#ifdef STOCH_INT
   endif
#endif
   return
#undef NX
#undef NY
#undef NXG
#undef NYG
#undef MPGF2P
   end subroutine stc_forcing
!-------------------------------------------------------------------------------
   subroutine stc_interval_uniform(weight,range)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
   real, intent(inout)    ::  weight  ! weighting
   real, intent(in)       ::  range   ! range
!
! local variables
!
   real                   ::  rnum
!
!  generate random number
!
   call random_number(rnum)
   weight=1.+range*(rnum-0.5)
!
   return
   end subroutine stc_interval_uniform
!
!-------------------------------------------------------------------------------
   subroutine stc_interval_gaussian(mu,sigma)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
   real, intent(inout)    ::  mu      ! mean
   real, intent(in)       ::  sigma   ! standard deviation
!
! local variables
!
   integer                ::  i
   real                   ::  sum,x,z,y
!
! random number on Gaussian distribution
!
   sum=0.
   do i = 1,12
     call random_number(x)
     sum=sum+x
   enddo
   z=sum-6.
   y=sigma*z
!
! sigma cut off (default=5)
!
   y=min(y,sigma*5.)   ! 5 sigma
!
! update
!
   mu=mu*(y+1.)
!
   return
   end subroutine stc_interval_gaussian
#endif
!-------------------------------------------------------------------------------
   end module stochastic
