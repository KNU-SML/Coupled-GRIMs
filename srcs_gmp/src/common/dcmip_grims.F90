#include <define.h>
module dcmip_grims
#ifdef DCMIP
!-------------------------------------------------------------------------------
   use paramodel, only : levs_
#ifdef DFS
   use dfsvar   , only : ib,jbw,mta,jlg
#else
   use paramodel, only : lnt22_,LONF2S,LATG2S
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer   , parameter                   ::  cc  = 1                        ,&
                                               cfv = 1
   real(8)   , parameter                   ::  p0  = 100000.0d0               ,&
                                               t0  = 300.0d0                  ,&
                                               rd  = 287.0d0                  ,&
                                               grav= 9.80616d0
   logical   , parameter                   ::  hybrid_eta=.true.
   integer                                 ::  icase
   real(8)                                 ::  ztop                           ,&
                                               dz                             ,&
                                               H                              ,&
                                               ptop
#ifndef DFS
   real(8), allocatable, dimension(:)      ::  gphf     ! orography    (m)
#endif
   real(8), allocatable, dimension(:)      ::  hyai  ,& ! hybrid a (interfacial)
                                               hybi  ,& ! hybrid b (interfacial)
                                               hyam  ,& ! hybrid a (full)
                                               hybm     ! hybrid b (full)
   real(8), allocatable, dimension(:,:)    ::  phis  ,& ! sfc geop.    (m^2/s^2)
                                               hsfc  ,& ! sfc gph      (m)
                                               ps    ,& ! sfc pressure (Pa)
                                               psfc     ! log(ps)      (cb)
   real(8), allocatable, dimension(:,:,:)  ::  rho   ,& ! density      (kg/m^3)
                                               zl    ,& ! geop. height (m)
                                               zi    ,& ! geop. height (m)
                                               pl    ,& ! pressure     (Pa)
                                               pi    ,& ! pressure     (Pa)
                                               u     ,& ! u-wind       (m/s)
                                               v     ,& ! v-wind       (m/s)
                                               us    ,& ! u-wind/cos   (m/s)
                                               vs    ,& ! v-wind/cos   (m/s)
                                               w     ,& ! veritcal wind(m/s)
                                               pwl   ,& ! pressure vel.(Pa/s)
                                               pwi   ,& ! pressure vel.(Pa/s)
                                               t     ,& ! temperature  (K)  
                                               tv    ,& ! virtual temp (K)  
                                               q0    ,& ! spec. humid. (kg/kg)
                                               q1    ,& ! tracer 1     (kg/kg)
                                               q2    ,& ! tracer 1     (kg/kg)
                                               q3    ,& ! tracer 1     (kg/kg)
                                               q4       ! tracer 1     (kg/kg)
!                     
contains
!-------------------------------------------------------------------------------
#ifdef HYBRID
   subroutine dcmip_grims_init(ak,bk)
#else
   subroutine dcmip_grims_init(si,sl)
#endif
!-------------------------------------------------------------------------------
!
! passing variables
!
#ifdef HYBRID
   real   , dimension(levs_+1), intent(out)   ::  ak,bk
#endif
   real   , dimension(levs_)                  ::  sl
   real   , dimension(levs_+1)                ::  si
!
! local vairables
!
   real   , parameter                         ::  cp   = 1004.5               ,&
                                                  nn   = 0.01                 ,&
                                                  kapa = 2./7.                ,&
                                                  p_eq = 100000.              ,&
                                                  t_eq = 300.                 ,&
                                                  laps = 0.0065               ,&
                                                  gg   = grav**2/(nn**2*cp)
   integer                                    ::  i,k,j
   real                                       ::  etatop
   real   , dimension(levs_)                  ::  etam
   real   , dimension(levs_+1)                ::  z,p,etai
#ifdef HYBRID
   real   , dimension(31)                     ::  hyai0,hybi0
   data hyai0/0.00225523952394724, 0.00503169186413288, 0.0101579474285245, &
              0.0185553170740604 , 0.0306691229343414 , 0.0458674766123295, &
              0.0633234828710556 , 0.0807014182209969 , 0.0949410423636436, &
              0.11169321089983   , 0.131401270627975  , 0.154586806893349 , &
              0.181863352656364  , 0.17459799349308   , 0.166050657629967 , &
              0.155995160341263  , 0.14416541159153   , 0.130248308181763 , &
              0.113875567913055  , 0.0946138575673103 , 0.0753444507718086, &
              0.0576589405536652 , 0.0427346378564835 , 0.0316426791250706, &
              0.0252212174236774 , 0.0191967375576496 , 0.0136180268600583, &
              0.00853108894079924, 0.00397881818935275, 0.                , &
              0. /
   data hybi0/0.                 , 0.                 , 0.                , &
              0.                 , 0.                 , 0.                , &
              0.                 , 0.                 , 0.                , &
              0.                 , 0.                 , 0.                , &
              0.                 , 0.0393548272550106 , 0.0856537595391273, &
              0.140122056007385  , 0.204201176762581  , 0.279586911201477 , &
              0.38274360895157   , 0.47261056303978   , 0.576988518238068 , &
              0.672786951065063  , 0.753628432750702  , 0.813710987567902 , &
              0.848494648933411  , 0.881127893924713  , 0.911346435546875 , &
              0.938901245594025  , 0.963559806346893  , 0.985112190246582 , &
              1. /
#endif
   namelist /dcmiplist/ icase
!-------------------------------------------------------------------------------
!
! allocate
!
#ifdef DFS
#define IJ   ib,jbw
#define IJK  ib,jbw,levs_
#define IJK1 ib,jbw,levs_+1
#else /* SPH */
#define IJ   LONF2S,LATG2S
#define IJK  LONF2S,levs_,LATG2S
#define IJK1 LONF2S,levs_+1,LATG2S
   if(.not.allocated(gphf)) allocate( gphf(lnt22_)        )
#endif /* DFS end */
   if(.not.allocated(hyai)) allocate( hyai(levs_+1)       )
   if(.not.allocated(hybi)) allocate( hybi(levs_+1)       )
   if(.not.allocated(hyam)) allocate( hyam(levs_)         )
   if(.not.allocated(hybm)) allocate( hybm(levs_)         )
   if(.not.allocated(phis)) allocate( phis(IJ  ) )
   if(.not.allocated(hsfc)) allocate( hsfc(IJ  ) )
   if(.not.allocated(ps)  ) allocate( ps  (IJ  ) )
   if(.not.allocated(psfc)) allocate( psfc(IJ  ) )
   if(.not.allocated(rho) ) allocate( rho (IJK ) )
   if(.not.allocated(zl)  ) allocate( zl  (IJK ) )
   if(.not.allocated(zi)  ) allocate( zi  (IJK1) )
   if(.not.allocated(pl)  ) allocate( pl  (IJK ) )
   if(.not.allocated(pi)  ) allocate( pi  (IJK1) )
   if(.not.allocated(u)   ) allocate( u   (IJK ) )
   if(.not.allocated(v)   ) allocate( v   (IJK ) )
   if(.not.allocated(us)  ) allocate( us  (IJK ) )
   if(.not.allocated(vs)  ) allocate( vs  (IJK ) )
   if(.not.allocated(w)   ) allocate( w   (IJK ) )
   if(.not.allocated(pwl) ) allocate( pwl (IJK ) )
   if(.not.allocated(pwi) ) allocate( pwi (IJK1) )
   if(.not.allocated(t)   ) allocate( t   (IJK ) )
   if(.not.allocated(tv)  ) allocate( tv  (IJK ) )
   if(.not.allocated(q0)  ) allocate( q0  (IJK ) )
   if(.not.allocated(q1)  ) allocate( q1  (IJK ) )
   if(.not.allocated(q2)  ) allocate( q2  (IJK ) )
   if(.not.allocated(q3)  ) allocate( q3  (IJK ) )
   if(.not.allocated(q4)  ) allocate( q4  (IJK ) )
#undef IJ
#undef IJK
#undef IJK1
!
! read dcmip.parm
!
   open(151,file='dcmip.parm')
   read(151,dcmiplist)
   close(151)
#ifdef DFS
#define LONF2S ib
#define LATG2S jbw
#define IJK i,j,k
#else
#define IJK i,k,j
#endif
!
! model top
!
   if (icase.eq.11 .or. icase.eq.12 .or. icase.eq.13 .or. icase.eq.200 ) then
     ztop=12000.0d0
   else if (icase.eq.410 .or. icase.eq.42 .or. icase.eq.43 .or. icase.eq.51) then
     ztop=44000.0d0
     ! interface
     do k = 1,levs_+1
       hyai(k)=hyai0(levs_+2-k)
       hybi(k)=hybi0(levs_+2-k)
       etai(k)=hyai(k)+hybi(k)
       do j = 1,LATG2S
         do i = 1,LONF2S
           pi(IJK)=etai(k)*p0
         enddo
       enddo
#ifdef HYBRID
       ! [t2b] ak,bk
       ak(levs_+2-k) = hyai(k)*p0
       bk(levs_+2-k) = hybi(k)
#endif
     enddo
     ! model
     do k = 1,levs_
       hyam(k)=0.5*(hyai(k)+hyai(k+1))
       hybm(k)=0.5*(hybi(k)+hybi(k+1))
       etam(k)=hyam(k)+hybm(k)
       do j = 1,LATG2S
         do i = 1,LONF2S
           pl(IJK)=etam(k)*p0
         enddo
       enddo
     enddo
     return
   endif
   H      = rd*t0/grav
   dz     = ztop/levs_
   if(icase.eq.11.or.icase.eq.12 .or.icase.eq.13) then
     ptop = p_eq*exp(-grav*ztop/(rd*t0))
   elseif(icase.eq.200) then
     ptop = p0*(1.-laps/t0*ztop)**(grav/(rd*laps))
!     ptop = 20544.8
   elseif(icase.eq.31) then
     ptop = p_eq*(gg/t_eq*exp(-nn**2*ztop/grav)+1-gg/t_eq)**(1./kapa)
   endif
   etatop = ptop/p0
!
! at half (interface) levels
!
   do k = 1,levs_+1
     z(k)=(k-1)*dz
     !
     if ( icase.eq.11 .or. icase.eq.12 .or. icase.eq.13 ) then
       p(k)=p_eq*exp(-grav*z(k)/(rd*t0))
     elseif ( icase.eq.200 ) then
       p(k)=p0*(1.-laps/t0*z(k))**(grav/(rd*laps))
     elseif ( icase.eq.31 ) then
       p(k)=p_eq*(gg/t_eq*exp(-nn**2*z(k)/grav)+1-gg/t_eq)**(1./kapa)
     endif
     !
     si(k)         = (p(k)-ptop)/(p0-ptop)
     etai(k)       = p(k)/p0
     hybi(k)       = ((etai(k)-etatop)/(1.-etatop))**cc
     hyai(k)       = etai(k)-hybi(k)
#ifdef HYBRID
     ! [t2b] ak,bk
     ak(levs_+2-k) = hyai(k)*p0
     bk(levs_+2-k) = hybi(k)
#endif
     do j = 1,LATG2S
       do i = 1,LONF2S
         pi(IJK)=p(k)
       enddo
     enddo
   enddo
!
! at full (model) levels
!
   do k = 1,levs_
     z(k)=(k-1)*dz+dz*0.5
     if ( icase.eq.11 .or. icase.eq.12 .or. icase.eq.13 .or. &
          icase.eq.21 .or. icase.eq.22 ) then
       p(k)=p_eq*exp(-grav*z(k)/(rd*t0))
     elseif ( icase.eq.20 ) then
       p(k)=p0*(1.-laps/t0*z(k))**(grav/(rd*laps))
     elseif ( icase.eq.31 ) then
       p(k)=p_eq*(gg/t_eq*exp(-nn**2*z(k)/grav)+1-gg/t_eq)**(1./kapa)
     endif
     !
     sl(k)         = (p(k)-ptop)/(p0-ptop)
     hyam(k)       = 0.5*(hyai(k)+hyai(k+1))
     hybm(k)       = 0.5*(hybi(k)+hybi(k+1))
     etam(k)       = hyam(k)+hybm(k)
     do j = 1,LATG2S
       do i = 1,LONF2S
         pl(IJK)=p(k)
       enddo
     enddo
   enddo
#undef IJK
#ifdef DFS
#undef LONF2S
#undef LATG2S
#endif
   !
   return
   end subroutine dcmip_grims_init
!
#endif /* DCMIP end */
end module dcmip_grims
