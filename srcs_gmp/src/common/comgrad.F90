#include <define.h>
   module comgrad
!-------------------------------------------------------------------------------
#ifndef RMP
   use paramodel, only      : LONF2S,LATG2S,LATGS,levs_
#ifdef SWRMDC
   use  aerparm, only       : nxc,ndn,imxae,jmxae
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   private                 ::  LONF2S,LATG2S,LATGS,levs_
#ifdef SWRMDC
   private                 ::  nxc,ndn,imxae,jmxae
#endif
!
   integer, parameter      ::  mcld=3,nseal=2,nbin=100,nlon=2,nlat=4
   real                    ::  rhcl(nbin,nlon,nlat,mcld,nseal)
   real                    ::  solc,rsin1,rcos1,rcos2,tsmin,tsmax,shmin,shmax
   real                    ::  raddt,fjd,r1,alf,dlt,dtlw
#ifdef DG
   real                    ::  dtacc
#endif
   real, allocatable       ::  albedr(:,:),slmskr(:,:)
   real, allocatable       ::  acoszer(:,:),coszdg(:,:)
   real, allocatable       ::  cvr  (:,:),cvtr (:,:)
   real, allocatable       ::  cvbr (:,:)
   real, allocatable       ::  alvbr(:,:),alnbr(:,:)
   real, allocatable       ::  alvdr(:,:),alndr(:,:)
   real, allocatable       ::  avecld(:,:),cldl(:,:)
   real, allocatable       ::  avecv(:,:),zonht(:,:)
   real, allocatable       ::  cldsig(:,:)
   real, allocatable       ::  alvsf(:,:,:),alnsf(:,:,:)
   real, allocatable       ::  alvwf(:,:,:),alnwf(:,:,:)
   real, allocatable       ::  facsf(:,:)  ,facwf(:,:)
#ifndef SWRMDC
   integer                 ::  nalaer
   real, allocatable       ::  paerr(:,:,:)
   real, allocatable       ::  paerf(:,:,:)
#else
   integer                 ::  nfalb,nfaer
   integer, allocatable    ::  idxcg(:,:,:),kprfg(:,:)
   real, allocatable       ::  cmixg(:,:,:),denng(:,:,:)
#endif
!
   integer                 ::  ier
#ifdef NIM
   real                    ::  rrs2
#endif
   integer                 ::  istrat,jo3,kalb,jsno,itimsw,itimlw
   integer                 ::  icfc,icwp,iswsrc(5),ibnd,ko3
   logical                 ::  runrad
!
   contains
!-------------------------------------------------------------------------------
   subroutine comgrad_init
!-------------------------------------------------------------------------------
   allocate(albedr(LONF2S,LATG2S),slmskr(LONF2S,LATG2S))
   allocate(acoszer(LONF2S,LATG2S),coszdg(LONF2S,LATG2S))
   allocate(cvr  (LONF2S,LATG2S),cvtr (LONF2S,LATG2S))
   allocate(cvbr (LONF2S,LATG2S))
   allocate(alvbr(LONF2S,LATG2S),alnbr(LONF2S,LATG2S))
   allocate(alvdr(LONF2S,LATG2S),alndr(LONF2S,LATG2S))
   allocate(avecld(levs_,LATGS),cldl(4,LATGS))
   allocate(avecv(3,LATGS),zonht(levs_,LATGS))
   allocate(cldsig(levs_,18))
   allocate(alvsf(LONF2S,LATG2S,4),alnsf(LONF2S,LATG2S,4))
   allocate(alvwf(LONF2S,LATG2S,4),alnwf(LONF2S,LATG2S,4))
   allocate(facsf(LONF2S,LATG2S)  ,facwf(LONF2S,LATG2S))
#ifndef SWRMDC
   allocate(paerr(LONF2S,5,LATG2S))
   allocate(paerf(LONF2S,LATG2S,5))
#else
   allocate(idxcg(nxc,LONF2S,LATG2S),kprfg(LONF2S,LATG2S))
   allocate(cmixg(nxc,LONF2S,LATG2S),denng(ndn,LONF2S,LATG2S))
#endif
#ifndef SWRMDC
   paerr=0.   !initialized here for SMS
#endif
!
   end subroutine comgrad_init
!-------------------------------------------------------------------------------
#endif /* not RMP */
   end module comgrad
