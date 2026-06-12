#include <define.h>   
   module rscompost
#ifdef RMP
!-------------------------------------------------------------------------------
   integer  ::  ksz,ksd,kst,ksq,kspsx,kspsy,ksu,ksv,ksps,kszs,kszsx,kszsy,     &
                ksqc,ksqr,ksqi,ksqs,ksqg,ksnccn,ksnc,ksnr,                     &
                kso3,kstke,kspw
   logical lpqr,lpqc,lpqi,lpqs,lpqg
   logical lpnccn,lpnc,lpnr,lpo3,lptke
#ifdef IDEAL
   integer  ::  kstt
#ifdef POST_SEC
   integer  ::  postsec
#endif
#endif
#ifdef NONHYD
   integer  ::  kspn,kstn,kswn
#endif
!-------------------------------------------------------------------------------
#endif /* RMP */
   end module rscompost
   
