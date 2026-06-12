#include <define.h>
   module rdparm
!-------------------------------------------------------------------------------
#if defined (RRTMGSW) || defined (RRTMGLW)
   use paramodel, only  : LONF2S,levs_,levh_
!-------------------------------------------------------------------------------
   private             :: LONF2S,levs_,levh_
#else
   use paramodel, only  : LONF2S,levs_
!-------------------------------------------------------------------------------
   private             :: LONF2S,levs_
#endif
!
!     parameter settings for the longwave and shortwave radiation code:
!          imax   =  no. points along the lat. circle used in calcs.
!          l      =  no. vertical levs_ (also layers) in model
!   note: the user normally will modify only the imax and l parameters
!          nblw   =  no. freq. bands for approx computations. see
!                      bandta for definition
!          nblx   =  no. freq bands for approx cts computations
!          nbly   =  no. freq. bands for exact cts computations. see
!                      bdcomb for definition
!          inlte  =  no. levs_ used for nlte calcs.
!          nnlte  =  index no. of freq. band in nlte calcs.
!          nb,ko2 are shortwave parameters; other quantities are derived
!                    from the above parameters.
   integer, parameter   ::  nblw=163,nblx=47,nbly=15
   integer, parameter   ::  nblm=nbly-1
   integer, parameter   ::  nb=12
   integer, parameter   ::  inlte=3,inltep=inlte+1,nnlte=56
   integer, parameter   ::  nb1=nb-1
   integer, parameter   ::  ko2=12
   integer, parameter   ::  ko21=ko2+1,ko2m=ko2-1
!
!  for 4 uv bands sw
!
   integer, parameter   ::  nbd=8,nvb=4,nrb=4,nk0=10,nae=6,nsrc=5
   integer              ::  NVECT, imax, imbx
!
   integer              ::  l
#if defined (RRTMGSW) || defined (RRTMGLW)
   integer              ::  lh
#endif
   integer              ::  lp1
   integer              ::  lp2
   integer              ::  lp3
   integer              ::  lm1
   integer              ::  lm2
   integer              ::  lm3
   integer              ::  ll
   integer              ::  llp1
   integer              ::  llp2
   integer              ::  llp3
   integer              ::  llm1
   integer              ::  llm2
   integer              ::  llm3
   integer              ::  lp1m
   integer              ::  lp1m1
   integer              ::  lp1v
   integer              ::  lp121
   integer              ::  ll3p
!
   contains
!-------------------------------------------------------------------------------
   subroutine rdparm_init
!-------------------------------------------------------------------------------
#if defined(REDUCE_GRID) || defined(SMP)
   NVECT=2
#else
   NVECT=LONF2S
#endif
   imax=NVECT
   imbx=imax
   l=levs_
#if defined (RRTMGSW) || defined (RRTMGLW)
   lh=levh_
#endif
   lp1=l+1
   lp2=l+2
   lp3=l+3
   lm1=l-1
   lm2=l-2
   lm3=l-3
   ll=2*l
   llp1=ll+1
   llp2=ll+2
   llp3=ll+3
   llm1=ll-1
   llm2=ll-2
   llm3=ll-3
   lp1m=lp1*lp1
   lp1m1=lp1m-1
   lp1v=lp1*(1+2*l/2)
   lp121=lp1*nbly
   ll3p=3*l+2
!
   end subroutine rdparm_init
!-------------------------------------------------------------------------------
   end module rdparm
