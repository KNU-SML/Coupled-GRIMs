#include <define.h>
   module comdyn
!-------------------------------------------------------------------------------
#if defined(RAW) || defined(HORA)
#ifdef DFS
   use dfsvar      , only : mt,jlg,levsp,levhp
#else
   use paramodel   , only : LNT22S,levs_,levh_
#endif
#endif /* RAW or HORA end */
!-------------------------------------------------------------------------------
#if defined(RAW) || defined(HORA)
#ifdef RAW
   real   , parameter                              ::  alpha = 0.53           ,&
                                                       nu    = 0.2
#else /* HORA */
   real   , parameter                              ::  beta  = 0.4
#endif
   real   , allocatable, target  , dimension(:,:)  ::  workt,workd,workz,workr
   real                                            ::  coef1,coef2,coef3,coef4
#endif /* RAW or HORA end */
!
contains
!-------------------------------------------------------------------------------
   subroutine comdyn_init
!-------------------------------------------------------------------------------
#if defined(RAW) || defined(HORA)
!
! allocation
!
#ifdef DFS
   allocate(   workt (1:mt*jlg,1:levsp) )
   allocate(   workd (1:mt*jlg,1:levsp) )
   allocate(   workz (1:mt*jlg,1:levsp) )
   allocate(   workr (1:mt*jlg,1:levhp) )
#else
   allocate(   workt (1:LNT22S,1:levs_) )
   allocate(   workd (1:LNT22S,1:levs_) )
   allocate(   workz (1:LNT22S,1:levs_) )
   allocate(   workr (1:LNT22S,1:levh_) )
#endif
!
! coefficient
!
#ifdef RAW
   coef1 =    0.5*nu*alpha
   coef2 =     1.-nu*alpha
   coef3 =        nu*(alpha-1.)
   coef4 = 1.+0.5*nu*(alpha-1.)
#else /* HORA */
   coef1 =      -0.5*beta
   coef2 =    3.*0.5*beta
   coef3 = 1.-3.*0.5*beta
   coef4 =       0.5*beta
#endif /* RAW end */

#endif /* RAW or HORA end */
   return
   end subroutine comdyn_init
!
#ifdef HORA
!-------------------------------------------------------------------------------
   subroutine hora_n2(tem,dim,zem,rm,t,d,z,r,lnt22,levs,levh)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                         ::  lnt22,levs,levh
   real   , dimension(lnt22,levs)  ::  tem,dim,zem,t,d,z
   real   , dimension(lnt22,levh)  ::  rm         ,r
!
   t(1:lnt22,1:levs)=tem(1:lnt22,1:levs)
   d(1:lnt22,1:levs)=dim(1:lnt22,1:levs)
   z(1:lnt22,1:levs)=zem(1:lnt22,1:levs)
   r(1:lnt22,1:levh)=rm(1:lnt22,1:levh)

   return
   end subroutine hora_n2
#endif /* HORA end */
!-------------------------------------------------------------------------------
   end module comdyn
