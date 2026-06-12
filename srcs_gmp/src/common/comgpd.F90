#include <define.h>
   module comgpd
!-------------------------------------------------------------------------------
#ifdef DGP
   use paramodel, only   :  ltstp_,slvark_,mlvark_,levs_,lpnt_
!-------------------------------------------------------------------------------
   private              ::  ltstp_,slvark_,mlvark_,levs_,lpnt_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  nstken=ltstp_
   integer              ::  nvrken=slvark_+mlvark_*levs_
   integer              ::  nptken=lpnt_
   integer, allocatable ::  igrd(nptken),jgrd(nptken)
   integer              ::  itnum,npoint,isave,isshrt,ilshrt,ikfreq
   real, allocatable    ::  svdata(nvrken,nptken,nstken)
!
   contains
!-------------------------------------------------------------------------------
   subroutine comgpd_init
!-------------------------------------------------------------------------------
   nstken=ltstp_
   nvrken=slvark_+mlvark_*levs_
   nptken=lpnt_
   allocate ( igrd(nptken),jgrd(nptken) )
   allocate ( svdata(nvrken,nptken,nstken) )
!
   end subroutine comgpd_init
!-------------------------------------------------------------------------------
#endif
   end module comgpd
