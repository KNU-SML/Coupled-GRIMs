#include <define.h>
   module rscombgt
!-------------------------------------------------------------------------------
#ifdef RMP
#ifdef SSI
   use paramodel, only  :  lngrd_,levs_,nt_,nq_,nu_,nv_
!-------------------------------------------------------------------------------
   private             ::  lngrd_,levs_,nt_,nq_,nu_,nv_
#ifdef A
!
!  include budget common
!
#endif
#ifdef T
   real, allocatable, dimension(:,:,:)  ::  wt
#endif
#ifdef Q
   real, allocatable, dimension(:,:,:)  ::  wq
#endif
#ifdef U
   real, allocatable, dimension(:,:,:)  ::  wu
#endif
#ifdef V
   real, allocatable, dimension(:,:,:)  ::  wv
#endif
#ifdef P
   real, allocatable, dimension(:,:,:)  ::  wp
#endif
#ifdef A
   real, allocatable, dimension(:,:)    ::  tmpbgt
#endif
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscombgt_init
!-------------------------------------------------------------------------------
#ifdef T
   allocate(   wt      (lngrd_,levs_,nt_+1)      )
#endif
#ifdef Q
   allocate(   wq      (lngrd_,levs_,nq_+1)      )
#endif
#ifdef U
   allocate(   wu      (lngrd_,levs_,nu_+1)      )
#endif
#ifdef V
   allocate(   wv      (lngrd_,levs_,nv_+1)      )
#endif
#ifdef P
   allocate(   wp      (lngrd_,1,np_+1)          )
#endif
#ifdef A
   allocate(   tmpbgt  (lngrd_,levs_)            )
#endif
!
   end subroutine rscombgt_init
!-------------------------------------------------------------------------------
#endif /* SSI end */
#endif /* RMP end */
   end module rscombgt
