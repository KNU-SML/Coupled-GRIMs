#include <define.h>
#include <machine.h>
   module rscomfc
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only : lnuv_,lnt22_,clngrd_,levs_,levh_
!-------------------------------------------------------------------------------
   private            :: lnuv_,lnt22_,clngrd_,levs_,levh_
#define LNUVS lnuv_
#define LNT22S lnt22_
#define CLNGRDS clngrd_ 
!
!  begin (rscomfc)
!    version with stacked transforms
!
#ifdef G2R
   integer,              dimension(4)    ::  idate
   real   , allocatable, dimension(:)    ::  si                               ,&
                                             sl                               ,&
                                             eps                              ,&
                                             epsi                             ,&
                                             gz                               ,&
                                             q
   real   , allocatable, dimension(:,:)  ::  ze                               ,&
                                             di                               ,&
                                             te                               ,&
                                             rq
#endif
#ifdef C2R
   integer,              dimension(4)    ::  idate
   real   , allocatable, dimension(:)    ::  si                               ,&
                                             sl                               ,&
                                             gz                               ,&
                                             q                                ,&
                                             flat                             ,&
                                             flon                             ,&
                                             fm2                              ,&
                                             fm2x                             ,&
                                             fm2y
   real   , allocatable, dimension(:,:)  ::  te                               ,&
                                             rq                               ,&
                                             uu                               ,&
                                             vv
#endif
   real                                  ::  filta
#ifdef CRA
   integer                               ::  ncpus                            ,&
                                             ncpus1
#endif
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscomfc_init
!-------------------------------------------------------------------------------
#ifdef G2R
   allocate(   si    (levs_+1)                                                ,&
               sl    (levs_)                                                  ,&
               eps   (LNUVS)                                                  ,&
               epsi  (LNUVS)                                                  ,&
               gz    (LNT22S)                                                 ,&
               q     (LNT22S)                                                  )
   allocate(   ze    (LNT22S,levs_)                                           ,&
               di    (LNT22S,levs_)                                           ,&
               te    (LNT22S,levs_)                                           ,&
               rq    (LNT22S,levh_)                                            )
#endif
#ifdef C2R
   allocate(   si    (levs_+1)                                                ,&
               sl    (levs_)                                                  ,&
               gz    (CLNGRDS)                                                ,&
               q     (CLNGRDS)                                                ,&
               flat  (CLNGRDS)                                                ,&
               flon  (CLNGRDS)                                                ,&
               fm2   (CLNGRDS)                                                ,&
               fm2x  (CLNGRDS)                                                ,&
               fm2y  (CLNGRDS)                                                 )
   allocate(   te    (CLNGRDS,levs_)                                          ,&
               rq    (CLNGRDS,levh_)                                          ,&
               uu    (CLNGRDS,levs_)                                          ,&
               vv    (CLNGRDS,levs_)                                           )
#endif
!
   end subroutine rscomfc_init
!-------------------------------------------------------------------------------
#endif /* RMP end */
   end module rscomfc
