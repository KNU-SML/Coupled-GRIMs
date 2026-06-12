#include <define.h>
   module rscomrinpg
!-------------------------------------------------------------------------------
   use paramodel, only : levs_, levh_, levp1_, lnt22_, lngrd_
!-------------------------------------------------------------------------------
   real   , allocatable, dimension(:)      ::  si                             ,&
                                               sl                             ,&
                                               gz                             ,&
                                               q                              ,&
                                               gm
   real   , allocatable, dimension(:,:)    ::  te                             ,&
                                               rq                             ,&
                                               uu                             ,&
                                               vv                             ,&
                                               di                             ,&
                                               ze                             ,&
                                               uln                            ,&
                                               vln
   real   , allocatable, dimension(:)      ::  rcsln                          ,&
                                               rsnln                          ,&
                                               xm                             ,&
                                               flat                           ,&
                                               flon                           ,&
                                               fm2                            ,&
                                               fm2x                           ,&
                                               fm2y
   real    , allocatable, dimension(:)     ::  ggz                            ,&
                                               ggzo                           ,&
                                               gq
   real    , allocatable, dimension(:,:)   ::  gte                            ,&
                                               grq                            ,&
                                               guu                            ,&
                                               gvv
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscomrinpg_init
!-------------------------------------------------------------------------------
   allocate(   si      (levp1_)                                               ,&
               sl      (levs_)                                                ,&
               gz      (lnt22_)                                               ,&
               q       (lnt22_)                                               ,&
               gm      (lnt22_)                                                )
   allocate(   te      (lnt22_,levs_)                                         ,&
               rq      (lnt22_,levh_)                                         ,&
               uu      (lnt22_,levs_)                                         ,&
               vv      (lnt22_,levs_)                                         ,&
               di      (lnt22_,levs_)                                         ,&
               ze      (lnt22_,levs_)                                         ,&
               uln     (lnt22_,levs_)                                         ,&
               vln     (lnt22_,levs_)                                          )
   allocate(   rcsln   (lngrd_)                                               ,&
               rsnln   (lngrd_)                                               ,&
               xm      (lngrd_)                                               ,&
               flat    (lngrd_)                                               ,&
               flon    (lngrd_)                                               ,&
               fm2     (lngrd_)                                               ,&
               fm2x    (lngrd_)                                               ,&
               fm2y    (lngrd_)                                                )
   allocate(   ggz     (lngrd_)                                               ,&
               ggzo    (lngrd_)                                               ,&
               gq      (lngrd_)                                                )
   allocate(   gte     (lngrd_,levs_)                                         ,&
               grq     (lngrd_,levh_)                                         ,&
               guu     (lngrd_,levs_)                                         ,&
               gvv     (lngrd_,levs_)                                          )
!
   end subroutine rscomrinpg_init
!-------------------------------------------------------------------------------
   end module rscomrinpg
