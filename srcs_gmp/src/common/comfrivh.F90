#include <define.h>
   module comfrivh
!-------------------------------------------------------------------------------
#ifdef RIVER
   use paramodel, only : LONF2S, LATG2S, io2_, jo2_
!-------------------------------------------------------------------------------
!
!  begin comfriv (kei)
!
   real, allocatable, dimension(:,:)    :: trunof                             ,&
                                           imap                               ,&
                                           gdriv                              ,&
                                           rflow                              ,&
                                           roff
!
   contains
!-------------------------------------------------------------------------------
   subroutine comfrivh_init
!-------------------------------------------------------------------------------
   allocate(   trunof  (LONF2S,LATG2S)                                        ,&
               imap    (io2_,jo2_)                                            ,&
               gdriv   (io2_,jo2_)                                            ,&
               rflow   (io2_,jo2_)                                            ,&
               roff    (io2_,jo2_) )
   end subroutine comfrivh_init
!-------------------------------------------------------------------------------
#endif
   end module comfrivh
