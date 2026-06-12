#include <define.h>
   module rscomltb
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only  : LNGRDS,bgf_,border_,levs_,levh_ 
   use rsparltb
!-------------------------------------------------------------------------------
   private             :: LNGRDS,bgf_,border_,levs_,levh_
!
!  begin comltb
!
   real   , allocatable, dimension(:)    ::  rltb                             ,&
                                             sltb
   real   , allocatable, dimension(:)    ::  blat                             ,&
                                             blon                             ,&
                                             bcsln                            ,&
                                             bsnln                            ,&
                                             qb                               ,&
                                             gzb                              ,&
                                             qltb
   real   , allocatable, dimension(:,:)  ::  wsltb                            ,&
                                             dwsltb                           ,&
                                             uub                              ,&
                                             vvb                              ,&
                                             teb                              ,&
                                             rqb                              ,&
                                             uultb                            ,&
                                             vvltb                            ,&
                                             teltb                            ,&
                                             rqltb
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscomltb_init
!-------------------------------------------------------------------------------
   allocate(   rltb    (LNGRDS)                                               ,&
               sltb    (LNGRDS)                                                )
   allocate(   blat    (lngrdb)                                               ,&
               blon    (lngrdb)                                               ,&
               bcsln   (lngrdb)                                               ,&
               bsnln   (lngrdb)                                               ,&
               qb      (lngrdb)                                               ,&
               gzb     (lngrdb)                                               ,&
               qltb    (lngrdb)                                                )
   allocate(   wsltb   (bgf_,1-border_:border_)                               ,&
               dwsltb  (bgf_,1-border_:border_)                               ,&
               uub     (lngrdb,levs_)                                         ,&
               vvb     (lngrdb,levs_)                                         ,&
               teb     (lngrdb,levs_)                                         ,&
               rqb     (lngrdb,levh_)                                         ,&
               uultb   (lngrdb,levs_)                                         ,&
               vvltb   (lngrdb,levs_)                                         ,&
               teltb   (lngrdb,levs_)                                         ,&
               rqltb   (lngrdb,levh_)                                          )
!
   end subroutine rscomltb_init
!-------------------------------------------------------------------------------
#endif /* RMP end */
   end module rscomltb
