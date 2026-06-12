#include <define.h>
   module rscomgpd
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only  :  levs_,rslvark_,rmlvark_,rlpnt_,rltstp_
!-------------------------------------------------------------------------------
   private             ::  levs_,rslvark_,rmlvark_,rlpnt_,rltstp_
!  
#ifdef RKN
   integer             ::  nvrken, nptken, nstken
   integer             ::  itnum,npoint,isave,isshrt,ilshrt,ikfreq,            &
                           imodk,irstkn
   real   , allocatable, dimension(:,:,:)   ::  svdata                        ,&
   integer, allocatable, dimension(:,:)     ::  igpd                          ,&
                                                jgpd                          ,&
                                                igpdr                         ,&
                                                jgpdr                         ,&
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscomgpd_init
!-------------------------------------------------------------------------------
   nvrken=rslvark_+rmlvark_*levs_
   nptken=rlpnt_
   nstken=rltstp_
   allocate(  svdata   (nvrken,nptken,nstken)  )
   allocate(  igpd     (nptken)                                               ,&
              jgpd     (nptken)                                               ,&
              igpdr    (nptken)                                               ,&
              jgpdr    (nptken)                )
!
   end subroutine subroutine rscomgpd_init
!-------------------------------------------------------------------------------
#endif /* RKN end */
#endif /* RMP end */
   end module rscomgpd
