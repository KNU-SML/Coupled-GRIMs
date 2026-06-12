#include <define.h>
   module compsplit
!-------------------------------------------------------------------------------
   use paramodel, only   :  LONF22S,LATG2S,levs_
!-------------------------------------------------------------------------------
   private              ::  LONF22S,LATG2S,levs_
!
   real, allocatable    ::  dudtm(:,:,:)                                      ,&
                            dvdtm(:,:,:)                                      ,&
                            dtdtm(:,:,:)                                      ,&
                            drdtm(:,:,:)                                      ,&
                            dpdtm(:,:)
   contains
!-------------------------------------------------------------------------------
   subroutine compsplit_init
!-------------------------------------------------------------------------------
   allocate(                dudtm(LONF22S,levs_,LATG2S)                       ,&
                            dvdtm(LONF22S,levs_,LATG2S)                       ,&
                            dtdtm(LONF22S,levs_,LATG2S)                       ,&
                            drdtm(LONF22S,levs_,LATG2S)                       ,&
                            dpdtm(LONF22S,      LATG2S)                        )
!
   end subroutine compsplit_init
!-------------------------------------------------------------------------------
   end module compsplit
