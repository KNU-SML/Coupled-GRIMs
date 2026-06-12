#include <define.h>
   module rscommap
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only : IGRD12S,JGRD12S
!-------------------------------------------------------------------------------
   private            :: IGRD12S,JGRD12S
!
!   begin commap
!
   real                                  ::  xm2m
   real   , allocatable, dimension(:,:)  ::  rlat                             ,&
                                             rlon                             ,&
                                             xm2                              ,&
                                             xm2p                             ,&
                                             xm2px                            ,&
                                             xm2py                            ,&
                                             gzdx                             ,&
                                             gzdy                             ,&
                                             corf                             ,&
#ifdef NONHYD
                                             corf2                            ,&
#endif
                                             xm
!
   contains
!-------------------------------------------------------------------------------
   subroutine rscommap_init
!-------------------------------------------------------------------------------
   allocate(   rlat  (IGRD12S,JGRD12S)                                        ,&
               rlon  (IGRD12S,JGRD12S)                                        ,&
               xm2   (IGRD12S,JGRD12S)                                        ,&
               xm2p  (IGRD12S,JGRD12S)                                        ,&
               xm2px (IGRD12S,JGRD12S)                                        ,&
               xm2py (IGRD12S,JGRD12S)                                        ,&
               gzdx  (IGRD12S,JGRD12S)                                        ,&
               gzdy  (IGRD12S,JGRD12S)                                        ,&
               corf  (IGRD12S,JGRD12S)                                        ,&
#ifdef NONHYD
               corf2 (IGRD12S,JGRD12S)                                        ,&
#endif
               xm    (IGRD12S,JGRD12S)                                         )
!
   end subroutine rscommap_init
!-------------------------------------------------------------------------------
#endif /* RMP end */
   end module rscommap
