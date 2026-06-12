#include <define.h>
   module sumfidx 
!-------------------------------------------------------------------------------
#ifdef RMP
#ifdef VECSUM
   use paramodel, only  :  IWAV1S,LEVHS,LEVSS,jgrd12_,jwav1_
!-------------------------------------------------------------------------------
   private             ::  IWAV1S,LEVHS,LEVSS,jgrd12_,jwav1_
!
   real   , allocatable, dimension(:,:,:)  ::  idxffca1                       ,&
                                               idxffca2
   integer, allocatable, dimension(:,:,:)  ::  idxffcb                        ,&
                                               idxffcc
!
   contains
!-------------------------------------------------------------------------------
   subroutine sumfidx_init
!-------------------------------------------------------------------------------
   allocate(  idxffca1(IWAV1S*2*jgrd12_*(LEVHS*2+LEVSS*2+2),jwav1_,7)         ,&
              idxffca2(IWAV1S*2*jgrd12_*(LEVHS*2+LEVSS*2+2),jwav1_,7)          )
   allocate(  idxffcb (IWAV1S*2*jgrd12_*(LEVHS*2+LEVSS*2+2),3,7)              ,&
              idxffcc (IWAV1S*2*jgrd12_*(LEVHS*2+LEVSS*2+2),jwav1_,7)          )
!
   end subroutine sumfidx_init
!-------------------------------------------------------------------------------
#endif
#endif
   end module sumfidx 
