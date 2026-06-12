#include <define.h>   
   module compspec
!-------------------------------------------------------------------------------
#if defined(MP) 
   use paramodel  , only: lln22p_,levsp_,levhp_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! local spectral coefficients
!
! n time step
!
   real   , allocatable, target  , dimension(:,:)  ::  zea
   real                , pointer , dimension(:,:)  ::  dia,tea,ulna,vlna
   real                , pointer , dimension(:)    ::  dpdphia,dpdlama,qa,qlapa
!
! n+1 time step
!
   real   , allocatable, target  , dimension(:,:)  ::  za
   real                , pointer , dimension(:,:)  ::  uua,vva,ya,xa,wa
#ifdef NISLQ
   real   , allocatable          , dimension(:,:)  ::  rqa,rta
#else
   real                , pointer , dimension(:,:)  ::  rqa,rta
#endif
!
   contains
!-------------------------------------------------------------------------------
   subroutine compspec_init
!-------------------------------------------------------------------------------
   use comkindex
!-------------------------------------------------------------------------------
!
! n time step
!
   allocate(   zea              (1:lln22p_,kszs:lotss       )                  )
               dia     =>zea    (1:lln22p_,ksds:ksds+levsp_-1)
               tea     =>zea    (1:lln22p_,ksts:ksts+levsp_-1)
#ifndef NISLQ
               rqa     =>zea    (1:lln22p_,ksrs:ksrs+levhp_-1)
#endif
               ulna    =>zea    (1:lln22p_,ksus:ksus+levsp_-1)
               vlna    =>zea    (1:lln22p_,ksvs:ksvs+levsp_-1)
               dpdphia =>zea    (1:lln22p_,kspphis)
               dpdlama =>zea    (1:lln22p_,ksplams)
               qa      =>zea    (1:lln22p_,ksps   )
               qlapa   =>zea    (1:lln22p_,ksplaps)
!
! n+1 time step
!
   allocate(   za               (1:lln22p_,kzs:lotts      )                    )
               uua     =>za     (1:lln22p_,kus:kus+levsp_-1)
               vva     =>za     (1:lln22p_,kvs:kvs+levsp_-1)
               ya      =>za     (1:lln22p_,kys:kys+levsp_-1)
#ifndef NISLQ
               rta     =>za     (1:lln22p_,krs:krs+levhp_-1)
#endif
               xa      =>za     (1:lln22p_,kxs:kxs+levsp_-1)
               wa      =>za     (1:lln22p_,kws:kws+levsp_-1)
#ifdef NISLQ
!
! for nislq
!
   allocate(   rqa              (1:lln22p_,1:levhp_)                          ,&
               rta              (1:lln22p_,1:levhp_)                           )
#endif
   end subroutine compspec_init 
!-------------------------------------------------------------------------------
#endif /* MP&~DFS end */
   end module compspec
