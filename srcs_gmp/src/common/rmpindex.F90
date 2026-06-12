#include <define.h>
   module rmpindex
!-------------------------------------------------------------------------------
#ifdef RMP
   use paramodel, only  :  igrd12_,jgrd12_,levs_,lnwav_,lngrd_ 
   use rsparltb, only   :  lngrdb
#ifdef MP
   use paramodel, only  :  igrd12p_,jgrd12p_,levsp_,llwavp_,                   &
                           lngrdp_,lnwavp_,npes_
!-------------------------------------------------------------------------------
   private             ::  igrd12p_,jgrd12p_,levsp_,llwavp_,                   &
                           lngrdp_,lnwavp_,npes_
#else
!-------------------------------------------------------------------------------
#endif
   private             ::  igrd12_,jgrd12_,levs_,lnwav_,lngrd_ 
   private             ::  lngrdb
!
   integer, allocatable  ::                                                    &
          igrdn(:,:,:),                                                        &
          jgrdn(:,:,:),                                                        &
          kgrdn(:,:,:),                                                        &
          lnwavn(:,:),                                                         &
          lnwav2n(:,:),                                                        &
          klnwavn(:,:),                                                        &
          rlngrdn(:),                                                          &
#ifdef MP                                        
          lngrdbn(:,:),                                                        &
          igrdnp(:,:,:),                                                       &
          jgrdnp(:,:,:),                                                       &
          kgrdnp(:,:,:),                                                       &
          lnwavnp(:,:),                                                        &
          klnwavnp(:,:),                                                       &
          rlngrdnp(:),                                                         &
          lngrdn2p(:,:),                                                       &
          llwavnp(:,:),                                                        &
          kllwavnp(:,:),                                                       &
          llwav2np(:,:),                                                       &
          lmype(:),                                                            &
          kmype(:),                                                            &
          imype(:),                                                            &
          jmype(:)                                    
#else
          lngrdbn(:,:)                                    
#endif     
!
contains
   subroutine rmpindex_init
   allocate(                                                                   &
          igrdn(igrd12_,jgrd12_,levs_),                                        &
          jgrdn(igrd12_,jgrd12_,levs_),                                        &
          kgrdn(igrd12_,jgrd12_,levs_),                                        &
          lnwavn(lnwav_,levs_)        ,                                        &
          lnwav2n(lnwav_,levs_)       ,                                        &
          klnwavn(lnwav_,levs_)       ,                                        &
          rlngrdn(lngrd_)             ,                                        &
#ifdef MP                                        
          lngrdbn(lngrdb,levs_)       ,                                        &
          igrdnp(igrd12p_,jgrd12p_,levsp_),                                    &
          jgrdnp(igrd12p_,jgrd12p_,levsp_),                                    &
          kgrdnp(igrd12p_,jgrd12p_,levsp_),                                    &
          lnwavnp(lnwavp_,levs_)          ,                                    &
          klnwavnp(lnwavp_,levs_)         ,                                    &
          rlngrdnp(lngrdp_)               ,                                    &
          lngrdn2p(igrd12p_,jgrd12p_)     ,                                    &
          llwavnp(llwavp_,levsp_)         ,                                    &
          kllwavnp(llwavp_,levsp_)        ,                                    &
          llwav2np(llwavp_,levsp_)        ,                                    &
          lmype(0:npes_-1)                ,                                    &
          kmype(0:npes_-1)                ,                                    &
          imype(0:npes_-1)                ,                                    &
          jmype(0:npes_-1)      )
#else
          lngrdbn(lngrdb,levs_)  )                                  
#endif     
!
   end subroutine rmpindex_init
!-------------------------------------------------------------------------------
#endif    /* RMP end */ 
   end module rmpindex
