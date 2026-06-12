#include <define.h>
   module combgt
!-------------------------------------------------------------------------------
   use paramodel, only : levs_, lnt22_
!-------------------------------------------------------------------------------
!
! common block for diagnostics
!
! the number of diagnostic terms + 2
!
#ifdef T
   parameter ( ntermt = 15 )
   parameter ( nsumt = ntermt - 1 )
   parameter ( nbgtt = ntermt - 2 )
   integer  :: levbgtt
#endif
#ifdef V
   parameter ( ntermv = 14 )
   parameter ( nsumv = ntermv - 1 )
   parameter ( nbgtv = ntermv - 2 )
   integer  :: levbgtv
#endif
#ifdef DA
   parameter ( nterm2 = 2 )
#endif
#ifdef DT
   parameter ( nglpa_start = 1 )
   parameter ( nglpa_end = 5 )
   parameter ( nglpb_start = 6 )
   parameter ( nglpb_end = 13 )
   parameter ( nimp_start = 3 )
   parameter ( nimp_end = 5 )
   parameter ( ndebug = 4 )
#endif
#ifdef DV
   parameter ( nglpa_start = 1 )
   parameter ( nglpa_end = 9 )
   parameter ( nglpb_start = 10 )
   parameter ( nglpb_end = 12 )
   parameter ( ndebug = 2 )
#endif
#ifdef T
   integer, allocatable  ::  levtem (:)
#endif
#ifdef V
   integer, allocatable  :: levvor (:)
#endif
!
! common blocks
!
#ifdef T
   real, allocatable  ::  btermst(:,:,:)
#endif
#ifdef V
   real, allocatable  ::  btermsv(:,:,:)
#endif
#ifdef A
   real, allocatable  ::  ulnr(:,:), vlnr(:,:)
#endif
#ifdef DT
   real, allocatable  ::  bterms2(:,:,:,:)
#endif
#ifdef DV
   real, allocatable  ::  bterms2(:,:,:,:)
#endif
   contains
!-------------------------------------------------------------------------------
   subroutine combgt_init
!-------------------------------------------------------------------------------
#ifdef T
   levbgtt = levs_
   allocate( levtem (levbgtt) )
   allocate( btermst(lnt22_, levbgtt, ntermt) )
#endif
#ifdef V
   levbgtv = levs_
   allocate( levvor (levbgtv) )
   allocate( btermsv(lnt22_, levbgtv, ntermv) )
#endif
#ifdef A
   allocate( ulnr(lnt22_,levs_), vlnr(lnt22_,levs_) )
#endif
#ifdef DT
   allocate( bterms2(lnt22_, levbgtt, nterm2, ndebug) )
#endif
#ifdef DV
   allocate( bterms2(lnt22_, levbgtv, nterm2, ndebug) )
#endif
   end subroutine combgt_init
!-------------------------------------------------------------------------------
   end module combgt
