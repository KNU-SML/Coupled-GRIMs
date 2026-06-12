#include <define.h>
   subroutine sfcfcnv(sfcftypin,sfcfcsin,idim,jdim,numsfcsin)
!-------------------------------------------------------------------------------
!
! convert input surface sfcftypein file to model sfcftype 
!
!    sfcftypin: character ['osu1','osu2','noa1','vic1']
!    sfcfcsin : real input sfc file record
!    idim,jdim    : integer dimension of sfcfsin array
!    numsfcsin: integer number of input sfc file record
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   character(len=4)  ::  sfcftypin
   integer           ::  idim,jdim,numsfcsin
   real              ::  sfcfcsin(idim,jdim,numsfcsin)
!-------------------------------------------------------------------------------
!
   if(sfcftypin.eq.'osu1') then
#ifdef OSULSM1
     print *,'no surface file type conversion performed'
#endif
#ifdef OSULSM2
     call osu1toosu2(sfcfcsin,idim,jdim,numsfcsin)
#endif
#ifdef NOALSM1
     call osu1tonoa1(sfcfcsin,idim,jdim,numsfcsin)
#endif
#ifdef VICLSM1
     call osu1tovic1(sfcfcsin,idim,jdim,numsfcsin)
#endif
#ifdef VICLSM2
     call osu1tovic2(sfcfcsin,idim,jdim,numsfcsin)
#endif
!
   elseif(sfcftypin.eq.'osu2') then
!
#ifdef OSULSM1
     call osu2toosu1
#endif
#ifdef OSULSM2
     print *,'no surface file type conversion performed'
#endif
#ifdef NOALSM1
     call osu2tonoa1(sfcfcsin,idim,jdim,numsfcsin)
#endif
#ifdef VICLSM1
     call osu2tovic1(sfcfcsin,idim,jdim,numsfcsin)
#endif
#ifdef VICLSM2
     call osu2tovic2(sfcfcsin,idim,jdim,numsfcsin)
#endif
!
   elseif(sfcftypin.eq.'noa1') then
! 
#ifdef OSULSM1
     call noa1toosu1(sfcfcsin,idim,jdim,numsfcsin)
#endif
#ifdef OSULSM2
     call noa1toosu2
#endif
#ifdef NOALSM1
     print *,'no surface file type conversion performed'
#endif
#ifdef VICLSM1
     call noa1tovic1(sfcfcsin,idim,jdim,numsfcsin)
#endif
#ifdef VICLSM2
     call noa1tovic2(sfcfcsin,idim,jdim,numsfcsin)
#endif
!
   elseif(sfcftypin.eq.'vic1')then
!
#ifdef OSULSM1
     call vic1toosu1
#endif
#ifdef OSULSM2
     call vic1toosu2
#endif
#ifdef NOALSM1
     call vic1tonoa1
#endif
#ifdef VICLSM1
     print *,'no surface file type conversion performed'
#endif
#ifdef VICLSM2
     call vic1tovic2(sfcfcsin,idim,jdim,numsfcsin)
#endif
!
   elseif(sfcftypin.eq.'vic2')then
!
#ifdef OSULSM1
     call vic2toosu1
#endif
#ifdef OSULSM2
     call vic2toosu2
#endif
#ifdef NOALSM1
     call vic2tonoa1
#endif
#ifdef VICLSM1
     call vic2tovic1
#endif
#ifdef VICLSM2
     print *,'no surface file type conversion performed'
#endif
!
   else
     print *,'illegal intput surface file type'
     call abort
   endif
!
   return
   end
!-------------------------------------------------------------------------------
