#include "define.h"
   subroutine rad_ozone_setup (idate, fhour, ko3p, ko3l)
!-------------------------------------------------------------------------------
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use comfcst, only : jo3, ko3, prdin, disin
   use comio
#ifdef MP
   use commpi
#endif
#if defined(DFS)
   use dfsvar, only : iope
#endif
!-------------------------------------------------------------------------------
   real, parameter      ::  blat=-85.0,dphi=10.0
   integer              ::  idate(4)
   data nskipo /-1/
   save nskipo
!-------------------------------------------------------------------------------
!fpp$ noconcur r
!
   call incdte(idate(4),idate(2),idate(3),idate(1),                            &
        jyy,jmm,jdd,jhh,nint(fhour))
   mon=jmm
   iday=jdd
!
   nskip = (mon-1)*3 + (iday-1) / 10
   if (iday .eq. 31) nskip = nskip -1
!
   if (nskip .ne. nskipo) then
     if( iope ) then
       rewind (ko3p)
       rewind (ko3l)
       do i = 1,nskip
         read (ko3p,*) prdin
         read (ko3l,*) disin
       enddo
       read (ko3p,*) prdin
       read (ko3l,*) disin
     endif
#ifdef MP
#ifdef RMP
     call rmpbcastr(prdin,jo3*ko3)
     call rmpbcastr(disin,jo3*ko3)
#else
     call mpbcastr(prdin,jo3*ko3)
     call mpbcastr(disin,jo3*ko3)
#endif
#endif
!
     nskipo = nskip
   endif
!
   return
   end subroutine rad_ozone_setup
