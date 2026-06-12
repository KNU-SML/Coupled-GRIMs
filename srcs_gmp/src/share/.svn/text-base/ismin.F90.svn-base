#include <define.h>
   function ismin(len,f,inc)
!-------------------------------------------------------------------------------
!
!     find index of the first occurence of min
!
!-------------------------------------------------------------------------------
#ifdef DG
   real                 ::  f(len)
#endif
!
   ismin=1
#ifdef DG
   rmin=f(1)
!
   do i = 1,len,inc
     if(f(i).le.rmin) rmin=f(i)
   enddo
!
   do i = 1,len,inc
     if(f(i).eq.rmin) then
       ismin=i
       return
     endif
   enddo
#endif
!
   return
   end function ismin
!-------------------------------------------------------------------------------
