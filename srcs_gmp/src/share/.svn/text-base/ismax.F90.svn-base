#include <define.h>
   function ismax(len,f,inc)
!-------------------------------------------------------------------------------
!
!     find index of the first occurence of max
!
!-------------------------------------------------------------------------------
#ifdef DG
   real                 ::  f(len)
!
   rmax=f(1)
!
   do i = 1,len,inc
     if(f(i).ge.rmax) rmax=f(i)
   enddo
!
   do i = 1,len,inc
     if(f(i).eq.rmax) then
       ismax=i
       return
     endif
   enddo
#else
   ismax=1
#endif
!
   return
   end function ismax
!-------------------------------------------------------------------------------
