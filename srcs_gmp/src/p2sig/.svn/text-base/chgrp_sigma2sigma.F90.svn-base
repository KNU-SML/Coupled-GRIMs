#include <define.h>
#define OLD
   subroutine chgrp_sigma2sigma(ps1,s1,var1,ps2,s2,var2,im,jm,km1,km2,         &
                               in,tensn,iuv)
!-------------------------------------------------------------------------------
!
! subprogram: chgrp_sigma2sigma
!
! abstract  : this subroutine transfer from one defined sigma coordinate (s1) to
!             another sigma coordinate (s2).
!
! program histroy log:
!   1993-02-06  henry juang            development
!   2000-03-09  songyou hong           cvs verion setup
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! input    ps1     primary ground pressure   alog(psfc)
!          ps2     secondary gound pressure alog(psfc(cb))
!          s1      primary sigma coordinate
!          var1    primary 3 dimensional variable
!          s2      secondary sigma coordinate (s2(im,jm,km2))
!          var2    secondary 3 dimensional variable
!          im      dimension in x
!          jm      dimension in y
!          km1     dimension in z for primary coordinate and variable
!          km2     dimension in z for secondary coordinate and variable
!          in      control index:   0 for initial transformation
!                                   1 for transformation as previous
!          tensn   factor of tension 0 for cubic spline
!                                   50 for linear interpolation
!          iuv     index for variable
!                 1 : u/v
!                 2 : temp
!                 3 : humidity
!
! local dimensions: st1 st2 sv1 sv2 q vinc
!
! common block /spl/ ovh sh iflag theta
!
! routines related: trispl valts
!
!-------------------------------------------------------------------------------
   use constant, only : qmin_
   save
!-------------------------------------------------------------------------------
   common /spl/ ovh(100), sh(100), iflag, jflag, theta
   real                 ::  s1(km1),s2(im,jm,km2)
   real                 ::  var1(im,jm,km1),var2(im,jm,km2)
   real                 ::  st1(100),st2(100),st0(100),ps1(im,jm),ps2(im,jm)
   real                 ::  sv1(100),sv2(100),q(100),vinc(100)
!
! theta is tension factor
!
   theta=tensn

   if( in .eq. 0 ) then
     in=1
     print *, ' initial values for transformation.'
     iflag = 0
     print *,'  iflag=',iflag,' tension factor=',theta
     do k = 1,km1
       kk=km1-k+1
       if( s1(kk) .gt. 0.0  .and.  s1(kk) .le. 1.0 ) then
         st0(k)=log(s1(kk))
       else
         print *, ' error in s1(k) values at k=',k,' sigma=',s1(kk)
         call abort
       endif
     enddo
!
     do k = 1,km2
       do j = 1,jm
         do i = 1,im
           kk=km2-k+1
           if( s2(i,j,kk) .gt. 0.0  .and.  s2(i,j,kk) .lt. 1.0 ) then
             st2(k)=log(s2(i,j,kk))
           else
             print *, ' error in s2(i,j,k) values at k=',k,'                   &
                        sigma=',s2(i,j,kk)
             call abort
           endif
         enddo
       enddo
     enddo
   endif
!
   do i = 1,im
     do j = 1,jm
!
       do k = 1,km1
         st1(k)=st0(k)+ps1(i,j)-ps2(i,j)
         sv1(k)=var1(i,j,km1-k+1)
       enddo
       do k = 1,km2
         kk=km2-k+1
         st2(k)=log(s2(i,j,kk))
       enddo
       vst=st2(1)
       do k = 1,km2-1
         vinc(k)=st2(k+1)-st2(k)
       enddo
!
       call trispl(km1,st1,q,sv1)
       call valts(q,vst,vinc,st1,sv1,km2,km1,sv2)
!
       kmm=km2/2 + 1
       do k = kmm,1,-1
         !
         ! extrapolate upper boundary of output
         !
         if( sv2(k) .eq. 99999.9 ) then
           do kk = 1,km1
             if(st1(kk).gt.st2(k)) then
               diff1=st1(kk)-st2(k)
               if(kk.gt.1) then
                 diff2=st1(kk-1)-st2(k)
               else
                 diff2=999.
               endif
#ifndef OLD
               sv2(k)=sv1(kk)+(sv1(kk)-sv2(kk+1))/                             &
                        (st1(kk)-st1(kk+1))*(st2(k)-st1(kk))
#else
               if(diff1.lt.diff2) then
                 sv2(k)=sv1(kk)*exp((st2(k)-st1(kk)))
               else
                 sv2(k)=sv1(kk-1)*exp((st2(k)-st1(kk-1)))
               endif
#endif
               exit
             endif
           enddo
           if (kk.gt.km1) then
             print*,'extp uppert error'
               print*,'st1=',exp(st1(1:km1))
               print*,'st2=',exp(st2(1:km2))
               call abort
           endif
         endif
       enddo ! k
!
       do k = kmm,km2
         !
         ! extrapolate lower boundary of output
         !
         if( sv2(k) .eq. 99999.9 ) then
           do kk = km1,1,-1
             if(st1(kk).lt.st2(k)) then
               diff1=abs(st1(kk)-st2(k))
               if(kk.lt.km1) then
                 diff2=abs(st1(kk+1)-st2(k))
               else
                 diff2=999.
               endif
#ifndef OLD
               sv2(k)=sv1(kk)+(sv1(kk)-sv2(kk-1))/                             &
                        (st1(kk)-st1(kk-1))*(st2(k)-st1(kk))
#else
               if(diff1.lt.diff2) then
                 sv2(k)=sv1(kk)
               else
                 sv2(k)=sv1(kk+1)
               endif
#endif
               exit
             endif
           enddo
           if (kk.lt.1) then
             sv2(k)=sv1(km1)
           endif
         endif
       enddo ! k
!
       if (iuv.eq.3) then
         do k = 1,km2
           var2(i,j,k)=max(qmin_,sv2(km2-k+1))
         enddo
       else
         do k = 1,km2
           var2(i,j,k)=sv2(km2-k+1)
         enddo
       endif

     enddo
   enddo
!
   return
   end subroutine chgrp_sigma2sigma
