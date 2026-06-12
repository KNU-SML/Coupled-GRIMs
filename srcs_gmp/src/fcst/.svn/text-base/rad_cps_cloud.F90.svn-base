#include <define.h>
   subroutine rad_cps_cloud
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [phys_main_solver] --- [rad_cps_cloud]
!                               |
!                               |--- [cps_cloudiness_slingo] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   end subroutine rad_cps_cloud
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine cps_cloudiness_slingo(ims2,imx2,clstp,                          &
                                     rn,kbot,ktop,cv,cvb,cvt,count)
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS
!-------------------------------------------------------------------------------
!
! abstract: computes convective cloud cover and cloud tops and bottoms
!   after the deep convection is invoked. cloud cover is interpolated
!   from a table relating cloud cover to precipitation rate.
!
! program history log:
!   1991-05-07  iredell                development
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call cps_cloudiness_slingo(cmean,lat,iistp,dt,rn,kbot,ktop,cv,cvb,cvt)
!
!   input argument list:
!     cmean    - real flag (ge 0 to accumulate, eq 99 to return values)
!     lat      - integer latitude index
!     iistp    - integer time step number
!     dt       - real time step in seconds
!     rn       - real (nx) convective rain in meters
!     kbot     - integer (nx) cloud bottom level
!     ktop     - integer (nx) cloud top level
!
!   output argument list:
!     cv       - real (nx,ny) convective cloud cover
!     cvb      - real (nx,ny) convective cloud base level
!     cvt      - real (nx,ny) convective cloud top level
!
!-------------------------------------------------------------------------------
   real                 ::  rn(imx2),cv(imx2),cvb(imx2),cvt(imx2)
   integer              ::  kbot(imx2),ktop(imx2)
!
!  local work variables and arrays
!
   integer              ::  nmd(ILOTS)
   real                 ::  pmd(ILOTS)
!
!  local save variables and arrays
!
   integer,parameter    :: ncc=9
   real                 :: cc(ncc),p(ncc),cvb0
   data cc/0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/
   data p/.14,.31,.70,1.6,3.4,7.7,17.,38.,85./
   data cvb0/100./
!-------------------------------------------------------------------------------
   im=ims2
!
!  initialize convective rain and range
!
   if(clstp.le.0..and.clstp.gt.-10.) then
     do i = 1,im
       cv(i)=0.
       cvb(i)=cvb0
       cvt(i)=0.
!         cvb(i)=0.
!         cvt(i)=0.
     enddo
   endif
!
!  accumulate convective rain and range
!
   if(clstp.gt.-99. .and. count.ne.0. ) then
     do i = 1,im
       if(rn(i).gt.0.) then
         cv(i)=cv(i)+rn(i)
         cvb(i)=min(cvb(i),float(kbot(i)))
         cvt(i)=max(cvt(i),float(ktop(i)))
!           cvt(i)=max(cvt(i),float(ktop(i)+1))
!           cvb(i)=cvb(i)+kbot(i)*rn(i)
!           cvt(i)=cvt(i)+(ktop(i))*rn(i)
!           cvt(i)=cvt(i)+(ktop(i)+1)*rn(i)
       endif
     enddo
   endif
!
!  convert precipitation rate into cloud fraction
!
   if (clstp.gt.0..or.(clstp.lt.0.and.clstp.gt.-10.)) then
     do i = 1,im
       if(cv(i).gt.0.) then
!           cvb(i)=nint(cvb(i)/cv(i))
!           cvt(i)=nint(cvt(i)/cv(i))
       else
         cvb(i)=cvb0
         cvt(i)=0.
       endif
       pmd(i)=cv(i)*(24.e+3/abs(clstp))
       nmd(i)=0
     enddo
     do n = 1,ncc
       do i = 1,im
         if(pmd(i).gt.p(n)) nmd(i)=n
       enddo
     enddo
     do i = 1,im
       if(nmd(i).eq.0) then
         cv(i)=0.
         cvb(i)=cvb0
         cvt(i)=0.
       elseif(nmd(i).eq.ncc) then
         cv(i)=cc(ncc)
       else
         cc1=cc(nmd(i))
         cc2=cc(nmd(i)+1)
         p1=p(nmd(i))
         p2=p(nmd(i)+1)
         cv(i)=cc1+(cc2-cc1)*(pmd(i)-p1)/(p2-p1)
       endif
     enddo
   endif
!
   return
   end subroutine cps_cloudiness_slingo
!-------------------------------------------------------------------------------
