#include <define.h>
   subroutine chgrp_pressure2sigma(ps,pi,sl,din,dou,ijmax,kg,kmax,indx)
!-------------------------------------------------------------------------------
!
! subprogram: chgrp_pressure2sigma
!
! abstract : this program transfers p coordinate to sigma coordinate 
!                                                                               
! input argument list
!       ps       - ln(Ps) ( ln(kPa) )
!       pi       - p-level values (hPa)
!       sl       - sigma values
!       din      - input (p-level)
!       ijmax    - number of horizontal grids
!       kg       - number of vertical grids (p-level)
!       kmax     - number of vertical grids (sigma-level)
!       indx     - variables
!                  linear or log linear interpolation depending variables
!                      linear to p for wind (indx=2 & 3)
!                      log linear to p for mass (indx=4 & 5)
!
! output argument
!       dou      - output (sigma level)
!
! local dimension : wp, wpm
!
! program history           
!     song-you hong        April  1986
!     young-hwa byun       29 Jan 2003                           
!
!-------------------------------------------------------------------------------
   use constant, only : akapa_
!-------------------------------------------------------------------------------
   integer              ::  ijmax,kg,kmax,indx
   real                 ::  ps(ijmax), pi(kg), sl(ijmax,kmax)
   real                 ::  din(ijmax,kg), dou(ijmax,kmax)
   real                 ::  wp(100), wpm(100)
!-------------------------------------------------------------------------------
!
! linear or log-linear to p 
!
   if (indx.le.3 .and. indx.ne.1) then
     do k = 1,kg
       wp(k) = pi(k)*0.1           ! kPa
     enddo
     print *, ' Linear interpolation ...'
   else
     do k = 1,kg
       wp(k) = alog(pi(k)*0.1)     ! ln(kPa)
     enddo
     print *, ' Log-Linear interpolation ...'
   endif
!
   do k = 1,kg-1
     wpm(k) = 0.5*(wp(k)+wp(k+1))
   enddo
!
! for temperature, theta is used to interpolate
!
   if (indx.eq.4) then
     do k = 1,kg
       do ij = 1,ijmax
         theta = din(ij,k)*(pi(k)*0.001)**(-akapa_)
         din(ij,k) = theta
       enddo
     enddo
   endif
!
! vertical interpolation
!
   ncnt1 = 0
   ncnt2 = 0
   ncnt3 = 0
!
   do ij = 1,ijmax
!
     do ks = 1,kmax
!
       psig = sl(ij,ks)*exp(ps(ij))
       ppsig = psig
       if (indx.ge.4) psig = alog(psig)
!
       do k = 1,kg-1
         if (psig.le.wp(k) .and. psig.gt.wp(k+1)) then
           nkp1 = k
           nkp2 = k+1
           if (psig.gt.wpm(k)) then
             nkp3 = k-1
           else
             nkp3 = k+2
           endif
           go to 500
         else
           if (psig.gt.wp(1)) go to 100        ! below lowest p-level
           if (psig.le.wp(kg)) go to 200       ! above p-top
         endif
       enddo
!
500    ncnt1 = ncnt1 + 1
!
       if (nkp1.le.nkp3) then
         nk1 = nkp1 
         nk2 = nkp2
         nk3 = nkp3
       else
         nk1 = nkp3
         nk2 = nkp1
         nk3 = nkp2
       endif
!
       if (nk1.ne.0 .and. nk3.ne.kg+1) then
         wup = 0.5*(din(ij,nk3)+din(ij,nk2))
         wdn = 0.5*(din(ij,nk2)+din(ij,nk1))
         wgt = (wpm(nk1) - psig)/(wpm(nk1) - wpm(nk2))
       elseif (nk1.eq.0) then
         wup = din(ij,nk3)
         wdn = din(ij,nk2)
         wgt = (wp(nk2) - psig)/(wp(nk2) - wp(nk3))
       elseif (nk3.eq.kg+1) then
         wup = din(ij,nk2)
         wdn = din(ij,nk1)
         wgt = (wp(nk1) - psig)/(wp(nk1) - wp(nk2))
       else
         print *, ' Also, something wrong ...? '
       endif
       wsig =  wgt*wup + (1.-wgt)*wdn 
       go to 300
!
! below lowest p-level
!
100    continue
       ncnt2 = ncnt2 + 1
       wsig = din(ij,1) + (din(ij,2)-din(ij,1))/(wp(2)-wp(1))*(psig-wp(1)) 
       go to 300
!
! above p-top 
200    continue
       ncnt3 = ncnt3 + 1
       wsig = din(ij,kg) +                                                     &
              (din(ij,kg)-din(ij,kg-1))/(wp(kg)-wp(kg-1))*(psig-wp(kg)) 
!
300    continue
!
       if (indx.eq.4) wsig = wsig*(ppsig*0.01)**(akapa_)
       dou(ij,ks) = wsig
!
     enddo
!
   enddo
!
   print *, 'chgrp_pressure2sigma ',                                           &
             ncnt1,ncnt2,ncnt3,ncnt1+ncnt2+ncnt3,ijmax*kmax
!
   return                                                                    
   end                                                                       
