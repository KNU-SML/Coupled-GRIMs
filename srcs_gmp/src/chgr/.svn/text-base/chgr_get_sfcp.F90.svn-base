#include <define.h>
   subroutine chgr_get_sfcp(t,q,ps,hold,hnew,psnew,si)
!-------------------------------------------------------------------------------
!
! subprogram:    chgr_get_sfcp       interpolate surface pressure.
!
! abstract: using the hydrostatic equation, the surface pressure
!           is interpolated from input orography to output orography.
!           below the input surface, the temperature lapse rate
!           is fixed at -6.5k/km.  above the input surface,
!           surface pressure is interpolated from height profiles
!           computed as in the mrf model.
!
! program history log:
!   1991-03-15  mark iredell
!   2000-01-01  song-you hong          cvs version
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call chgr_get_sfcp(t,q0,h0,h1,q1,si,sl)
!   input argument list:
!     t        - temperature
!     q0       - old ln(psfc)
!     h0       - old orography
!     h1       - new orography
!     si       - edge sigma values
!     sl       - full sigma values
!
!   output argument list:
!     q1       - new ln(psfc)
!
!   subprograms called:
!     sph_matrix_init     - compute mrf temperature to height matrix
!
!-------------------------------------------------------------------------------
   use constant, only : g_,rd_
   use parmchgr, only : kdimi
   use paramter, only : idim
   save
!-------------------------------------------------------------------------------
   real                 ::  t(idim*2,kdimi),ps(idim*2),psnew(idim*2)
   real                 ::  q(idim*2,kdimi)
   real                 ::  hold(idim*2),hnew(idim*2),si(idim*2,kdimi+1)
   real                 ::  tau(idim*2,kdimi)
   real                 ::  rlsig(idim*2,kdimi+1)
   real                 ::  h(idim*2,kdimi+1)
   real                 ::  a(idim*2)
   real                 ::  absa(idim*2)
!-------------------------------------------------------------------------------
!
! compute interface heights.
!
!  print *,' ------ into newpsfc ---'
   rog = rd_/g_
   do i = 1,idim*2
     do k = 1,kdimi
       rlsig(i,k) = - alog(si(i,k))
       tau(i,k) = (1.+0.61*q(i,k))*t(i,k) * rog
     enddo
     rlsig(i,kdimi+1) = - alog ( si(i,kdimi+1))
   enddo
!
   do i = 1,idim*2
     psnew(i)= - ps(i)
     h(i,1)=hold(i)
   enddo
!
   do k = 1,kdimi
     do i = 1,idim*2
       h(i,k+1)=h(i,k)+tau(i,k)*(rlsig(i,k+1)-rlsig(i,k))
     enddo
   enddo
!
! loop over layers, testing to see if new sfc pressure is
! in layer, and obtaining it, if so.
!
   eps=1.e-1
   do k = 1,kdimi
!
! compute lapse rate
!
     kp=min(k+1,kdimi)
     km=max(1,k-1)
     kppz=min(k+2,kdimi)
     kpz=kppz-1
     kmmz=max(1,k-2)
     kmz=kmmz+1
     do iq = 1,idim*2
       a(iq)=2.0e0*(tau(iq,kp)-tau(iq,km))/                                    &
               (rlsig(iq,kpz)+rlsig(iq,kppz)-rlsig(iq,kmz)-rlsig(iq,kmmz))
       absa(iq)=abs(a(iq))
     enddo
     do i = 1,idim*2
       if(hnew(i).ge.h(i,k).and.hnew(i).le.h(i,k+1).and.                       &
               absa(i).gt.eps) psnew(i)=rlsig(i,k)                             &
               -ps(i)+(sqrt(tau(i,k)**2+                                       &
               2.0e0*a(i)*(hnew(i)-h(i,k))) - tau(i,k) )/a(i)
       if(hnew(i).ge.h(i,k).and.hnew(i).le.h(i,k+1).and.                       &
               absa(i).le.eps) psnew(i)=rlsig(i,k)                             &
               -ps(i)+(hnew(i)-h(i,k))/tau(i,k)
     enddo
   enddo
!
! do points which fall below first layer (use fixed lapse rate
!  of 6.5 deg per km.
!
   gamma=6.5e-3
   gascon=rd_
   g=g_
   c=gascon*gamma/g
!
   do iq = 1,idim*2
     a(iq)=tau(iq,1)*                                                          &
            (1.e0-exp(c*(rlsig(iq,1)-rlsig(iq,2))))/(rlsig(iq,1)-rlsig(iq,2))
   enddo
!
   do i = 1,idim*2
     if( hnew(i).le.h(i,1) ) psnew(i)=rlsig(i,1)                               &
        -ps(i)+(sqrt(tau(i,1)**2+                                              &
        2.e0*a(i)*(hnew(i)-h(i,1))) - tau(i,1))/a(i)
   enddo
!
   icount = 0
   do i = 1,idim*2
     psnew(i) = - psnew(i)
!    if( psnew(i).eq.ps(i) ) icount = icount + 1
     if( psnew(i).eq.ps(i) ) then
       icount = icount + 1
!    else
!      print *,'  gz gzold ps psnew ',hold(i),hnew(i),ps(i),psnew(i)
     endif
   enddo
!
!  write(6,300)icount
   300  format(' number of surface pressure points not updated=',i6)
!
   return
   end
!-------------------------------------------------------------------------------
