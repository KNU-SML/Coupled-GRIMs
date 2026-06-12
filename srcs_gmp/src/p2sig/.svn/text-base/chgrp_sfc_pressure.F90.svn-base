!
   subroutine chgrp_sfc_pressure(t,q,kdimprs,ps,hold,hnew,psnew,si,sl)
!-------------------------------------------------------------------------------
!
! subprogram:    newps       interpolate surface pressure.
!
! abstract: using the hydrostatic equation, the surface pressure
!           is interpolated from input orography to output orography.
!           below the input surface, the temperature lapse rate
!           is fixed at -6.5k/km.  above the input surface,
!           surface pressure is interpolated from height profiles
!           computed as in the mrf model.
!
! program history log:
!   1991-03-15  mark iredell           development
!   2000-03-09  songyou hong           cvs verion setup
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call newps(t,kdimprs,q0,h0,h1,q1,si)
!   input argument list:
!     t        - temperature
!     q        - specific humidity
!     q0       - old ln(psfc)
!     h0       - old orography
!     h1       - new orography
!     si       - edge sigma values
!
!   output argument list:
!     q1       - new ln(psfc)
!
!   subprograms called:
!     sph_matrix_init     - compute mrf temperature to height matrix
!
!-------------------------------------------------------------------------------
   use constant, only : g_,rd_
   use paramter
   use parmchgr
   use comchgr, only  : idimt, jdimhf
!-------------------------------------------------------------------------------
   save
!-------------------------------------------------------------------------------
   real                 ::  t(idimt,kdimprs),ps(idimt),psnew(idimt)
   real                 ::  q(idimt,kdimprs)
   real                 ::  hold(idimt,kdimprs),hnew(idimt),si(kdimprs+1)
   real                 ::  sl(kdimprs)
   real                 ::  tau(idimt,kdimprs)
   real                 ::  rlsig(kdimprs+1)
   real                 ::  rllsig(kdimprs)
   real                 ::  h(idimt,kdimprs+1)
   real                 ::  a(idimt)
   real                 ::  absa(idimt)
!-------------------------------------------------------------------------------
!--------
!-------- compute interface heights.
!--------
!        print *,' ------ into newpsfc ---'
!
   rog = rd_/g_
   do k = 1,kdimprs
     rlsig(k) = - alog(si(k))
     rllsig(k) = - alog(sl(k))
     do i = 1,idimt
       tau(i,k) = (1.+0.61*q(i,k))*t(i,k) * rog
     enddo
   enddo
!
   rlsig(kdimprs+1) = - alog (si(kdimprs+1))
   do i = 1,idimt
     psnew(i)= - ps(i)
     h(i,1)=hold(i,1)
   enddo
   do k = 1,kdimprs
     do i = 1,idimt
       h(i,k+1)=hold(i,k)+tau(i,k)*(rlsig(k+1)-rllsig(k))
     enddo
   enddo
!
!--------
!-------- loop over layers, testing to see if new sfc pressure is
!-------- in layer, and obtaining it, if so.
!--------
!
   eps=1.e-1
   do k = 1,kdimprs
     !----------
     !---------- compute lapse rate
     !----------
     kp=min(k+1,kdimprs)
     km=max(1,k-1)
     kppz=min(k+2,kdimprs)
     kpz=kppz-1
     kmmz=max(1,k-2)
     kmz=kmmz+1
     do iq = 1,idimt
       a(iq)=2.0e0*(tau(iq,kp)-tau(iq,km))/                                    &
            (rlsig(kpz)+rlsig(kppz)-rlsig(kmz)-rlsig(kmmz))
       absa(iq)=abs(a(iq))
     enddo
     do i = 1,idimt
       if(hnew(i).ge.h(i,k).and.hnew(i).le.h(i,k+1).and.                       &
          absa(i).gt.eps) psnew(i)=rlsig(k)                                    &
                                  -ps(i)+(sqrt(tau(i,k)**2+                    &
          2.0e0*a(i)*(hnew(i)-h(i,k))) - tau(i,k) )/a(i)
       if(hnew(i).ge.h(i,k).and.hnew(i).le.h(i,k+1).and.                       &
          absa(i).le.eps) psnew(i)=rlsig(k)                                    &
                                  -ps(i)+(hnew(i)-h(i,k))/tau(i,k)
     enddo
   enddo
!
!--------
!-------- do points which fall below first layer (use fixed lapse rate
!--------  of 6.5 deg per km.
!--------
!
   gamma=6.5e-3
   gascon=rd_
   g=g_
   c=gascon*gamma/g
!
   do iq = 1,idimt
     a(iq)=tau(iq,1)*(1.e0-exp(c*(rlsig(1)-rlsig(2))))/(rlsig(1)-rlsig(2))
   enddo
!
   do i = 1,idimt
     if( hnew(i).le.h(i,1) ) psnew(i)=rlsig(1)                                 &
                                     -ps(i)+(sqrt(tau(i,1)**2+                 &
       2.e0*a(i)*(hnew(i)-h(i,1))) - tau(i,1))/a(i)
   enddo
!
   icount = 0
   do i = 1,idimt
     psnew(i) = - psnew(i)
!cc      if( psnew(i).eq.ps(i) ) icount = icount + 1
     if( psnew(i).eq.ps(i) ) then
       icount = icount + 1
     else
!        print *,'  gz gzold ps psnew ',hold(i),hnew(i),ps(i),psnew(i)
     endif
   enddo
!
!  write(6,300)icount
300 format(' number of surface pressure points not updated=',i6)
!
   return
   end
