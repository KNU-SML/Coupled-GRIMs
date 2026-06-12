#include <define.h>
   subroutine dyn_sigma_setup(ci, si, del, sl, cl, rpi)
!-------------------------------------------------------------------------------
!
! subprogram: dyn_sigma_setup
!
!   sigma file contains si and sl, so
!   input: si sl
!   output: ci cl del rpi
!
!-------------------------------------------------------------------------------
   use constant, only : cp_,rd_
   use parmchgr, only : kdimi
   use comchgr, only : kdimpi, kdimmi
   save
!-------------------------------------------------------------------------------
   real                 ::  ci(kdimpi),si(kdimpi),del(kdimi),                  &
                            sl(kdimi),cl(kdimi),                               &
                            rpi(kdimmi)
   real                 ::  rk,rk1,rkinv
!-------------------------------------------------------------------------------
!
   rk  = rd_/cp_
   rk1 = rk + 1.
   rkinv=1./rk
!
   do k = 1,kdimpi
     ci(k)=1. - si(k)
   enddo
!
   do k = 1,kdimi
     del(k)=ci(k+1)-ci(k)
   enddo
!
   do k = 1,kdimi
     cl(k)=1. - sl(k)
   enddo
!
!     compute pi ratios for temp. matrix.
!
   do le = 1,kdimmi
     rpi(le) = (sl(le+1)/sl(le))
   enddo
!
   do le = 1,kdimmi
     rpi(le) = rpi(le)**rk
   enddo
!
   do le = 1,kdimpi
     print 100, le, ci(le), si(le)
100  format (1h , 'level=', i2, 2x, 'ci=', f6.3, 2x, 'si=', f6.3)
   enddo
!
   print 200
   200 format (1h0)
!
   do le = 1,kdimi
     print 300, le, cl(le), sl(le), del(le)
300  format (1h , 'layer=', i2, 2x, 'cl=', f6.3, 2x, 'sl=', f6.3, 2x,          &
                   'del=', f6.3)
   enddo
!
   print 400, (rpi(le), le=1,kdimmi)
   400 format (1h0, 'rpi=', (18(1x,f6.3)) )
!
   return
   end
!-------------------------------------------------------------------------------
