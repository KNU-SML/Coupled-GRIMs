#include <define.h>
   subroutine newsig2(si,sl,del)
!-------------------------------------------------------------------------------
   use paramodel, only : kdim=>levs_, kdimp=>levp1_, kdimm=>levm1_
   use constant, only : cp_,rd_
!-------------------------------------------------------------------------------
!
! this routine sets coordinates for levels and use phillips method
!
   real                 ::  si(kdimp),sl(kdim),del(kdim)
   real                 ::  rk,rk1,rkinv
   real, allocatable    ::  ci(:),cl(:),rpi(:)
   real, dimension(100) ::  delmdl
!-------------------------------------------------------------------------------
   namelist /modlsig/ delmdl
   delmdl = 0.
   read(1,modlsig)
!
   allocate (ci(kdimp),cl(kdim),rpi(kdimm))
!
   ci(1) = 0.
   do k = 1,kdim
     del(k)=delmdl(k)
     ci(k+1)=ci(k)+delmdl(k)
   enddo
   ci(kdimp)=1.
!
   rk  = rd_/cp_
   rk1 = rk + 1.
   rkinv=1./rk
!
   levs=kdim
!
   do li = 1,kdimp
     si(li) = 1. - ci(li)
   enddo
!
   do le = 1,kdim
     dif = si(le)**rk1 - si(le+1)**rk1
     dif = dif / (rk1*(si(le)-si(le+1)))
     sl(le) = dif**rkinv
     cl(le) = 1. - sl(le)
   enddo
!
!     compute pi ratios for temp. matrix.
!
   do le = 1,kdimm
     rpi(le) = (sl(le+1)/sl(le))
   enddo
   do le = 1,kdimm
     rpi(le) = rpi(le)**rk
   enddo
!
   do le = 1,kdimp
     write(6,"(1h , 'level=', i2, 2x, 'ci=', f6.3, 2x, 'si=', f6.3)")          &
           le, ci(le), si(le)
   enddo
!
   write(6,'(1h0)')
   do le=1,kdim
     write(6,"(1h , 'layer=', i2, 2x, 'cl=', f6.3, 2x, 'sl=', f6.3, 2x,        &
       'del=', f6.3)") le, cl(le), sl(le), del(le)
   enddo
!
   write(6,"(1h0, 'rpi=', (18(1x,f6.3)))")(rpi(le), le=1,kdimm)
!
   deallocate (ci,cl,rpi)
!
   return
   end
