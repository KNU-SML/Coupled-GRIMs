#include <define.h>
   subroutine osu2tovic1(sfcfcsin,idim,jdim,numsfcsin)
!-------------------------------------------------------------------------------
!
! fill vic1 type surface file records from osu2 type
!
!-------------------------------------------------------------------------------
   use varsfc, only : lsoil_,nsoil_
   use comsfc
!-------------------------------------------------------------------------------
   integer        ::  idim,jdim,numsfcsin
   real           ::  sfcfcsin(idim,jdim,numsfcsin)
   integer        ::  ind,i,j,k,l
   real           ::  undef
!-------------------------------------------------------------------------------
#ifdef VICLSM1
!
! tsea
!
   ind=1
   do j = 1,jdim
     do i = 1,idim
       tsea(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! smc
!
   do k = 1,2
     do j = 1,jdim
       do i = 1,idim
         smc(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+2
!
! snow
!
   do j = 1,jdim
     do i = 1,idim
       snoweq(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! stc
!
   do k = 1,2
     do j = 1,jdim
       do i = 1,idim
         stc(i,j,k)=sfcfcsin(i,j,ind+k-1)
           if(k.eq.1) then
             stc(i,j,k)=tsea(i,j)
           endif
         enddo
      enddo
   enddo
   ind=ind+2
!
! tg3
!
   do j = 1,jdim
     do i = 1,idim
       tg3(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! z0cm
!
   do j = 1,jdim
     do i = 1,idim
       z0cm(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! cv
!
   do j = 1,jdim
     do i = 1,idim
       cv(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! cvb
!
   do j = 1,jdim
     do i = 1,idim
       cvb(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! cvt
!
   do j = 1,jdim
     do i = 1,idim
       cvt(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! albedo    !osu2 kalbd=4
!
   do k = 1,4
     do j = 1,jdim
       do i = 1,idim
         albedo(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   end do
   ind=ind+4
!
! slmsk
!
   do j = 1,jdim
     do i = 1,idim
       slmsk(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! vegetation cover
!
   do j = 1,jdim
     do i = 1,idim
       vfrac(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! canopy
!
   do j = 1,jdim
     do i = 1,idim
       canopy(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! f10m
!
   do j = 1,jdim
     do i = 1,idim
       f10m(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! vegitation type
!
   do j = 1,jdim
     do i = 1,idim
       vtype(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! soil type (in VIC no soil type)
!
   ind=ind+1
!
! albedo fraction type
!
   do k = 1,2
     do j = 1,jdim
       do i = 1,idim
         facalf(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+2
!
! ustar
!
   do j = 1,jdim
     do i = 1,idim
       uustar(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! ffmm
!
   do j = 1,jdim
     do i = 1,idim
       ffmm(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! ffhh
!
   do j = 1,jdim
     do i = 1,idim
       ffhh(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
!
   if(ind.ne.numsfcsin) then
     print *,'counting error in osu2tovic1'
     call abort
   endif
!
!  fill 3-rd soil layer moisture and temperature
!
   do k = 3,lsoil_
     do j = 1,jdim
       do i = 1,idim
         smc(i,j,k)=smc(i,j,2)
       enddo
     enddo
   enddo
!
   do k = 3,nsoil_
     do j = 1,jdim
       do i = 1,idim
         stc(i,j,k)=stc(i,j,2)
       enddo
     enddo
   enddo
!
!  initialize srflag from tsea (t850 is not available)
!
   print *,'initialize srflag-use tsea as surrogate for t850'
   do j = 1,jdim
     do i = 1,idim
       srflag(i,j)=0.
       if(tsea(i,j).le.273.16) srflag(i,j)=1.
     enddo
   enddo
!
!  fill prcp by zero
!
   do j = 1,jdim
     do i = 1,idim
       prcp  (i,j)=0.
     enddo
   enddo
!
!  fill vroot, binf, Ds, Dsm, Ws,cef, expt, kest, dph, bub, qrt, 
!       bkd, sld, wcr, wpw, smr, smx, dphn, smxn, expn, bubn, 
!       alpn, betn, gamn, flai, silz, snwz, sic, csno, rsno,
!       tsf, tpk, sfw, pkw, lstsn
!
!  These fields will be replaced with climatology field in sfc program.
!
   undef=1.e30
   do j = 1,jdim
     do i = 1,idim
       binf(i,j) = undef
       Ds  (i,j) = undef
       Dsm (i,j) = undef
       Ws  (i,j) = undef
       cef (i,j) = undef
       silz(i,j) = undef
       snwz(i,j) = undef
       csno(i,j) = undef
       rsno(i,j) = undef
       tsf (i,j) = undef
       tpk (i,j) = undef
       sfw (i,j) = undef
       pkw (i,j) = undef
       lstsn(i,j)= undef
     enddo
   enddo
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         vroot(i,j,k)= undef
         expt(i,j,k) = undef
         kest(i,j,k) = undef
         dph (i,j,k) = undef
         bub (i,j,k) = undef
         qrt (i,j,k) = undef
         bkd (i,j,k) = undef
         sld (i,j,k) = undef
         wcr (i,j,k) = undef
         wpw (i,j,k) = undef
         smr (i,j,k) = undef
         smx (i,j,k) = undef
         sic (i,j,k) = undef
       enddo
     enddo
   enddo
!
   do k = 1,nsoil_
     do j = 1,jdim
       do i = 1,idim
         dphn(i,j,k) = undef
         smxn(i,j,k) = undef
         expn(i,j,k) = undef
         bubn(i,j,k) = undef
         alpn(i,j,k) = undef
         betn(i,j,k) = undef
         gamn(i,j,k) = undef
       enddo
     enddo
   enddo
!
   do j = 1,jdim
     do i = 1,idim
       hml0(i,j)  = undef
       kprfi(i,j) = undef
     enddo
   enddo
!
   do k = 1,5
     do j = 1,jdim
       do i = 1,idim
         paer(i,j,k) = undef
         idxci(i,j,k) = undef
         cmixi(i,j,k) = undef
       end do
     enddo
   enddo
!
   do k = 1,2
     do j = 1,jdim
       do i = 1,idim
         denni(i,j,k) = undef
       end do
     enddo
   enddo
!
#endif
   return
   end
!-------------------------------------------------------------------------------
