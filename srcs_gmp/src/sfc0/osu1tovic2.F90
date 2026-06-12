#include <define.h>
   subroutine osu1tovic2(sfcfcsin,idim,jdim,numsfcsin)
!------------------------------------------------------------------------------
!
! fill vic1 type surface file records from osu1 type
!
!------------------------------------------------------------------------------
   use varsfc, only : kslmb_,nslmb_,lsoil_,msub_,nsoil_,lalbd_
   use comsfc
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer               ::  idim,jdim,numsfcsin
   real                  ::  sfcfcsin(idim,jdim,numsfcsin)
!
   integer               ::  ind,i,j,k,l,nv
   real                  ::  alog30,undef
!------------------------------------------------------------------------------
#ifdef VICLSM2
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
         if(k.eq.1) then
           stc(i,j,k) = tsea(i,j)  ! the first node is ground surface
         else
           stc(i,j,k)=sfcfcsin(i,j,ind+k-1-1)
           stc(i,j,k+1)=sfcfcsin(i,j,ind+k-1)
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
! albedo
!
   do j = 1,jdim
     do i = 1,idim
       albedo(i,j,1)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
! slmsk
!
   do j = 1,jdim
     do i = 1,idim
       slmsk(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
!
   ind=ind+1
#ifdef DBG
   open(166,file='CHECK_ini_temp.dat',status='unknown')
   do j = 1,jdim
     do i = 1,idim
       if(slmsk(i,j).eq.1)then
         if(abs(stc(i,j,1)-stc(i,j,2)).gt.5.0) then
           write(166,1166) i,j,stc(i,j,1),stc(i,j,2),stc(i,j,3)
1166       format('cell id=',2i5,1x,'temp=',3f15.3)
           if(stc(i,j,1).gt.stc(i,j,2)) then
             stc(i,j,2) = stc(i,j,1) - 5.0
           else
             stc(i,j,2) = stc(i,j,1) + 5.0
           endif
         endif
       endif
     enddo
   enddo
   close(166)
#endif
!
! plantr to be discarded
!
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
!
   if(ind.ne.numsfcsin) then
     print *,'counting error in osu1toosu2'
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
   do k = 4,nsoil_
     do j = 1,jdim
       do i = 1,idim
         stc(i,j,k)=stc(i,j,3)
       enddo
     enddo
   enddo
!
!  fill subgrid canopy water
!
   do j = 1,jdim
     do i = 1,idim
       if(msub_.gt.1) then
         do nv = 1,msub_
           mcnp(i,j,nv) = canopy(i,j)
         end do
       endif
     enddo
   enddo
!
!  uustar, ffmm, ffhh, canopy are filled with reasonable constant
!
   alog30=log(30.)
   do j = 1,jdim
     do i = 1,idim
       uustar(i,j)=1.
       ffmm(i,j)=alog30
       ffhh(i,j)=alog30
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
!  fill vfrac vegtyp, albdeo, albedo fraction with 1.e30.  
!       vroot, binf, Ds, Dsm, Ws,cef, expt, kst, dph, bub, qrt,
!       bkd, sld, wcr, wpw, smr, smx, dphn, smxn, expn, bubn, 
!       alpn, betn, gamn, silz, snwz,nveg,flai,msno,msmc,msic,
!       mstc, csno, rsno,tsf, tpk, sfw, pkw, lstsn
!
!  These fields will be replaced with climatology field 
!                        by sfc_merge_field_driver call
!
   undef=1.e30
   do j = 1,jdim
     do i = 1,idim
       vfrac(i,j)=undef
       vtype(i,j)=undef
       binf(i,j) = undef
       Ds  (i,j) = undef
       Dsm (i,j) = undef
       Ws  (i,j) = undef
       cef (i,j) = undef
       silz(i,j) = undef
       snwz(i,j) = undef
       nveg(i,j) = undef
     enddo
   enddo
!
   do nv = 1,msub_
     do j = 1,jdim
       do i = 1,idim
         mvfr(i,j,nv)=undef
         mvty(i,j,nv)=undef
       enddo
     enddo
   enddo
!
   do l = 1,2
     do j = 1,jdim
       do i = 1,idim
         facalf(i,j,l)=undef
       enddo
     enddo
   enddo
!
! modified by hoon from 1 to 2
!
   do l = 2,lalbd_
     do j = 1,jdim
       do i = 1,idim
         albedo(i,j,l)=undef
       enddo
     enddo
   enddo
!
   do nv = 1,msub_
     do j = 1,jdim
       do i = 1,idim
         csno(i,j,nv) = undef
         rsno(i,j,nv) = undef
         tsf (i,j,nv) = undef
         tpk (i,j,nv) = undef
         sfw (i,j,nv) = undef
         pkw (i,j,nv) = undef
         lstsn(i,j,nv)= undef
       enddo
     enddo
   enddo
!
   do k = 1,kslmb_
     do j = 1,jdim
       do i = 1,idim
         vroot(i,j,k) = undef
         msmc (i,j,k) = undef
         msic (i,j,k) = undef
       enddo
     enddo
   enddo
!
   do k = 1,nslmb_
     do j = 1,jdim
       do i = 1,idim
         mstc(i,j,k) = undef
       enddo
     enddo
   enddo
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         expt(i,j,k) = undef
         kst (i,j,k) = undef
         dph (i,j,k) = undef
         bub (i,j,k) = undef
         qrt (i,j,k) = undef
         bkd (i,j,k) = undef
         sld (i,j,k) = undef
         wcr (i,j,k) = undef
         wpw (i,j,k) = undef
         smr (i,j,k) = undef
         smx (i,j,k) = undef
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
!------------------------------------------------------------------------------
