#include <define.h>
   subroutine noa1tovic2(sfcfcsin,idim,jdim,numsfcsin)
!-------------------------------------------------------------------------------
   use comsfc
   use varsfc, only : lsoil_,msub_,nsoil_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! fill vic2 type surface file records from noa1 type
!
   integer             ::  idim,jdim,numsfcsin
   real                ::  sfcfcsin(idim,jdim,numsfcsin)
!
   integer             ::  ind,ij,k,l,nv
   real                ::  undef
#ifdef VICLSM2
   integer, parameter  ::  kslmb_=lsoil_*msub_, nslmb_=nsoil_*msub_
!-------------------------------------------------------------------------------
!
! tsea
!
   ind=1
   do j = 1,jdim
     do i = 1,idim
       tsea(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! smc
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         smc(ij,k)=sfcfcsin(ij,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+4
!
! snow
!
   do j = 1,jdim
     do i = 1,idim
       snoweq(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! stc
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         stc(ij,k)=sfcfcsin(ij,ind+k-1)
         if(k.eq.1) then
           stc(ij,k)=tsea(ij)
         else if (k.eq.4) then
           stc(ij,k)=sfcfcsin(ij,ind+k-1-1)
           stc(ij,k+1)=sfcfcsin(ij,ind+k-1)
         else
           stc(ij,k)=sfcfcsin(ij,ind+k-1-1)
         endif
       enddo
     enddo
   enddo
   ind=ind+4
!
! tg3
!
   do j = 1,jdim
     do i = 1,idim
       tg3(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! z0cm
!
   do j = 1,jdim
     do i = 1,idim
       z0cm(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! cv
!
   do j = 1,jdim
     do i = 1,idim
       cv(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! cvb
!
   do j = 1,jdim
     do i = 1,idim
       cvb(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! cvt
!
   do j = 1,jdim
     do i = 1,idim
       cvt(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! albedo    !osu2 kalbd=4
!
   do k = 1,4
     do j = 1,jdim
       do i = 1,idim
         albedo(ij,k)=sfcfcsin(ij,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+4
!
! slmsk
!
   do j = 1,jdim
     do i = 1,idim
       slmsk(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! vegetation cover
!
   do j = 1,jdim
     do i = 1,idim
       vfrac(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! canopy
!
   do j = 1,jdim
     do i = 1,idim
       canopy(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! f10m
!
   do j = 1,jdim
     do i = 1,idim
       f10m(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! vegitation type
!
   do j = 1,jdim
     do i = 1,idim
       vtype(ij)=sfcfcsin(ij,ind)
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
         facalf(ij,k)=sfcfcsin(ij,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+2
!
! ustar
!
   do j = 1,jdim
     do i = 1,idim
       uustar(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! ffmm
!
   do j = 1,jdim
     do i = 1,idim
       ffmm(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
! ffhh
!
   do j = 1,jdim
     do i = 1,idim
       ffhh(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
!  prcp
!
   do j = 1,jdim
     do i = 1,idim
       prcp(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
   ind=ind+1
!
!  srflag
!
   do j = 1,jdim
     do i = 1,idim
       srflag(ij)=sfcfcsin(ij,ind)
     enddo
   enddo
!
   ind = ind + 9
!
   if(ind.ne.numsfcsin) then
     print *,'counting error in noa1tovic1'
     call abort
   endif
!
!  fill 3-rd soil layer moisture and temperature
!
   do k = 5,nsoil_
     do j = 1,jdim
       do i = 1,idim
         stc(ij,k)=stc(ij,5)
       enddo
     enddo
   enddo
!
!  fill subgrid canopy water
!
   do nv = 2, msub_
     do j = 1,jdim
       do i = 1,idim
         mcnp(ij,nv) = canopy(ij)
       end do
     enddo
   enddo
!
!  fill vroot, binf, Ds, Dsm, Ws,cef, expt, kst, dph, bub, qrt, 
!       bkd, sld, wcr, wpw, smr, smx, dphn, smxn, expn, bubn, 
!       alpn, betn, gamn, silz, snwz,nveg,flai,msno,msmc,msic,
!       mstc, csno, rsno,tsf, tpk, sfw, pkw, lstsn
!
!  These fields will be replaced with climatology field in sfc program.
!
   undef=1.e30
   do nv = 1,msub_
     do j = 1,jdim
       do i = 1,idim
         mvfr(ij,nv)=undef
         mvty(ij,nv)=undef
       enddo
     enddo
   enddo
!
   do j = 1,jdim
     do i = 1,idim
       binf(ij) = undef
       Ds  (ij) = undef
       Dsm (ij) = undef
       Ws  (ij) = undef
       cef (ij) = undef
       silz(ij) = undef
       snwz(ij) = undef
       nveg(ij) = undef
     enddo
   enddo
!
   do nv = 1,msub_
     do j = 1,jdim
       do i = 1,idim
         csno(ij,nv) = undef
         rsno(ij,nv) = undef
         tsf (ij,nv) = undef
         tpk (ij,nv) = undef
         sfw (ij,nv) = undef
         pkw (ij,nv) = undef
         lstsn(ij,nv)= undef
       end do
     enddo
   enddo
!
   do k = 1,kslmb_
     do j = 1,jdim
       do i = 1,idim
         vroot(ij,k) = undef
         msmc (ij,k) = undef
         msic (ij,k) = undef
       enddo
     enddo
   enddo
!
   do k = 1,nslmb_
     do j = 1,jdim
       do i = 1,idim
         mstc(ij,k) = undef
       enddo
     enddo
   enddo
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         expt(ij,k) = undef
         kst (ij,k) = undef
         dph (ij,k) = undef
         bub (ij,k) = undef
         qrt (ij,k) = undef
         bkd (ij,k) = undef
         sld (ij,k) = undef
         wcr (ij,k) = undef
         wpw (ij,k) = undef
         smr (ij,k) = undef
         smx (ij,k) = undef
       end do
     enddo
   enddo
!
   do k = 1,nsoil_
     do j = 1,jdim
       do i = 1,idim
         dphn(ij,k) = undef
         smxn(ij,k) = undef
         expn(ij,k) = undef
         bubn(ij,k) = undef
         alpn(ij,k) = undef
         betn(ij,k) = undef
         gamn(ij,k) = undef
       end do
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
         paer(i,j,k)  = undef
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
