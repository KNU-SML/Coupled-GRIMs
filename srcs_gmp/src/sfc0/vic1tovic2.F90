#include <define.h>
   subroutine vic1tovic2(sfcfcsin,idim,jdim,numsfcsin)
!-------------------------------------------------------------------------------
   use varsfc, only : kslmb_,nslmb_,lsoil_,msub_,nsoil_,lalbd_
   use comsfc
!-------------------------------------------------------------------------------
!
! fill vic2 type surface file records from noa1 type
!
   integer        ::  idim,jdim,numsfcsin
   real           ::  sfcfcsin(idim,jdim,numsfcsin)
!
   integer        ::  ind,i,j,k,l,nv
   real           ::  undef
!-------------------------------------------------------------------------------
!
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
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         smc(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
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
   do k = 1,nsoil_
     do j = 1,jdim
       do i = 1,idim
         stc(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+nsoil_
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
! albedo    !vic1 kalbd=4
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
   ind=ind+1
!
!  prcp
!
   do j = 1,jdim
     do i = 1,idim
       prcp(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
!  srflag
!
   do j = 1,jdim
     do i = 1,idim
       srflag(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
!  binf
!
   do j = 1,jdim
     do i = 1,idim
       binf(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
!  Ds
!
   do j = 1,jdim
     do i = 1,idim
       Ds(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
!  Dsm
!
   do j = 1,jdim
     do i = 1,idim
       Dsm(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
!  Ws
!
   do j = 1,jdim
     do i = 1,idim
       Ws(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
!  cef
!
   do j = 1,jdim
     do i = 1,idim
       cef(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind=ind+1
!
!  expt
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         expt(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
!
!  kst
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         kst(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
!
!  dph
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         dph(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
!
!  bub
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         bub(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
!
!  qrt
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         qrt(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
!
!  bkd
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         bkd(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
!
!  sld
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         sld(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
!
!  wcr
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         wcr(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
!
!  wpw
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         wpw(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
!
!  smr
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         smr(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
!
!  smx
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         smx(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
!
!  dphn
!
   do k = 1,nsoil_
     do j = 1,jdim
       do i = 1,idim
         dphn(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+nsoil_
!
!  smxn
!
   do k = 1,nsoil_
     do j = 1,jdim
       do i = 1,idim
         smxn(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+nsoil_
!
!  expn
!
   do k = 1,nsoil_
     do j = 1,jdim
       do i = 1,idim
         expn(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+nsoil_
!
!  bubn
!
   do k = 1,nsoil_
     do j = 1,jdim
       do i = 1,idim
         bubn(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+nsoil_
!
!  alpn
!
   do k = 1,nsoil_
     do j = 1,jdim
       do i = 1,idim
         alpn(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+nsoil_
!
!  flai
!
   ind = ind + 1
!
!  silz
!
   do j = 1,jdim
     do i = 1,idim
       silz(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind = ind + 1
!
!  snwz
!
   do j = 1,jdim
     do i = 1,idim
       snwz(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind = ind + 1
!
!  sic
!
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         msic(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+lsoil_
!
!  csno
!
   do j = 1,jdim
     do i = 1,idim
       csno(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind = ind + 1
!
!  rsno
!
   do j = 1,jdim
     do i = 1,idim
       rsno(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind = ind + 1
!
!  tsf
!
   do j = 1,jdim
     do i = 1,idim
       tsf(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind = ind + 1
!
!  tpk
!
   do j = 1,jdim
     do i = 1,idim
       tpk(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind = ind + 1
!
!  sfw
!
   do j = 1,jdim
     do i = 1,idim
       sfw(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind = ind + 1
!
!  pkw
!
   do j = 1,jdim
     do i = 1,idim
       pkw(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind = ind + 1
!
!  lstsn
!
   do j = 1,jdim
     do i = 1,idim
       lstsn(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
   ind = ind + 1
!
   if(ind.ne.numsfcsin) then
     print *,'counting error in vic1tovic2'
     call abort
   endif
!
!  fill subgrid variables
!
   if(msub_.gt.1) then
     do nv = 1, msub_
       do j = 1,jdim
         do i = 1,idim
           mcnp(i,j,nv) = canopy(i,j)
           msno(i,j,nv) = msno(i,j)
           csno(i,j,nv) = csno(i,j)
           rsno(i,j,nv) = rsno(i,j)
           tsf(i,j,nv)  = tsf(i,j)
           tpk(i,j,nv)  = tpk(i,j)
           sfw(i,j,nv)  = sfw(i,j)
           pkw(i,j,nv)  = pkw(i,j)
           lstsn(i,j,nv)= lstsn(i,j)
         enddo
       enddo
     enddo
   endif
!
!  fill vroot, binf, nveg,flai,msno,msmc,msic,mstc
!
!  These fields will be replaced with climatology field in sfc program.
!
   undef=1.e30
   do nv = 1,msub_
     do j = 1,jdim
       do i = 1,idim
            mvfr(i,j,nv)=undef
            mvty(i,j,nv)=undef
       enddo
     enddo
   enddo
!
   do j = 1,jdim
     do i = 1,idim
       nveg(i,j) = undef
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

