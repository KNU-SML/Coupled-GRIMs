#include <define.h>
   subroutine osu2tonoa1(sfcfcsin,idim,jdim,numsfcsin)
!------------------------------------------------------------------------------
!
! fill noa1 type surface file records from osu1 type
!
!------------------------------------------------------------------------------
   use varsfc, only : lalbd_,lsoil_
   use comsfc
!------------------------------------------------------------------------------
   integer        ::  idim,jdim,numsfcsin
   real           ::  sfcfcsin(idim,jdim,numsfcsin)
!
   integer        ::  ind,i,j,k,l,ijdim
   real           ::  alog30,undef
!------------------------------------------------------------------------------
#ifdef NOALSM1
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
   do k = 1,lalbd_
     do j = 1,jdim
       do i = 1,idim
         albedo(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
     ind=ind+1
   enddo
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
! soil type
!
   do j = 1,jdim
     do i = 1,idim
       stype(i,j)=sfcfcsin(i,j,ind)
     enddo
   enddo
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
     ind=ind+1
   enddo
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
     print *,'counting error in osu2tonoa1'
     call abort
   endif
!
!  fill 3-rd and 4-th soil layer moisture and temperature
!
   do k = 3,lsoil_
     do j = 1,jdim
       do i = 1,idim
         smc(i,j,k)=smc(i,j,2)
         stc(i,j,k)=stc(i,j,2)
       enddo
     enddo
   enddo
!
!  initialize snwdph from snoweq
!
   do j = 1,jdim
     do i = 1,idim
       snwdph(i,j)=snoweq(i,j)*5.0
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
!  initialize lsoil layer slc from smc and stc
!
   ijdim=idim*jdim
   call getslc(slmsk,stype,stc,smc,ijdim,lsoil_,slc)
!
!  fill prcp by zero
!
   do j = 1,jdim
     do i = 1,idim
       prcp  (i,j)=0.
     enddo
   enddo
!
!  fill vfrac vegtyp, stype, albdeo, albedo fraction with 1.e30.  
!  These fields will be replaced with climatology field in sfc program.
!
   undef=1.e30
   do j = 1,jdim
     do i = 1,idim
       shdmin(i,j)=undef
       shdmax(i,j)=undef
       snoalb(i,j)=undef
       slope (i,j)=undef
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
