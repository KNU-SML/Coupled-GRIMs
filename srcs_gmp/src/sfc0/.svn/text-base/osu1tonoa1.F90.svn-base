#include <define.h>
   subroutine osu1tonoa1(sfcfcsin,idim,jdim,numsfcsin)
!-------------------------------------------------------------------------------
!
! fill noa1 type surface file records from osu1 type
!
!-------------------------------------------------------------------------------
   use comsfc
   use varsfc, only : lsoil_,lalbd_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  idim,jdim,numsfcsin,i,j
   integer              ::  ind,k,l
   real                 ::  sfcfcsin(idim,jdim,numsfcsin)
   real                 ::  alog30,undef
!
   integer              ::  ijdim
!-------------------------------------------------------------------------------
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
   ind=ind+1
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
!  uustar, ffmm, ffhh, canopy are filled with reasonable constant
!
   alog30=log(30.)
   do j = 1,jdim
     do i = 1,idim
       uustar(i,j)=1.
       ffmm(i,j)=alog30
       ffhh(i,j)=alog30
       prcp(i,j)=0.
     enddo
   enddo
!
!  fill vfrac vegtyp, stype, albdeo, albedo fraction with 1.e30.  
!  These fields will be replaced with climatology field in sfc program.
!
   undef=1.e30
   do j = 1,jdim
     do i = 1,idim
       vfrac(i,j)=undef
       vtype(i,j)=undef
       stype(i,j)=undef
       do l=1,2
         facalf(i,j,l)=undef
       enddo
       do l=1,lalbd_
         albedo(i,j,l)=undef
       enddo
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
