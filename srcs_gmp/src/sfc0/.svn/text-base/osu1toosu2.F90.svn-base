#include "define.h"
   subroutine osu1toosu2(sfcfcsin,idim,jdim,numsfcsin)
!-------------------------------------------------------------------------------
!
! fill osu2 type surface file records from osu1 type
!
!-------------------------------------------------------------------------------
   use varsfc, only : lalbd_
   use comsfc
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  idim,jdim,numsfcsin
   real                 ::  sfcfcsin(idim,jdim,numsfcsin)
!
   integer              ::  i,j
   integer              ::  ind,k,l
   real                 ::  alog30,undef
!-------------------------------------------------------------------------------
#ifdef OSULSM2
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
   do k=1,2
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
   do k=1,2
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
   if (ind.ne.numsfcsin) then
     print *,'counting error in osu1toosu2'
     call abort
   endif
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
       enddo
     enddo
   enddo
!
   do k = 1,2
     do j = 1,jdim
       do i = 1,idim
         denni(i,j,k) = undef
       enddo
     enddo
   enddo
!
#endif
   return
   end
!-------------------------------------------------------------------------------
