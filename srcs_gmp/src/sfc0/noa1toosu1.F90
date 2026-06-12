#include <define.h>
   subroutine noa1toosu1(sfcfcsin,idim,jdim,numsfcsin)
!-------------------------------------------------------------------------------
   use comsfc
   use varsfc, only : lsoil_,lalbd_
!
! fill osu1 type surface file records from noa1 type
!
   integer   ::  idim,jdim,numsfcsin
   real      ::  sfcfcsin(idim,jdim,numsfcsin)
!
   integer   ::  ind,i,j,k,l
   real      ::  undef
!-------------------------------------------------------------------------------
#ifdef OSULSM1
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
   ind=ind+4
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
   do k = 1,lsoil_
     do j = 1,jdim
       do i = 1,idim
         stc(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
   ind=ind+4
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
! albedo    !osu1 kalbd=1
!
   do k = 1,lalbd_
     do j = 1,jdim
       do i = 1,idim
         albedo(i,j,k)=sfcfcsin(i,j,ind+k-1)
       enddo
     enddo
   enddo
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
! vegetation cover to be discarded
!
   ind=ind+1
!  do j=1,jdim
!    do i=1,idim
!         vfrac(i,j)=sfcfcsin(i,j,ind)
!    enddo
!  enddo
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
   ind=ind+37
!
   if(ind.ne.numsfcsin) then
      print *,'counting error in noa1toosu1'
      call abort
   endif
!
!  fill plantr with 1.e30.
!
!  These fields will be replaced with climatology field in sfc program.
!
   undef=1.e30
   do j = 1,jdim
     do i = 1,idim
       plantr(i,j) = undef
     enddo
   enddo
!
#endif
   return
   end
!-------------------------------------------------------------------------------
