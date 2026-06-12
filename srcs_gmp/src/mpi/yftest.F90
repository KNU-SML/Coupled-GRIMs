#include <define.h>
   program yftest
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_,latg_,levs_,lonfp_,latgp_,levsr_,levsp_,        &
                         jcap_, para_init
   use commpi, only : commpi_init, mype, master
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real, allocatable    ::  a(:,:,:)
   real, allocatable    ::  bx(:,:,:)
   real, allocatable    ::  by(:,:,:)
   real, allocatable    ::  c(:,:,:)
   real*8 stim,etim
   real*4, allocatable  ::  dum(:,:)
   integer              ::  k
!
   call para_init
   call commpi_init
   call mpinit(stim)
   call mpdimset(jcap_,levs_,lonf_,latg_)
!
   allocate( a(lonf_,latg_,levs_)   )
   allocate( bx(lonfp_,latg_,levsr_)   )
   allocate( by(lonf_,latgp_,levsp_)   )
   allocate( c(lonf_,latg_,levs_)   )
   allocate( dum(lonf_,latg_)   )
!
   if (mype.eq.master) then
     read(11) a(:,:,1)
     a(:,:,10)=a(:,:,1)
     a(:,:,3)=a(:,:,1)
     do k = 1,levs_
       dum(:,:)=a(:,:,k)
       write(100) dum(:,:)
     enddo
   endif
!
! to Y partial
!
   call mpgf2yk(a,lonf_,latg_,levs_,1,by,latgp_,levsp_)
   if (mype.eq.master) then
     write(6,*)'mpgf2yk end'
     call flush(6)
   endif
   call mpyk2f(by,lonf_,latgp_,levsp_,1,c,latg_,levs_)
   if (mype.eq.master) then
     do k = 1,levs_
       dum(:,:)=c(:,:,k)
       write(100) dum(:,:)
     enddo
   endif
!
! to X partial
!
   call mpgf2xk(c,lonf_,latg_,levs_,1,bx,lonfp_,levsr_)
   if (mype.eq.master) then
     write(6,*)'mpgf2xk end'
     call flush(6)
   endif
   call mpxk2f(bx,lonfp_,latg_,levsr_,1,a,lonf_,levs_)
   if (mype.eq.master) then
     do k = 1,levs_
       dum(:,:)=a(:,:,k)
       write(100) dum(:,:)
     enddo
   endif
   call mpfine(etim)
!
   deallocate ( a, bx, by, c, dum )
!
   end program yftest
