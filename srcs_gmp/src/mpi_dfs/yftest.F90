#include "define.h"
   program yftest
!-------------------------------------------------------------------------------
!
! program: yftest
!
! program history log:
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!
!-------------------------------------------------------------------------------
   use paramodel
   use dfsvar, only : iba,jbwa,ib,jbw,get_dfs_dim,iope,                        &
                      levs,levsp,kler
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real, allocatable, dimension(:,:,:)  ::  a,bx,by,c
   real*4, allocatable, dimension(:,:)  ::  dum
   real*8                               ::  stim,etim
   integer k
!
   call para_init
   call mpinit(stim)
   call get_dfs_dim(jcap_,levs_,ngases_,nwater_,lonf_,latg_)
!
   allocate(a(iba,jbwa,levs),c(iba,jbwa,levs))
   allocate(bx(ib,jbwa,kler),by(iba,jbw,levsp))
   allocate(dum(iba,jbwa))

   if (iope) then
     read(11) a(:,:,1)
     a(:,:,10)=a(:,:,1)
     a(:,:,2)=a(:,:,1)
     do k = 1,levs
       dum(:,:)=a(:,:,k)
       write(100) dum(:,:)
     enddo
   endif
!
   call mpgf2yk(a,iba,jbwa,levs,1,by,jbw,levsp)
!
   if (iope) then
     write(6,*)'mpgf2yk end'
     call flush(6)
   endif
!
   call mpyk2f(by,iba,jbw,levsp,1,c,jbwa,levs)
!
   if (iope) then
     do k = 1,levs
       dum(:,:)=c(:,:,k)
       write(100) dum(:,:)
     enddo
   endif
!
   call mpgf2xk(c,iba,jbwa,levs,1,bx,ib,kler)
!
   if (iope) then
     write(6,*)'mpgf2xk end'
     call flush(6)
   endif
!
   call mpxk2f(bx,ib,jbwa,kler,1,a,iba,levs)
!
   if (iope) then
     do k = 1,levs
       dum(:,:)=a(:,:,k)
       write(100) dum(:,:)
     enddo
   endif
!
   call mpfine(etim)
!
   deallocate (a, bx, by, c, dum)
!
   end program yftest
