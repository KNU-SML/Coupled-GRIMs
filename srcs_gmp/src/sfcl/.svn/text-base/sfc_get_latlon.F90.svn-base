#include <define.h>
   subroutine sfc_get_latlon(idim,jdim,rlat,rlon)
!-------------------------------------------------------------------------------
   use constant, only     :  pi_
#ifdef RMP
   use module_trans, only :  rmp_trans2output_grid
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer               ::  idim,jdim
   real                  ::  rlat(idim,jdim),rlon(idim,jdim)
   real                  ::  rdelx,rdely,dlamda0,dx
   real                  ::  gaul(jdim)
!
   integer i,j
!
!  compute latitude and longitude
!
#ifdef RMP
   call rmp_setup_rsm_grid(rlat,rlon,rdelx,rdely,dlamda0)
   call rmp_trans2output_grid(rlat,1)
   call rmp_trans2output_grid(rlon,1)
!
   do j = 1,jdim
     do i = 1,idim
       rlat(i,j) = rlat(i,j) * 180. / pi_
       rlon(i,j) = rlon(i,j) * 180. / pi_
     enddo
   enddo
#else
!
!  compute gaussian latitude and longitude
!
   dx=360./float(idim)
   call gaulat(gaul,jdim)
   do j = 1,jdim
     do i = 1,idim
       rlat(i,j)=90.-gaul(j)
       rlon(i,j)=float(i-1)*dx
       if(rlon(i,j).gt.180.) then
         rlon(i,j)=rlon(i,j)-360.
       endif
     enddo
   enddo
#endif
!
   return
   end subroutine sfc_get_latlon
!-------------------------------------------------------------------------------
