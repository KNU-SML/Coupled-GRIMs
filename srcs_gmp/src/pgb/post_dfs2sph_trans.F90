#include "define.h"
   subroutine post_dfs2sph_trans(sph,dfs,jcap,jl,km,ib,jbw,coslat,add)
!-------------------------------------------------------------------------------
!#define BIN_DBG
   use constant, only          : pi_
   use module_file_write, only : file_write_bin
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ... 
!
!    [post_dfs2sph_trans] * --- [post_swap_lat] * 
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                                   ::  jcap,jl,km,ib,jbw
   real, dimension(2*jcap+1,jl,km)           ::  dfs
   real, dimension((jcap+1)*(jcap+2),km)     ::  sph
   real, dimension(jbw,3)                    ::  coslat
   real, dimension(km)                       ::  add
!
! local variables
!
   integer                                   ::  k,idir
   real                                      ::  grid(ib*jbw,km)
!-------------------------------------------------------------------------------
   idir=101
   call dfs_fft_driver(-1,grid,ib,jbw,km,dfs,2*jcap+1,jl,km,                   &
                       jl,km,km,1,coslat,+1 )
   do k = 1,km
     sph(:,k)=0.
     if (add(k).ne.0.) then
       grid(:,k)=grid(:,k)+add(k)
     endif
   enddo
!
! change grid arrangement to gaussian grid
!
   call post_swap_lat(grid,ib,jbw,km)
   do k = 1,km
     call post_trans_wave_grid(idir,grid(1,k),sph(1,k),0,0.,ib,jbw,jcap,0)
   enddo      ! k
!
#ifdef BIN_DBG
   do k = 1,km
     call post_trans_wave_grid(-idir,grid(1,k),sph(1,k),0,0.,ib,jbw,jcap,0)
   enddo
   call post_swap_lat(grid,ib,jbw,km)
   call file_write_bin(221,grid,ib,jbw,km,0)
#endif
!
   return
   end subroutine post_dfs2sph_trans
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_swap_lat(a,im,jm,km)
!-------------------------------------------------------------------------------
   real      ::  a(im,jm,km)
   integer   ::  im,jm,km
!
   real      ::  dum(im)
   integer   ::  j,k
!-------------------------------------------------------------------------------
!
   do k = 1,km
     do j = 1,jm/2
       dum(:)=a(:,j,k)
       a(:,j,k)=a(:,jm-j+1,k)
       a(:,jm-j+1,k)=dum(:)
     enddo
   enddo
!   
   return
   end subroutine post_swap_lat
!-------------------------------------------------------------------------------
