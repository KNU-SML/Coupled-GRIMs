#include <define.h>
   module module_trans
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [module_trans] --- [dyn_trans2output_grid] *
!                   --- [dyn_trans2model_grid] *
!                   --- [rmp_trans2output_grid] *
!                   --- [rmp_trans2model_grid] *
!
!-------------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine dyn_trans2output_grid(a,km)
!-------------------------------------------------------------------------------
!
! abstract: separate northern and southern latitudes of gaussian grid.
!
! program history log:
!   1988-04-12  joseph sela            initial mrf
!   2000-01-01  hann-ming henry juang  mpi, reduced grid
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call dyn_trans2output_grid (a,lot)
!   input argument list:
!     a        - array of gaussian grid with paired
!                northern and southern latitudes.
!                input array is overwritten by output array.
!
!   output argument list:
!     a        - array of gaussian grid with separated
!                northern and southern latitudes.
!
!-------------------------------------------------------------------------------
   use paramodel, only : latg2_,lonf2_,lonf_,latg_
#ifdef REDUCE_GRID
   use comreduce, only : lonfd
   use module_sph_reduced_grid, only : sph_reduced_grid_interp
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                     ::  km
   real                        ::  a(lonf_,latg_,km)
!
! local variables
!
   real                       ::  b(lonf_,km)
   integer                    ::  jdone(latg_)
   integer                    ::  i,j,k,jsave,jput,jget
!
#ifdef REDUCE_GRID
   do k = 1 ,km
     call sph_reduced_grid_interp(a(1,1,k),lonfd,lonf2_,latg2_)
   enddo
#endif
   do k = 1,latg_
     jdone(k)=0
   enddo
   jsave=2
!
 3 continue
!
   if(jsave.le.latg2_) then
     jget=2*jsave-1
   else
     jget=2*(latg_+1-jsave)
   endif
   if(jget.eq.jsave) then
     jdone(jsave)=1
     go to 35
   endif
   do k = 1,km
     do i = 1,lonf_
       b(i,k)=a(i,jsave,k)
       a(i,jsave,k)=a(i,jget,k)
     enddo
   enddo
   jdone(jsave)=1
!
10 continue
!
   jput=jget
   if(jput.le.latg2_) then
     jget=2*jput-1
   else
     jget=2*(latg_+1-jput)
   endif
!
   if(jget.eq.jsave) go to 20
!
   do k = 1,km
     do i = 1,lonf_
       a(i,jput,k)=a(i,jget,k)
     enddo
   enddo
   jdone(jput)=1
!
   go to 10
!
20 continue
!
   do k = 1,km
     do i = 1,lonf_
       a(i,jput,k)=b(i,k)
     enddo
   enddo
   jdone(jput)=1
!
35 continue
!
   do j = jsave,latg_
     if(jdone(j).eq.0) then
       jsave=j
       go to 3
     endif
   enddo
!
   return
   end subroutine dyn_trans2output_grid
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine dyn_trans2model_grid(a,km)
!-------------------------------------------------------------------------------
!
! abstract: pair northern and southern latitudes of gaussian grid.
!
! usage:    call dyn_trans2model_grid (a,lot)
!   input argument list:
!     a        - array of gaussian grid with separated
!                northern and southern latitudes.
!                input array is overwritten by output array.
!
!   output argument list:
!     a        - array of gaussian grid with paired
!                northern and southern latitudes.
!
!-------------------------------------------------------------------------------
   use paramodel, only : latg_,lonf_
#ifdef REDUCE_GRID
   use comreduce
   use module_sph_reduced_grid, only : sph_reduced_grid_pick
#endif
!-------------------------------------------------------------------------------
   real                 ::  a(lonf_,latg_,km)
!
! local
!
   real                 ::  b(lonf_,km)
   integer              ::  jdone(latg_)
!
   do j = 1,latg_
     jdone(j)=0
   enddo
   jsave=2
!
 3 continue
!
   if(mod(jsave,2).eq.0) then
     jget=latg_+1-jsave/2
   else
     jget=(jsave+1)/2
   endif
   if(jget.eq.jsave) then
     jdone(jsave)=1
!
     go to 35
!
   endif
   do k = 1,km
     do i = 1,lonf_
       b(i,k)=a(i,jsave,k)
       a(i,jsave,k)=a(i,jget,k)
     enddo
   enddo
   jdone(jsave)=1
!
10 continue
!
   jput=jget
   if(mod(jput,2).eq.0) then
     jget=latg_+1-jput/2
   else
     jget=(jput+1)/2
   endif
!
   if(jget.eq.jsave) go to 20
!
   do k = 1,km
     do i = 1,lonf_
       a(i,jput,k)=a(i,jget,k)
     enddo
   enddo
   jdone(jput)=1
!
   go to 10
!
20 continue
!
   do k = 1,km
     do i = 1,lonf_
       a(i,jput,k)=b(i,k)
     enddo
   enddo
   jdone(jput)=1
!
35 continue
!
   do k = jsave,latg_
     if(jdone(k).eq.0) then
       jsave=k
!
       go to 3
!
     endif
   enddo
#ifdef REDUCE_GRID
!
   lonf2=lonf_*2
   latg2=latg_/2
   do k =1,km
     call sph_reduced_grid_pick(a(1,1,k),lonfd,lonf2,latg2)
   enddo
#endif
!
   return
   end subroutine dyn_trans2model_grid
!-------------------------------------------------------------------------------
#ifdef RMP
!
!-------------------------------------------------------------------------------
   subroutine rmp_trans2output_grid(a,km)
!-------------------------------------------------------------------------------
   use paramodel, only : igrd1_,igrd12_,jgrd1_,jgrd12_
!-------------------------------------------------------------------------------
!
!  abstract:  transform model grid (i,j) to regular grid (2i,j/2) with
!             north and south grid combined together as global model for
!             regional model computation.
!             the model grid is
!                (1,1),(2,1),.......(i,1),(1,j  ),(2,j  ),...,(i,j  )
!                (1,2),(2,2),......,(i,2),(1,j-1),(2,j-1),...,(i,j-1)
!                 ...................................................
!                (1,j/2),(2,j/2),..,(i,j/2),(1,j-j/2+1),...,(i,j-j/2+1)
!             the regular grid is
!                (1,1),(2,1),.......(i,1)
!                (1,2),(2,2),.......(i,2)
!                ........................
!                (1,j),(2,j),.......(i,j)
!
!  program history log:
!   1992-01-01  hann-ming henry juang  initial development
!
!  usage:    call rmp_trans2output_grid(a,k)
!    input argument list:
!      a 	- model grid with dimension (i,j,k)
!      k 	- the third dimension of a
!
!    output argument list:
!      a 	- regular output grid with dimension (i,j,k)
!
!    input files: none
!
!    output files: none
!
!-------------------------------------------------------------------------------
!
! arrage j row for output (regional regular grid)
!
   real  ::  a(igrd1_,jgrd1_,km)
   real  ::  tmp(igrd1_,jgrd1_)
!
   do k = 1,km
     do j = 1,jgrd12_
       jj = 2*j-1
       js = j
       jn = jgrd1_+1 - j
       do i = 1,igrd1_
         tmp(i,js) = a(i,jj  ,k)
         tmp(i,jn) = a(i,jj+1,k)
       enddo
     enddo
     do j = 1,jgrd1_
       do i = 1,igrd1_
         a(i,j,k) = tmp(i,j)
       enddo
     enddo
   enddo
!
   return
   end subroutine rmp_trans2output_grid
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rmp_trans2model_grid(a,km)
!-------------------------------------------------------------------------------
   use paramodel, only : igrd1_,igrd12_,jgrd1_,jgrd12_
!-------------------------------------------------------------------------------
!
!  abstract:  transform regulat grid (i,j) to model grid (2i,j/2) with
!             north and south grid combined together as global model for
!             regional model computation.
!             the regular grid is
!                (1,1),(2,1),.......(i,1)
!                (1,2),(2,2),.......(i,2)
!                ........................
!                (1,j),(2,j),.......(i,j)
!             the model grid is
!                (1,1),(2,1),.......(i,1),(1,j  ),(2,j  ),...,(i,j  )
!                (1,2),(2,2),......,(i,2),(1,j-1),(2,j-1),...,(i,j-1)
!                ...................................................
!                (1,j/2),(2,j/2),..,(i,j/2),(1,j-j/2+1),...,(i,j-j/2+1)
!
!  program history log:
!   1992-01-01  hann-ming henry juang  initial development
!
!  usage:    call rmp_trans2model_grid(a,k)
!    input argument list:
!      a 	- regular grid with dimension (i,j,k)
!      k 	- the third dimension of a
!
!    output argument list:
!      a 	- model grid with dimension (i,j,k)
!
!    input files: none
!
!    output files: none
!
!-------------------------------------------------------------------------------
!
!  arrage j row for model integration
!
   real  ::  a(igrd1_,jgrd1_,km)
   real  ::  tmp(igrd1_,jgrd1_)
!
   do k = 1,km
     do j = 1,jgrd12_
       jj = 2*j-1
       js = j
       jn = jgrd1_+1 - j
       do i = 1,igrd1_
            tmp(i,jj)   = a(i,js,k)
            tmp(i,jj+1) = a(i,jn,k)
       enddo
     enddo
     do j = 1,jgrd1_
       do i = 1,igrd1_
         a(i,j,k) = tmp(i,j)
       enddo
     enddo
   enddo
!
   return
   end subroutine rmp_trans2model_grid
!-------------------------------------------------------------------------------
#endif
   end module module_trans
