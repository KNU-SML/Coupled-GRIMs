#include <define.h>
   module module_sph_reduced_grid
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [module_sph_reduced_grid]
!      |
!      |--- [sph_reduced_grid_init] *
!      |--- [sph_reduced_grid_pick] *
!      |--- [sph_reduced_grid_interp] *
!
!-------------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_reduced_grid_init(qtt,lnt2,jcap,qttcut,lcapd,lonfd)
!-------------------------------------------------------------------------------
!
! purpose: global reduce grid initial routine to compute lcapd and lonfd
!          the maximal value of qtt and the accuray of the digit (ndigit)
!          are used to determine the resolution lcapd and grid point
!          lonfd for each latitude. The lonfd is determined by the
!          checks of factors only 2, 3 and 5 with at least one 2, and
!          it is larger than lcapd*3+1.
!
!   2000-01-01  hann-ming henry juang  initial development
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! input:
!   qtt   coefficient for spectral transform at given lat
!   lnt2   dimension of qtt at given lat
!   jcap   wave dimension
!   qttcut   minimax value of accuracy for qtt
! output:
!   lcapd   wave resolution for reduced grid
!   lonfd   number of reduced grid point for latitude
!
!-------------------------------------------------------------------------------
   real                 ::   qtt(lnt2)
!
   ind=0
   mwave=0
   do m = 0,jcap
     need=0
     do n = m,jcap
       ind=ind+1
       ireal=2*ind-1
       imagi=2*ind
       if(abs(qtt(ireal)).ge.qttcut.and.                                       &
          abs(qtt(imagi)).ge.qttcut) then
         need=1
       endif
     enddo
     mwave=mwave+need
   enddo
   lcapd=mwave
!
   lonfi=3*(lcapd-1)+1
   lonff=lonfi+mod(lonfi,2)
   do ii = 0,100,2
     lonf=lonff+ii
     lonfo=lonf
     jtime=nint(log(float(lonf))/log(2.))
#ifdef DCRFT
     jtime=min(jtime,25)
#endif
     do j = 1,jtime
       if( mod(lonf,2).eq.0 ) then
         lonf=lonf/2
         if( lonf.eq.1 ) go to 200
       endif
     enddo
     ktime=nint(log(float(lonf))/log(3.))
#ifdef DCRFT
     ktime=min(ktime,2)
#endif
     do k = 1,ktime
       if( mod(lonf,3).eq.0 ) then
         lonf=lonf/3
         if( lonf.eq.1 ) go to 200
       endif
     enddo
     ltime=nint(log(float(lonf))/log(5.))
#ifdef DCRFT
     ltime=min(ltime,1)
#endif
     do l = 1,ltime
       if( mod(lonf,5).eq.0 ) then
         lonf=lonf/5
         if( lonf.eq.1 ) go to 200
       endif
     enddo
#ifdef DCRFT
     if( mod(lonf,7).eq.0 ) then
       lonf=lonf/7
       if( lonf.eq.1 ) go to 200
     endif
     if( mod(lonf,11).eq.0 ) then
       lonf=lonf/11
       if( lonf.eq.1 ) go to 200
     endif
#endif
   enddo
200 lonfd=lonfo
!
   return
   end subroutine sph_reduced_grid_init
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_reduced_grid_pick(a,lonfd,lonf2,latg2)
!-------------------------------------------------------------------------------
!
! subprogram: sph_reduced_grid_pick
!
! abstract: pick the reduce-grid value from the nearest regular grid
!   then fill the tailing points to be the same as the last
!   reduce-grid point.
!
! program history log:
!   2000        hann-ming henry juang 
!
!------------------------------------------------------------------------------
   real                 ::  a(lonf2,latg2)
   integer              ::  lonfd(latg2)
   real,allocatable     ::  tmp(:)
!
   allocate(tmp(lonf2))
!
   lonf=lonf2/2
   dg=360./float(lonf)
!
   do j = 1,latg2
     dr=360./float(lonfd(j))
     do i = 1,lonfd(j)
       ii=nint((i-1.)*dr/dg + 1.0)
       tmp(i)=a(ii,j)
       tmp(i+lonfd(j))=a(ii+lonf,j)
     enddo
     ii=lonfd(j)*2
     do i = 1,ii
       a(i,j)=tmp(i)
     enddo
     do i = ii+1,lonf2
       a(i,j)=tmp(ii)
     enddo
   enddo
   deallocate(tmp)
!
   return
   end subroutine sph_reduced_grid_pick
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_reduced_grid_interp(a,lonfd,lonf2,latg2)
!-------------------------------------------------------------------------------
!
! subprogram: sph_reduced_grid_interp
!
! abstract: return the reduce-grid value to the nearest regular grid
!   then, the restored value is propogated eastward to the regular
!   grid which has un-specified value.
!
! program history log:
!   2000        hann-ming henry juang 
!
!------------------------------------------------------------------------------
   real                 ::  a(lonf2,latg2)
   integer              ::  lonfd(latg2)
   real,allocatable     ::  tmp(:)
!
   allocate(tmp(lonf2))
!
   lonf=lonf2/2
   dg=360./float(lonf)
   do j = 1,latg2
     do i = 1,lonf2
       tmp(i)=-999.999
     enddo
     dr=360./float(lonfd(j))
     do i = 1,lonfd(j)
       ii=nint((i-1.)*dr/dg + 1.0)
       tmp(ii)=a(i,j)
       tmp(ii+lonf)=a(i+lonfd(j),j)
     enddo
     do i = 1,lonf
       if( tmp(i).ne.-999.999 ) rep=tmp(i)
       tmp(i)=rep
     enddo
     do i = 1,lonf
       if( tmp(i+lonf).ne.-999.999 ) rep=tmp(i+lonf)
       tmp(i+lonf)=rep
     enddo
     do i = 1,lonf2
       a(i,j)=tmp(i)
     enddo
   enddo
   deallocate(tmp)
!
   return
   end subroutine sph_reduced_grid_interp
!-------------------------------------------------------------------------------
   end module module_sph_reduced_grid
