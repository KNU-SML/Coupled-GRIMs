#include <define.h>
   program mtn_driver
!-------------------------------------------------------------------------------
!
! main program:  mtn_driver  terrain maker for global spectral model
!
! abstract: this program creates terrain-related files
!   computed from the usgs 30 sec terrain dataset.
!   the model physics grid parameters and spectral truncation
!   are read by this program as input.
!   the files produced are respectively:
!     1) sea-land mask on model physics grid used by asfc program
!     2) gridded orography on model physics grid used by asfc program
!     3) mountain variance plus 10 statics (4oa,oc,4ol) for kagwd 
!     4) spectral orography in spectral domain used by chgr program
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          usgs topo, kagwd options 
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:
!
!   input files:
!     unit5      - physics longitudes (im), physics latitudes (jm),
!                  spectral truncation (nm) and rhomboidal flag (nr)
!                  respectively read in free format.
!     unit11     - packed usgs terrain dataset in lib/con
!                  see hong (2000, ncep office note, new orography data set)
!
!   output files:
!     unit51     - sea-land mask (im,jm)
!     unit52     - gridded orography (im,jm)
!     unit53     - mountain variance (im,jm)
!     unit54     - spectral orography ((nm+1)*((nr+1)*nm+2))
!     unit55     - unfiltered gridded orography (im,jm)
!
!   subprograms called:
!     unique:
!     mtn_solver               - main subprogram
!     mtn_trans_wave_grid      - spherical transform
!     mtn_gaussian_lat         - compute gaussian latitudes
!     mtn_sph_lat              - compute equally-spaced latitudes
!     mtn_legendre             - compute legendre polynomials
!     fftfax     - fft (library call can be substituted)
!     rfftmlt    - fft (library call can be substituted)
!     library:
!     gbytes     - unpack bits
!     qpass      - fft utility
!     rpass      - fft utility
!
!-------------------------------------------------------------------------------
   use paramodel, only : para_init
#ifndef SMP_RA2SFC
   use paramodel, only : im=>lonf_, jm=>latg_, nm=>jcap_
!
   integer, parameter  ::  nr=0
#endif
   integer, parameter  ::  nv=1
!
#ifdef NCO_TAG
   call w3tagb('clim_mtn',2001,0000,0000,'np51   ')
!
#endif
   call para_init
#ifdef SMP
#ifdef SMP_RA2SFC
   im = 192
   jm = 94
   nm = 62
   nr = 0
#endif
#endif
   print*, im,jm,nm,nr,nv
!
   if(nr.ne.0) then
     print *,'rohmboidal truncation not allowed'
     stop
   endif
!
   nw=(nm+1)*((nr+1)*nm+2)
   ! call mtn_solver(im,jm,nm,nr,nw,nv)
#ifdef SMP
#ifdef SMP_RA2SFC
   call mtn_solver(im,jm,nm,nr,nw)
#else
   call mtn_solver_scm(im,jm,nm,nr,nw)
#endif
#else
   call mtn_solver(im,jm,nm,nr,nw)
#endif
#ifdef NCO_TAG
!
   call w3tage('clim_mtn')
!
#endif
   stop
   end program mtn_driver
!-------------------------------------------------------------------------------
