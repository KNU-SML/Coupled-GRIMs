#include <define.h>
   program post_sgb_main
!-------------------------------------------------------------------------------
! main program documentation block
!
! main program:  post_sgb_main     transform sigma to sigma grib
!
! abstract: program transforms sigma input to sigma grib output.
!   the output consists of data on a regular lat/lon grid.
!   geopotential height, wind components, relative humidity,
!   temperature and vertical velocity are output on model sigma layer.
!   also output are post_sgb_2d_diag fields consisting of
!   precipitable water, three lower level relative humidities,
!   lower level potential temperature and wind components,
!   surface temperature, pressure, omega and relative humidity,
!   and tropopause temperature, pressure, wind components and shear.
!   first namsgb namelist is read to determine output format.
!   then a sigma (grid or spectral) file is read from unit 11 and
!   the program produces and writes a sigma grib1 file to unit 51.
!   then a sigma file is read from unit 12 and
!   the program produces and writes a sigma grib1 file to unit 52.
!   the program continues until an empty input file is encountered.
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! namelists:
!   namsgb:      parameters determining output format
!     io         number of longitude points
!     jo         number of latitude points
!     ncpus      number of parallel processes (default: 1)
!     ids(255)   decimal scaling of packed data
!                (default: set by subprogram idsdef)
!     icen       forecast center identifier (default: 7)
!     icen2      forecast sub-center identifier (default: 0)
!     igen       model generating code (default: from sigma file)
!
! input files:
!   unit   11-?  sigma file(s)
!
! output files:
!   unit   51-?  sigma grib1 file(s)
!
! subprograms called:
!   idsdef                       set defaults for decimal scaling
!   funct_svp_init               compute saturated vapor pressure table
!   funct_dew_point_temp_init    compute dewpoint temerature table
!   funct_pot_temp_init          compute equivalent potential temperature table
!   funct_moist_adiabat_init     compute moist adiabat table
!   post_read_sigma              read a sigma file header
!   post_sgb_solver              transform one sigma file to sigma grib
!
!-------------------------------------------------------------------------------
   use paramodel, only : para_init,lonf_,latg_,levs_,jcap_
#ifdef DFS
   use dfsvar   , only : get_dfs_dim
   use paramodel, only : ngases_,nwater_
#endif
!-------------------------------------------------------------------------------
   integer, parameter  ::  levmax=100
   real                ::  si(levmax+1),sl(levmax)
   integer             ::  idate(4)
   integer             ::  ids(255)
   namelist/namsgb/ ids,icen,icen2,igen,fhs,fhe,fhinc,keid
!
   data ncpus/1/
   data keid/0/
   data ids/255*0/
   data icen/7/,icen2/1/,igen/0/
   data fhs/0./,fhe/48./,fhinc/6./
!-------------------------------------------------------------------------------
   call para_init
#ifdef DFS
   call get_dfs_dim(jcap_,levs_,ngases_,nwater_,lonf_,latg_)
#endif

#ifdef NCO_TAG
   call w3tagb('clim_sgb',2001,0000,0000,'np51   ')
#endif
!
!  set defaults and read namelist
!
   call idsdef(2,ids)
   read(*,namsgb,end=5)
5 continue
   call funct_svp_init
   call funct_dew_point_temp_init
   call funct_pot_temp_init
   call funct_moist_adiabat_init
   nsig=11
   nsgb=50
!
   fh=fhs
   do while(fh.le.fhe)
!
!  read sigma header record
!
     call post_read_sigma(nsig,fh,fhour,idate,si,sl,iret,igases,iwater)
     if(iret.ne.0) then
       print *,'read sigma file error'
       call abort
     endif
!
!  transform to sigma grib and attempt to read next sigma header
!
     iymdh=idate(4)*1000000+idate(2)*10000+idate(3)*100+idate(1)
     print *,' posting date ',iymdh,'+',nint(fhour),                           &
             '   sigma spectral t',jcap_,' l',levs_,                           &
             '   sigma grid ',lonf_,'x',latg_,'x',levs_
     call post_sgb_solver(fhour,idate,nsig,si,sl,                              &
               nsgb,ncpus,ids,icen,icen2,igen,igases,iwater,keid)
     close(nsig)
     close(nsgb)
     fh=fh+nint(fhinc)
   enddo
!
#ifdef NCO_TAG
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   call w3tage('clim_sgb')
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
   stop
   end program post_sgb_main
!-------------------------------------------------------------------------------
