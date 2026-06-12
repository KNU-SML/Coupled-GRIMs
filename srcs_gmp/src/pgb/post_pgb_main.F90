#include <define.h>
   program post_pgb_main
!-------------------------------------------------------------------------------
! main program documentation block
!
! main program:  post_pgb_main     transform sigma to pressure grib
!
! abstract: program transforms sigma input to pressure grib output.
!   the output consists of data on a regular lat/lon grid.
!   geopotential height, wind components, relative humidity,
!   temperature and vertical velocity are output on mandatory pressures.
!   also output are post_pgb_2d_diag fields consisting of
!   precipitable water, three lower level relative humidities,
!   lower level potential temperature and wind components,
!   surface temperature, pressure, omega and relative humidity,
!   and tropopause temperature, pressure, wind components and shear.
!   first nampgb namelist is read to determine output format.
!   then a sigma (grid or spectral) file is read from unit 11 and
!   the program produces and writes a pressure grib1 file to unit 51.
!   then a sigma file is read from unit 12 and
!   the program produces and writes a pressure grib1 file to unit 52.
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
!   nampgb:      parameters determining output format
!     io_        number of longitude points
!     jo_        number of latitude points
!     ko_        number of pressure levels
!     po(ko_)    pressures in mb (default: mandatory levels)
!     ncpus      number of parallel processes (default: 1)
!     mxbit      maximum number of bits to pack data (default: 16)
!     ids(255)   decimal scaling of packed data
!                (default: set by subprogram idsdef)
!     pot(255)   highest pressure in mb to output data
!                as a function of parameter indicator
!                (default: 300 for rh, 100 for omega, 0 otherwise)
!     icen       forecast center identifier (default: 7)
!     icen2      forecast sub-center identifier (default: 0)
!     igen       model generating code (default: from sigma file)
!
! input files:
!   unit   11-?  sigma file(s)
!
! output files:
!   unit   51-?  pressure grib1 file(s)
!
! subprograms called:
!   idsdef                       set defaults for decimal scaling
!   funct_svp_init               compute saturated vapor pressure table
!   funct_dew_point_temp_init    compute dewpoint temerature table
!   funct_pot_temp_init          compute equivalent potential temperature table
!   funct_moist_adiabat_init     compute moist adiabat table
!   post_read_sigma              read a sigma file header
!   post_pgb_solver              transform one sigma file to pressure grib
!
!-------------------------------------------------------------------------------
   use paramodel, only : para_init                                            ,&
                         ngases_,nwater_,ntotal_,io_,jo_,ko_,levs_,jcap_      ,&
                         lonf_,latg_,lonf2_,latg2_,siodd_
#ifdef DFS
   use dfsvar   , only : get_dfs_dim
   use paramodel, only : lonf_,latg_
#else
   use paramodel, only : iod_,iodd_,johf_,twoj1_
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer, parameter            ::  levmax=100                               ,&
                                     komax=100                                ,&
                                     mxbit=16
   integer                       ::  icen,icen2,igen,keid,ncpus,nsig,npgb     ,&
                                     iret,igases,iwater,iymdh,mgrid           ,&
                                     io,jo,io2,io22,johf,nctop
   real                          ::  fh,fhs,fhe,fhinc,fhour
   integer, dimension(4)         ::  idate
   integer, dimension(255)       ::  ids
   real   , dimension(levmax)    ::  sl
   real   , dimension(levmax+1)  ::  si
   real   , dimension(komax)     ::  po
   real   , dimension(255)       ::  pot
   namelist/nampgb/ ids,pot,icen,icen2,igen,fhs,fhe,fhinc,keid,mgrid
   namelist /postpres/ po
!
   data ncpus/1/
   data keid/0/
   data mgrid/0/
   data ids/255*0/,pot/255*0./
   data icen/7/,icen2/1/,igen/0/
   data fhs/0./,fhe/120./,fhinc/24./
!-------------------------------------------------------------------------------
   call para_init
#ifdef DFS
   call get_dfs_dim(jcap_,levs_,ngases_,nwater_,lonf_,latg_)
#endif
#ifdef NCO_TAG
   call w3tagb('clim_pgb',2001,0000,0000,'np51   ')
#endif
!
   po = 0.
   read(1,postpres)
!
!  set defaults and read namelist
!
   call idsdef(2,ids)
!
! pdot
!
   pot(33)=po(ko_)
   pot(34)=po(ko_)
   pot( 7)=po(ko_)
   pot(11)=po(ko_)
   pot(52)=100.
   pot(39)=100.
   pot(152)=100.
   pot(153)=100.
   pot(154)=po(ko_)
   read(*,nampgb,end=5)
5 continue
   call funct_svp_init
   call funct_dew_point_temp_init
   call funct_pot_temp_init
   call funct_moist_adiabat_init
   nsig=11
   npgb=50
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
! grid information
!
#ifdef DFS
     if (mgrid.eq.1) then
       io=lonf_
       jo=latg_
     else
       io=io_
       jo=jo_
     endif
     io2=io
     io22=io
     johf=jo
     nctop=0  ! not used
#else
     if(mgrid.eq.1) then
       io=lonf_
       jo=latg_
       io2=lonf2_
       io22=siodd_
       johf=latg2_
     else
       io=io_
       jo=jo_
       io2=iod_
       io22=iodd_
       johf=johf_
     endif
     nctop=twoj1_
#endif
!
!  transform to pressure grib and attempt to read next sigma header
!
     iymdh=idate(4)*1000000+idate(2)*10000+idate(3)*100+idate(1)
     print *,' posting date ',iymdh,'+',nint(fhour),                           &
             '   sigma spectral t',jcap_,' l',levs_,                           &
             '   pressure grid ',io,'x',jo,'x',ko_
!
     call post_pgb_solver(fhour,idate,nsig,si,sl,po,mgrid,                     &
            npgb,ncpus,ids,pot,icen,icen2,igen,igases,iwater,keid,             &
            ngases_,nwater_,ntotal_,jcap_,ko_,levs_,io,jo,io2,io22,johf,nctop)
     close(nsig)
     close(npgb)
     fh=fh+nint(fhinc)
   enddo
!
#ifdef NCO_TAG
   call w3tage('clim_pgb')
!
#endif
   stop
   end program post_pgb_main
