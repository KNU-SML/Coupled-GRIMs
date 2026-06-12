#include "define.h"
   program global_model_program
!-------------------------------------------------------------------------------
!
! program: global_model_program
!
! abstract: make global forecast 
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatble with namelist input
!
! references : 
!   hong et al. (2013, apjas): global/regional integrated model system (grims)
!   park et al. (2013, mwr): dfs dynamical core
!   byun and hong (2007, j. climate): single-column model (scm)
!   kanamitsu et al. (2002, bams): ncep dynamical seasonal forecast system 2000
!
! input files:
!   unit   11    sigma file (analysis or at time t-dt)
!   unit   12    sigma file (at time t if not analysis)
!   unit   14    surface file
!   unit   15    co2 constants (dependent on vertical resolution)
!   unit   24    mountain variance (dependent on horizontal resolution)
!   unit   43    cloud tuning
!
! output files:
!   unit   51    sigma file (at time t-dt)
!   unit   52    sigma file (at time t)
!   unit   53    surface file
!   unit   61    initial zonal diagnostics
!   unit   63    flux diagnostics
!   unit   64    final zonal diagnostics
!   unit   67    grid point diagnostics
!
!-------------------------------------------------------------------------------
#ifndef NIM
   use paramodel, only : LONF2S,LATG2S
   use comsfc          ! init_comsfc
#ifdef MP
   use commpi, only    : mype,master,mpi_comm_private
#endif
   use comfibm
   use comcon
   use comgpd
#ifdef LFM
   use comlfm
#endif
!soojin_couple
#ifdef AOMG
   use coupling
   use couple_oasis
#endif /* AOMG */
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#ifdef INTWALL
   include 'mpif.h'
   real*8   ::  stwall,etwall
#endif
   integer  :: info 
   real*8   ::  strwtime,endwtime
#ifndef MP
   integer  ::  clock_count,clock_rate,clock_max
!
   call system_clock(clock_count,clock_rate,clock_max)
   strwtime=real(clock_count,kind=8) / real(clock_rate,kind=8)
#endif
!
! initialize model configuration
!
   call gmp_init_array(strwtime)
!
#ifdef NCO_TAG
#ifdef MP
   if( mype.eq.0 ) then
#endif
     call w3tagb('clim_fcst',2001,0000,0000,'np51')
#ifdef MP
   endif
#endif
#endif /* NCO_TAG end */
!
! initialize surface variables
!
   call init_comsfc(LONF2S,LATG2S)
!
! initial call to global spectral model for gmp_start_setup etc.
!
   call gmp_start
!soojin_couple
#ifdef AOMG
! oasis initilize : define exchange variables
   call coupling_init
#endif /* AOMG */
!!!!!!!!!!!!!!!!!!!!!
!
! start time loop
!
#ifdef INTWALL
   stwall=mpi_wtime()
#endif
   do jdt = limlow,maxstp
!
#ifdef TIMEF
     tstep0=timef()*0.001
#endif
!
! call one time step for gmp related works
!    (can be forward-step or centered-step)
!
     call gmp_integrate
!
#ifdef TIMEF
     tstep=timef()*0.001-tstep0
#ifndef NOPRINT
     if(mype.eq.master) 
1    write(6,*)' time step, time per step: ',jdt,tstep
#endif
#endif
   enddo
#ifdef INTWALL
   etwall=mpi_wtime()
#endif
!
#ifdef MP
#ifdef DBG
   if(mype.eq.master) write(6,*) ' time loop ended '
#endif
#endif


!
#ifndef BARO_TEST
!  end  time  loop
!
   call gmp_end
!
#endif

!soojin_couple
#ifdef AOMG
!   write(*,*) "mype1:", info, mype
   call cpl_finalize
#endif
!
#ifdef MP
!   write(*,*) "mype2:", info, mype
   call mpfine(endwtime)
#else
#ifdef T3E
   endwtime=secondr()
#else
   call system_clock(clock_count,clock_rate,clock_max)
   endwtime=real(clock_count,kind=8) / real(clock_rate,kind=8)
#endif
#endif /* MP end */
!
! wall-clock time
!
#ifdef MP
   if (mype.eq.master) then
#endif
#ifdef INTWALL
   write(6,*) ' Integration Wall-clock-time is ',etwall-stwall
#endif
   write(6,*) ' Wall-clock-time is ',endwtime-strwtime
   call flush(6)
#ifdef MP
   endif
#endif
!
#ifdef LINUX_PGI
   call exit
#else
   stop
#endif


#endif /* only for GRIMs */
   end program global_model_program
