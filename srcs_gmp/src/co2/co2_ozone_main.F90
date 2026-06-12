   program co2_ozone_main
!-------------------------------------------------------------------------------
!
! main program documentation block
!
! main program:  co2_ozone_main
!
! abstract:
!  program co2_ozone_main = consolidation of a number of dan schwarzkopf,gfdl
!                     codes to produce a file of co2 hgt data
!                     for any vertical coordinate (read by subroutine
!                     rad_co2_read in the gfdl radiation codes)-k.a.c. jun89.
!nov89--updated (nov 89) for latest gfdl lw radiation.....k.a.c.
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
! namelists:
!
! input files:
!   unit
!
! output files:
!   unit   66
!
! subprograms called:
!     co2o3_step1
!     co2o3_step2
!     co2o3_step3
!     co2o3_step4
!
!-------------------------------------------------------------------------------
   use paramodel, only : para_init, l=>levs_,lp1=>levp1_,lp2=>levp2_
   use comco2, only    : comco2_init
!-------------------------------------------------------------------------------
   real,allocatable,target  ::  sgtemp(:,:),co2d1d(:,:),co2d2d(:,:,:)
   real,allocatable,target  ::  co2iq2(:,:,:),co2iq3(:,:,:),co2iq5(:,:,:)
   real,allocatable,target  ::  t41(:,:),t42(:),t43(:,:),t44(:)
   real,allocatable,target  ::  t20(:,:,:),t21(:,:,:)
   real,allocatable,target  ::  t22(:,:,:),t23(:,:,:)
   real,allocatable,target  ::  sglvnu(:),siglnu(:)
!
   real,pointer      ::  stemp(:),gtemp(:)
   real,pointer      ::  cdtm51(:),co2m51(:),c2dm51(:)
   real,pointer      ::  cdtm58(:),co2m58(:),c2dm58(:)
   real,pointer      ::  cdt51(:,:),co251(:,:),c2d51(:,:)
   real,pointer      ::  cdt58(:,:),co258(:,:),c2d58(:,:)
   real,pointer      ::  cdt31(:),co231(:),c2d31(:)
   real,pointer      ::  cdt38(:),co238(:),c2d38(:)
   real,pointer      ::  cdt71(:),co271(:),c2d71(:)
   real,pointer      ::  cdt78(:),co278(:),c2d78(:)
   real,pointer      ::  co211(:),co218(:)
!
   call para_init
   call comco2_init
#ifdef NCO_TAG
   call w3tagb('clim_co2',2001,0000,0000,'np51')
#endif
!
   allocate(sgtemp(lp1,2),co2d1d(l,6),co2d2d(lp1,lp1,6) )
   allocate(co2iq2(lp1,lp1,6),co2iq3(lp1,lp1,6),co2iq5(lp1,lp1,6) )
   allocate(t41(lp2,2),t42(lp1),t43(lp2,2),t44(lp1) )
   allocate(t20(lp1,lp1,3),t21(lp1,lp1,3) )
   allocate(t22(lp1,lp1,3),t23(lp1,lp1,3) )
   allocate(sglvnu(lp1),siglnu(l) )
   cdt31=>co2iq2(:,1,1); co231=>co2iq2(:,1,2); c2d31=>co2iq2(:,1,3)
   cdt38=>co2iq2(:,1,4); co238=>co2iq2(:,1,5); c2d38=>co2iq2(:,1,6)
   cdt71=>co2iq3(:,1,1); co271=>co2iq3(:,1,2); c2d71=>co2iq3(:,1,3)
   cdt78=>co2iq3(:,1,4); co278=>co2iq3(:,1,5); c2d78=>co2iq3(:,1,6)
   co211=>co2iq5(:,1,2); co218=>co2iq5(:,1,5)
   stemp=>sgtemp(:,1); gtemp=>sgtemp(:,2)
   cdtm51=>co2d1d(:,1); co2m51=>co2d1d(:,2); c2dm51=>co2d1d(:,3)
   cdtm58=>co2d1d(:,4); co2m58=>co2d1d(:,5); c2dm58=>co2d1d(:,6)
   cdt51=>co2d2d(:,:,1); co251=>co2d2d(:,:,2); c2d51=>co2d2d(:,:,3)
   cdt58=>co2d2d(:,:,4); co258=>co2d2d(:,:,5); c2d58=>co2d2d(:,:,6)
!
!===>  get sgtemp and output which used to be on units 41,42,43,44....
!
   lread = 0
   call co2o3_step1(sgtemp,t41,t42,t43,t44,sglvnu,siglnu,lread)
!
!===>  interpolate desired co2 data from the detailed(109,109) grid..
!         ir=1,iq=1 is for common /co2bd3/ in radiation code...
!           for the consolidated 490-850 cm-1 band...
!
   ico2tp=61
   ir = 1
!
!===>  ratio is the actual co2 mixing ratio in units of 330 ppmv.....
!        e.g. sinv in co2_tau_interp has been obtained using rco2=330 ppmv
!
   ratio = 1.0
   read(5,*,end=599) co2ppm
   ratio = co2ppm / 330.
   599   continue
   nmethd = 2
   call co2o3_step2(ico2tp,t41,t42,t22,ratio,ir,nmethd)
   ir = 1
   nmethd = 1
   call co2o3_step2(ico2tp,t41,t42,t20,ratio,ir,nmethd)
   ir = 1
   nmethd = 2
   call co2o3_step2(ico2tp,t43,t44,t23,ratio,ir,nmethd)
   ir = 1
   nmethd = 1
   call co2o3_step2(ico2tp,t43,t44,t21,ratio,ir,nmethd)
!
!===>    fill up the co2d1d array
!       the following gets co2 transmission functions and
!         their derivatives for tau(i,i+1),i=1,levs,
!         where the values are not obtained by quadrature but are the
!         actual transmissivities,etc,between a pair of pressures. these
!         are used only for nearby layer calculations including h2o..
!
!
   iq = 1
   call co2o3_step3(t20,t21,co2d1d,iq)
!
!
!===>    fill up the co2d2d array
!    the following gets co2 transmission functions and their derivatives
!        from 109-level line-by-line calculations made using the 1982
!        mcclatchy tape (12511 lines),consolidated,interpolated
!        to the mrf vertical coordinate,and re-consolidated to a
!        200 cm-1 bandwidth. the interpolation method is described in
!        schwarzkopf and fels (j.g.r.,1985).
!
   call co2o3_step4(t22,t23,co2d2d,iq)
!
!nov89
!===>  interpolate desired co2 data from the detailed(109,109) grid..
!         ir=2,iq=2 is for common /co2bd2/ in radiation code...
!           for the consolidated 490-670 cm-1 band...
!
   ico2tp=62
   ir = 2
   nmethd = 2
   call co2o3_step2(ico2tp,t41,t42,t22,ratio,ir,nmethd)
   call co2o3_step2(ico2tp,t43,t44,t23,ratio,ir,nmethd)
   iq = 2
   call co2o3_step4(t22,t23,co2iq2,iq)
!
!===>  interpolate desired co2 data from the detailed(109,109) grid..
!         ir=3,iq=3 is for common /co2bd4/ in radiation code...
!           for the consolidated 670-850 cm-1 band...
!
   ico2tp=63
   ir = 3
   nmethd = 2
   call co2o3_step2(ico2tp,t41,t42,t22,ratio,ir,nmethd)
   call co2o3_step2(ico2tp,t43,t44,t23,ratio,ir,nmethd)
   iq = 3
   call co2o3_step4(t22,t23,co2iq3,iq)
!
!---      following code not working and not needed yet
!===>  interpolate desired co2 data from the detailed(109,109) grid..
!         ir=4,iq=5 is for common /co2bd5/ in radiation code...
!           for the 4.3 micron band....this for rco2=300 ppmv...
! not used yet      ico2tp=65
! not used yet      ir = 4
! not used yet      ratio = 0.909091
! not used yet      nmethd = 2
! not used yet      call co2o3_step2(ico2tp,t41,t42,t22,ratio,ir,nmethd)
! not used yet      call co2o3_step2(ico2tp,t43,t44,t23,ratio,ir,nmethd)
! not used yet      iq = 5
! not used yet      call co2o3_step4(t22,t23,co2iq5,iq)
!nov89
!...     write data to disk..
!            ...since these codes are compiled with autodbl,the co2 data
!               is converted to single precision in a later job step..
!
   rewind 66
   write(66) stemp
   write(66) gtemp
   write(66) cdtm51
   write(66) co2m51
   write(66) c2dm51
   write(66) cdtm58
   write(66) co2m58
   write(66) c2dm58
   write(66) cdt51
   write(66) co251
   write(66) c2d51
   write(66) cdt58
   write(66) co258
   write(66) c2d58
   write(66) cdt31
   write(66) co231
   write(66) c2d31
   write(66) cdt38
   write(66) co238
   write(66) c2d38
   write(66) cdt71
   write(66) co271
   write(66) c2d71
   write(66) cdt78
   write(66) co278
   write(66) c2d78
   write(66) co2ppm*1.e-6
!
! not used yet      write(66) co211
! not used yet      write(66) co218
!
   rewind 66
#ifdef NCO_TAG
   call w3tage('clim_co2')
#endif
   deallocate(sgtemp,co2d1d,co2d2d )
   deallocate(co2iq2,co2iq3,co2iq5 )
   deallocate(t41,t42,t43,t44 )
   deallocate(t20,t21 )
   deallocate(t22,t23 )
   deallocate(sglvnu,siglnu )
!
   stop
   end program co2_ozone_main
!-------------------------------------------------------------------------------
