










MODULE trcbc
   !!======================================================================
   !!                     ***  MODULE  trcbc  ***
   !! TOP :  module for passive tracer boundary conditions
   !!=====================================================================
   !! History :  3.5 !  2014 (M. Vichi, T. Lovato)  Original
   !!            3.6 !  2015 (T . Lovato) Revision and BDY support
   !!            4.0 !  2016 (T . Lovato) Include application of sbc and cbc
   !!----------------------------------------------------------------------
   !!   trc_bc       :  Apply tracer Boundary Conditions
   !!----------------------------------------------------------------------
   USE par_trc       !  passive tracers parameters
   USE oce_trc       !  shared variables between ocean and passive tracers
   USE trc           !  passive tracers common variables
   USE iom           !  I/O manager
   USE lib_mpp       !  MPP library
   USE fldread       !  read input fields
   USE bdy_oce,  ONLY: ln_bdy, nb_bdy , idx_bdy, ln_coords_file, rn_time_dmp, rn_time_dmp_out

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_bc         ! called in trcstp.F90 or within TOP modules
   PUBLIC   trc_bc_ini     ! called in trcini.F90 

   INTEGER  , SAVE, PUBLIC                             :: nb_trcobc    ! number of tracers with open BC
   INTEGER  , SAVE, PUBLIC                             :: nb_trcsbc    ! number of tracers with surface BC
   INTEGER  , SAVE, PUBLIC                             :: nb_trccbc    ! number of tracers with coastal BC
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_indobc ! index of tracer with OBC data
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_indsbc ! index of tracer with SBC data
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_indcbc ! index of tracer with CBC data
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trsfac    ! multiplicative factor for SBC tracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trcsbc    ! structure of data input SBC (file informations, fields read)
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trcfac    ! multiplicative factor for CBC tracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trccbc    ! structure of data input CBC (file informations, fields read)
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trofac    ! multiplicative factor for OBCtracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:), TARGET  :: sf_trcobc

   !!----------------------------------------------------------------------
   !!   Dummy module                              NO 3D passive tracer data
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_bc_ini( ntrc, Kmm )        ! Empty routine
      INTEGER, INTENT(IN) :: ntrc                           ! number of tracers
      INTEGER, INTENT(in) :: Kmm                            ! time level index
      WRITE(*,*) 'trc_bc_ini: You should not have seen this print! error?', ntrc, Kmm
   END SUBROUTINE trc_bc_ini
   SUBROUTINE trc_bc( kt, Kmm, Krhs )        ! Empty routine
      INTEGER, INTENT(in) :: kt, Kmm, Krhs ! time level indices
      WRITE(*,*) 'trc_bc: You should not have seen this print! error?', kt, Kmm, Krhs 
   END SUBROUTINE trc_bc

   !!======================================================================
END MODULE trcbc
