










MODULE iom
   !!======================================================================
   !!                    ***  MODULE  iom ***
   !! Input/Output manager :  Library to read input files
   !!======================================================================
   !! History :  2.0  ! 2005-12  (J. Belier) Original code
   !!            2.0  ! 2006-02  (S. Masson) Adaptation to NEMO
   !!            3.0  ! 2007-07  (D. Storkey) Changes to iom_gettime
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie and G. Reffray)  add C1D case
   !!            3.6  ! 2014-15  DIMG format removed
   !!            3.6  ! 2015-15  (J. Harle) Added procedure to read REAL attributes
   !!            4.0  ! 2017-11  (M. Andrejczuk) Extend IOM interface to write any 3D fields
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   iom_open       : open a file read only
   !!   iom_close      : close a file or all files opened by iom
   !!   iom_get        : read a field (interfaced to several routines)
   !!   iom_varid      : get the id of a variable in a file
   !!   iom_rstput     : write a field in a restart file (interfaced to several routines)
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE domutl          !
   USE flo_oce         ! floats module declarations
   USE lbclnk          ! lateal boundary condition / mpp exchanges
   USE iom_def         ! iom variables definitions
   USE iom_nf90        ! NetCDF format with native NetCDF library
   USE in_out_manager  ! I/O manager
   USE lib_mpp           ! MPP library
   USE sbc_oce  , ONLY :   nn_fsbc, ght_abl, ghw_abl, e3t_abl, e3w_abl, jpka, jpkam1
   USE icb_oce  , ONLY :   nclasses, class_num       !  !: iceberg classes
   USE ice      , ONLY :   jpl
   USE phycst          ! physical constants
   USE dianam          ! build name of file
   USE ioipsl, ONLY :  ju2ymds    ! for calendar
   USE crs             ! Grid coarsening
   USE lib_fortran
   USE diu_bulk, ONLY : ln_diurnal_only, ln_diurnal
   USE iom_nf90
   USE netcdf

   IMPLICIT NONE
   PUBLIC   !   must be public to be able to access iom_def through iom

   LOGICAL, PUBLIC, PARAMETER ::   lk_iomput = .FALSE.       !: iom_put flag
   PUBLIC iom_init, iom_init_closedef, iom_swap, iom_open, iom_close, iom_setkt, iom_varid, iom_get, iom_get_var
   PUBLIC iom_chkatt, iom_getatt, iom_putatt, iom_getszuld, iom_rstput, iom_delay_rst, iom_put
   PUBLIC iom_use, iom_context_finalize, iom_update_file_name, iom_miss_val
   PUBLIC iom_xios_setid

   PRIVATE iom_rp0d_sp, iom_rp1d_sp, iom_rp2d_sp, iom_rp3d_sp
   PRIVATE iom_rp0d_dp, iom_rp1d_dp, iom_rp2d_dp, iom_rp3d_dp
   PRIVATE iom_get_123d
   PRIVATE iom_g0d_sp, iom_g1d_sp, iom_g2d_sp, iom_g3d_sp
   PRIVATE iom_g0d_dp, iom_g1d_dp, iom_g2d_dp, iom_g3d_dp
   PRIVATE iom_p1d_sp, iom_p2d_sp, iom_p3d_sp, iom_p4d_sp
   PRIVATE iom_p1d_dp, iom_p2d_dp, iom_p3d_dp, iom_p4d_dp
   PRIVATE set_xios_context
   PRIVATE iom_set_rstw_active

   INTERFACE iom_get
      MODULE PROCEDURE iom_g0d_sp, iom_g1d_sp, iom_g2d_sp, iom_g3d_sp
      MODULE PROCEDURE iom_g0d_dp, iom_g1d_dp, iom_g2d_dp, iom_g3d_dp
   END INTERFACE
   INTERFACE iom_getatt
      MODULE PROCEDURE iom_g0d_iatt, iom_g1d_iatt, iom_g0d_ratt, iom_g1d_ratt, iom_g0d_catt
   END INTERFACE
   INTERFACE iom_putatt
      MODULE PROCEDURE iom_p0d_iatt, iom_p1d_iatt, iom_p0d_ratt, iom_p1d_ratt, iom_p0d_catt
   END INTERFACE
   INTERFACE iom_rstput
      MODULE PROCEDURE iom_rp0d_sp, iom_rp1d_sp, iom_rp2d_sp, iom_rp3d_sp
      MODULE PROCEDURE iom_rp0d_dp, iom_rp1d_dp, iom_rp2d_dp, iom_rp3d_dp
   END INTERFACE
   INTERFACE iom_put
      MODULE PROCEDURE iom_p0d_sp, iom_p1d_sp, iom_p2d_sp, iom_p3d_sp, iom_p4d_sp
      MODULE PROCEDURE iom_p0d_dp, iom_p1d_dp, iom_p2d_dp, iom_p3d_dp, iom_p4d_dp
   END INTERFACE iom_put

   !! * Substitutions




   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: iom.F90 15033 2021-06-21 10:24:45Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE iom_init( cdname, kdid, ld_closedef )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE   ***
      !!
      !! ** Purpose :
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=*),           INTENT(in)  :: cdname
      INTEGER         , OPTIONAL, INTENT(in)  :: kdid
      LOGICAL         , OPTIONAL, INTENT(in)  :: ld_closedef
      !
   END SUBROUTINE iom_init

   SUBROUTINE iom_init_closedef(cdname)
      !!----------------------------------------------------------------------
      !!            ***  SUBROUTINE iom_init_closedef  ***
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose : Closure of context definition
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=*), OPTIONAL, INTENT(IN) :: cdname
      IF( .FALSE. )   WRITE(numout,*) 'iom_init_closedef: should not see this'   ! useless statement to avoid compilation warnings

   END SUBROUTINE iom_init_closedef

   SUBROUTINE iom_set_vars_active(idnum)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_set_vars_active  ***
      !!
      !! ** Purpose :  define filename in XIOS context for reading file,
      !!               enable variables present in a file for reading with XIOS
      !!               id of the file is assumed to be rrestart.
      !!---------------------------------------------------------------------
      INTEGER, INTENT(IN) :: idnum

   END SUBROUTINE iom_set_vars_active

   SUBROUTINE iom_set_rstw_file(cdrst_file)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE iom_set_rstw_file   ***
      !!
      !! ** Purpose :  define file name in XIOS context for writing restart
      !!---------------------------------------------------------------------
      CHARACTER(len=*) :: cdrst_file
   END SUBROUTINE iom_set_rstw_file


   SUBROUTINE iom_set_rstw_active(sdfield, rd0, rs0, rd1, rs1, rd2, rs2, rd3, rs3)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE iom_set_rstw_active   ***
      !!
      !! ** Purpose :  define file name in XIOS context for writing restart
      !!               enable variables present in restart file for writing
      !!---------------------------------------------------------------------
!sets enabled = .TRUE. for each field in restart file
      CHARACTER(len = *), INTENT(IN)                     :: sdfield
      REAL(dp), OPTIONAL, INTENT(IN)                     :: rd0
      REAL(sp), OPTIONAL, INTENT(IN)                     :: rs0
      REAL(dp), OPTIONAL, INTENT(IN), DIMENSION(:)       :: rd1
      REAL(sp), OPTIONAL, INTENT(IN), DIMENSION(:)       :: rs1
      REAL(dp), OPTIONAL, INTENT(IN), DIMENSION(:, :)    :: rd2
      REAL(sp), OPTIONAL, INTENT(IN), DIMENSION(:, :)    :: rs2
      REAL(dp), OPTIONAL, INTENT(IN), DIMENSION(:, :, :) :: rd3
      REAL(sp), OPTIONAL, INTENT(IN), DIMENSION(:, :, :) :: rs3
   END SUBROUTINE iom_set_rstw_active

   FUNCTION iom_axis(idlev) result(axis_ref)
      !!---------------------------------------------------------------------
      !!                   ***  FUNCTION  iom_axis  ***
      !!
      !! ** Purpose : Used for grid definition when XIOS is used to read/write
      !!              restart. Returns axis corresponding to the number of levels
      !!              given as an input variable. Axes are defined in routine
      !!              iom_set_rst_context
      !!---------------------------------------------------------------------
      INTEGER, INTENT(IN) :: idlev
      CHARACTER(len=lc)   :: axis_ref
      CHARACTER(len=12)   :: str
      IF(idlev == jpk) THEN
         axis_ref="nav_lev"
      ELSEIF(idlev == jpka) THEN
         axis_ref="nav_hgt"
      ELSEIF(idlev == jpl) THEN
         axis_ref="numcat"
      ELSE
         write(str, *) idlev
         CALL ctl_stop( 'iom_axis', 'Definition for axis with '//TRIM(ADJUSTL(str))//' levels missing')
      ENDIF
   END FUNCTION iom_axis

   FUNCTION iom_xios_setid(cdname) result(kid)
     !!---------------------------------------------------------------------
      !!                   ***  FUNCTION    ***
      !!
      !! ** Purpose : this function returns first available id to keep information about file
      !!              sets filename in iom_file structure and sets name
      !!              of XIOS context depending on cdcomp
      !!              corresponds to iom_nf90_open
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   ) :: cdname      ! File name
      INTEGER                         :: kid      ! identifier of the opened file
      INTEGER                         :: jl

      kid = 0
      DO jl = jpmax_files, 1, -1
         IF( iom_file(jl)%nfid == 0 )   kid = jl
      ENDDO

      iom_file(kid)%name   = TRIM(cdname)
      iom_file(kid)%nfid   = 1
      iom_file(kid)%nvars  = 0
      iom_file(kid)%irec   = -1

   END FUNCTION iom_xios_setid

   SUBROUTINE iom_set_rst_context(ld_rstr)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_set_rst_context  ***
      !!
      !! ** Purpose : Define domain, axis and grid for restart (read/write)
      !!              context
      !!
      !!---------------------------------------------------------------------
      LOGICAL, INTENT(IN)               :: ld_rstr
      INTEGER :: ji
   END SUBROUTINE iom_set_rst_context


   SUBROUTINE set_xios_context(kdid, cdcont)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_set_rst_context  ***
      !!
      !! ** Purpose : set correct XIOS context based on kdid
      !!
      !!---------------------------------------------------------------------
      INTEGER,           INTENT(IN)     :: kdid           ! Identifier of the file
      CHARACTER(LEN=lc), INTENT(OUT)    :: cdcont         ! name of the context for XIOS read/write

      cdcont = "NONE"

      IF(lrxios) THEN
         IF(kdid == numror) THEN
            cdcont = cr_ocerst_cxt
         ELSEIF(kdid == numrir) THEN
            cdcont = cr_icerst_cxt
         ELSEIF(kdid == numrar) THEN
            cdcont = cr_ablrst_cxt
         ELSEIF(kdid == numrtr) THEN
            cdcont = cr_toprst_cxt
         ELSEIF(kdid == numrsr) THEN
            cdcont = cr_sedrst_cxt
         ENDIF
      ENDIF

      IF(lwxios) THEN
         IF(kdid == numrow) THEN
            cdcont = cw_ocerst_cxt
         ELSEIF(kdid == numriw) THEN
            cdcont = cw_icerst_cxt
         ELSEIF(kdid == numraw) THEN
            cdcont = cw_ablrst_cxt
         ELSEIF(kdid == numrtw) THEN
            cdcont = cw_toprst_cxt
         ELSEIF(kdid == numrsw) THEN
            cdcont = cw_sedrst_cxt
         ENDIF
      ENDIF
   END SUBROUTINE set_xios_context


   SUBROUTINE iom_swap( cdname )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_swap  ***
      !!
      !! ** Purpose :  swap context between different agrif grid for xmlio_server
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdname
      CHARACTER(len=256)           :: clname, cltmpn
      INTEGER                      :: iln
      !
   END SUBROUTINE iom_swap


   SUBROUTINE iom_open( cdname, kiomid, ldwrt, ldstop, ldiof, kdlev, cdcomp )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_open  ***
      !!
      !! ** Purpose :  open an input file (return 0 if not found)
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   )           ::   cdname   ! File name
      INTEGER         , INTENT(  out)           ::   kiomid   ! iom identifier of the opened file
      LOGICAL         , INTENT(in   ), OPTIONAL ::   ldwrt    ! open in write modeb          (default = .FALSE.)
      LOGICAL         , INTENT(in   ), OPTIONAL ::   ldstop   ! stop if open to read a non-existing file (default = .TRUE.)
      LOGICAL         , INTENT(in   ), OPTIONAL ::   ldiof    ! Interp On the Fly, needed for AGRIF (default = .FALSE.)
      INTEGER         , INTENT(in   ), OPTIONAL ::   kdlev    ! number of vertical levels
      CHARACTER(len=3), INTENT(in   ), OPTIONAL ::   cdcomp   ! name of component calling iom_nf90_open
      !
      CHARACTER(LEN=256)    ::   clname    ! the name of the file based on cdname [[+clcpu]+clcpu]
      CHARACTER(LEN=256)    ::   cltmpn    ! tempory name to store clname (in writting mode)
      CHARACTER(LEN=10)     ::   clsuffix  ! ".nc"
      CHARACTER(LEN=15)     ::   clcpu     ! the cpu number (max jpmax_digits digits)
      CHARACTER(LEN=256)    ::   clinfo    ! info character
      LOGICAL               ::   llok      ! check the existence
      LOGICAL               ::   llwrt     ! local definition of ldwrt
      LOGICAL               ::   llstop    ! local definition of ldstop
      LOGICAL               ::   lliof     ! local definition of ldiof
      INTEGER               ::   icnt      ! counter for digits in clcpu (max = jpmax_digits)
      INTEGER               ::   iln, ils  ! lengths of character
      INTEGER               ::   istop     !
      ! local number of points for x,y dimensions
      ! position of first local point for x,y dimensions
      ! position of last local point for x,y dimensions
      ! start halo size for x,y dimensions
      ! end halo size for x,y dimensions
      !---------------------------------------------------------------------
      ! Initializations and control
      ! =============
      kiomid = -1
      clinfo = '                    iom_open ~~~  '
      istop = nstop
      ! if iom_open is called for the first time: initialize iom_file(:)%nfid to 0
      ! (could be done when defining iom_file in f95 but not in f90)
      IF( Agrif_Root() ) THEN
         IF( iom_open_init == 0 ) THEN
            iom_file(:)%nfid = 0
            iom_open_init = 1
         ENDIF
      ENDIF
      ! do we read or write the file?
      IF( PRESENT(ldwrt) ) THEN   ;   llwrt = ldwrt
      ELSE                        ;   llwrt = .FALSE.
      ENDIF
      ! do we call ctl_stop if we try to open a non-existing file in read mode?
      IF( PRESENT(ldstop) ) THEN   ;   llstop = ldstop
      ELSE                         ;   llstop = .TRUE.
      ENDIF
      ! are we using interpolation on the fly?
      IF( PRESENT(ldiof) ) THEN   ;   lliof = ldiof
      ELSE                        ;   lliof = .FALSE.
      ENDIF
      ! create the file name by added, if needed, TRIM(Agrif_CFixed()) and TRIM(clsuffix)
      ! =============
      clname   = TRIM(cdname)
      IF ( .NOT. Agrif_Root() .AND. .NOT. lliof ) THEN
         iln    = INDEX(clname,'/', BACK=.TRUE.)
         cltmpn = clname(1:iln)
         clname = clname(iln+1:LEN_TRIM(clname))
         clname = TRIM(cltmpn)//TRIM(Agrif_CFixed())//'_'//TRIM(clname)
      ENDIF
      ! which suffix should we use?
      clsuffix = '.nc'
      ! Add the suffix if needed
      iln = LEN_TRIM(clname)
      ils = LEN_TRIM(clsuffix)
      IF( iln <= ils .OR. INDEX( TRIM(clname), TRIM(clsuffix), back = .TRUE. ) /= iln - ils + 1 )   &
         &   clname = TRIM(clname)//TRIM(clsuffix)
      cltmpn = clname   ! store this name
      ! try to find if the file to be opened already exist
      ! =============
      INQUIRE( FILE = clname, EXIST = llok )
      IF( .NOT.llok ) THEN
         ! we try to add the cpu number to the name
         WRITE(clcpu,*) narea-1

         clcpu  = TRIM(ADJUSTL(clcpu))
         iln = INDEX(clname,TRIM(clsuffix), back = .TRUE.)
         clname = clname(1:iln-1)//'_'//TRIM(clcpu)//TRIM(clsuffix)
         icnt = 0
         INQUIRE( FILE = clname, EXIST = llok )
         ! we try different formats for the cpu number by adding 0
         DO WHILE( .NOT.llok .AND. icnt < jpmax_digits )
            clcpu  = "0"//TRIM(clcpu)
            clname = clname(1:iln-1)//'_'//TRIM(clcpu)//TRIM(clsuffix)
            INQUIRE( FILE = clname, EXIST = llok )
            icnt = icnt + 1
         END DO
      ELSE
         lxios_sini = .TRUE.
      ENDIF
      ! Open the NetCDF file
      ! =============
      ! do we have some free file identifier?
      IF( MINVAL(iom_file(:)%nfid) /= 0 )   &
         &   CALL ctl_stop( TRIM(clinfo), 'No more free file identifier', 'increase jpmax_files in iom_def' )
      ! if no file was found...
      IF( .NOT. llok ) THEN
         IF( .NOT. llwrt ) THEN   ! we are in read mode
            IF( llstop ) THEN   ;   CALL ctl_stop( TRIM(clinfo), 'File '//TRIM(cltmpn)//'* not found' )
            ELSE                ;   istop = nstop + 1   ! make sure that istop /= nstop so we don't open the file
            ENDIF
         ELSE                     ! we are in write mode so we
            clname = cltmpn       ! get back the file name without the cpu number
         ENDIF
      ELSE
         IF( llwrt .AND. .NOT. ln_clobber ) THEN   ! we stop as we want to write in a new file
            CALL ctl_stop( TRIM(clinfo), 'We want to write in a new file but '//TRIM(clname)//' already exists...' )
            istop = nstop + 1                      ! make sure that istop /= nstop so we don't open the file
         ELSEIF( llwrt ) THEN     ! the file exists and we are in write mode with permission to
            clname = cltmpn       ! overwrite so get back the file name without the cpu number
         ENDIF
      ENDIF
      IF( istop == nstop ) THEN   ! no error within this routine
         CALL iom_nf90_open( clname, kiomid, llwrt, llok, kdlev = kdlev, cdcomp = cdcomp )
      ENDIF
      !
   END SUBROUTINE iom_open


   SUBROUTINE iom_close( kiomid )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_close  ***
      !!
      !! ** Purpose : close an input file, or all files opened by iom
      !!--------------------------------------------------------------------
      INTEGER, INTENT(inout), OPTIONAL ::   kiomid   ! iom identifier of the file to be closed
      !                                              ! return 0 when file is properly closed
      !                                              ! No argument: all files opened by iom are closed

      INTEGER ::   jf         ! dummy loop indices
      INTEGER ::   i_s, i_e   ! temporary integer
      CHARACTER(LEN=100)    ::   clinfo    ! info character
      !---------------------------------------------------------------------
      !
      IF( iom_open_init == 0 )   RETURN   ! avoid to use iom_file(jf)%nfid that us not yet initialized
      !
      clinfo = '                    iom_close ~~~  '
      IF( PRESENT(kiomid) ) THEN
         i_s = kiomid
         i_e = kiomid
      ELSE
         i_s = 1
         i_e = jpmax_files
      ENDIF

      IF( i_s > 0 ) THEN
         DO jf = i_s, i_e
            IF( iom_file(jf)%nfid > 0 ) THEN
               CALL iom_nf90_close( jf )
               iom_file(jf)%nfid       = 0          ! free the id
               IF( PRESENT(kiomid) )   kiomid = 0   ! return 0 as id to specify that the file was closed
               IF(lwp) WRITE(numout,*) TRIM(clinfo)//' close file: '//TRIM(iom_file(jf)%name)//' ok'
            ELSEIF( PRESENT(kiomid) ) THEN
               WRITE(ctmp1,*) '--->',  kiomid
               CALL ctl_stop( TRIM(clinfo)//' Invalid file identifier', ctmp1 )
            ENDIF
         END DO
      ENDIF
      !
   END SUBROUTINE iom_close


   FUNCTION iom_varid ( kiomid, cdvar, kdimsz, kndims, lduld, ldstop )
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  iom_varid  ***
      !!
      !! ** Purpose : get the id of a variable in a file (return 0 if not found)
      !!-----------------------------------------------------------------------
      INTEGER              , INTENT(in   )           ::   kiomid   ! file Identifier
      CHARACTER(len=*)     , INTENT(in   )           ::   cdvar    ! name of the variable
      INTEGER, DIMENSION(:), INTENT(  out), OPTIONAL ::   kdimsz   ! size of each dimension
      INTEGER              , INTENT(  out), OPTIONAL ::   kndims   ! number of dimensions
      LOGICAL              , INTENT(  out), OPTIONAL ::   lduld    ! true if the last dimension is unlimited (time)
      LOGICAL              , INTENT(in   ), OPTIONAL ::   ldstop   ! stop if looking for non-existing variable (default = .TRUE.)
      !
      INTEGER                        ::   iom_varid, iiv, i_nvd
      LOGICAL                        ::   ll_fnd
      CHARACTER(LEN=100)             ::   clinfo                   ! info character
      LOGICAL                        ::   llstop                   ! local definition of ldstop
      !!-----------------------------------------------------------------------
      iom_varid = 0                         ! default definition
      ! do we call ctl_stop if we look for non-existing variable?
      IF( PRESENT(ldstop) ) THEN   ;   llstop = ldstop
      ELSE                         ;   llstop = .TRUE.
      ENDIF
      !
      IF( kiomid > 0 ) THEN
         clinfo = 'iom_varid, file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(cdvar)
         IF( iom_file(kiomid)%nfid == 0 ) THEN
            CALL ctl_stop( TRIM(clinfo), 'the file is not open' )
         ELSE
            ll_fnd  = .FALSE.
            iiv = 0
            !
            DO WHILE ( .NOT.ll_fnd .AND. iiv < iom_file(kiomid)%nvars )
               iiv = iiv + 1
               ll_fnd  = ( TRIM(cdvar) == TRIM(iom_file(kiomid)%cn_var(iiv)) )
            END DO
            !
            IF( .NOT.ll_fnd ) THEN
               iiv = iiv + 1
               IF( iiv <= jpmax_vars ) THEN
                  iom_varid = iom_nf90_varid( kiomid, cdvar, iiv, kdimsz, kndims, lduld )
               ELSE
                  CALL ctl_stop( TRIM(clinfo), 'Too many variables in the file '//iom_file(kiomid)%name,   &
                        &                      'increase the parameter jpmax_vars')
               ENDIF
               IF( llstop .AND. iom_varid == -1 )   CALL ctl_stop( TRIM(clinfo)//' not found' )
            ELSE
               iom_varid = iiv
               IF( PRESENT(kdimsz) ) THEN
                  i_nvd = iom_file(kiomid)%ndims(iiv)
                  IF( i_nvd <= size(kdimsz) ) THEN
                     kdimsz(1:i_nvd) = iom_file(kiomid)%dimsz(1:i_nvd,iiv)
                  ELSE
                     WRITE(ctmp1,*) i_nvd, size(kdimsz)
                     CALL ctl_stop( TRIM(clinfo), 'error in kdimsz size'//TRIM(ctmp1) )
                  ENDIF
               ENDIF
               IF( PRESENT(kndims) )  kndims = iom_file(kiomid)%ndims(iiv)
               IF( PRESENT( lduld) )  lduld  = iom_file(kiomid)%luld( iiv)
            ENDIF
         ENDIF
      ENDIF
      !
   END FUNCTION iom_varid


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_get
   !!----------------------------------------------------------------------
   SUBROUTINE iom_g0d_sp( kiomid, cdvar, pvar, ktime )
      INTEGER         , INTENT(in   )                 ::   kiomid    ! Identifier of the file
      CHARACTER(len=*), INTENT(in   )                 ::   cdvar     ! Name of the variable
      REAL(sp)        , INTENT(  out)                 ::   pvar      ! read field
      REAL(dp)                                        ::   ztmp_pvar ! tmp var to read field
      INTEGER         , INTENT(in   ),     OPTIONAL   ::   ktime     ! record number
      !
      INTEGER                                         ::   idvar     ! variable id
      INTEGER                                         ::   idmspc    ! number of spatial dimensions
      INTEGER         , DIMENSION(1)                  ::   itime     ! record number
      CHARACTER(LEN=100)                              ::   clinfo    ! info character
      CHARACTER(LEN=100)                              ::   clname    ! file name
      CHARACTER(LEN=1)                                ::   cldmspc   !
      CHARACTER(LEN=lc)                               ::   context
      !
      CALL set_xios_context(kiomid, context)

      IF(context == "NONE") THEN  ! read data using default library
         itime = 1
         IF( PRESENT(ktime) ) itime = ktime
         !
         clname = iom_file(kiomid)%name
         clinfo = '          iom_g0d, file: '//TRIM(clname)//', var: '//TRIM(cdvar)
         !
         IF( kiomid > 0 ) THEN
            idvar = iom_varid( kiomid, cdvar )
            IF( iom_file(kiomid)%nfid > 0 .AND. idvar > 0 ) THEN
               idmspc = iom_file ( kiomid )%ndims( idvar )
               IF( iom_file(kiomid)%luld(idvar) )  idmspc = idmspc - 1
               WRITE(cldmspc , fmt='(i1)') idmspc
               IF( idmspc > 0 )  CALL ctl_stop( TRIM(clinfo), 'When reading to a 0D array, we do not accept data', &
                                    &                         'with 1 or more spatial dimensions: '//cldmspc//' were found.' , &
                                    &                         'Use ncwa -a to suppress the unnecessary dimensions' )
               CALL iom_nf90_get( kiomid, idvar, ztmp_pvar, itime )
               pvar = ztmp_pvar
            ENDIF
         ENDIF
      ELSE
         WRITE(ctmp1,*) 'Can not use XIOS in iom_g0d, file: '//TRIM(clname)//', var:'//TRIM(cdvar)
         CALL ctl_stop( 'iom_g0d', ctmp1 )
      ENDIF
   END SUBROUTINE iom_g0d_sp

   SUBROUTINE iom_g0d_dp( kiomid, cdvar, pvar, ktime )
      INTEGER         , INTENT(in   )                 ::   kiomid    ! Identifier of the file
      CHARACTER(len=*), INTENT(in   )                 ::   cdvar     ! Name of the variable
      REAL(dp)        , INTENT(  out)                 ::   pvar      ! read field
      INTEGER         , INTENT(in   ),     OPTIONAL   ::   ktime     ! record number
      !
      INTEGER                                         ::   idvar     ! variable id
      INTEGER                                         ::   idmspc    ! number of spatial dimensions
      INTEGER         , DIMENSION(1)                  ::   itime     ! record number
      CHARACTER(LEN=100)                              ::   clinfo    ! info character
      CHARACTER(LEN=100)                              ::   clname    ! file name
      CHARACTER(LEN=1)                                ::   cldmspc   !
      CHARACTER(LEN=lc)                               ::   context
      !
      CALL set_xios_context(kiomid, context)

      IF(context == "NONE") THEN  ! read data using default library
         itime = 1
         IF( PRESENT(ktime) ) itime = ktime
         !
         clname = iom_file(kiomid)%name
         clinfo = '          iom_g0d, file: '//TRIM(clname)//', var: '//TRIM(cdvar)
         !
         IF( kiomid > 0 ) THEN
            idvar = iom_varid( kiomid, cdvar )
            IF( iom_file(kiomid)%nfid > 0 .AND. idvar > 0 ) THEN
               idmspc = iom_file ( kiomid )%ndims( idvar )
               IF( iom_file(kiomid)%luld(idvar) )  idmspc = idmspc - 1
               WRITE(cldmspc , fmt='(i1)') idmspc
               IF( idmspc > 0 )  CALL ctl_stop( TRIM(clinfo), 'When reading to a 0D array, we do not accept data', &
                                    &                         'with 1 or more spatial dimensions: '//cldmspc//' were found.' , &
                                    &                         'Use ncwa -a to suppress the unnecessary dimensions' )
               CALL iom_nf90_get( kiomid, idvar, pvar, itime )
            ENDIF
         ENDIF
      ELSE
         WRITE(ctmp1,*) 'Can not use XIOS in iom_g0d, file: '//TRIM(clname)//', var:'//TRIM(cdvar)
         CALL ctl_stop( 'iom_g0d', ctmp1 )
      ENDIF
   END SUBROUTINE iom_g0d_dp

   SUBROUTINE iom_g1d_sp( kiomid, kdom, cdvar, pvar, ktime, kstart, kcount )
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(sp)        , INTENT(  out), DIMENSION(:)           ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      INTEGER         , INTENT(in   ), DIMENSION(1), OPTIONAL ::   kstart    ! start axis position of the reading
      INTEGER         , INTENT(in   ), DIMENSION(1), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom, cdvar, pvsp1d = pvar,   &
               &                                            ktime = ktime, kstart = kstart, kcount = kcount )
      ENDIF
   END SUBROUTINE iom_g1d_sp


   SUBROUTINE iom_g1d_dp( kiomid, kdom, cdvar, pvar, ktime, kstart, kcount )
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(dp)        , INTENT(  out), DIMENSION(:)           ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      INTEGER         , INTENT(in   ), DIMENSION(1), OPTIONAL ::   kstart    ! start axis position of the reading
      INTEGER         , INTENT(in   ), DIMENSION(1), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom, cdvar, pvdp1d = pvar,   &
            &                                               ktime = ktime, kstart = kstart, kcount = kcount)
      ENDIF
   END SUBROUTINE iom_g1d_dp

   SUBROUTINE iom_g2d_sp( kiomid, kdom, cdvar, pvar, ktime, cd_type, psgn, kfill, kstart, kcount)
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(sp)        , INTENT(  out), DIMENSION(:,:)         ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      CHARACTER(len=1), INTENT(in   )              , OPTIONAL ::   cd_type   ! nature of grid-points (T, U, V, F, W)
      REAL(sp)        , INTENT(in   )              , OPTIONAL ::   psgn      ! -1.(1.): (not) change sign across the north fold
      INTEGER         , INTENT(in   )              , OPTIONAL ::   kfill     ! value of kfillmode in lbc_lbk
      INTEGER         , INTENT(in   ), DIMENSION(2), OPTIONAL ::   kstart    ! start axis position of the reading
      INTEGER         , INTENT(in   ), DIMENSION(2), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom, cdvar, pvsp2d = pvar,                 &
            &                                               cd_type = cd_type, psgn_sp = psgn, kfill = kfill,   &
            &                                               ktime = ktime, kstart = kstart, kcount = kcount )
      ENDIF
   END SUBROUTINE iom_g2d_sp

   SUBROUTINE iom_g2d_dp( kiomid, kdom, cdvar, pvar, ktime, cd_type, psgn, kfill, kstart, kcount)
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(dp)        , INTENT(  out), DIMENSION(:,:)         ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      CHARACTER(len=1), INTENT(in   )              , OPTIONAL ::   cd_type   ! nature of grid-points (T, U, V, F, W)
      REAL(dp)        , INTENT(in   )              , OPTIONAL ::   psgn      ! -1.(1.): (not) change sign across the north fold
      INTEGER         , INTENT(in   )              , OPTIONAL ::   kfill     ! value of kfillmode in lbc_lbk
      INTEGER         , INTENT(in   ), DIMENSION(2), OPTIONAL ::   kstart    ! start axis position of the reading
      INTEGER         , INTENT(in   ), DIMENSION(2), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom, cdvar, pvdp2d = pvar,                 &
            &                                               cd_type = cd_type, psgn_dp = psgn, kfill = kfill,   &
            &                                               ktime = ktime, kstart = kstart, kcount = kcount )
      ENDIF
   END SUBROUTINE iom_g2d_dp

   SUBROUTINE iom_g3d_sp( kiomid, kdom, cdvar, pvar, ktime, cd_type, psgn, kfill, kstart, kcount )
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(sp)        , INTENT(  out), DIMENSION(:,:,:)       ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      CHARACTER(len=1), INTENT(in   )              , OPTIONAL ::   cd_type   ! nature of grid-points (T, U, V, F, W)
      REAL(sp)        , INTENT(in   )              , OPTIONAL ::   psgn      ! -1.(1.) : (not) change sign across the north fold
      INTEGER         , INTENT(in   )              , OPTIONAL ::   kfill     ! value of kfillmode in lbc_lbk
      INTEGER         , INTENT(in   ), DIMENSION(3), OPTIONAL ::   kstart    ! start axis position of the reading
      INTEGER         , INTENT(in   ), DIMENSION(3), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom, cdvar, pvsp3d = pvar,                 &
            &                                               cd_type = cd_type, psgn_sp = psgn, kfill = kfill,   &
            &                                               ktime = ktime, kstart = kstart, kcount = kcount )
      ENDIF
   END SUBROUTINE iom_g3d_sp

   SUBROUTINE iom_g3d_dp( kiomid, kdom, cdvar, pvar, ktime, cd_type, psgn, kfill, kstart, kcount )
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(dp)        , INTENT(  out), DIMENSION(:,:,:)       ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      CHARACTER(len=1), INTENT(in   )              , OPTIONAL ::   cd_type   ! nature of grid-points (T, U, V, F, W)
      REAL(dp)        , INTENT(in   )              , OPTIONAL ::   psgn      ! -1.(1.) : (not) change sign across the north fold
      INTEGER         , INTENT(in   )              , OPTIONAL ::   kfill     ! value of kfillmode in lbc_lbk
      INTEGER         , INTENT(in   ), DIMENSION(3), OPTIONAL ::   kstart    ! start axis position of the reading
      INTEGER         , INTENT(in   ), DIMENSION(3), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom, cdvar, pvdp3d = pvar,                 &
            &                                               cd_type = cd_type, psgn_dp = psgn, kfill = kfill,   &
            &                                               ktime = ktime, kstart = kstart, kcount = kcount )
      ENDIF
   END SUBROUTINE iom_g3d_dp

   !!----------------------------------------------------------------------

   SUBROUTINE iom_get_123d( kiomid , kdom, cdvar, pvsp1d, pvsp2d, pvsp3d, pvdp1d, pvdp2d, pvdp3d,   &
         &                  ktime , cd_type, psgn_sp, psgn_dp, kfill, kstart, kcount )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_get_123d  ***
      !!
      !! ** Purpose : read a 1D/2D/3D variable
      !!
      !! ** Method : read ONE record at each CALL
      !!-----------------------------------------------------------------------
      INTEGER                    , INTENT(in   )           ::   kiomid    ! Identifier of the file
      INTEGER                    , INTENT(in   )           ::   kdom      ! Type of domain to be read
      CHARACTER(len=*)           , INTENT(in   )           ::   cdvar     ! Name of the variable
      REAL(sp), DIMENSION(:)     , INTENT(  out), OPTIONAL ::   pvsp1d    ! read field (1D case), single precision
      REAL(sp), DIMENSION(:,:)   , INTENT(  out), OPTIONAL ::   pvsp2d    ! read field (2D case), single precision
      REAL(sp), DIMENSION(:,:,:) , INTENT(  out), OPTIONAL ::   pvsp3d    ! read field (3D case), single precision
      REAL(dp), DIMENSION(:)     , INTENT(  out), OPTIONAL ::   pvdp1d    ! read field (1D case), double precision
      REAL(dp), DIMENSION(:,:)   , INTENT(  out), OPTIONAL ::   pvdp2d    ! read field (2D case), double precision
      REAL(dp), DIMENSION(:,:,:) , INTENT(  out), OPTIONAL ::   pvdp3d    ! read field (3D case), double precision
      INTEGER                    , INTENT(in   ), OPTIONAL ::   ktime     ! record number
      CHARACTER(len=1)           , INTENT(in   ), OPTIONAL ::   cd_type   ! nature of grid-points (T, U, V, F, W)
      REAL(sp)                   , INTENT(in   ), OPTIONAL ::   psgn_sp   ! -1.(1.) : (not) change sign across the north fold
      REAL(dp)                   , INTENT(in   ), OPTIONAL ::   psgn_dp   ! -1.(1.) : (not) change sign across the north fold
      INTEGER                    , INTENT(in   ), OPTIONAL ::   kfill     ! value of kfillmode in lbc_lbk
      INTEGER , DIMENSION(:)     , INTENT(in   ), OPTIONAL ::   kstart    ! start position of the reading in each axis
      INTEGER , DIMENSION(:)     , INTENT(in   ), OPTIONAL ::   kcount    ! number of points to be read in each axis
      !
      LOGICAL                        ::   llok        ! true if ok!
      INTEGER                        ::   jl          ! loop on number of dimension
      INTEGER                        ::   idom        ! type of domain
      INTEGER                        ::   idvar       ! id of the variable
      INTEGER                        ::   inbdim      ! number of dimensions of the variable
      INTEGER                        ::   idmspc      ! number of spatial dimensions
      INTEGER                        ::   itime       ! record number
      INTEGER                        ::   istop       ! temporary value of nstop
      INTEGER                        ::   ix1, ix2, iy1, iy2   ! subdomain indexes
      INTEGER                        ::   ji, jj      ! loop counters
      INTEGER                        ::   irankpv     !
      INTEGER                        ::   ind1, ind2  ! substring index
      INTEGER, DIMENSION(jpmax_dims) ::   istart      ! starting point to read for each axis
      INTEGER, DIMENSION(jpmax_dims) ::   icnt        ! number of value to read along each axis
      INTEGER, DIMENSION(jpmax_dims) ::   idimsz      ! size of the dimensions of the variable
      INTEGER, DIMENSION(jpmax_dims) ::   ishape      ! size of the dimensions of the variable
      REAL(sp)                       ::   zscf_sp, zofs_sp  ! sacle_factor and add_offset, single precision
      REAL(dp)                       ::   zscf_dp, zofs_dp  ! sacle_factor and add_offset, double precision
      REAL(sp)                       ::   zsgn_sp     ! local value of psgn, single precision
      REAL(dp)                       ::   zsgn_dp     ! local value of psgn, double precision
      INTEGER                        ::   itmp        ! temporary integer
      CHARACTER(LEN=256)             ::   clinfo      ! info character
      CHARACTER(LEN=256)             ::   clname      ! file name
      CHARACTER(LEN=1)               ::   clrankpv, cldmspc      !
      CHARACTER(LEN=1)               ::   cl_type     ! local value of cd_type
      LOGICAL                        ::   ll_only3rd  ! T => if kstart, kcount present then *only* use values for 3rd spatial dimension.
      LOGICAL                        ::   llis1d, llis2d, llis3d
      LOGICAL                        ::   llsp        ! use single precision
      INTEGER                        ::   inlev       ! number of levels for 3D data
      !---------------------------------------------------------------------
      CHARACTER(LEN=lc)                               ::   context
      !
      CALL set_xios_context(kiomid, context)
      !
      llsp = PRESENT(pvsp1d) .OR. PRESENT(pvsp2d) .OR. PRESENT(pvsp3d)
      IF( llsp ) THEN 
         llis1d = PRESENT(pvsp1d)   ;   IF( llis1d )   ishape(1:1) = SHAPE(pvsp1d)
         llis2d = PRESENT(pvsp2d)   ;   IF( llis2d )   ishape(1:2) = SHAPE(pvsp2d)
         llis3d = PRESENT(pvsp3d)   ;   IF( llis3d )   ishape(1:3) = SHAPE(pvsp3d)
      ELSE
         llis1d = PRESENT(pvdp1d)   ;   IF( llis1d )   ishape(1:1) = SHAPE(pvdp1d)
         llis2d = PRESENT(pvdp2d)   ;   IF( llis2d )   ishape(1:2) = SHAPE(pvdp2d)
         llis3d = PRESENT(pvdp3d)   ;   IF( llis3d )   ishape(1:3) = SHAPE(pvdp3d)
      ENDIF
      inlev = -1
      IF( llis3d )   inlev = ishape(3)
      !
      idom = kdom
      istop = nstop
      !
      IF(context == "NONE") THEN
         clname = iom_file(kiomid)%name   !   esier to read
         clinfo = '          iom_get_123d, file: '//TRIM(clname)//', var: '//TRIM(cdvar)
         ! check kcount and kstart optionals parameters...
         IF( PRESENT(kcount) .AND. .NOT. PRESENT(kstart) ) CALL ctl_stop(TRIM(clinfo), 'kcount present needs kstart present')
         IF( PRESENT(kstart) .AND. .NOT. PRESENT(kcount) ) CALL ctl_stop(TRIM(clinfo), 'kstart present needs kcount present')
         IF( PRESENT(kstart) .AND. idom /= jpdom_unknown .AND. idom /= jpdom_auto_xy ) &
            &          CALL ctl_stop(TRIM(clinfo), 'kstart present needs idom = jpdom_unknown or idom = jpdom_auto_xy')
         IF( idom == jpdom_auto_xy .AND. .NOT. PRESENT(kstart) ) &
            &          CALL ctl_stop(TRIM(clinfo), 'idom = jpdom_auto_xy requires kstart to be present')
         !
         ! Search for the variable in the data base (eventually actualize data)
         !
         idvar = iom_varid( kiomid, cdvar )
         IF( idvar > 0 ) THEN
            !
            idimsz(:) = iom_file(kiomid)%dimsz(:, idvar)      ! to write iom_file(kiomid)%dimsz in a shorter way
            inbdim = iom_file(kiomid)%ndims(idvar)            ! number of dimensions in the file
            idmspc = inbdim                                   ! number of spatial dimensions in the file
            IF( iom_file(kiomid)%luld(idvar) )   idmspc = inbdim - 1
            IF( idmspc > 3 )   CALL ctl_stop(TRIM(clinfo), 'the file has more than 3 spatial dimensions this case is not coded...')
            !
            ! Identify the domain in case of jpdom_auto definition
            IF( idom == jpdom_auto .OR. idom == jpdom_auto_xy ) THEN
               idom = jpdom_global   ! default
               ! else: if the file name finishes with _xxxx.nc with xxxx any number
               ind1 = INDEX( clname, '_', back = .TRUE. ) + 1
               ind2 = INDEX( clname, '.', back = .TRUE. ) - 1
               IF( ind2 > ind1 ) THEN   ;   IF( VERIFY( clname(ind1:ind2), '0123456789' ) == 0 )   idom = jpdom_local   ;   ENDIF
            ENDIF
            !
            ! check the consistency between input array and data rank in the file
            !
            ! initializations
            itime = 1
            IF( PRESENT(ktime) ) itime = ktime
            !
            irankpv = 1 * COUNT( (/ llis1d /) ) + 2 * COUNT( (/ llis2d /) ) + 3 * COUNT( (/ llis3d /) )
            WRITE(clrankpv, fmt='(i1)') irankpv
            WRITE(cldmspc , fmt='(i1)') idmspc
            !
            IF(     idmspc <  irankpv ) THEN                     ! it seems we want to read more than we can...
               IF(     llis3d .AND. idmspc == 2 ) THEN           !   3D input array from 2D spatial data in the file:
                  llok = inlev == 1                              !     -> 3rd dimension must be equal to 1
               ELSEIF( llis3d .AND. idmspc == 1 ) THEN           !   3D input array from 1D spatial data in the file:
                  llok = inlev == 1 .AND. ishape(2) == 1         !     -> 2nd and 3rd dimensions must be equal to 1
               ELSEIF( llis3d .AND. idmspc == 2 ) THEN           !   2D input array from 1D spatial data in the file:
                  llok = ishape(2) == 1                          !     -> 2nd dimension must be equal to 1
               ELSE
                  llok = .FALSE.
               ENDIF
               IF( .NOT. llok )   CALL ctl_stop( TRIM(clinfo), 'The file has only '//cldmspc//' spatial dimension',   &
                  &                                            '=> cannot read a true '//clrankpv//'D array from this file...' )
            ELSEIF( idmspc == irankpv ) THEN
               IF( llis1d .AND. idom /= jpdom_unknown )   &
                  &   CALL ctl_stop( TRIM(clinfo), 'case not coded...You must use jpdom_unknown' )
            ELSEIF( idmspc >  irankpv ) THEN                     ! it seems we want to read less than we should...
                  IF( llis2d .AND. itime == 1 .AND. idimsz(3) == 1 .AND. idmspc == 3 ) THEN
                     CALL ctl_warn( TRIM(clinfo), '2D array input but 3 spatial dimensions in the file...'              ,   &
                           &         'As the size of the z dimension is 1 and as we try to read the first record, ',   &
                           &         'we accept this case, even if there is a possible mix-up between z and time dimension' )
                     idmspc = idmspc - 1
                  !!GS: possibility to read 3D ABL atmopsheric forcing and use 1st level to force BULK simulation
                  !ELSE
                  !   CALL ctl_stop( TRIM(clinfo), 'To keep iom lisibility, when reading a '//clrankpv//'D array,',   &
                  !      &                         'we do not accept data with '//cldmspc//' spatial dimensions'  ,   &
                  !      &                         'Use ncwa -a to suppress the unnecessary dimensions' )
                  ENDIF
            ENDIF
            !
            ! definition of istart and icnt
            !
            icnt  (:) = 1              ! default definition (simple way to deal with special cases listed above)
            istart(:) = 1              ! default definition (simple way to deal with special cases listed above)
            istart(idmspc+1) = itime   ! temporal dimenstion
            !
            IF( idom == jpdom_unknown ) THEN
               IF( PRESENT(kstart) .AND. idom /= jpdom_auto_xy ) THEN
                  istart(1:idmspc) = kstart(1:idmspc)
                  icnt  (1:idmspc) = kcount(1:idmspc)
               ELSE
                  icnt  (1:idmspc) = idimsz(1:idmspc)
               ENDIF
            ELSE   !   not a 1D array as pv(sd)p1d requires jpdom_unknown
               ! we do not read the overlap and the extra-halos -> from Nis0 to Nie0 and from Njs0 to Nje0
               IF( idom == jpdom_global )   istart(1:2) = (/ mig0(Nis0), mjg0(Njs0) /)
               icnt(1:2) = (/ Ni_0, Nj_0 /)
               IF( llis3d ) THEN
                  IF( idom == jpdom_auto_xy ) THEN
                     istart(3) = kstart(3)
                     icnt  (3) = kcount(3)
                  ELSE
                     icnt  (3) = inlev
                  ENDIF
               ENDIF
            ENDIF
            !
            ! check that istart and icnt can be used with this file
            !-
            DO jl = 1, jpmax_dims
               itmp = istart(jl)+icnt(jl)-1
               IF( itmp > idimsz(jl) .AND. idimsz(jl) /= 0 ) THEN
                  WRITE( ctmp1, FMT="('(istart(', i1, ') + icnt(', i1, ') - 1) = ', i5)" ) jl, jl, itmp
                  WRITE( ctmp2, FMT="(' is larger than idimsz(', i1,') = ', i5)"         ) jl, idimsz(jl)
                  CALL ctl_stop( TRIM(clinfo), 'start and count too big regarding to the size of the data, ', ctmp1, ctmp2 )
               ENDIF
            END DO
            !
            ! check that icnt matches the input array
            !-
            IF( idom == jpdom_unknown ) THEN
               ctmp1 = 'd'
            ELSE                    ! we must redefine ishape as we don't read the full array
               IF( llis2d ) THEN
                  IF( llsp ) THEN   ;   ishape(1:2) = SHAPE(pvsp2d(Nis0:Nie0,Njs0:Nje0  ))
                  ELSE              ;   ishape(1:2) = SHAPE(pvdp2d(Nis0:Nie0,Njs0:Nje0  ))
                  ENDIF
                  ctmp1 = 'd(Nis0:Nie0,Njs0:Nje0)'
               ENDIF
               IF( llis3d ) THEN
                  IF( llsp ) THEN   ;   ishape(1:3) = SHAPE(pvsp3d(Nis0:Nie0,Njs0:Nje0,:))
                  ELSE              ;   ishape(1:3) = SHAPE(pvdp3d(Nis0:Nie0,Njs0:Nje0,:))
                  ENDIF
                  ctmp1 = 'd(Nis0:Nie0,Njs0:Nje0,:)'
               ENDIF
            ENDIF
            DO jl = 1, irankpv
               WRITE( ctmp2, FMT="(', ', i1,'): ', i5,' /= icnt(', i1,'):', i5)" ) jl, ishape(jl), jl, icnt(jl)
               IF( llsp ) THEN
                  IF( ishape(jl) /= icnt(jl) )   CALL ctl_stop( TRIM(clinfo), 'size(pvsp'//clrankpv//TRIM(ctmp1)//TRIM(ctmp2) )
               ELSE
                  IF( ishape(jl) /= icnt(jl) )   CALL ctl_stop( TRIM(clinfo), 'size(pvdp'//clrankpv//TRIM(ctmp1)//TRIM(ctmp2) )
               ENDIF
            END DO

         ENDIF

         ! read the data
         !-
         IF( idvar > 0 .AND. istop == nstop ) THEN   ! no additional errors until this point...
            !
            ! find the right index of the array to be read
            IF( idom /= jpdom_unknown ) THEN   ;   ix1 = Nis0   ;   ix2 = Nie0      ;   iy1 = Njs0   ;   iy2 = Nje0
            ELSE                               ;   ix1 = 1      ;   ix2 = icnt(1)   ;   iy1 = 1      ;   iy2 = icnt(2)
            ENDIF

            CALL iom_nf90_get( kiomid, idvar, inbdim, istart, icnt, ix1, ix2, iy1, iy2,   &
               &               pvsp1d, pvsp2d, pvsp3d, pvdp1d, pvdp2d, pvdp3d )

            IF( istop == nstop ) THEN   ! no additional errors until this point...
               IF(lwp) WRITE(numout,"(10x,' read ',a,' (rec: ',i6,') in ',a,' ok')") TRIM(cdvar), itime, TRIM(iom_file(kiomid)%name)

               cl_type = 'T'
               IF( PRESENT(cd_type) )   cl_type = cd_type
               IF( llsp ) THEN
                  zsgn_sp = 1._sp
                  IF( PRESENT(psgn_sp) )   zsgn_sp = psgn_sp
                  IF(     llis2d .AND. idom /= jpdom_unknown .AND. cl_type /= 'Z' ) THEN
                     CALL lbc_lnk( 'iom', pvsp2d, cl_type, zsgn_sp, kfillmode = kfill )
                  ELSEIF( llis3d .AND. idom /= jpdom_unknown .AND. cl_type /= 'Z' ) THEN
                     CALL lbc_lnk( 'iom', pvsp3d, cl_type, zsgn_sp, kfillmode = kfill )
                  ENDIF
               ELSE
                  zsgn_dp = 1._dp
                  IF( PRESENT(psgn_dp) )   zsgn_dp = psgn_dp
                  IF(     llis2d .AND. idom /= jpdom_unknown .AND. cl_type /= 'Z' ) THEN
                     CALL lbc_lnk( 'iom', pvdp2d, cl_type, zsgn_dp, kfillmode = kfill )
                  ELSEIF( llis3d .AND. idom /= jpdom_unknown .AND. cl_type /= 'Z' ) THEN
                     CALL lbc_lnk( 'iom', pvdp3d, cl_type, zsgn_dp, kfillmode = kfill )
                  ENDIF
               ENDIF
               !--- overlap areas and extra hallows (mpp)
               !
            ELSE
               ! return if istop == nstop is false
               RETURN
            ENDIF
         ELSE
            ! return if statment idvar > 0 .AND. istop == nstop is false
            RETURN
         ENDIF
         !
      ELSE        ! read using XIOS. Only if key_xios is defined
         istop = istop + 1
         clinfo = 'Can not use XIOS in iom_get_123d, file: '//TRIM(clname)//', var:'//TRIM(cdvar)
      ENDIF

      !--- Apply scale_factor and offset
      IF( llsp ) THEN
         zscf_sp = iom_file(kiomid)%scf(idvar)      ! scale factor
         zofs_sp = iom_file(kiomid)%ofs(idvar)      ! offset
         IF(     llis1d ) THEN
            IF( zscf_sp /= 1._sp )   pvsp1d(:    ) = pvsp1d(:    ) * zscf_sp
            IF( zofs_sp /= 0._sp )   pvsp1d(:    ) = pvsp1d(:    ) + zofs_sp
         ELSEIF( llis2d ) THEN
            IF( zscf_sp /= 1._sp )   pvsp2d(:,:  ) = pvsp2d(:,:  ) * zscf_sp
            IF( zofs_sp /= 0._sp )   pvsp2d(:,:  ) = pvsp2d(:,:  ) + zofs_sp
         ELSEIF( llis3d ) THEN
            IF( zscf_sp /= 1._sp )   pvsp3d(:,:,:) = pvsp3d(:,:,:) * zscf_sp
            IF( zofs_sp /= 0._sp )   pvsp3d(:,:,:) = pvsp3d(:,:,:) + zofs_sp
         ENDIF
      ELSE
         zscf_dp = iom_file(kiomid)%scf(idvar)      ! scale factor
         zofs_dp = iom_file(kiomid)%ofs(idvar)      ! offset
         IF(     llis1d ) THEN
            IF( zscf_dp /= 1._dp )   pvdp1d(:    ) = pvdp1d(:    ) * zscf_dp
            IF( zofs_dp /= 0._dp )   pvdp1d(:    ) = pvdp1d(:    ) + zofs_dp
         ELSEIF( llis2d ) THEN
            IF( zscf_dp /= 1._dp )   pvdp2d(:,:  ) = pvdp2d(:,:  ) * zscf_dp
            IF( zofs_dp /= 0._dp )   pvdp2d(:,:  ) = pvdp2d(:,:  ) + zofs_dp
         ELSEIF( llis3d ) THEN
            IF( zscf_dp /= 1._dp )   pvdp3d(:,:,:) = pvdp3d(:,:,:) * zscf_dp
            IF( zofs_dp /= 0._dp )   pvdp3d(:,:,:) = pvdp3d(:,:,:) + zofs_dp
         ENDIF
      ENDIF
      !
   END SUBROUTINE iom_get_123d

   SUBROUTINE iom_get_var( cdname, z2d)
      CHARACTER(LEN=*), INTENT(in ) ::   cdname
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d
      IF( .FALSE. )   WRITE(numout,*) cdname, z2d ! useless test to avoid compilation warnings
   END SUBROUTINE iom_get_var


   FUNCTION iom_getszuld ( kiomid )
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  iom_getszuld  ***
      !!
      !! ** Purpose : get the size of the unlimited dimension in a file
      !!              (return -1 if not found)
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kiomid   ! file Identifier
      !
      INTEGER                ::   iom_getszuld
      !!-----------------------------------------------------------------------
      iom_getszuld = -1
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%iduld > 0 )   iom_getszuld = iom_file(kiomid)%lenuld
      ENDIF
   END FUNCTION iom_getszuld


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_chkatt
   !!----------------------------------------------------------------------
   SUBROUTINE iom_chkatt( kiomid, cdatt, llok, ksize, cdvar )
      INTEGER         , INTENT(in   )                 ::   kiomid    ! Identifier of the file
      CHARACTER(len=*), INTENT(in   )                 ::   cdatt     ! Name of the attribute
      LOGICAL         , INTENT(  out)                 ::   llok      ! Error code
      INTEGER         , INTENT(  out), OPTIONAL       ::   ksize     ! Size of the attribute array
      CHARACTER(len=*), INTENT(in   ), OPTIONAL       ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_chkatt( kiomid, cdatt, llok, ksize=ksize, cdvar=cdvar )
      ENDIF
      !
   END SUBROUTINE iom_chkatt

   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_getatt
   !!----------------------------------------------------------------------
   SUBROUTINE iom_g0d_iatt( kiomid, cdatt, katt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      INTEGER               , INTENT(  out)           ::   katt0d    ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt,  katt0d =  katt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g0d_iatt

   SUBROUTINE iom_g1d_iatt( kiomid, cdatt, katt1d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      INTEGER, DIMENSION(:) , INTENT(  out)           ::   katt1d    ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt,  katt1d =  katt1d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g1d_iatt

   SUBROUTINE iom_g0d_ratt( kiomid, cdatt, patt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      REAL(wp)              , INTENT(  out)           ::   patt0d    ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt,  patt0d =  patt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g0d_ratt

   SUBROUTINE iom_g1d_ratt( kiomid, cdatt, patt1d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      REAL(wp), DIMENSION(:), INTENT(  out)           ::   patt1d    ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt,  patt1d =  patt1d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g1d_ratt

   SUBROUTINE iom_g0d_catt( kiomid, cdatt, cdatt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      CHARACTER(len=*)      , INTENT(  out)           ::   cdatt0d   ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt, cdatt0d = cdatt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g0d_catt


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_putatt
   !!----------------------------------------------------------------------
   SUBROUTINE iom_p0d_iatt( kiomid, cdatt, katt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      INTEGER               , INTENT(in   )           ::   katt0d    ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt,  katt0d =  katt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p0d_iatt

   SUBROUTINE iom_p1d_iatt( kiomid, cdatt, katt1d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      INTEGER, DIMENSION(:) , INTENT(in   )           ::   katt1d    ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt,  katt1d =  katt1d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p1d_iatt

   SUBROUTINE iom_p0d_ratt( kiomid, cdatt, patt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      REAL(wp)              , INTENT(in   )           ::   patt0d    ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt,  patt0d =  patt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p0d_ratt

   SUBROUTINE iom_p1d_ratt( kiomid, cdatt, patt1d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      REAL(wp), DIMENSION(:), INTENT(in   )           ::   patt1d    ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt,  patt1d =  patt1d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p1d_ratt

   SUBROUTINE iom_p0d_catt( kiomid, cdatt, cdatt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt0d   ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt, cdatt0d = cdatt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p0d_catt


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_rstput
   !!----------------------------------------------------------------------
   SUBROUTINE iom_rp0d_sp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(sp)        , INTENT(in)                         ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx                ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvsp0d = pvar )
      !
   END SUBROUTINE iom_rp0d_sp

   SUBROUTINE iom_rp0d_dp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(dp)        , INTENT(in)                         ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx                ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvdp0d = pvar )
      !
   END SUBROUTINE iom_rp0d_dp


   SUBROUTINE iom_rp1d_sp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(sp)        , INTENT(in), DIMENSION(          :) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx                ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvsp1d = pvar )
      !
   END SUBROUTINE iom_rp1d_sp

   SUBROUTINE iom_rp1d_dp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(dp)        , INTENT(in), DIMENSION(          :) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx                ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvdp1d = pvar )
      !
   END SUBROUTINE iom_rp1d_dp

   
   SUBROUTINE iom_rp2d_sp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(sp)        , INTENT(in), DIMENSION(:,    :    ) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL            :: llx
      INTEGER            :: ivid   ! variable id
      CHARACTER(LEN=lc)  :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvsp2d = pvar )
      !
   END SUBROUTINE iom_rp2d_sp

   SUBROUTINE iom_rp2d_dp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(dp)        , INTENT(in), DIMENSION(:,    :    ) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvdp2d = pvar )
      !
   END SUBROUTINE iom_rp2d_dp


   SUBROUTINE iom_rp3d_sp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(sp)        , INTENT(in),       DIMENSION(:,:,:) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx                 ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvsp3d = pvar )
      !
   END SUBROUTINE iom_rp3d_sp

   SUBROUTINE iom_rp3d_dp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(dp)        , INTENT(in),       DIMENSION(:,:,:) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx                 ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvdp3d = pvar )
      !
   END SUBROUTINE iom_rp3d_dp

   SUBROUTINE iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvsp0d, pvsp1d, pvsp2d, pvsp3d, pvdp0d, pvdp1d, pvdp2d, pvdp3d )
      INTEGER                    , INTENT(in)           ::   kt       ! ocean time-step
      INTEGER                    , INTENT(in)           ::   kwrite   ! writing time-step
      INTEGER                    , INTENT(in)           ::   kiomid   ! Identifier of the file
      CHARACTER(len=*)           , INTENT(in)           ::   cdvar    ! time axis name
      INTEGER                    , INTENT(in), OPTIONAL ::   ktype    ! variable external type
      REAL(sp)                   , INTENT(in), OPTIONAL ::   pvsp0d    ! read field (0D case), single precision
      REAL(sp), DIMENSION(:)     , INTENT(in), OPTIONAL ::   pvsp1d    ! read field (1D case), single precision
      REAL(sp), DIMENSION(:,:)   , INTENT(in), OPTIONAL ::   pvsp2d    ! read field (2D case), single precision
      REAL(sp), DIMENSION(:,:,:) , INTENT(in), OPTIONAL ::   pvsp3d    ! read field (3D case), single precision
      REAL(dp)                   , INTENT(in), OPTIONAL ::   pvdp0d    ! read field (0D case), double precision
      REAL(dp), DIMENSION(:)     , INTENT(in), OPTIONAL ::   pvdp1d    ! read field (1D case), double precision
      REAL(dp), DIMENSION(:,:)   , INTENT(in), OPTIONAL ::   pvdp2d    ! read field (2D case), double precision
      REAL(dp), DIMENSION(:,:,:) , INTENT(in), OPTIONAL ::   pvdp3d    ! read field (3D case), double precision
      !
      LOGICAL           :: llx                 ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL set_xios_context(kiomid, context)

      llx = .NOT. (context == "NONE")

      IF( llx ) THEN
      ELSE
         IF( kiomid > 0 ) THEN
            IF( iom_file(kiomid)%nfid > 0 ) THEN
               ivid = iom_varid( kiomid, cdvar, ldstop = .FALSE. )
               CALL iom_nf90_rstput( kt, kwrite, kiomid, cdvar, ivid, ktype, pvsp0d, pvsp1d, pvsp2d, pvsp3d,   &
                  &                                                          pvdp0d, pvdp1d, pvdp2d, pvdp3d )
            ENDIF
         ENDIF
      ENDIF
   END SUBROUTINE iom_rp0123d


  SUBROUTINE iom_delay_rst( cdaction, cdcpnt, kncid )
      !!---------------------------------------------------------------------
      !!   Routine iom_delay_rst: used read/write restart related to mpp_delay
      !!
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   ) ::   cdaction        !
      CHARACTER(len=*), INTENT(in   ) ::   cdcpnt
      INTEGER         , INTENT(in   ) ::   kncid
      !
      INTEGER  :: ji
      INTEGER  :: indim
      LOGICAL  :: llattexist
      REAL(wp), ALLOCATABLE, DIMENSION(:) ::   zreal1d
      !!---------------------------------------------------------------------
      !
      !                                      ===================================
      IF( TRIM(cdaction) == 'READ' ) THEN   ! read restart related to mpp_delay !
         !                                   ===================================
         DO ji = 1, nbdelay
            IF ( c_delaycpnt(ji) == cdcpnt ) THEN
               CALL iom_chkatt( kncid, 'DELAY_'//c_delaylist(ji), llattexist, indim )
               IF( llattexist )  THEN
                  ALLOCATE( todelay(ji)%z1d(indim) )
                  CALL iom_getatt( kncid, 'DELAY_'//c_delaylist(ji), todelay(ji)%z1d(:) )
                  ndelayid(ji) = 0   ! set to 0 to specify that the value was read in the restart
               ENDIF
           ENDIF
         END DO
         !                                   ====================================
      ELSE                                  ! write restart related to mpp_delay !
         !                                   ====================================
         DO ji = 1, nbdelay   ! save only ocean delayed global communication variables
            IF ( c_delaycpnt(ji) == cdcpnt ) THEN
               IF( ASSOCIATED(todelay(ji)%z1d) ) THEN
                  CALL mpp_delay_rcv(ji)   ! make sure %z1d is received
                  CALL iom_putatt( kncid, 'DELAY_'//c_delaylist(ji), todelay(ji)%z1d(:) )
               ENDIF
            ENDIF
         END DO
         !
      ENDIF

   END SUBROUTINE iom_delay_rst



   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_put
   !!----------------------------------------------------------------------
   SUBROUTINE iom_p0d_sp( cdname, pfield0d )
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      REAL(sp)        , INTENT(in) ::   pfield0d
      !!      REAL(wp)        , DIMENSION(jpi,jpj) ::   zz     ! masson
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield0d   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_p0d_sp

   SUBROUTINE iom_p0d_dp( cdname, pfield0d )
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      REAL(dp)        , INTENT(in) ::   pfield0d
!!      REAL(wp)        , DIMENSION(jpi,jpj) ::   zz     ! masson
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield0d   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_p0d_dp


   SUBROUTINE iom_p1d_sp( cdname, pfield1d )
      CHARACTER(LEN=*)          , INTENT(in) ::   cdname
      REAL(sp),     DIMENSION(:), INTENT(in) ::   pfield1d
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield1d   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_p1d_sp

   SUBROUTINE iom_p1d_dp( cdname, pfield1d )
      CHARACTER(LEN=*)          , INTENT(in) ::   cdname
      REAL(dp),     DIMENSION(:), INTENT(in) ::   pfield1d
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield1d   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_p1d_dp

   SUBROUTINE iom_p2d_sp( cdname, pfield2d )
      CHARACTER(LEN=*)            , INTENT(in) ::   cdname
      REAL(sp),     DIMENSION(:,:), INTENT(in) ::   pfield2d
      IF( iom_use(cdname) ) THEN
         WRITE(numout,*) pfield2d   ! iom_use(cdname) = .F. -> useless test to avoid compilation warnings
      ENDIF
   END SUBROUTINE iom_p2d_sp

   SUBROUTINE iom_p2d_dp( cdname, pfield2d )
      CHARACTER(LEN=*)            , INTENT(in) ::   cdname
      REAL(dp),     DIMENSION(:,:), INTENT(in) ::   pfield2d
      IF( iom_use(cdname) ) THEN
         WRITE(numout,*) pfield2d   ! iom_use(cdname) = .F. -> useless test to avoid compilation warnings
      ENDIF
   END SUBROUTINE iom_p2d_dp

   SUBROUTINE iom_p3d_sp( cdname, pfield3d )
      CHARACTER(LEN=*)                , INTENT(in) ::   cdname
      REAL(sp),       DIMENSION(:,:,:), INTENT(in) ::   pfield3d
      IF( iom_use(cdname) ) THEN
         WRITE(numout,*) pfield3d   ! iom_use(cdname) = .F. -> useless test to avoid compilation warnings
      ENDIF
   END SUBROUTINE iom_p3d_sp

   SUBROUTINE iom_p3d_dp( cdname, pfield3d )
      CHARACTER(LEN=*)                , INTENT(in) ::   cdname
      REAL(dp),       DIMENSION(:,:,:), INTENT(in) ::   pfield3d
      IF( iom_use(cdname) ) THEN
         WRITE(numout,*) pfield3d   ! iom_use(cdname) = .F. -> useless test to avoid compilation warnings
      ENDIF
   END SUBROUTINE iom_p3d_dp

   SUBROUTINE iom_p4d_sp( cdname, pfield4d )
      CHARACTER(LEN=*)                , INTENT(in) ::   cdname
      REAL(sp),       DIMENSION(:,:,:,:), INTENT(in) ::   pfield4d
      IF( iom_use(cdname) ) THEN
         WRITE(numout,*) pfield4d   ! iom_use(cdname) = .F. -> useless test to avoid compilation warnings
      ENDIF
   END SUBROUTINE iom_p4d_sp

   SUBROUTINE iom_p4d_dp( cdname, pfield4d )
      CHARACTER(LEN=*)                , INTENT(in) ::   cdname
      REAL(dp),       DIMENSION(:,:,:,:), INTENT(in) ::   pfield4d
      IF( iom_use(cdname) ) THEN
         WRITE(numout,*) pfield4d   ! iom_use(cdname) = .F. -> useless test to avoid compilation warnings
      ENDIF
   END SUBROUTINE iom_p4d_dp

   !!----------------------------------------------------------------------
   !!   NOT 'key_xios'                               a few dummy routines
   !!----------------------------------------------------------------------
   SUBROUTINE iom_setkt( kt, cdname )
      INTEGER         , INTENT(in)::   kt
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      IF( .FALSE. )   WRITE(numout,*) kt, cdname   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_setkt

   SUBROUTINE iom_context_finalize( cdname )
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      IF( .FALSE. )   WRITE(numout,*)  cdname   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_context_finalize

   SUBROUTINE iom_update_file_name( cdid )
      CHARACTER(LEN=*), INTENT(in) ::   cdid
      IF( .FALSE. )   WRITE(numout,*)  cdid   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_update_file_name


   LOGICAL FUNCTION iom_use( cdname )
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      iom_use = .FALSE.
   END FUNCTION iom_use

   SUBROUTINE iom_miss_val( cdname, pmiss_val )
      CHARACTER(LEN=*), INTENT(in ) ::   cdname
      REAL(wp)        , INTENT(out) ::   pmiss_val
      REAL(dp)                      ::   ztmp_pmiss_val
      IF( .FALSE. )   WRITE(numout,*) cdname, pmiss_val   ! useless test to avoid compilation warnings
      IF( .FALSE. )   pmiss_val = 0._wp                   ! useless assignment to avoid compilation warnings
   END SUBROUTINE iom_miss_val

   !!======================================================================
END MODULE iom
