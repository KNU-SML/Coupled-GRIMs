!
   module comcmf
!-------------------------------------------------------------------------------
!
! Common block for moist convective mass flux procedure
!
! program history log:
!   2000-03-09  songyou hong           cvs verion setup
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   real     ::  cp      ,&     ! specific heat of dry air
                hlat    ,&     ! latent heat of vaporization
                grav    ,&     ! gravitational constant       
                c0      ,&     ! rain water autoconversion coefficient
                betamn  ,&     ! minimum overshoot parameter
                rhlat   ,&     ! reciprocal of hlat
                rcp     ,&     ! reciprocal of cp
                rgrav   ,&     ! reciprocal of grav
                cmftau  ,&     ! characteristic adjustment time scale
                rhoh2o  ,&     ! density of liquid water (STP)
                rgas    ,&     ! gas constant for dry air
                dzmin   ,&     ! minimum convective depth for precipitation
                tiny    ,&     ! arbitrary small num used in transport estimates
                eps     ,&     ! convergence criteria (machine dependent)
                tpmax   ,&     ! maximum acceptable t perturbation (degrees C)
                shpmax         ! maximum acceptable q perturbation (g/g)           
   integer  ::  iloc    ,&     ! longitude location for diagnostics
                jloc    ,&     ! latitude  location for diagnostics
                nsloc   ,&     ! nstep for which to produce diagnostics
                limcnv         ! top interface level limit for convection
   logical  ::  rlxclm         ! logical to relax column versus cloud triplet
!
   end module comcmf 
