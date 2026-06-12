!
   module comadj
!-------------------------------------------------------------------------------
!
! Convective adjustment
!
! program history log:
!   2000-03-09  songyou hong           cvs verion setup
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
   real     ::  cappa      ,&   ! R/cp
                cpair      ,&   ! Specific heat of dry air
                epsilo     ,&   ! Ratio of h2o to dry air molecular weights 
                gravit     ,&   ! Gravitational acceleration 
                latvap     ,&   ! Latent heat of vaporization
                rhoh2o     ,&   ! Density of liquid water (STP)
                cldcp      ,&   ! Latvap/cpair (L/cp)
                clrh2o     ,&   ! Ratio of latvap to water vapor gas const
                latice     ,&   ! Latent heat of fusion
                rair       ,&   ! Gas constant for dry air
                rh2o       ,&   ! Gas constant for water vapor
                zvir       ,&   ! rh2o/rair - 1
                t0              ! Reference temperature for t-prime computations
   integer  ::  nlvdry          ! Number of levels to apply dry adjustment
!
   end module comadj 
