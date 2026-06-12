! ------------------------ code history --------------------------------
! source file:       vic_soileb.h
! purpose:           variable common block for solving soil temperature
! PRGMMR:            Ji Chen
! ORG:               ECPC/CRD/SIO/UCSD 
! DATE:              November 2003
! ----------------------------------------------------------------------

! ------------------- input variables ----------------------------------
   real  ::  TL             ! low node soil temperature (k)
   real  ::  TU             ! upper node soil temperature (k)
   real  ::  soilm          ! current soil moisture (m3/m3)
   real  ::  soilmx         ! maximum soil moisture (m3/m3)
   real  ::  bubble         ! bubbling pressure of soil (cm)
   real  ::  expt           ! parameter for Ksat with soil moisture (N/A)
   real  ::  ice0           ! soil ice at previous time step (m3/m3)
   real  ::  gamma          ! coef for computing soiltemp profile
   real  ::  fprime         ! coef for computing soiltemp profile
   real  ::  AA, BB, CC, DD, EE  ! coef for computing soiltemp profile
! ----------------------------------------------------------------------

! ------------------- output variables ---------------------------------
   real  ::  ice            ! current soil ice content
! ----------------------------------------------------------------------

   common /vicsoileb/TL,   TU, soilm, soilmx,                                  &
                 bubble, expt,  ice0,  gamma,                                  &
                 fprime,   AA,    BB,     CC,                                  &
                     DD,   EE,   ice
