#include <define.h>
!
module coupling
#ifdef AOMG
! OASIS3-MCT variables for atmos-ocean coupling
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_, latg_
!   use mod_oasis
!-------------------------------------------------------------------------------
!
!  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(6,37)   ! real
   INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307) ! double

   integer :: comp_id                   ! component indentification 
   integer :: localComm                 ! local MPI communicator and Initialized
   integer :: il_paral_size
   integer :: part_id
   integer :: var_nodims(2)
   integer :: var_type
   integer :: var_actual_shape(4)
   character(len=5) :: comp_name 
   character(len=3) :: cplmodel
   integer, dimension(:), allocatable :: il_paral
!!ice coupled
   integer :: var_id(15)                !-- array size = number of coupling variables

   character(len=6)  :: var_name1
   character(len=7)  :: var_name2
   character(len=7)  :: var_name3
   character(len=7)  :: var_name4
   character(len=8)  :: var_name5
   character(len=8)  :: var_name6
   character(len=9)  :: var_name7
   character(len=9)  :: var_name8
   character(len=6)  :: var_name9
   character(len=7)  :: var_name10
   character(len=7)  :: var_name11
   character(len=6)  :: var_name12
   character(len=7)  :: var_name13
   character(len=7)  :: var_name14
   character(len=7)  :: var_name15

!!oce only
!   integer :: var_id(6)                !-- array size = number of coupling variables

!   character(len=6)  :: var_name1
!   character(len=8)  :: var_name2
!   character(len=8)  :: var_name3
!   character(len=6)  :: var_name4
!   character(len=7)  :: var_name5
!   character(len=6)  :: var_name6

!
#endif /* AOMG */
end module coupling
