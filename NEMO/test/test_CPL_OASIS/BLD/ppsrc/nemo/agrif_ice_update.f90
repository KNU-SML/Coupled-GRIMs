











MODULE agrif_ice_update
   !!=====================================================================================
   !!                       ***  MODULE agrif_ice_update ***
   !! Nesting module :  update surface ocean boundary condition over ice from a child grid
   !!=====================================================================================
   !! History :  2.0   !  04-2008  (F. Dupont)               initial version
   !!            3.4   !  08-2012  (R. Benshila, C. Herbaut) update and EVP
   !!            4.0   !  2018     (C. Rousset)              SI3 compatibility
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Empty module                                             no sea-ice
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE agrif_ice_update_empty
      WRITE(*,*)  'agrif_ice_update : You should not have seen this print! error?'
   END SUBROUTINE agrif_ice_update_empty

   !!======================================================================
END MODULE agrif_ice_update
