










MODULE agrif_ice_interp
   !!=====================================================================================
   !!                       ***  MODULE agrif_ice_interp ***
   !! Nesting module :  interp surface ice boundary condition from a parent grid
   !!=====================================================================================
   !! History :  2.0   !  04-2008  (F. Dupont)               initial version
   !!            3.4   !  09-2012  (R. Benshila, C. Herbaut) update and EVP
   !!            4.0   !  2018     (C. Rousset)              SI3 compatibility
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Empty module                                             no sea-ice
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE agrif_ice_interp_empty
      WRITE(*,*)  'agrif_ice_interp : You should not have seen this print! error?'
   END SUBROUTINE agrif_ice_interp_empty

   !!======================================================================
END MODULE agrif_ice_interp
