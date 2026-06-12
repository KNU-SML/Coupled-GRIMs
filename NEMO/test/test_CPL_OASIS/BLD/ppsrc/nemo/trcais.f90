










MODULE trcais
   !!======================================================================
   !!                         ***  MODULE trcais  ***
   !!  Module for passive tracers in Antarctic ice sheet
   !!  delivered by iceberg and ice shelf freshwater fluxes
   !!======================================================================
   !! History :  1.0  ! 2020    (R. Person, O. Aumont, C. Ethe),
   !!======================================================================
   !!----------------------------------------------------------------------
   !!   Dummy module                              NO 3D passive tracer data
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ais_ini   ! Empty routine
   END SUBROUTINE trc_ais_ini
   SUBROUTINE trc_ais( kt, Kmm, Krhs )        ! Empty routine
      INTEGER, INTENT(in) :: kt, Kmm, Krhs ! time level indices
      WRITE(*,*) 'trc_ais: You should not have seen this print! error?', kt, Kmm, Krhs
   END SUBROUTINE trc_ais

   !!======================================================================
END MODULE trcais
