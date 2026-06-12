










MODULE vremap
!$AGRIF_DO_NOT_TREAT
   !!======================================================================
   !!                       ***  MODULE  vremap  ***
   !! Ocean physics:  Vertical remapping routines
   !!
   !!======================================================================
   !! History : 4.0  !  2019-09  (Jérôme Chanut)  Original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!
   !!----------------------------------------------------------------------
   USE par_oce
   USE ppr_1d   ! D. Engwirda piecewise polynomial reconstruction library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   reconstructandremap, remap_linear

   !! * Substitutions




   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: vremap 11573 2019-09-19 09:18:03Z jchanut $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE reconstructandremap(ptin, phin, ptout, phout, kjpk_in, kjpk_out, kn_var)
      !!----------------------------------------------------------------------
      !!                    *** ROUTINE  reconstructandremap ***
      !!
      !! ** Purpose :   Conservative remapping of a vertical column 
      !!                from one set of layers to an other one.
      !!
      !! ** Method  :   Uses D. Engwirda Piecewise Polynomial Reconstruction library.
      !!                https://github.com/dengwirda/PPR
      !!                
      !!
      !! References :   Engwirda, Darren & Kelley, Maxwell. (2015). A WENO-type 
      !!                slope-limiter for a family of piecewise polynomial methods. 
      !!                https://arxiv.org/abs/1606.08188
      !!-----------------------------------------------------------------------
      INTEGER , INTENT(in   )                      ::   kjpk_in    ! Number of input levels
      INTEGER , INTENT(in   )                      ::   kjpk_out   ! Number of output levels
      INTEGER , INTENT(in   )                      ::   kn_var     ! Number of variables
      REAL(wp), INTENT(in   ), DIMENSION(kjpk_in)  ::   phin       ! Input thicknesses
      REAL(wp), INTENT(in   ), DIMENSION(kjpk_out) ::   phout      ! Output thicknesses
      REAL(wp), INTENT(in   ), DIMENSION(kjpk_in , kn_var) ::   ptin       ! Input data
      REAL(wp), INTENT(inout), DIMENSION(kjpk_out, kn_var) ::   ptout      ! Remapped data
      !
      INTEGER, PARAMETER :: ndof = 1
      INTEGER  :: jk, jn
      REAL(dp) ::  zwin(kjpk_in+1) ,  ztin(ndof, kn_var, kjpk_in)    ! rmap1d uses dp
      REAL(dp) :: zwout(kjpk_out+1), ztout(ndof, kn_var, kjpk_out)   ! rmap1d uses dp
      TYPE(rmap_work) :: work
      TYPE(rmap_opts) :: opts
      TYPE(rcon_ends) :: bc_l(kn_var)
      TYPE(rcon_ends) :: bc_r(kn_var)
      !!--------------------------------------------------------------------
     
      ! Set interfaces and input data:
      zwin(1) = 0._wp
      DO jk = 2, kjpk_in + 1
         zwin(jk) = zwin(jk-1) + phin(jk-1) 
      END DO
      
      DO jn = 1, kn_var 
         DO jk = 1, kjpk_in
            ztin(ndof, jn, jk) =  ptin(jk, jn)
         END DO
      END DO

      zwout(1) = 0._wp
      DO jk = 2, kjpk_out + 1
         zwout(jk) = zwout(jk-1) + phout(jk-1) 
      END DO

      ! specify methods
!      opts%edge_meth = p1e_method     ! 1st-order edge interp.
!      opts%cell_meth = pcm_method
!      opts%cell_meth = plm_method     ! PLM method in cells
      opts%edge_meth = p3e_method     ! 3rd-order edge interp.
      opts%cell_meth = ppm_method     ! PPM method in cells    
!      opts%edge_meth = p5e_method     ! 5th-order edge interp.
!      opts%cell_meth = pqm_method     ! PQM method in cells

      ! limiter
!      opts%cell_lims = null_limit     ! no lim.
!      opts%cell_lims = weno_limit
      opts%cell_lims = mono_limit     ! monotone limiter   
 
      ! set boundary conditions
      bc_l%bcopt = bcon_loose         ! "loose" = extrapolate
      bc_r%bcopt = bcon_loose
!      bc_l%bcopt = bcon_slope        
!      bc_r%bcopt = bcon_slope

      ! init. method workspace
      CALL work%init(kjpk_in+1, kn_var, opts)

      ! remap
      CALL rmap1d(kjpk_in+1, kjpk_out+1, kn_var, ndof, &
      &           zwin, zwout, ztin, ztout,            &
      &           bc_l, bc_r, work, opts)

      ! clear method workspace
      CALL work%free()

      DO jn = 1, kn_var 
         DO jk = 1, kjpk_out
            ptout(jk, jn) = ztout(1, jn, jk)
         END DO
      END DO
            
   END SUBROUTINE reconstructandremap

   SUBROUTINE remap_linear(ptin, pzin, ptout, pzout, kjpk_in, kjpk_out, kn_var)
      !!----------------------------------------------------------------------
      !!                    *** ROUTINE  remap_linear ***
      !!
      !! ** Purpose :   Linear interpolation based on input/ouputs depths
      !!
      !!-----------------------------------------------------------------------
      INTEGER , INTENT(in   )                      ::   kjpk_in    ! Number of input levels
      INTEGER , INTENT(in   )                      ::   kjpk_out   ! Number of output levels
      INTEGER , INTENT(in   )                      ::   kn_var     ! Number of variables
      REAL(wp), INTENT(in   ), DIMENSION(kjpk_in)  ::   pzin       ! Input depths
      REAL(wp), INTENT(in   ), DIMENSION(kjpk_out) ::   pzout      ! Output depths
      REAL(wp), INTENT(in   ), DIMENSION(kjpk_in , kn_var) ::   ptin       ! Input data
      REAL(wp), INTENT(inout), DIMENSION(kjpk_out, kn_var) ::   ptout      ! Interpolated data
      !
      INTEGER  :: jkin, jkout, jn
      !!--------------------------------------------------------------------
      !      
      DO jkout = 1, kjpk_out !  Loop over destination grid
         !
         IF     ( pzout(jkout) <=  pzin(  1    ) ) THEN ! Surface extrapolation	
            DO jn = 1, kn_var 
! linear
!               ptout(jkout,jn) = ptin(1 ,jn) + &
!                               & (pzout(jkout) - pzin(1)) / (pzin(2)    - pzin(1)) &
!                               &                          * (ptin(2,jn) - ptin(1,jn))
               ptout(jkout,jn) = ptin(1,jn)
            END DO
         ELSEIF ( pzout(jkout) >= pzin(kjpk_in) ) THEN ! Bottom extrapolation 
            DO jn = 1, kn_var 
! linear
!               ptout(jkout,jn) = ptin(kjpk_in ,jn) + &
!                               & (pzout(jkout) - pzin(kjpk_in)) / (pzin(kjpk_in)    - pzin(kjpk_in-1)) &
!                               &                                * (ptin(kjpk_in,jn) - ptin(kjpk_in-1,jn))
               ptout(jkout,jn) = ptin(kjpk_in ,jn)
            END DO
         ELSEIF ( ( pzout(jkout) > pzin(1) ).AND.( pzout(jkout) < pzin(kjpk_in) )) THEN
            DO jkin = 1, kjpk_in - 1 !  Loop over source grid
               IF ( pzout(jkout) < pzin(jkin+1) ) THEN
                  DO jn = 1, kn_var
                     ptout(jkout,jn) =  ptin(jkin,jn) + &
                                     & (pzout(jkout) - pzin(jkin)) / (pzin(jkin+1)    - pzin(jkin)) &
                                     &                             * (ptin(jkin+1,jn) - ptin(jkin,jn))
                  END DO  
                  EXIT
               ENDIF  
            END DO
         ENDIF
         !
      END DO

   END SUBROUTINE remap_linear

   !!======================================================================
!$AGRIF_END_DO_NOT_TREAT
END MODULE vremap
