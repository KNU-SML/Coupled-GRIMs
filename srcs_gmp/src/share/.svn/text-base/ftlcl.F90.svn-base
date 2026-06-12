!
   function ftlcl(t,tdpd)                                                    
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: ftlcl        compute lcl temperature.                             
!                                                                               
! abstract: compute temperature at the lifting condensation level               
!   from temperature and dewpoint depression.  the formula used is              
!   a polynomial taken from phillips scv_moist_adiabat routine which 
!   empirically approximates the original exact implicit relationship.     
!   (this kind of approximation is customary (inman, 1969), but                 
!   the original source for this particular one is not yet known. -mi)          
!   its accuracy is about 0.03 kelvin for a dewpoint depression of 30.          
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   1991-05-07  iredell                made into inlinable function
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:  tlcl=ftlcl(t,tdpd)                                                    
!                                                                               
!   input argument list:                                                        
!     t        - real temperature in kelvin                                     
!     tdpd     - real dewpoint depression in kelvin                             
!                                                                               
!   output argument list:                                                       
!     ftlcl    - real temperature at the lcl in kelvin                          
!                                                                               
!-------------------------------------------------------------------------------
   real, parameter   ::  clcl1= 0.954442e+0,clcl2= 0.967772e-3,                &
                         clcl3=-0.710321e-3,clcl4=-0.270742e-5                          
!
   ftlcl=t-tdpd*(clcl1+clcl2*t+tdpd*(clcl3+clcl4*t))                         
!
   return                                                                    
   end function ftlcl                                                                       
!-------------------------------------------------------------------------------
