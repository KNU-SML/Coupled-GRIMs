   module blockdata_kappa                                                        
!-------------------------------------------------------------------------------
!
! subprogram: blockdata_kappa       set default coefficients for p**kappa                
!                                                                               
! abstract: set default coefficients for p**kappa.  subprogram gpkap was        
!   used to compute a rational weighted chebyshev approximation to              
!   pressure raised to the kappa power for use by the function fpkap.           
!   the numerator is of order 2 and the denominator is of order 4.              
!   the pressure range is 40-110 kpa and kappa is 287.05/1004.6.                
!   imsl subprogram ratch is used to compute coefficients.                      
!   see subprograms gpkap and fpkap for more details.                           
!                                                                               
! program history log:                                                          
!   1982-12-29  phillips    
!   1994-12-30  iredell                                                           
!                                                                               
! common blocks:                                                                
!   compkap  - coefficients for function fpkap                                  
!                                                                               
!-------------------------------------------------------------------------------
   common/compkap/ cn0,cn1,cn2,cd0,cd1,cd2,cd3,cd4                           
   data cn0,cn1,cn2                                                            &
    /   3.13198449e-1,5.78544829e-2, 8.35491871e-4/                          
   data cd0,cd1,cd2,cd3,cd4                                                    &
    /1.,8.15968401e-2,5.72839518e-4,-4.86959812e-7,5.24459889e-10/           
!
   end module blockdata_kappa 
