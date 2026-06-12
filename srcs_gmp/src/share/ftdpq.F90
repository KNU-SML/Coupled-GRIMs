!
   function ftdpq(pv)                                                        
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: ftdpq        compute saturation vapor pressure                    
!                                                                               
! abstract: compute dewpoint temperature from vapor pressure.                   
!   a quadratic interpolation is done between values in a lookup table          
!   computed in funct_dew_point_temp_init. 
!   see documentation for ftdpxg for details.                 
!   input values outside table range are reset to table extrema.                
!   the interpolation accuracy is better than 0.00001 kelvin                    
!   for dewpoint temperatures greater than 250 kelvin,                          
!   but decreases to 0.002 kelvin for a dewpoint around 230 kelvin.             
!   on the cray, ftdpq is about 60 times faster than exact calculation.         
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   1991-05-07  iredell                made into inlinable function
!   1994-12-30  iredell                exact computation
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:   tdp=ftdpq(pv)                                                        
!                                                                               
!   input argument list:                                                        
!     pv       - real vapor pressure in kilopascals (cb)                        
!                                                                               
!   output argument list:                                                       
!     ftdpq    - real dewpoint temperature in kelvin                            
!                                                                               
! common blocks:                                                                
!   comtdp   - scaling parameters and table computed in funct_dew_point_temp_init. 
!                                                                               
!-------------------------------------------------------------------------------
   integer, parameter   :: nx=5001                                                        
   real                 ::  tbtdp(nx)                                                       
   common/comtdp/ c1xtdp,c2xtdp,tbtdp                                        
!
   xj=min(max(c1xtdp+c2xtdp*pv,1.),float(nx))                                
   jx=min(max(nint(xj),2),nx-1)                                              
   dxj=xj-jx                                                                 
   fj1=tbtdp(jx-1)                                                           
   fj2=tbtdp(jx)                                                             
   fj3=tbtdp(jx+1)                                                           
   ftdpq=(((fj3+fj1)/2-fj2)*dxj+(fj3-fj1)/2)*dxj+fj2                         
!
   return                                                                    
   end function ftdpq                                                                      
!-------------------------------------------------------------------------------
