!
   function fpvsq0(t)                                                         
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: fpvsq        compute saturation vapor pressure                    
!                                                                               
! abstract: compute saturation vapor pressure from the temperature.             
!   a quadratic interpolation is done between values in a lookup table          
!   computed in funct_svp_init. see documentation for fpvsx for details. 
!   input values outside table range are reset to table extrema.                
!   the interpolation accuracy is almost 9 decimal places.                      
!   on the cray, fpvsq is about 3 times faster than exact calculation.          
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   1991-05-07  iredell                made into inlinable function
!   1994-12-30  iredell                exact computation
!   1998-12-30  hong                   for water phase
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:   pvs=fpvsq(t)                                                         
!                                                                               
!   input argument list:                                                        
!     t        - real temperature in kelvin                                     
!                                                                               
!   output argument list:                                                       
!     fpvsq    - real saturation vapor pressure in kilopascals (cb)             
!                                                                               
! common blocks:                                                                
!   compvs   - scaling parameters and table computed in funct_svp_init.   
!                                                                               
!-------------------------------------------------------------------------------
   integer, parameter   ::  nx=7501
   real                 ::  tbpvs(nx),tbpvs0(nx)
   common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
   common/compvs/ c1xpvs,c2xpvs,tbpvs
!
   xj=min(max(c1xpvs0+c2xpvs0*t,1.),float(nx))
   jx=min(max(nint(xj),2),nx-1)
   dxj=xj-jx
   fj1=tbpvs0(jx-1)
   fj2=tbpvs0(jx)
   fj3=tbpvs0(jx+1)
   fpvsq0=(((fj3+fj1)/2-fj2)*dxj+(fj3-fj1)/2)*dxj+fj2
!
   return
   end function fpvsq0
!-------------------------------------------------------------------------------
