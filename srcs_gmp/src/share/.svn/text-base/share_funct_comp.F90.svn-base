#include "define.h"
   subroutine share_funct_comp
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [share_funct_comp] : use functions in share_funct_init
!      |
!      |--- [function: fpvs] *
!      |--- [function: fpvs0] *
!      |--- [function: fpvs_pa] *
!      |--- [function: fpvs0_pa] *
!      |--- [function: ftdp] *
!      |--- [function: fthe] *
!      |--- [function: ftma] *
!      |--- [function: fpkap] *
!
!-------------------------------------------------------------------------------
   end subroutine share_funct_comp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function fpvs(t)                                                          
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: fpvs         compute saturation vapor pressure                    
!                                                                               
! abstract: compute saturation vapor pressure from the temperature.             
!   a linear interpolation is done between values in a lookup table             
!   computed in funct_svp_init. see documentation for fpvsx for details. 
!   input values outside table range are reset to table extrema.                
!   the interpolation accuracy is almost 6 decimal places.                      
!   on the cray, fpvs is about 4 times faster than exact calculation.           
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   1991-05-07  iredell             made into inlinable function                  
!   1994-12-30  iredell             expand table                                  
!   1996-02-19  hong                ice effect
!   2009-10-01  jung-eun kim        f90 format with standard physics modules
!   2012-01-01  hong                pascal for saturation vapor pressure
!                                                                               
! usage:   pvs=fpvs(t)                                                          
!                                                                               
!   input argument list:                                                        
!     t        - real temperature in kelvin                                     
!                                                                               
!   output argument list:                                                       
!     fpvs     - real saturation vapor pressure in kilopascals (cb)             
!                                                                               
! common blocks:                                                                
!   compvs   - scaling parameters and table computed in funct_svp_init. 
!                                                                               
!-------------------------------------------------------------------------------
   integer, parameter   ::  nx=7501                                                        
   real                 ::  tbpvs(nx)                                                       
   common/compvs/ c1xpvs,c2xpvs,tbpvs                                        
!
   xj=min(max(c1xpvs+c2xpvs*t,1.),float(nx))                                 
   jx=min(xj,nx-1.)                                                          
   fpvs=tbpvs(jx)+(xj-jx)*(tbpvs(jx+1)-tbpvs(jx))                            
!
   return                                                                    
   end function fpvs                                                                      
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function fpvs0(t)                                                          
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: fpvs         compute saturation vapor pressure                    
!                                                                               
! abstract: compute saturation vapor pressure from the temperature.             
!   a linear interpolation is done between values in a lookup table             
!   computed in funct_svp_init. see documentation for fpvsx for details.
!   input values outside table range are reset to table extrema.                
!   the interpolation accuracy is almost 6 decimal places.                      
!   on the cray, fpvs is about 4 times faster than exact calculation.           
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   91-05-07  iredell             made into inlinable function                  
!   94-12-30  iredell             expand table                                  
!   96-02-19  hong                ice effect
!                                                                               
! usage:   pvs=fpvs(t)                                                          
!                                                                               
!   input argument list:                                                        
!     t        - real temperature in kelvin                                     
!                                                                               
!   output argument list:                                                       
!     fpvs     - real saturation vapor pressure in kilopascals (cb)             
!                                                                               
! common blocks:                                                                
!   compvs   - scaling parameters and table computed in funct_svp_init.  
!                                                                               
!-------------------------------------------------------------------------------
   integer, parameter   ::  nx=7501                                                        
   real                 ::  tbpvs0(nx)
   common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
!
   xj1=min(max(c1xpvs0+c2xpvs0*t,1.),float(nx))
   jx1=min(xj1,nx-1.)
   fpvs0=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
!
   return
   end function fpvs0
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function fpvs_pa(t)                                                          
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: fpvs         compute saturation vapor pressure                    
!                                                                               
! abstract: compute saturation vapor pressure from the temperature.             
!   a linear interpolation is done between values in a lookup table             
!   computed in funct_svp_init. see documentation for fpvsx for details. 
!   input values outside table range are reset to table extrema.                
!   the interpolation accuracy is almost 6 decimal places.                      
!   on the cray, fpvs is about 4 times faster than exact calculation.           
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   1991-05-07  iredell             made into inlinable function                  
!   1994-12-30  iredell             expand table                                  
!   1996-02-19  hong                ice effect
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!                                                                               
! usage:   pvs=fpvs(t)                                                          
!                                                                               
!   input argument list:                                                        
!     t        - real temperature in kelvin                                     
!                                                                               
!   output argument list:                                                       
!     fpvs     - real saturation vapor pressure (pa)             
!                                                                               
! common blocks:                                                                
!   compvs   - scaling parameters and table computed in funct_svp_init. 
!                                                                               
!-------------------------------------------------------------------------------
   integer, parameter   ::  nx=7501                                                        
   real                 ::  tbpvs(nx)                                                       
   common/compvs/ c1xpvs,c2xpvs,tbpvs                                        
!
   xj=min(max(c1xpvs+c2xpvs*t,1.),float(nx))                                 
   jx=min(xj,nx-1.)                                                          
   fpvs_pa=tbpvs(jx)+(xj-jx)*(tbpvs(jx+1)-tbpvs(jx))                            
   fpvs_pa=fpvs_pa * 1.e3
!
   return                                                                    
   end function fpvs_pa                                                                      
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function fpvs0_pa(t)                                                          
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: fpvs_pa         compute saturation vapor pressure                    
!                                                                               
! abstract: compute saturation vapor pressure from the temperature.             
!   a linear interpolation is done between values in a lookup table             
!   computed in funct_svp_init. see documentation for fpvsx for details.
!   input values outside table range are reset to table extrema.                
!   the interpolation accuracy is almost 6 decimal places.                      
!   on the cray, fpvs is about 4 times faster than exact calculation.           
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   91-05-07  iredell             made into inlinable function                  
!   94-12-30  iredell             expand table                                  
!   96-02-19  hong                ice effect
!                                                                               
! usage:   pvs=fpvs(t)                                                          
!                                                                               
!   input argument list:                                                        
!     t        - real temperature in kelvin                                     
!                                                                               
!   output argument list:                                                       
!     fpvs     - real saturation vapor pressure (pa)             
!                                                                               
! common blocks:                                                                
!   compvs   - scaling parameters and table computed in funct_svp_init.  
!                                                                               
!-------------------------------------------------------------------------------
   integer, parameter   ::  nx=7501                                                        
   real                 ::  tbpvs0(nx)
   common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
!
   xj1=min(max(c1xpvs0+c2xpvs0*t,1.),float(nx))
   jx1=min(xj1,nx-1.)
   fpvs0_pa=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
   fpvs0_pa=fpvs_pa0 * 1.e3
!
   return
   end function fpvs0_pa
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function ftdp(pv)                                                         
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: ftdp         compute saturation vapor pressure                    
!                                                                               
! abstract: compute dewpoint temperature from vapor pressure.                   
!   a linear interpolation is done between values in a lookup table             
!   computed in funct_dew_point_temp_init. 
!   see documentation for ftdpxg for details.                 
!   input values outside table range are reset to table extrema.                
!   the interpolation accuracy is better than 0.0005 kelvin                     
!   for dewpoint temperatures greater than 250 kelvin,                          
!   but decreases to 0.02 kelvin for a dewpoint around 230 kelvin.              
!   on the cray, ftdp is about 75 times faster than exact calculation.          
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   91-05-07  iredell             made into inlinable function                  
!   94-12-30  iredell             expand table                                  
!   96-02-19  hong                ice effect
!                                                                               
! usage:   tdp=ftdp(pv)                                                         
!                                                                               
!   input argument list:                                                        
!     pv       - real vapor pressure in kilopascals (cb)                        
!                                                                               
!   output argument list:                                                       
!     ftdp     - real dewpoint temperature in kelvin                            
!                                                                               
! common blocks:                                                                
!   comtdp   - scaling parameters and table computed in funct_dew_point_temp_init. 
!                                                                               
!-------------------------------------------------------------------------------
   integer, parameter   ::  nx=5001                                                        
   real                 ::  tbtdp(nx)                                                       
   common/comtdp/ c1xtdp,c2xtdp,tbtdp                                        
!
   xj=min(max(c1xtdp+c2xtdp*pv,1.),float(nx))                                
   jx=min(xj,nx-1.)                                                          
   ftdp=tbtdp(jx)+(xj-jx)*(tbtdp(jx+1)-tbtdp(jx))                            
!
   return                                                                    
   end function ftdp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function fthe(t,pk)                                                       
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: fthe         compute saturation vapor pressure                    
!                                                                               
! abstract: compute equivalent potential temperature at the lcl                 
!   from temperature and pressure over 100 kpa to the kappa power.              
!   a bilinear interpolation is done between values in a lookup table           
!   computed in funct_pot_temp_init. 
!   see documentation for fthex for details.                  
!   input values outside table range are reset to table extrema,                
!   except zero is returned for too cold or high lcls.                          
!   the interpolation accuracy is better than 0.01 kelvin.                      
!   on the cray, fthe is almost 6 times faster than exact calculation.          
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   91-05-07  iredell             made into inlinable function                  
!   94-12-30  iredell             expand table                                  
!   96-02-19  hong                ice effect
!                                                                               
! usage:   the=fthe(pv)                                                         
!                                                                               
!   input argument list:                                                        
!     t        - real lcl temperature in kelvin                                 
!     pk       - real lcl pressure over 100 kpa to the kappa power              
!                                                                               
!   output argument list:                                                       
!     fthe     - real equivalent potential temperature in kelvin                
!                                                                               
! common blocks:                                                                
!   comthe   - scaling parameters and table computed in funct_pot_temp_init.
!                                                                               
!-------------------------------------------------------------------------------
   integer, parameter   ::  nx=241,ny=151                                                  
   real                 ::  tbthe(nx,ny)                                                    
   common/comthe/ c1xthe,c2xthe,c1ythe,c2ythe,tbthe                          
!
   xj=min(c1xthe+c2xthe*t,float(nx))                                         
   yj=min(c1ythe+c2ythe*pk,float(ny))                                        
!
   if(xj.ge.1..and.yj.ge.1.) then                                            
     jx=min(xj,nx-1.)                                                        
     jy=min(yj,ny-1.)                                                        
     ftx1=tbthe(jx,jy)+(xj-jx)*(tbthe(jx+1,jy)-tbthe(jx,jy))                 
     ftx2=tbthe(jx,jy+1)+(xj-jx)*(tbthe(jx+1,jy+1)-tbthe(jx,jy+1))           
     fthe=ftx1+(yj-jy)*(ftx2-ftx1)                                           
   else                                                                      
     fthe=0.                                                                 
   endif                                                                     
!
   return                                                                    
   end function fthe
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function ftma(the,pk,qma)                                                 
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: ftma         compute moist adiabat temperature                    
!                                                                               
! abstract: compute temperature and specific humidity of a parcel               
!   lifted up a moist adiabat from equivalent potential temperature             
!   at the lcl and pressure over 100 kpa to the kappa power.                    
!   bilinear interpolations are done between values in a lookup table           
!   computed in funct_moist_adiabat_init. 
!   see documentation for ftmaxg for details.                 
!   input values outside table range are reset to table extrema.                
!   the interpolation accuracy is better than 0.01 kelvin                       
!   and 5.e-6 kg/kg for temperature and humidity, respectively.                 
!   on the cray, ftma is about 35 times faster than exact calculation.          
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   91-05-07  iredell             made into inlinable function                  
!   94-12-30  iredell             expand table                                  
!   96-02-19  hong                ice effect
!                                                                               
! usage:   tma=ftma(the,pk,qma)                                                 
!                                                                               
!   input argument list:                                                        
!     the      - real equivalent potential temperature in kelvin                
!     pk       - real pressure over 100 kpa to the kappa power                  
!                                                                               
!   output argument list:                                                       
!     ftma     - real parcel temperature in kelvin                              
!     qma      - real parcel specific humidity in kg/kg                         
!                                                                               
! common blocks:                                                                
!   comma    - scaling parameters and table computed in funct_moist_adiabat_init. 
!                                                                               
!-------------------------------------------------------------------------------
   integer, parameter   ::  nx=151,ny=121                                                  
   real                 ::  tbtma(nx,ny),tbqma(nx,ny)                                       
   common/comma/ c1xma,c2xma,c1yma,c2yma,tbtma,tbqma                         
!
   xj=min(max(c1xma+c2xma*the,1.),float(nx))                                 
   yj=min(max(c1yma+c2yma*pk,1.),float(ny))                                  
   jx=min(xj,nx-1.)                                                          
   jy=min(yj,ny-1.)                                                          
   ftx1=tbtma(jx,jy)+(xj-jx)*(tbtma(jx+1,jy)-tbtma(jx,jy))                   
   ftx2=tbtma(jx,jy+1)+(xj-jx)*(tbtma(jx+1,jy+1)-tbtma(jx,jy+1))             
   ftma=ftx1+(yj-jy)*(ftx2-ftx1)                                             
   qx1=tbqma(jx,jy)+(xj-jx)*(tbqma(jx+1,jy)-tbqma(jx,jy))                    
   qx2=tbqma(jx,jy+1)+(xj-jx)*(tbqma(jx+1,jy+1)-tbqma(jx,jy+1))              
   qma=qx1+(yj-jy)*(qx2-qx1)                                                 
!
   return                                                                    
   end function ftma
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function fpkap(p)                                                         
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: fpkap        raise surface pressure to the kappa power.           
!                                                                               
! abstract: raise surface pressure over 100 kpa to the kappa power              
!   using a rational weighted chebyshev approximation.                          
!   the numerator is of order 2 and the denominator is of order 4.              
!   the pressure range is 40-110 kpa and kappa is defined in fpkapx.            
!   the coeffiecients are set by calling gpkap or including bdpkap.             
!   the accuracy of this approximation is almost 8 decimal places.              
!   on the cray, fpkap is over 10 times faster than exact calculation.          
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   91-05-07  iredell             made into inlinable function                  
!   94-12-30  iredell             standardized kappa,                           
!   96-02-19  hong                ice effect
!                                 increased range and accuracy                  
!                                                                               
! usage:  pkap=fpkap(p)                                                         
!                                                                               
!   input argument list:                                                        
!     p        - real surface pressure in kilopascals (cb)                      
!                p should be in the range 40. to 110.                           
!                                                                               
!   output argument list:                                                       
!     fpkap    - real p/100 to the kappa power                                  
!                                                                               
! common blocks:                                                                
!   compkap  - coefficients for function fpkap                                  
!                                                                               
!-------------------------------------------------------------------------------
   common/compkap/ cn0,cn1,cn2,cd0,cd1,cd2,cd3,cd4                           
!
   fpkap=(cn0+p*(cn1+p*cn2))/(cd0+p*(cd1+p*(cd2+p*(cd3+p*cd4))))             
!
   return                                                                    
   end function fpkap                                                                       
!-------------------------------------------------------------------------------
