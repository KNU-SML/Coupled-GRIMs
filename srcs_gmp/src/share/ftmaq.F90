!
   function ftmaq(the,pk,qma)                                                
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: ftmaq        compute moist adiabat temperature                    
!                                                                               
! abstract: compute temperature and specific humidity of a parcel               
!   lifted up a moist adiabat from equivalent potential temperature             
!   at the lcl and pressure over 100 kpa to the kappa power.                    
!   biquadratic interpolations are done between values in a lookup table        
!   computed in funct_moist_adiabat_init. 
!   see documentation for ftmaxg for details.                 
!   input values outside table range are reset to table extrema.                
!   the interpolation accuracy is better than 0.0005 kelvin                     
!   and 1.e-7 kg/kg for temperature and humidity, respectively.                 
!   on the cray, ftmaq is about 25 times faster than exact calculation.         
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   1991-05-07  iredell                made into inlinable function
!   1994-12-30  iredell                exact computation
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:   tma=ftmaq(the,pk,qma)                                                
!                                                                               
!   input argument list:                                                        
!     the      - real equivalent potential temperature in kelvin                
!     pk       - real pressure over 100 kpa to the kappa power                  
!                                                                               
!   output argument list:                                                       
!     ftmaq    - real parcel temperature in kelvin                              
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
   jx=min(max(nint(xj),2),nx-1)                                              
   jy=min(max(nint(yj),2),ny-1)                                              
   dxj=xj-jx                                                                 
   dyj=yj-jy                                                                 
   ft11=tbtma(jx-1,jy-1)                                                     
   ft12=tbtma(jx-1,jy)                                                       
   ft13=tbtma(jx-1,jy+1)                                                     
   ft21=tbtma(jx,jy-1)                                                       
   ft22=tbtma(jx,jy)                                                         
   ft23=tbtma(jx,jy+1)                                                       
   ft31=tbtma(jx+1,jy-1)                                                     
   ft32=tbtma(jx+1,jy)                                                       
   ft33=tbtma(jx+1,jy+1)                                                     
   ftx1=(((ft31+ft11)/2-ft21)*dxj+(ft31-ft11)/2)*dxj+ft21                    
   ftx2=(((ft32+ft12)/2-ft22)*dxj+(ft32-ft12)/2)*dxj+ft22                    
   ftx3=(((ft33+ft13)/2-ft23)*dxj+(ft33-ft13)/2)*dxj+ft23                    
   ftmaq=(((ftx3+ftx1)/2-ftx2)*dyj+(ftx3-ftx1)/2)*dyj+ftx2                   
   q11=tbqma(jx-1,jy-1)                                                      
   q12=tbqma(jx-1,jy)                                                        
   q13=tbqma(jx-1,jy+1)                                                      
   q21=tbqma(jx,jy-1)                                                        
   q22=tbqma(jx,jy)                                                          
   q23=tbqma(jx,jy+1)                                                        
   q31=tbqma(jx+1,jy-1)                                                      
   q32=tbqma(jx+1,jy)                                                        
   q33=tbqma(jx+1,jy+1)                                                      
   qx1=(((q31+q11)/2-q21)*dxj+(q31-q11)/2)*dxj+q21                           
   qx2=(((q32+q12)/2-q22)*dxj+(q32-q12)/2)*dxj+q22                           
   qx3=(((q33+q13)/2-q23)*dxj+(q33-q13)/2)*dxj+q23                           
   qma=(((qx3+qx1)/2-qx2)*dyj+(qx3-qx1)/2)*dyj+qx2                           
!
   return                                                                    
   end function ftmaq                                                                      
!-------------------------------------------------------------------------------
