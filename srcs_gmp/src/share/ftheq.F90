!
   function ftheq(t,pk)                                                      
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: ftheq        compute saturation vapor pressure                    
!                                                                               
! abstract: compute equivalent potential temperature at the lcl                 
!   from temperature and pressure over 100 kpa to the kappa power.              
!   a biquadratic interpolation is done between values in a lookup table        
!   computed in funct_pot_temp_init. 
!   see documentation for fthex for details.                  
!   input values outside table range are reset to table extrema,                
!   except zero is returned for too cold or high lcls.                          
!   the interpolation accuracy is better than 0.0002 kelvin.                    
!   on the cray, ftheq is almost 3 times faster than exact calculation.         
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   1991-05-07  iredell                made into inlinable function
!   1994-12-30  iredell                exact computation
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:   the=ftheq(pv)                                                        
!                                                                               
!   input argument list:                                                        
!     t        - real lcl temperature in kelvin                                 
!     pk       - real lcl pressure over 100 kpa to the kappa power              
!                                                                               
!   output argument list:                                                       
!     ftheq    - real equivalent potential temperature in kelvin                
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
     jx=min(max(nint(xj),2),nx-1)                                            
     jy=min(max(nint(yj),2),ny-1)                                            
     dxj=xj-jx                                                               
     dyj=yj-jy                                                               
     ft11=tbthe(jx-1,jy-1)                                                   
     ft12=tbthe(jx-1,jy)                                                     
     ft13=tbthe(jx-1,jy+1)                                                   
     ft21=tbthe(jx,jy-1)                                                     
     ft22=tbthe(jx,jy)                                                       
     ft23=tbthe(jx,jy+1)                                                     
     ft31=tbthe(jx+1,jy-1)                                                   
     ft32=tbthe(jx+1,jy)                                                     
     ft33=tbthe(jx+1,jy+1)                                                   
     ftx1=(((ft31+ft11)/2-ft21)*dxj+(ft31-ft11)/2)*dxj+ft21                  
     ftx2=(((ft32+ft12)/2-ft22)*dxj+(ft32-ft12)/2)*dxj+ft22                  
     ftx3=(((ft33+ft13)/2-ft23)*dxj+(ft33-ft13)/2)*dxj+ft23                  
     ftheq=(((ftx3+ftx1)/2-ftx2)*dyj+(ftx3-ftx1)/2)*dyj+ftx2                 
   else                                                                      
     ftheq=0.                                                                
   endif                                                                     
!
   return                                                                    
   end function ftheq
!-------------------------------------------------------------------------------
