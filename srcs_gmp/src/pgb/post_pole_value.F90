#include <define.h>
   subroutine post_pole_value
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [post_pole_value]
!           |
!           |-- [post_pole_vector] *
!           |-- [post_pole_extrap] *
!
!-------------------------------------------------------------------------------
   end subroutine post_pole_value
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_pole_vector(iptv,mpf)                    
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: post_pole_vector         sets default pole vector flags  
!                                                                               
! abstract: sets field identifier defaults for various parameters.              
!   a flag of 0 means scalar, 1 means vector, and 2 means flag.                 
!   these identifiers are used in interpolation.                                
!                                                                               
! program history log:                                                          
!   1993-10-21  iredell                                                           
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!                                                                               
! usage:    call post_pole_vector(iptv,mpf)               
!   input arguments:                                                            
!     iptv         paramter table version (only 1 or 2 is recognized)           
!   output arguments:                                                           
!     mpf          integer (255) field parameter identifiers                    
!                                                                               
!-------------------------------------------------------------------------------
   dimension mpf(0:255)                                                        
!-------------------------------------------------------------------------------
!
   do i = 1,255                                                                
     mpf(i)=0                                                                  
   enddo                                                                     
!
   if(iptv.eq.1.or.iptv.eq.2) then                                           
     do i = 33,34                                                              
!       mpf(033:034)=1                                                          
       mpf(i)=1                                                                
     enddo                                                                   
     do i = 49,50                                                              
!       mpf(049:050)=1                                                          
       mpf(i)=1                                                                
     enddo                                                                   
     do i = 95,96                                                              
!       mpf(095:096)=1                                                          
       mpf(i)=1                                                                
     enddo                                                                   
     do i = 124,125                                                            
!       mpf(124:125)=1                                                          
       mpf(i)=1                                                                
     enddo                                                                   
     do i = 181,182                                                            
!       mpf(181:182)=1                                                          
       mpf(i)=1                                                                
     enddo                                                                   
     do i = 183,184                                                            
!       mpf(183:184)=1                                                          
       mpf(i)=1                                                                
     enddo                                                                   
     do i = 247,248                                                            
!       mpf(247:248)=1                                                          
       mpf(i)=1                                                                
     enddo                                                                   
     mpf(081)=2                                                              
     mpf(091)=2                                                              
     mpf(140)=2                                                              
     mpf(141)=2                                                              
     mpf(142)=2                                                              
     mpf(143)=2                                                              
     mpf(173)=2                                                              
     mpf(174)=2                                                              
     mpf(175)=2                                                              
     mpf(209)=2                                                              
   endif                                                                     
! 
   return                                                                    
   end subroutine post_pole_vector                                                                     
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_pole_extrap(mp,im,fnx,fsx,fn,fs)          
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:    post_pole_extrap      extrapolate a field to the poles. 
!                                                                               
! abstract: a global horizontal field is exterpolated to the poles.             
!   polar scalars are the average of the closest latitude circle values.        
!   polar vector components are taken from the wavenumber 1 component           
!   extracted from the values on the closest latitude circle.                   
!   polar flags are copied from the closest prime meridian value.               
!                                                                               
! program history log:                                                          
!    1992-10-31  iredell                development
!    2000-03-09  songyou hong           cvs verion setup
!    2009-10-01  jung-eun kim           f90 format with standard physics modules
!    2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call post_pole_extrap(mp,im,fnx,fsx,fn,fs)               
!   input argument list:                                                        
!     mp       - integer field parameter identifier                             
!                (0 for scalar, 1 for vector, 2 for flag)                       
!     im       - integer number of longitudes                                   
!     fnx      - real (im) field values on the closest latitude circle          
!                to the north pole                                              
!     fsx      - real (im) field values on the closest latitude circle          
!                to the south pole                                              
!                                                                               
!   output argument list:                                                       
!     fn       - real (im) field values extrapolated to the north pole          
!     fs       - real (im) field values extrapolated to the south pole          
!                                                                               
!-------------------------------------------------------------------------------
   real fnx(im),fsx(im),fn(im),fs(im)                                        
!-------------------------------------------------------------------------------
! 
!  get polar values for scalars or vectors                                      
!
   pi=acos(-1.)                                                              
   if(mp.eq.0) then                                                          
!
!  full scalar                                                                  
!
     fnp=0.                                                                  
     fsp=0.                                                                  
     do i = 1,im                                                          
       fnp=fnp+fnx(i)                                                        
       fsp=fsp+fsx(i)                                                        
     enddo    
     fnp=fnp/im
     fsp=fsp/im
     do i = 1,im                                                          
       fn(i)=fnp                                                             
       fs(i)=fsp                                                             
     enddo 
   elseif(mp.eq.1) then                                                      
!
!  full vector                                                                  
!
     fnpc=0.                                                                 
     fnps=0.                                                                 
     fspc=0.                                                                 
     fsps=0.                                                                 
     do i = 1,im                                                          
       ci=cos(2*pi*(i-1)/im)                                                 
       si=sin(2*pi*(i-1)/im)                                                 
       fnpc=fnpc+ci*fnx(i)                                                   
       fnps=fnps+si*fnx(i)                                                   
       fspc=fspc+ci*fsx(i)                                                   
       fsps=fsps+si*fsx(i)                                                   
     enddo 
     fnpc=2*fnpc/im                                                          
     fnps=2*fnps/im                                                          
     fspc=2*fspc/im                                                          
     fsps=2*fsps/im                                                          
     do i = 1,im                                                          
       ci=cos(2*pi*(i-1)/im)                                                 
       si=sin(2*pi*(i-1)/im)                                                 
       fn(i)=fnpc*ci+fnps*si                                                 
       fs(i)=fspc*ci+fsps*si                                                 
     enddo 
   elseif(mp.eq.2) then                                                      
!
!  full flag                                                                    
!
     do i = 1,im                                                          
       fn(i)=fnx(1)                                                          
       fs(i)=fsx(1)                                                          
     enddo 
   endif                                                                     
!
   return                                                                    
   end subroutine post_pole_extrap                                                                      
!-------------------------------------------------------------------------------
