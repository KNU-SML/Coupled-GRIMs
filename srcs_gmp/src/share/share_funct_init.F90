#include "define.h"
   subroutine share_funct_init
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [share_funct_init]
!      |
!      |--- [function: ftdpx] *           : initializing
!      |--- [funct_dew_point_temp_init] * : making tables
!      |
!      |--- [function: ftmax] *
!      |--- [funct_moist_adiabat_init] *
!      |
!      |--- [function: fthex] *
!      |--- [funct_pot_temp_init] *
!      |
!      |--- [function: fpvsx] *   <--- ice effect
!      |--- [function: fpvsx0] *  <--- no ice effect
!      |--- [funct_svp_init] *
!      |
!      |--- [function: fpkapx] *   <--- ice effect
!
!-------------------------------------------------------------------------------
   end subroutine share_funct_init
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function ftdpx(pv)                                                        
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: ftdpx        compute saturation vapor pressure                    
!                                                                               
! abstract: exactly compute dewpoint temperature from vapor pressure.           
!   an approximate dewpoint temperature for function ftdpxg                     
!   is obtained using ftdp so funct_dew_point_temp_init must be already called.
!   see documentation for ftdpxg for details.                                   
!                                                                               
! program history log:                                                          
!   1991-05-07  iredell                made into inlinable function                  
!   1994-12-30  iredell                exact computation                             
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:   tdp=ftdpx(pv)                                                        
!                                                                               
!   input argument list:                                                        
!     pv       - real vapor pressure in kilopascals (cb)                        
!                                                                               
!   output argument list:                                                       
!     ftdpx    - real dewpoint temperature in kelvin                            
!                                                                               
! subprograms called:                                                           
!   (ftdp)   - inlinable function to compute dewpoint temperature               
!   (ftdpxg) - inlinable function to compute dewpoint temperature               
!                                                                               
!-------------------------------------------------------------------------------
   tg=ftdp(pv)                                                               
   ftdpx=ftdpxg(tg,pv)                                                       
! 
   return                                                                    
   end function ftdpx                                                                       
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine funct_dew_point_temp_init   
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: funct_dew_point_temp_init   compute dewpoint temperature table 
!                                                                               
! abstract: compute dewpoint temperature table as a function of                 
!   vapor pressure for inlinable function ftdp.                                 
!   exact dewpoint temperatures are calculated in subprogram ftdpxg.            
!   the current implementation computes a table with a length                   
!   of 5001 for vapor pressures ranging from 0.001 to 10.001 kilopascals        
!   giving a dewpoint temperature range of 208.0 to 319.0 kelvin.               
!                                                                               
! program history log:                                                          
!   91-05-07  iredell                                                           
!   94-12-30  iredell             expand table                                  
!                                                                               
! usage:  call funct_dew_point_temp_init                     
!                                                                               
! subprograms called:                                                           
!   (ftdpxg) - inlinable function to compute dewpoint temperature               
!                                                                               
! common blocks:                                                                
!   comtdp   - scaling parameters and table for function ftdp.                  
!                                                                               
!-------------------------------------------------------------------------------
   integer, parameter  ::  nx=5001                                             
   real                ::  tbtdp(nx)                                           
   common/comtdp/ c1xtdp,c2xtdp,tbtdp                                        
!
   xmin= 0.001                                                               
   xmax=10.001                                                               
   xinc=(xmax-xmin)/(nx-1)                                                   
   c1xtdp=1.-xmin/xinc                                                       
   c2xtdp=1./xinc                                                            
   t=208.0                                                                   
!
   do jx = 1,nx                                                                
     x=xmin+(jx-1)*xinc                                                      
     pv=x                                                                    
     t=ftdpxg(t,pv)                                                          
     tbtdp(jx)=t                                                             
   enddo                                                                     
!
   return                                                                    
   end subroutine funct_dew_point_temp_init   
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function ftmax(the,pk,qma)                                                
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: ftmax        compute moist adiabat temperature                    
!                                                                               
! abstract: exactly compute temperature and humidity of a parcel                
!   lifted up a moist adiabat from equivalent potential temperature             
!   at the lcl and pressure over 100 kpa to the kappa power.                    
!   an approximate parcel temperature for function ftmaxg                       
!   is obtained using ftma so funct_moist_adiabat_init must be already called. 
!   see documentation for ftmaxg for details.                                   
!                                                                               
! program history log:                                                          
!   91-05-07  iredell             made into inlinable function                  
!   94-12-30  iredell             exact computation                             
!                                                                               
! usage:   tma=ftmax(the,pk,qma)                                                
!                                                                               
!   input argument list:                                                        
!     the      - real equivalent potential temperature in kelvin                
!     pk       - real pressure over 100 kpa to the kappa power                  
!                                                                               
!   output argument list:                                                       
!     ftmax    - real parcel temperature in kelvin                              
!     qma      - real parcel specific humidity in kg/kg                         
!                                                                               
! subprograms called:                                                           
!   (ftma)   - inlinable function to compute parcel temperature                 
!   (ftmaxg) - inlinable function to compute parcel temperature                 
!                                                                               
!-------------------------------------------------------------------------------
   tg=ftma(the,pk,qg)                                                        
   ftmax=ftmaxg(tg,the,pk,qma)                                               
! 
   return                                                                    
   end function ftmax                                                                       
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine funct_moist_adiabat_init           
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: funct_moist_adiabat_init  compute moist adiabat tables  
!                                                                               
! abstract: compute temperature and specific humidity tables                    
!   as a function of equivalent potential temperature and                       
!   pressure over 100 kpa to the kappa power for function ftma.                 
!   exact parcel temperatures are calculated in subprogram ftmaxg.              
!   the current implementation computes a table with a first dimension          
!   of 151 for equivalent potential temperatures ranging from 200 to 500        
!   kelvin and a second dimension of 121 for pressure over 100 kpa              
!   to the kappa power ranging from 0.01**rocp to 1.10**rocp.                   
!                                                                               
! program history log:                                                          
!   91-05-07  iredell                                                           
!   94-12-30  iredell             expand table                                  
!                                                                               
! usage:  call funct_moist_adiabat_init              
!                                                                               
! subprograms called:                                                           
!   (ftmaxg) - inlinable function to compute parcel temperature                 
!                                                                               
! common blocks:                                                                
!   comma    - scaling parameters and table for function ftma.                  
!                                                                               
!-------------------------------------------------------------------------------
   use constant, only : cp_,rd_
!-------------------------------------------------------------------------------
   integer, parameter   ::  nx=151,ny=121                                                  
   real, parameter      ::  cp=cp_,rd=rd_                                                  
   real, parameter      ::  rocp=rd/cp                                                     
   real                 ::  tbtma(nx,ny),tbqma(nx,ny)                                       
   common/comma/ c1xma,c2xma,c1yma,c2yma,tbtma,tbqma                         
!
   xmin=200.                                                                 
   xmax=500.                                                                 
   xinc=(xmax-xmin)/(nx-1)                                                   
   c1xma=1.-xmin/xinc                                                        
   c2xma=1./xinc                                                             
   ymin=0.01**rocp                                                           
   ymax=1.10**rocp                                                           
   yinc=(ymax-ymin)/(ny-1)                                                   
   c1yma=1.-ymin/yinc                                                        
   c2yma=1./yinc                                                             
!
   do jy = 1,ny                                                                
     y=ymin+(jy-1)*yinc                                                      
     pk=y                                                                    
     t=xmin*y                                                                
     do jx = 1,nx                                                              
       x=xmin+(jx-1)*xinc                                                    
       the=x                                                                 
       t=ftmaxg(t,the,pk,q)                                                  
       tbtma(jx,jy)=t                                                        
       tbqma(jx,jy)=q                                                        
     enddo                                                                   
   enddo                                                                     
!
   return                                                                    
   end subroutine funct_moist_adiabat_init           
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function fthex(t,pk)                                                      
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: fthex        compute saturation vapor pressure                    
!                                                                               
! abstract: exactly compute equivalent potential temperature at the lcl         
!   from temperature and pressure over 100 kpa to the kappa power.              
!   equivalent potential temperature is constant for a saturated parcel         
!   rising adiabatically up a moist adiabat when the heat and mass              
!   of the condensed water are neglected. the formula for                       
!   equivalent potential temperature (derived in holton) is                     
!       the=t*(pd**(-rocp))*exp(el*eps*pv/(cp*t*pd))                            
!   where t is the temperature, pv is the saturated vapor pressure,             
!   pd is the dry pressure p-pv, el is the temperature dependent                
!   latent heat of condensation hvap+dldt*(t-ttp), and other values             
!   are physical constants defined in parameter statements in the code.         
!   zero is returned if the input values make saturation impossible.            
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   91-05-07  iredell             made into inlinable function                  
!   94-12-30  iredell             exact computation                             
!                                                                               
! usage:   the=fthex(t,pk)                                                      
!                                                                               
!   input argument list:                                                        
!     t        - real lcl temperature in kelvin                                 
!     pk       - real lcl pressure over 100 kpa to the kappa power              
!                                                                               
!   output argument list:                                                       
!     fthex    - real equivalent potential temperature in kelvin                
!                                                                               
!-------------------------------------------------------------------------------
   use constant, only : cliq_,cp_,cvap_,hvap_,psat_,rd_,rv_,ttp_
!-------------------------------------------------------------------------------
   real, parameter  ::  cp=cp_,rd=rd_,rv=rv_,                                  &
                        ttp=ttp_,hvap=hvap_,psat=psat_,                        &
                        cliq=cliq_,cvap=cvap_          
   real, parameter  ::  psatk=psat*1.e-3                                               
   real, parameter  ::  rocp=rd/cp,cpor=cp/rd,psatb=psatk*1.e-2,eps=rd/rv,     &
                        dldt=cvap-cliq,xa=-dldt/rv,xb=xa+hvap/(rv*ttp)                 
!
   p=pk**cpor                                                                
   tr=ttp/t                                                                  
   pv=psatb*(tr**xa)*exp(xb*(1.-tr))                                         
   pd=p-pv                                                                   
!
   if(pd.gt.0.) then                                                         
     el=hvap+dldt*(t-ttp)                                                    
!    expo=el*eps*pv/(cp*t*pd)                                                
!    fthex=t*pd**(-rocp)*exp(expo)                                           
     expo=el*eps*pv/(cp*t*pd)                                                
     expo = min(expo,100.0)                                                  
     fthex=t*pd**(-rocp)*exp(expo)                                           
   else                                                                      
     fthex=0.                                                                
   endif                                                                     
!
   return                                                                    
   end function fthex                                                                      
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine funct_pot_temp_init   
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: funct_pot_temp_init compute equivalent potential temperature table
!                                                                               
! abstract: compute equivalent potential temperature table                      
!   as a function of lcl temperature and pressure over 100 kpa                  
!   to the kappa power for function fthe.                                       
!   equivalent potential temperatures are calculated in subprogram fthex        
!   the current implementation computes a table with a first dimension          
!   of 241 for temperatures ranging from 183.16 to 303.16 kelvin                
!   and a second dimension of 151 for pressure over 100 kpa                     
!   to the kappa power ranging from 0.04**rocp to 1.10**rocp.                   
!                                                                               
! program history log:                                                          
!   91-05-07  iredell                                                           
!   94-12-30  iredell             expand table                                  
!                                                                               
! usage:  call funct_pot_temp_init                 
!                                                                               
! subprograms called:                                                           
!   (fthex)  - inlinable function to compute equiv. pot. temperature            
!                                                                               
! common blocks:                                                                
!   comthe   - scaling parameters and table for function fthe.                  
!                                                                               
!-------------------------------------------------------------------------------
   use constant, only : cp_,rd_,ttp_
!-------------------------------------------------------------------------------
   real, parameter    ::  cp=cp_,rd=rd_,ttp=ttp_                                         
   real, parameter    ::  rocp=rd/cp                                                     
   integer, parameter ::  nx=241,ny=151                                                  
   real               ::  tbthe(nx,ny)                                                    
   common/comthe/ c1xthe,c2xthe,c1ythe,c2ythe,tbthe                          
!
   xmin=ttp-90.                                                              
   xmax=ttp+30.                                                              
   xinc=(xmax-xmin)/(nx-1)                                                   
   c1xthe=1.-xmin/xinc                                                       
   c2xthe=1./xinc                                                            
   ymin=0.04**rocp                                                           
   ymax=1.10**rocp                                                           
   yinc=(ymax-ymin)/(ny-1)                                                   
   c1ythe=1.-ymin/yinc                                                       
   c2ythe=1./yinc                                                            
!
   do jy = 1,ny                                                                
     y=ymin+(jy-1)*yinc                                                      
     pk=y                                                                    
     do jx = 1,nx                                                              
       x=xmin+(jx-1)*xinc                                                    
       t=x                                                                   
       tbthe(jx,jy)=fthex(t,pk)                                              
     enddo                                                                   
   enddo                                                                     
!
   return                                                                    
   end subroutine funct_pot_temp_init   
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function fpvsx(t)                                                         
!-------------------------------------------------------------------------------
   use constant, only : cliq_,cice_,cp_,cvap_,hvap_,hsub_,psat_,rd_,rv_,ttp_
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: fpvsx        compute saturation vapor pressure                    
!                                                                               
! abstract: exactly compute saturation vapor pressure from temperature.         
!   the water model assumes a perfect gas, constant specific heats              
!   for gas and liquid, and neglects the volume of the liquid.                  
!   the model does account for the variation of the latent heat                 
!   of condensation with temperature.  the ice option is not included.          
!   the clausius-clapeyron equation is integrated from the triple point         
!   to get the formula                                                          
!       pvs=psatk*(tr**xa)*exp(xb*(1.-tr))                                      
!   where tr is ttp/t and other values are physical constants                   
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   1991-05-07  iredell             made into inlinable function                  
!   1994-12-30  iredell             exact computation                             
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!                                                                               
! usage:   pvs=fpvsx(t)                                                         
!                                                                               
!   input argument list:                                                        
!     t        - real temperature in kelvin                                     
!                                                                               
!   output argument list:                                                       
!     fpvsx    - real saturation vapor pressure in kilopascals (cb)             
!                                                                               
!-------------------------------------------------------------------------------
   real, parameter   ::  cp=cp_,rd=rd_,rv=rv_,                                 &
                         ttp=ttp_,hvap=hvap_,psat=psat_,                       &
                         cliq=cliq_,cvap=cvap_,cice=cice_,hsub=hsub_
   real, parameter   ::  psatk=psat*1.e-3
   real, parameter   ::  dldt=cvap-cliq,xa=-dldt/rv,xb=xa+hvap/(rv*ttp)
   real, parameter   ::  dldti=cvap-cice,xai=-dldti/rv,xbi=xai+hsub/(rv*ttp)
!
   tr=ttp/t
   if(t.ge.ttp) then
     fpvsx=psatk*(tr**xa)*exp(xb*(1.-tr))
   else
     fpvsx=psatk*(tr**xai)*exp(xbi*(1.-tr))
   endif
!
   return
   end function fpvsx
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function fpvsx0(t)                                                         
!-------------------------------------------------------------------------------
   use constant, only : cp_,cvap_,cliq_,hvap_,psat_,rd_,rv_,ttp_
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: fpvsx        compute saturation vapor pressure                    
!                                                                               
! abstract: exactly compute saturation vapor pressure from temperature.         
!   the water model assumes a perfect gas, constant specific heats              
!   for gas and liquid, and neglects the volume of the liquid.                  
!   the model does account for the variation of the latent heat                 
!   of condensation with temperature.  the ice option is not included.          
!   the clausius-clapeyron equation is integrated from the triple point         
!   to get the formula                                                          
!       pvs=psatk*(tr**xa)*exp(xb*(1.-tr))                                      
!   where tr is ttp/t and other values are physical constants                   
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   91-05-07  iredell             made into inlinable function                  
!   94-12-30  iredell             exact computation                             
!                                                                               
! usage:   pvs=fpvsx(t)                                                         
!                                                                               
!   input argument list:                                                        
!     t        - real temperature in kelvin                                     
!                                                                               
!   output argument list:                                                       
!     fpvsx    - real saturation vapor pressure in kilopascals (cb)             
!                                                                               
!-------------------------------------------------------------------------------
   real, parameter   ::  cp=cp_,rd=rd_,rv=rv_,                                 &
                         ttp=ttp_,hvap=hvap_,psat=psat_,                       &
                         cliq=cliq_,cvap=cvap_
   real, parameter   ::  psatk=psat*1.e-3
   real, parameter   ::  dldt=cvap-cliq,xa=-dldt/rv,xb=xa+hvap/(rv*ttp)
!
   tr=ttp/t
   fpvsx0=psatk*(tr**xa)*exp(xb*(1.-tr))
!
   return
   end function fpvsx0
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine funct_svp_init                    
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: funct_svp_init    compute saturation vapor pressure table
!                                                                               
! abstract: compute saturation vapor pressure table as a function of            
!   temperature for the table lookup function fpvs.                             
!   exact saturation vapor pressures are calculated in subprogram fpvsx.        
!   the current implementation computes a table with a length                   
!   of 7501 for temperatures ranging from 180. to 330. kelvin.                  
!                                                                               
! program history log:                                                          
!   91-05-07  iredell                                                           
!   94-12-30  iredell             expand table                                  
!   96-02-19  hong                ice effect
!                                                                               
! usage:  call funct_svp_init                              
!                                                                               
! subprograms called:                                                           
!   (fpvsx)  - inlinable function to compute saturation vapor pressure          
!                                                                               
! common blocks:                                                                
!   compvs   - scaling parameters and table for function fpvs.                  
!                                                                               
!-------------------------------------------------------------------------------
   integer, parameter   ::  nx=7501                                                        
   real                 ::  tbpvs(nx),tbpvs0(nx)
   common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
   common/compvs/ c1xpvs,c2xpvs,tbpvs
!
   xmin=180.0
   xmax=330.0
   xinc=(xmax-xmin)/(nx-1)
   c1xpvs=1.-xmin/xinc
   c2xpvs=1./xinc
   c1xpvs0=1.-xmin/xinc
   c2xpvs0=1./xinc
!
   do jx = 1,nx
     x=xmin+(jx-1)*xinc
     t=x
     tbpvs(jx)=fpvsx(t)
     tbpvs0(jx)=fpvsx0(t)
   enddo
!
   return                                                                    
   end subroutine funct_svp_init                    
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function fpkapx(p)                                                        
!-------------------------------------------------------------------------------
   use constant, only : cp_,rd_
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: fpkapx       raise surface pressure to the kappa power.           
!                                                                               
! abstract: raise surface pressure over 100 kpa to the kappa power.             
!   kappa is equal to rd/cp where rd and cp are physical constants.             
!   this function should be expanded inline in the calling routine.             
!                                                                               
! program history log:                                                          
!   94-12-30  iredell             made into inlinable function                  
!                                                                               
! usage:  pkap=fpkapx(p)                                                        
!                                                                               
!   input argument list:                                                        
!     p        - real surface pressure in kilopascals (cb)                      
!                                                                               
!   output argument list:                                                       
!     fpkapx   - real p/100 to the kappa power                                  
!                                                                               
!-------------------------------------------------------------------------------
   real, parameter   ::  cp=cp_,rd=rd_                                                  
   real, parameter   ::  rocp=rd/cp                                                     
!
   fpkapx=(p/100.)**rocp                                                     
!
   return                                                                    
   end function fpkapx                                                                       
!-------------------------------------------------------------------------------
