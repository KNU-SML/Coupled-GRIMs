!
   function ftdpxg(tg,pv)                                                    
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: ftdpxg       compute saturation vapor pressure                    
!                                                                               
! abstract: exactly compute dewpoint temperature from vapor pressure.           
!   a guess dewpoint temperature must be provided.                              
!   the water model assumes a perfect gas, constant specific heats              
!   for gas and liquid, and neglects the volume of the liquid.                  
!   the model does account for the variation of the latent heat                 
!   of condensation with temperature.  the ice option is not included.          
!   the clausius-clapeyron equation is integrated from the triple point         
!   to get the formula                                                          
!       pvs=psatk*(tr**xa)*exp(xb*(1.-tr))                                      
!   where tr is ttp/t and other values are physical constants                   
!   the formula is inverted by iterating newtonian approximations               
!   for each pvs until t is found to within 1.e-6 kelvin.                       
!   this function can be expanded inline in the calling routine.                
!                                                                               
! program history log:                                                          
!   1982-12-30   n phillips            initial development
!   1991-05-07  iredell                made into inlinable function
!   1994-12-30  iredell                exact computation
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:   tdp=ftdpxg(tg,pv)                                                    
!                                                                               
!   input argument list:                                                        
!     tg       - real guess dewpoint temperature in kelvin                      
!     pv       - real vapor pressure in kilopascals (cb)                        
!                                                                               
!   output argument list:                                                       
!     ftdpxg   - real dewpoint temperature in kelvin                            
!                                                                               
!-------------------------------------------------------------------------------
   use constant, only : cliq_,cp_,cvap_,hvap_,psat_,rd_,rv_,ttp_
!-------------------------------------------------------------------------------
   real, parameter  ::  cp=cp_,rd=rd_,rv=rv_,                                  &
                        ttp=ttp_,hvap=hvap_,psat=psat_,                        &
                        cliq=cliq_,cvap=cvap_                                         
   real, parameter  ::  psatk=psat*1.e-3                                              
   real, parameter  ::  dldt=cvap-cliq,xa=-dldt/rv,xb=xa+hvap/(rv*ttp)                 
   real, parameter  ::  terrm=1.e-4                                                    
!
   t=tg                                                                      
   tr=ttp/t                                                                  
   pvt=psatk*(tr**xa)*exp(xb*(1.-tr))                                        
   el=hvap+dldt*(t-ttp)                                                      
   dpvt=el*pvt/(rv*t**2)                                                     
   terr=(pvt-pv)/dpvt                                                        
   t=t-terr                                                                  
!
   do while(abs(terr).gt.terrm)                                               
     tr=ttp/t                                                                
     pvt=psatk*(tr**xa)*exp(xb*(1.-tr))                                      
     el=hvap+dldt*(t-ttp)                                                    
     dpvt=el*pvt/(rv*t**2)                                                   
     terr=(pvt-pv)/dpvt                                                      
     t=t-terr                                                                
   enddo                                                                     
!
   ftdpxg=t                                                                  
!
   return                                                                    
   end function ftdpxg                                                                       
!-------------------------------------------------------------------------------
