!
   function ftmaxg(tg,the,pk,qma)                                            
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: ftmaxg       compute moist adiabat temperature                    
!                                                                               
! abstract: exactly compute temperature and humidity of a parcel                
!   lifted up a moist adiabat from equivalent potential temperature             
!   at the lcl and pressure over 100 kpa to the kappa power.                    
!   a guess parcel temperature must be provided.                                
!   equivalent potential temperature is constant for a saturated parcel         
!   rising adiabatically up a moist adiabat when the heat and mass              
!   of the condensed water are neglected. the formula for                       
!   equivalent potential temperature (derived in holton) is                     
!       the=t*(pd**(-rocp))*exp(el*eps*pv/(cp*t*pd))                            
!   where t is the temperature, pv is the saturated vapor pressure,             
!   pd is the dry pressure p-pv, el is the temperature dependent                
!   latent heat of condensation hvap+dldt*(t-ttp), and other values             
!   are physical constants defined in parameter statements in the code.         
!   the formula is inverted by iterating newtonian approximations               
!   for each the and p until t is found to within 1.e-4 kelvin.                 
!   the specific humidity is then computed from pv and pd.                      
!   this function can be expanded inline in the calling routine.                
!                                                                               
! program history log:                                                          
!   1991-05-07  iredell                made into inlinable function
!   1994-12-30  iredell                exact computation
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:   tma=ftmaxg(tg,the,pk,qma)                                            
!                                                                               
!   input argument list:                                                        
!     tg       - real guess parcel temperature in kelvin                        
!     the      - real equivalent potential temperature in kelvin                
!     pk       - real pressure over 100 kpa to the kappa power                  
!                                                                               
!   output argument list:                                                       
!     ftmaxg   - real parcel temperature in kelvin                              
!     qma      - real parcel specific humidity in kg/kg                         
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
   real, parameter  ::  terrm=1.e-4                                                    
!
   t=tg                                                                      
   p=pk**cpor                                                                
   tr=ttp/t                                                                  
   pv=psatb*(tr**xa)*exp(xb*(1.-tr))                                         
   pd=p-pv                                                                   
   el=hvap+dldt*(t-ttp)                                                      
   expo=el*eps*pv/(cp*t*pd)                                                  
   thet=t*pd**(-rocp)*exp(expo)                                              
   dthet=thet/t*(1.+expo*(dldt*t/el+el*p/(rv*t*pd)))                         
   terr=(thet-the)/dthet                                                     
   t=t-terr                                                                  
!
   do while(abs(terr).gt.terrm)                                               
     tr=ttp/t                                                                
     pv=psatb*(tr**xa)*exp(xb*(1.-tr))                                       
     pd=p-pv                                                                 
     el=hvap+dldt*(t-ttp)                                                    
     expo=el*eps*pv/(cp*t*pd)                                                
     thet=t*pd**(-rocp)*exp(expo)                                            
     dthet=thet/t*(1.+expo*(dldt*t/el+el*p/(rv*t*pd)))                       
     terr=(thet-the)/dthet                                                   
     t=t-terr                                                                
   enddo                                                                     
!
   ftmaxg=t                                                                  
   tr=ttp/t                                                                  
   pv=psatb*(tr**xa)*exp(xb*(1.-tr))                                         
   pd=p-pv                                                                   
   qma=eps*pv/(pd+eps*pv)                                                    
!
   return                                                                    
   end function ftmaxg
!-------------------------------------------------------------------------------
