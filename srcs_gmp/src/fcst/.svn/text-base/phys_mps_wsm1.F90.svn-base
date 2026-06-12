#include <define.h>
   subroutine phys_mps_wsm1(dt,t,q,prsl,delprsi,rn,lat,                        &
                       cp,g,rd,rv,hvap,                                        &
#ifdef SAS_DIAG
                       dlt,dlq,dlh,                                            &
#endif
#ifdef RASV2
                       clw,                                                    &
#endif
                       ids,ide, jds,jde, kds,kde,                              &
                       ims,ime, jms,jme, kms,kme,                              &
                       its,ite, jts,jte, kts,kte)
!-------------------------------------------------------------------------------
!
! subroutine:    phys_mps_wsm1     
!                                                                               
! abstract: calculates grid-scale condensation for one leap-frog                
!   timestep, produces rain, and adjusts temperature and specific               
!   humidity by wet-bulb process.  evaporation of part or all of the            
!   rain may occur as it traverses unsaturated layers on the way down.          
!   for conditionally unstable layers, a convective adjustment procedure        
!   is applied to adjust to a uniform theta-e.                                  
!                                                                               
! program history log:
!   1994-04-15  hua-lu pan                                                        
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call phys_mps_wsm1(im,im2,km,dt,ps,t,q,sl,del,slk,rain,lat)                
!                                                                               
!   input argument list:                                                        
!     im       - integer number of points                                       
!     im2      - real first dimension of t and q                                
!     km       - integer number of levels                                       
!     dt       - real time step in seconds                                      
!     ps       - real (im) surface pressure in kilopascals (cb)                 
!     t        - real (im2,km) current temperature in kelvin                    
!     q        - real (im2,km) current specific humidity in kg/kg               
!     sl       - real (km) sigma values                                         
!     del      - real (km) sigma layer thickness                                
!     slk      - real (km) sigma values ** kappa                                
!     lat      - integer latitude number                                        
!                                                                               
!   output argument list:                                                       
!     q        - real (im2,km) adjusted specific humidity in kg/kg              
!     t        - real (im2,km) adjusted temperature in kelvin                   
!     rn       - real (im) large-scale rain in meters                           
!                                                                               
! subprograms called:                                                           
!   fpvs     - function to compute saturation vapor pressure                    
!   fpkap    - function to compute p raised to the factor kappa                 
!   fthe     - function to compute theta-e                                      
!   ftma     - function to compute temperature and moisture along a             
!              moist adiabat                                                    
!                                                                               
! remarks: the precipitation reaching the ground should be halved               
!   before it is used to increment raintot, the running total.  this             
!   prevents double-counting of raintot, which is incremented every              
!   half leap-frog timestep.                                                    
!          the evaporation rate for falling precip is calculated                
!   according to a method devised by e.kessler, in which a mean drop            
!   surface area is obtained from the rainwater content.                        
!          ice is not considered.                                               
!          functions fpvs,fpkap,fthe,ftma are inlined by fpp.                   
!                                                                               
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  lat,                                               &
                            ids,ide, jds,jde, kds,kde,                         &
                            ims,ime, jms,jme, kms,kme,                         &
                            its,ite, jts,jte, kts,kte
   real                 ::  g,rd,rv,cp,hvap,dt
   real                 ::  delprsi(ims:ime,kms:kme),prsl(ims:ime,kms:kme)
   real                 ::  q(ims:ime,kms:kme),t(ims:ime,kms:kme),rn(ims:ime)
#ifdef RASV2
   real                 ::  clw(ims:ime,kms:kme)
#endif
#ifdef SAS_DIAG
   real                 ::  dlt(ims:ime,kms:kme), dlq(ims:ime,kms:kme),        &
                            dlh(ims:ime,kms:kme)
#endif
!
   integer              ::  i,k
   real                 ::  elocp,el2orc,eps,epsm1
   real                 ::  qs,es,dpovg,qevap,qcond,rnevap,fpvs,fpvs0
!
! constants
!
   elocp  = hvap/cp
   el2orc = hvap*hvap/(rv*cp)
   eps    = rd/rv
   epsm1  = rd/rv-1.
!
#ifdef SAS_DIAG
   do k = kte,kts,-1                                                             
     do i = its,ite
       dlt(i,k) = 0.0
       dlq(i,k) = 0.0
       dlh(i,k) = 0.0
     enddo
   enddo
#endif
!
!  condense rain to wetbulb temperature if supersaturated                       
!  or evaporate rain using kessler parameterization.                            
!
   do i = its,ite
     rn(i)=0.                                                                
   enddo
!
   do k = kte,kts,-1                                                             
     do i = its,ite
       dpovg=delprsi(i,k)/g
#ifdef ICE
       es=fpvs(t(i,k))                                                       
#else
       es=fpvs0(t(i,k))                                                       
#endif
       qs=eps*es/(prsl(i,k)+epsm1*es)                                      
       qcond=(q(i,k)-qs)/(1.+el2orc*qs/t(i,k)**2)                            
#ifdef RASV2
       rn(i)=rn(i)+clw(i,k)
#endif
       if(qcond.gt.0.) then                                                  
         q(i,k)=q(i,k)-qcond                                                 
         t(i,k)=t(i,k)+qcond*elocp                                           
         rn(i)=rn(i)+qcond*dpovg                                             
#ifdef SAS_DIAG
         dlq(i,k) = dlq(i,k) - qcond
         dlt(i,k) = dlt(i,k) + qcond*elocp
#endif
         elseif(rn(i).gt.0.) then                                              
         qevap=-qcond*(1.-exp(-0.32*sqrt(2.*dt*rn(i))))                      
         rnevap=min(qevap*dpovg,rn(i))                                       
         q(i,k)=q(i,k)+rnevap/dpovg                                          
         t(i,k)=t(i,k)-rnevap/dpovg*elocp                                    
         rn(i)=rn(i)-rnevap                                                  
#ifdef SAS_DIAG
         dlq(i,k) = dlq(i,k) + rnevap/dpovg
         dlt(i,k) = dlt(i,k) - rnevap/dpovg*elocp
#endif
       endif                                                                 
#ifdef SAS_DIAG
       dlh(i,k) = cp*dlt(i,k) + hvap*dlq(i,k)
#endif
     enddo                                                                   
   enddo                                                                     
!
   return                                                                    
   end subroutine phys_mps_wsm1
!-------------------------------------------------------------------------------
