#include <define.h>
   subroutine co2_sigma_interp(pstar,pd,gtemp,t41,t42,t43,t44,sglvnu,siglnu,   &
                               lread)     
!-------------------------------------------------------------------------------
   use paramodel, only : kd=>levs_,kp=>levp1_,km=>levm1_,kp2=>levp2_
!-------------------------------------------------------------------------------
   real     ::  q(kd),qmh(kp),pd(kp2),plm(kp),gtemp(kp),pdt(kp2)                
   real     ::  ci(kp),sglvnu(kp),del(kd),siglnu(kd),cl(kd),rpi(km)             
   real     ::  t41(kp2,2),t42(kp),                                            &
                t43(kp2,2),t44(kp)                                              
   integer  ::  idate(4)                                                        
!
!     for sigma models,q=sigma,qmh=0.5(q(i)+q(i+1),                             
!     pd=q*pss,plm=qmh*pss.pss=surface pressure(spec.)                          
!                                                                               
!.....   get sigma structure                                                
!
   call co2_new_sigma(ci,sglvnu,del,siglnu,cl,rpi)                                  
!
   do k = 1,kd                                                              
     q(k) = siglnu(kd+1-k)                                                  
   enddo
!
   pss=    1013250.                                                          
   qmh(1)=0.                                                                 
   qmh(kp)=1.                                                                
!
   do k = 2,kd                                                               
     qmh(k)=0.5*(q(k-1)+q(k))                                                  
   enddo
!
   pd(1)=0.                                                                  
   pd(kp2)=pss                                                               
   do k = 2,kp                                                               
     pd(k)=q(k-1)*pss                                                          
   enddo
!
   plm(1)=0.                                                                 
   do k = 1,km                                                               
     plm(k+1)=0.5*(pd(k+1)+pd(k+2))                                            
   enddo
!
   plm(kp)=pss                                                               
   do k = 1,kd                                                               
     gtemp(k)=pd(k+1)**0.2*(1.+pd(k+1)/30000.)**0.8/1013250.                   
   enddo 
!
   gtemp(kp)=0.                                                              
   write (6,100) (gtemp(k),k=1,kd)                                           
   write (6,100) (pd(k),k=1,kp2)                                             
   write (6,100) (plm(k),k=1,kp)                                             
!
!***tapes 41,42 are output to the co2 interpolation program (ps=1013mb)         
!  the following puts p-data into mb                                            
!
   do i = 1,kp                                                              
     pd(i)=pd(i)*1.0e-3                                                        
     plm(i)=plm(i)*1.0e-3                                                      
   enddo
!
   pd(kp2)=pd(kp2)*1.0e-3                                                    
!cc         write (41,101) (pd(k),k=1,kp2)                                      
!cc         write (41,101) (plm(k),k=1,kp)                                      
!cc         write (42,101) (plm(k),k=1,kp)                                      
   do k = 1,kp2                                                            
     t41(k,1) = pd(k)                                                         
   enddo
!
   do k = 1,kp                                                             
     t41(k,2) = plm(k)                                                        
     t42(k) = plm(k)                                                          
   enddo
!
!***store as pdt,so that right pd is returned to ptz                            
!
   do i = 1,kp2                                                             
     pdt(i)=pd(i)                                                              
   enddo
!
!***second pass: pss=810mb,gtemp not computed                                   
!
   pss=0.8*1013250.                                                          
   qmh(1)=0.                                                                 
   qmh(kp)=1.                                                                
   do k = 2,kd                                                             
     qmh(k)=0.5*(q(k-1)+q(k))                                                  
   enddo
!
   pd(1)=0.                                                                  
   pd(kp2)=pss                                                               
   do k = 2,kp                                                             
     pd(k)=q(k-1)*pss                                                          
   enddo
!
   plm(1)=0.                                                                 
   do k = 1,km                                                             
     plm(k+1)=0.5*(pd(k+1)+pd(k+2))                                            
   enddo
!
   plm(kp)=pss                                                               
   write (6,100) (pd(k),k=1,kp2)                                             
   write (6,100) (plm(k),k=1,kp)                                             
!
!***tapes 43,44 are output to the co2 interpolation program(ps=810 mb)          
!  the following puts p-data into mb                                            
!
   do i = 1,kp                                                             
     pd(i)=pd(i)*1.0e-3                                                        
     plm(i)=plm(i)*1.0e-3                                                      
   enddo
!
   pd(kp2)=pd(kp2)*1.0e-3                                                    
!cc       write (43,101) (pd(k),k=1,kp2)                                        
!cc       write (43,101) (plm(k),k=1,kp)                                        
!cc       write (44,101) (plm(k),k=1,kp)                                        
   do k = 1,kp2                                                            
     t43(k,1) = pd(k)                                                         
   enddo
!
   do k = 1,kp                                                             
     t43(k,2) = plm(k)                                                        
     t44(k) = plm(k)                                                          
   enddo
!
!***restore pd                                                                  
!
   do i = 1,kp2                                                            
     pd(i)=pdt(i)                                                              
   enddo
!
   100   format (1x,5e20.13)                                                       
   101   format (5e16.9)                                                           
!
   return                                                                    
   end subroutine co2_sigma_interp
!-------------------------------------------------------------------------------
