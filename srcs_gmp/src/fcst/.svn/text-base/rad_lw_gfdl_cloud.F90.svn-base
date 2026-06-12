#include <define.h>
   subroutine rad_lw_gfdl_cloud(ipts,cldfac,camt,nclds,kbtm,ktop)           
!-------------------------------------------------------------------------------
!                                                                               
! subroutine: rad_cloud_gfdl
!
! abstract:
!     subroutine clo88 computes cloud transmission functions for the            
!  longwave code,using code written by bert katz (301-763-8161).                
!  and modified by dan schwarzkopf in december,1988.                            
!
! program history log:
!   1988-05-06  kenneth campana        development
!   2000-01-01  song-you hong          physcis options
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel
   use hcon
   use rdparm
!-------------------------------------------------------------------------------
   integer              ::  nclds(imax),ktop(imbx,lp1),kbtm(imbx,lp1)                       
   real                 ::  camt(imbx,lp1),cldfac(imbx,lp1,lp1)                             
   real                 ::  cldrow(lp1)                                                     
!  dimension cldfip(lp1,lp1)                                                
   real                ::  cldipt(lp1,lp1,NVECT)                                         
!-------------------------------------------------------------------------------
!
   do iq = 1,ipts,2
     itop=iq+(2-1)
     if(itop.gt.ipts) itop=ipts
     jtop=itop-iq+1
!
! hoon
!
     do ip = 1,jtop                                                           
!
       ir=iq+ip-1                                                                
       if (nclds(ir).eq.0) then                                                  
         do j = 1,lp1                                                           
           do i = 1,lp1                                                           
             cldipt(i,j,ip)=1.                                                       
           enddo
         enddo
       endif                                                                     
!
       if (nclds(ir).ge.1) then                                                  
         xcld=1.-camt(ir,2)                                                    
         k1=ktop(ir,2)+1                                                      
         k2=kbtm(ir,2)                                                        
         do j = 1,lp1                                                         
           cldrow(j)=1.                                                      
         enddo
         do j = 1,k2                                                          
           cldrow(j)=xcld                                                    
         enddo
         kb=max(k1,k2+1)                                                       
         do k = kb,lp1                                                        
           do kp = 1,lp1                                                        
             cldipt(kp,k,ip)=cldrow(kp)                                       
           enddo
         enddo
         do j = 1,lp1                                                         
           cldrow(j)=1.                                                      
         enddo
         do j = k1,lp1                                                        
           cldrow(j)=xcld                                                    
         enddo
         kt=min(k1-1,k2)                                                       
         do k = 1,kt                                                          
           do kp = 1,lp1                                                        
             cldipt(kp,k,ip)=cldrow(kp)                                        
           enddo
         enddo
!
         if(k2+1.le.k1-1) then                                                 
           do j = k2+1,k1-1                                                   
             do i = 1,lp1                                                       
               cldipt(i,j,ip)=1.                                               
             enddo
           enddo
         else if(k1.le.k2) then                                                
           do j = k1,k2                                                       
             do i = 1,lp1                                                       
               cldipt(i,j,ip)=xcld                                             
             enddo
           enddo
         endif                                                                 
!
       endif ! nclds(ir).ge.1                                                                    
!
       if (nclds(ir).ge.2) then                                                  
!
         do nc = 2,nclds(ir)                                                    
           xcld=1.-camt(ir,nc+1)                                                 
           k1=ktop(ir,nc+1)+1                                                   
           k2=kbtm(ir,nc+1)                                                     
           do j = 1,lp1                                                         
             cldrow(j)=1.                                                      
           enddo
           do j = 1,k2                                                          
             cldrow(j)=xcld                                                    
           enddo
           kb=max(k1,k2+1)                                                       
           do k = kb,lp1                                                        
             do kp = 1,lp1                                                        
               cldipt(kp,k,ip)=cldipt(kp,k,ip)*cldrow(kp)                       
!              cldfip(kp,k)=cldrow(kp)
             enddo
           enddo
           do j = 1,lp1                                                         
             cldrow(j)=1.                                                      
           enddo
           do j = k1,lp1                                                        
             cldrow(j)=xcld                                                    
           enddo
           kt=min(k1-1,k2)                                                       
           do k = 1,kt                                                          
             do kp = 1,lp1                                                        
             cldipt(kp,k,ip)=cldipt(kp,k,ip)*cldrow(kp)                        
!             cldfip(kp,k)=cldrow(kp)                                           
             enddo
           enddo
!
!         if(k2+1.le.k1-1) then                                                 
!           do 51 j=k2+1,k1-1                                                   
!           do 51 i=1,lp1                                                       
!               cldipt(i,j,ip)=1.                                               
!51         continue                                                            
!
           if(k1.le.k2) then                                                     
             do j = k1,k2                                                       
               do i = 1,lp1                                                       
                 cldipt(i,j,ip)=cldipt(i,j,ip)*xcld                              
               enddo
             enddo
           endif                                                                 
!
!           do 65 j=1,lp1                                                       
!           do 65 i=1,lp1                                                       
!         cldipt(i,j,ip)=cldipt(i,j,ip)*cldfip(i,j)                             
!65       continue                                                              
         enddo
!
       endif ! nclds(ir).ge.2
!
     enddo ! ip
!
     do j = 1,lp1                                                             
       do i = 1,lp1                                                             
         do ip = 1,jtop                                                           
           ir=iq+ip-1                                                                
           cldfac(ir,i,j)=cldipt(i,j,ip)                                             
         enddo
       enddo
     enddo
!
   enddo ! iq
!
   return                                                                    
   end subroutine rad_lw_gfdl_cloud
