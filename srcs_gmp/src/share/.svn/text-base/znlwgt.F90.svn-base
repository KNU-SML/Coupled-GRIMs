!
   subroutine znlwgt(slmsk,snoweq,wgt,lat,                                     &
                     znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)        
!-------------------------------------------------------------------------------
   logical lsmsk(idimt,6)                                                    
   real   ::  znlsl(6,6,nrcznl),weis(6,6)                                     
!                                                                               
   real   ::  slmsk(1),snoweq(1)                                              
!                                                                               
   jdimhf=jdim/2                                                             
!                                                                               
!  zonal average weight                                                         
!                                                                               
   nzl1=jdim/6+1                                                             
   nzl1p=nzl1+1                                                              
   nzl2=jdim/3+1                                                             
   nzl2p=nzl2+1                                                              
!                                                                               
   n=idim+1                                                                  
!                                                                               
   if(lat.eq.1) then                                                         
     do j = 1,6                                                           
       do i = 1,6                                                           
         weis(i,j) = 0.0                                                           
       enddo
     enddo
   endif                                                                     
!                                                                               
!  bare/snow sea/land/ice averages                                              
!                                                                               
   do i = 1,idimt                                                       
     lsmsk(i,2)=(slmsk (i).eq.1.).and.                                         &
                 (snoweq(i).le.1.e-3)                                           
     lsmsk(i,3)=(slmsk (i).eq.1.).and.                                         &
                 (snoweq(i).gt.1.e-3)                                           
     lsmsk(i,4)=(slmsk (i).eq.2.).and.                                         &
                 (snoweq(i).le.1.e-3)                                           
     lsmsk(i,5)=(slmsk (i).eq.2.).and.                                         &
                 (snoweq(i).gt.1.e-3)                                           
     lsmsk(i,6)= slmsk (i).eq.0.0                                              
   enddo
!                                                                               
   if(lat.ge.1.and.lat.le.nzl1) then                                         
!                                                                               
!  northern and southern polar region                                           
!                                                                               
     weis(2,1)=weis(2,1)+         float(idim)*wgt                              
     weis(6,1)=weis(6,1)+         float(idim)*wgt                              
     do l = 2,6                                                               
       isum = 0                                                                  
       jsum = 0                                                                  
       do i = 1,idim                                                        
         if(lsmsk(i,l)) isum = isum + 1                                            
         if(lsmsk(i+idim,l)) jsum = jsum + 1                                       
       enddo
       weis(2,l) = weis(2,l) + isum * wgt                                        
       weis(6,l) = weis(6,l) + jsum * wgt                                        
     enddo
   endif                                                                     
!                                                                               
!  northern and southern middle latitudes                                       
!                                                                               
   if(lat.ge.nzl1p.and.lat.le.nzl2) then                                     
     weis(3,1)=weis(3,1)+         float(idim)*wgt                              
     weis(5,1)=weis(5,1)+         float(idim)*wgt                              
     do l = 2,6                                                               
       isum = 0                                                                  
       jsum = 0                                                                  
       do i = 1,idim                                                        
         if(lsmsk(i,l)) isum = isum + 1                                            
         if(lsmsk(i+idim,l)) jsum = jsum + 1                                       
       enddo
       weis(3,l) = weis(3,l) + isum * wgt                                        
       weis(5,l) = weis(5,l) + jsum * wgt                                        
     enddo
   endif                                                                     
!                                                                               
   if(lat.ge.nzl2p) then                                                     
     weis(4,1)=weis(4,1)+         float(idimt)*wgt                             
     do l = 2,6                                                               
       isum = 0                                                                  
       do i = 1,idimt                                                       
         if(lsmsk(i,l)) isum = isum + 1                                            
       enddo
       weis(4,l) = weis(4,l) + isum * wgt                                        
     enddo
   endif                                                                     
!                                                                               
   if(lat.eq.jdimhf) then                                                    
     znlsl(1,1,nrcznl) = 0.                                                    
     do j = 2,6                                                               
       znlsl(1,1,nrcznl)=znlsl(1,1,nrcznl)+weis(j,1)                             
     enddo
!
!*** normalize latitude band weights with global sum                            
!
     do j = 2,6                                                               
       znlsl(j,1,nrcznl)=weis(j,1)/znlsl(1,1,nrcznl)                             
     enddo
!                                                                               
     do j = 2,6                                                               
       do l = 2,6                                                               
       enddo
       if(weis(j,1) .ne. 0.0) then                                               
         do l = 2,6                                                             
           znlsl(j,l,nrcznl)=weis(j,l)/weis(j,1)                                   
         enddo
       else                                                                      
         znlsl(j,l,nrcznl)=999.0                                                 
       endif                                                                    
     enddo
     znlsl(1,1,nrcznl)=1.0                                                     
!
!*** compute global coverages by latitude band weighting                        
!
     do l = 2,6                                                             
       znlsl(1,l,nrcznl) = 0.0                                                   
       do j = 2,6                                                             
         znlsl(1,l,nrcznl) = znlsl(1,l,nrcznl)                                 &
                    + znlsl(j,l,nrcznl)*znlsl(j,1,nrcznl)                  
       enddo
     enddo
     do l = 1,6                                                             
       do j = 1,6                                                             
         znlsl(j,l,nrcznl) = znlsl(j,l,nrcznl)*100.0                               
       enddo
     enddo
     do l = 1,6                                                               
       weis(1,l)=0.0                                                             
       do j = 2,6                                                               
         weis(1,l)=weis(1,l)+weis(j,l)                                             
       enddo
       do j = 1,6                                                               
         if(weis(j,l).ne.0.0) weis(j,l)=1.0/weis(j,l)                              
       enddo
     enddo
   endif                                                                     
!                                                                               
   return                                                                    
   end subroutine znlwgt
!-------------------------------------------------------------------------------
