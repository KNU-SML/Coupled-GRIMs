!
   subroutine znlavs(f,lat,wgt,ind,                                            &
                     znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)        
!-------------------------------------------------------------------------------
   logical lsmsk(idimt,6)                                                    
   real   ::  znlsl(6,6,nrcznl),weis(6,6)                                     
!                                                                               
   real   ::  f(1)                                                            
!                                                                               
   if(ind.gt.nrcznl) return                                                  
!                                                                               
   jdimhf=jdim/2                                                             
!                                                                               
   n=idim+1                                                                  
!                                                                               
   if(lat.eq.1) then                                                         
     do j = 1,6                                                           
       do i = 1,6                                                           
         znlsl(i,j,ind) = 0.0                                                      
       enddo
     enddo
   endif                                                                     
!                                                                               
   if(lat.le.nzl1) then                                                      
     do i = 1,idim                                                        
       znlsl(2,1,ind) = znlsl(2,1,ind) + f(i)      * wgt                         
       znlsl(6,1,ind) = znlsl(6,1,ind) + f(i+idim) * wgt                         
     enddo
     do l = 2,6                                                               
       do i = 1,idim                                                        
         if(lsmsk(i,l)) then                                                       
           znlsl(2,l,ind) = znlsl(2,l,ind) + f(i) * wgt                            
         endif                                                                     
         if(lsmsk(i+idim,l)) then                                                  
           znlsl(6,l,ind) = znlsl(6,l,ind) + f(i+idim) * wgt                       
         endif                                                                     
       enddo
     enddo
   endif                                                                     
!                                                                               
   if(lat.gt.nzl1.and.lat.le.nzl2) then                                      
     do i = 1,idim                                                        
       znlsl(3,1,ind) = znlsl(3,1,ind) + f(i)      * wgt                         
       znlsl(5,1,ind) = znlsl(5,1,ind) + f(i+idim) * wgt                         
     enddo
     do l = 2,6                                                               
       do i = 1,idim                                                        
         if(lsmsk(i,l)) then                                                       
           znlsl(3,l,ind) = znlsl(3,l,ind) + f(i) * wgt                            
         endif                                                                     
         if(lsmsk(i+idim,l)) then                                                  
           znlsl(5,l,ind) = znlsl(5,l,ind) + f(i+idim) * wgt                       
         endif                                                                     
       enddo
     enddo
   endif                                                                     
!                                                                               
   if(lat.gt.nzl2) then                                                      
     do i = 1,idimt                                                       
       znlsl(4,1,ind) = znlsl(4,1,ind) + f(i)      * wgt                         
     enddo
     do l = 2,6                                                               
       do i = 1,idimt                                                       
         if(lsmsk(i,l)) then                                                       
           znlsl(4,l,ind) = znlsl(4,l,ind) + f(i) * wgt                            
         endif                                                                     
       enddo
     enddo
   endif                                                                     
!                                                                               
   if(lat.eq.jdimhf) then                                                    
     do l = 1,6                                                               
       znlsl(1,l,ind)=0.0                                                        
       do j = 2,6                                                               
         znlsl(1,l,ind)=znlsl(1,l,ind)+znlsl(j,l,ind)                              
       enddo
       do j = 1,6                                                               
         znlsl(j,l,ind)=znlsl(j,l,ind)*weis(j,l)                                   
       enddo
     enddo
   endif                                                                     
!                                                                               
   return                                                                    
   end subroutine znlavs                                                                      
!-------------------------------------------------------------------------------
