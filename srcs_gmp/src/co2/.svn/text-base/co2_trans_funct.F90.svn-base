!
   subroutine co2_trans_funct(itap1,ratstd,ratsm,ratio)
!-------------------------------------------------------------------------------
!      ratstd=value of higher std co2 concentration                             
!      ratsm=value of lower std co2 concentration                               
!      ratio=actual co2 concentration                                           
!      the 3 above quantities are in units of 330 ppmv.                         
!-------------------------------------------------------------------------------
   use comco2, only : p1, p2, trnslo, ia, ja, n
   use comco2, only : pa
   use comco2, only : transa
!-------------------------------------------------------------------------------
   real  ::  trns1(109,109),trns2(109,109)                                   
!
!   read in tfs of lower std co2 concentration                                  
!
   read (itap1,100) ((trns1(i,j),i=1,109),j=1,109)                           
   100  format (4f20.14)                                                          
   itap2=itap1+1                                                             
!
!   read in tfs of higher std co2 concentration                                 
!
   read (itap2,100) ((transa(i,j),i=1,109),j=1,109)                          
   call co2_tau_interp(ratstd)                                                       
!
   do 401 i = 1,109                                                            
     do 401 j = 1,i                                                              
       if (j.eq.i) go to 401                                                     
!
!  using higher co2 concentration,compute 1st guess co2 tfs for                 
!  actual co2 concentration.                                                    
!
       p2=(ratio+ratstd)*pa(i)/(2.*ratstd)  +                                &
            (ratstd-ratio)*pa(j)/(2.*ratstd)                                    
       p1=(ratstd-ratio)*pa(i)/(2.*ratstd)  +                                &
            (ratio+ratstd)*pa(j)/(2.*ratstd)                                    
       call co2_interp                                                            
       trnspr=trnslo                                                          
!
!  using higher co2 concentration,compute 1st guess co2 tfs for                 
!  lower std co2 concentration                                                  
!
       p2=(ratsm+ratstd)*pa(i)/(2.*ratstd)  +                                & 
             (ratstd-ratsm)*pa(j)/(2.*ratstd)                                       
       p1=(ratstd-ratsm)*pa(i)/(2.*ratstd)  +                                & 
             (ratsm+ratstd)*pa(j)/(2.*ratstd)                                       
       call co2_interp                                                               
       trnspm=trnslo                                                             
!
!  compute tfs for co2 concentration given by (ratio).                          
!   store temporarily in (trns2)                                                
!
       trns2(j,i)=trnspr+(ratstd-ratio)*(trns1(j,i)-                         &
                    trnspm)/(ratstd-ratsm)                                                   
       trns2(i,j)=trns2(j,i)                                                     
!
! we now can overwrite (trns1) and store in (trns1) the 1st guess               
!  co2 tfs for lower std co2 concentration                                      
!
       trns1(j,i)=trnslo                                                         
       trns1(i,j)=trnslo                                                         
401   continue                                                                  
!
!  set diagonal values of co2 tfs to unity                                      
!
   do i = 1,109                                                            
     trns1(i,i)=1.0                                                            
     trns2(i,i)=1.0                                                            
   enddo
!
!  now output the computed co2 tfs for (ratio) co2 conc. in (transa)            
!
   do i = 1,109                                                            
     do j = 1,109                                                            
       transa(j,i)=trns2(j,i)                                                    
     enddo
   enddo
!
   return                                                                    
   end subroutine co2_trans_funct
!-------------------------------------------------------------------------------
