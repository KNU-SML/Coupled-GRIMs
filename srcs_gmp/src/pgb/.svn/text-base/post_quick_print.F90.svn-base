!
   subroutine post_quick_print(data,imax,jmax,fact)                  
!-------------------------------------------------------------------------------
   real  ::  data(imax*jmax)                                                 
!-------------------------------------------------------------------------------
   ilast=0                                                                   
   i1=1                                                                      
   i2=80                                                                     
   1112 continue                                                                  
!
   if(i2.ge.imax) then                                                       
     ilast=1                                                                 
     i2=imax                                                                 
   endif                                                                     
!
   write(6,*) ' '                                                            
!
   do j = 1,jmax                                                               
     write(6,1111) (nint(data(imax*(j-1)+i)*fact),i=i1,i2)                   
   enddo                                                                     
!
   if(ilast.eq.1) return                                                     
     i1=i1+80                                                                  
     i2=i1+79                                                                  
   if(i2.ge.imax) then                                                       
     ilast=1                                                                 
     i2=imax                                                                 
   endif                                                                     
!
   go to 1112                                                                
   1111 format(80i1)                                                              
!
   return                                                                    
   end                                                                       
