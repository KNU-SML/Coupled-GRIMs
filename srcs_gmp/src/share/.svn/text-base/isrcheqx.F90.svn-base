#include <define.h>
   function isrcheqx(n,x,incx,target)                                        
!-------------------------------------------------------------------------------
   real   ::  x(*)                                                            
   real   ::  target                                                               
!
   j=1                                                                       
   isrcheqx=0                                                                
!
   if(n.le.0) return                                                         
   if(incx.lt.0) j=1-(n-1)*incx                                              
!
   do i = 1,n                                                                  
     if(x(j).eq.target) then                                                 
       isrcheqx=j                                                            
       return                                                                
     endif                                                                   
     j=j+incx                                                                
   enddo                                                                     
!
   return                                                                    
   end function isrcheqx
!-------------------------------------------------------------------------------
