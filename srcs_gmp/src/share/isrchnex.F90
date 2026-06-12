#include <define.h>
   function isrchnex(n,x,incx,target)                                        
!-------------------------------------------------------------------------------
   integer   ::  x(*), target                                                      
!
   j=1                                                                       
   isrchnex=0                                                                
!
   if(n.le.0) return                                                         
   if(incx.lt.0) j=1-(n-1)*incx                                              
!
   do i = 1,n                                                                  
     if(x(j).ne.target) then                                                 
       isrchnex=i                                                            
       return                                                                
     endif                                                                   
     j=j+incx                                                                
   enddo                                                                     
!
   return                                                                    
   end function isrchnex
!-------------------------------------------------------------------------------
