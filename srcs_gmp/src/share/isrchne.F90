!
   function isrchne(n,a,is,val)                                              
!-------------------------------------------------------------------------------
   integer  ::  val                                                               
   integer  ::  a(n)                                                              
!
   isrchne=n                                                                 
!
   do nn = is,n                                                                
     if( a(nn).ne.val ) then                                                   
       isrchne=nn                                                              
       return                                                                  
     endif                                                                     
   enddo                                                                     
!
   return                                                                    
   end function isrchne                                                                      
!-------------------------------------------------------------------------------
