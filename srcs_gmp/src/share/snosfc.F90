!
   subroutine snosfc(snoanl,tsfanl,tsfsmx,ijdim)                             
!-------------------------------------------------------------------------------
   real  ::   snoanl(ijdim)                                                   
   real  ::   tsfanl(ijdim)                                                   
!                                                                               
!     write(6,*) 'set snow temp to tsfsmx if greater'                           
!
   kount=0                                                                   
   do ij = 1,ijdim                                                             
     if(snoanl(ij).gt.0.) then                                               
       if(tsfanl(ij).gt.tsfsmx) tsfanl(ij)=tsfsmx                            
       kount=kount+1                                                         
     endif                                                                   
   enddo                                                                     
!
   if(kount.gt.0) then                                                       
     per=float(kount)/float(ijdim)*100.                                      
!    write(6,*) 'snow sfc.  tsf set to ',tsfsmx,' at ',                        &
!                  kount, ' points ',per,'percent'                             
   endif                                                                     
!
   return                                                                    
   end subroutine snosfc
!-------------------------------------------------------------------------------
