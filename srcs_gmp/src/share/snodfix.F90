!
   subroutine snodfix(snoanl,snofcs,ijdim)                                   
!-------------------------------------------------------------------------------
   real  ::  snoanl(ijdim),snofcs(ijdim)                                     
!
   do ij = 1,ijdim                                                             
     if(snoanl(ij).gt.0..and.snofcs(ij).gt.0.) snoanl(ij)=snofcs(ij)         
   enddo                                                                     
!
   return                                                                    
   end subroutine snodfix
!-------------------------------------------------------------------------------
