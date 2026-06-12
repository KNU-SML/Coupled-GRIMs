!
   subroutine setzro(fld,eps,ijdim)                                          
!-------------------------------------------------------------------------------
   real  ::  fld(ijdim)                                                      
!
   do ij = 1,ijdim                                                             
     if(abs(fld(ij)).lt.eps) fld(ij)=0.                                        
   enddo                                                                     
!
   return                                                                    
   end subroutine setzro
!-------------------------------------------------------------------------------
