!
   subroutine file_write_byte(lu,lc,c)   
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:    file_write_byte       write data out by bytes                            
!                                                                               
! abstract: efficiently write unformatted a characeter array.                   
!                                                                               
! program history log:                                                          
!   1991-10-31  mark iredell                                                      
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!                                                                               
! usage:    call file_write_byte(lu,lc,c)              
!                                                                               
!   input argument list:                                                        
!     lu       - integer unit to which to write                                 
!     lc       - integer number of characters or bytes to write                 
!     c        - characeter (lc) data to write                                  
!                                                                               
!-------------------------------------------------------------------------------
   character  ::  c(lc)                                                           
!
   write(lu) c                                                               
!
   return                                                                    
   end subroutine file_write_byte
!-------------------------------------------------------------------------------
