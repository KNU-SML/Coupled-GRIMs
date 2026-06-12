   subroutine post_byte2char(l,c,z)                                                     
!-------------------------------------------------------------------------------
!
! subprogram:    post_byte2char     convert byte to hexadecimal character pair         
!                                                                               
! abstract: converts an array of bytes to its hexadecimal representation        
!   (2 characters per byte) for diagnostic purposes.                            
!                                                                               
! program history log:                                                          
!   1991-10-31  mark iredell                                                      
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!                                                                               
! usage:    call post_byte2char(l,c,z)                                                     
!                                                                               
!   input argument list:                                                        
!     l        - integer number of bytes to represent                           
!     c        - character (l) byte data to convert                             
!                                                                               
!   output argument list:                                                       
!     z        - character (2*l) hexadecimal representation                     
!                                                                               
!-------------------------------------------------------------------------------
   character  ::  c(l)*1,z(l)*2                                                   
!-------------------------------------------------------------------------------
!
   do i = 1,l                                                                  
     write(z(i),'(z2)') ichar(c(i))                                          
   enddo                                                                     
!
   return                                                                    
   end subroutine post_byte2char
