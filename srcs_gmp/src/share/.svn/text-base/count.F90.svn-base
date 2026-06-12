!
   subroutine count(slimsk,sno,ijmax)                                        
!-------------------------------------------------------------------------------
   real   ::  slimsk(1),sno(1)                                                
!                                                                               
!  count number of points for the four surface conditions                       
!                                                                               
   l0=0                                                                      
   l1=0                                                                      
   l2=0                                                                      
   l3=0                                                                      
   l4=0                                                                      
!
   do ij = 1,ijmax                                                         
     if(slimsk(ij).eq.0.) l1=l1+1                                              
     if(slimsk(ij).eq.1. .and. sno(ij).le.0.) l0=l0+1                          
     if(slimsk(ij).eq.2. .and. sno(ij).le.0.) l2=l2+1                          
     if(slimsk(ij).eq.1. .and. sno(ij).gt.0.) l3=l3+1                          
     if(slimsk(ij).eq.2. .and. sno(ij).gt.0.) l4=l4+1                          
   enddo
!
   l5=l0+l3                                                                  
   l6=l2+l4                                                                  
   l7=l1+l6                                                                  
   l8=l1+l5+l6                                                               
   rl0=float(l0)/float(l8)*100.                                              
   rl3=float(l3)/float(l8)*100.                                              
   rl1=float(l1)/float(l8)*100.                                              
   rl2=float(l2)/float(l8)*100.                                              
   rl4=float(l4)/float(l8)*100.                                              
   rl5=float(l5)/float(l8)*100.                                              
   rl6=float(l6)/float(l8)*100.                                              
   rl7=float(l7)/float(l8)*100.                                              
   print *,'1) no. of not snow-covered land points   ',l0,' ',rl0,' '        
   print *,'2) no. of snow covered land points       ',l3,' ',rl3,' '        
   print *,'3) no. of open sea points                ',l1,' ',rl1,' '        
   print *,'4) no. of not snow-covered seaice points ',l2,' ',rl2,' '        
   print *,'5) no. of snow covered sea ice points    ',l4,' ',rl4,' '        
   print *,' '                                                               
   print *,'6) no. of land points                    ',l5,' ',rl5,' '        
   print *,'7) no. sea points (including sea ice)    ',l7,' ',rl7,' '        
   print *,'   (no. of sea ice points)          (',l6,')',' ',rl6,' '        
   print *,' '                                                               
   print *,'9) no. of total grid points               ',l8                   
   print *,' '                                                               
   print *,' '                                                               
!                                                                               
   return                                                                    
   end subroutine count                                                                       
!-------------------------------------------------------------------------------
