#include <define.h>
   subroutine sph_comp_index
!-------------------------------------------------------------------------------
!
! subprogram: sph_comp_index          
!
!-------------------------------------------------------------------------------
   use comchgr, only   : mwavep, indxnn, indxmm
!-------------------------------------------------------------------------------
!                                                                               
   l=0                                                                       
   do m = 1,mwavep                                                          
     nend=mwavep-m+1                                                           
     do nn = 1,nend                                                           
       n=nn+m-1          
       l=l+2            
       indx=(mwavep*(n-m)-(n-m)*(n-m-1)/2+m)*2-1  
       indxnn(l-1)=indx                          
       indxnn(l  )=indx+1                       
     enddo
   enddo
!                                                                               
   l=0                                                                       
   do nn = 1,mwavep                                                         
     lln=mwavep-nn+1                                                           
     do ll = 1,lln                                
       n=ll+nn-1                              
       m=ll                                  
       indx=(m*mwavep-(mwavep-n)-(m-1)*m/2)*2-1  
       l=l+2                                    
       indxmm(l-1)=indx                        
       indxmm(l  )=indx+1                     
     enddo
   enddo
!                                                                               
   return                                                                    
   end subroutine sph_comp_index
!-------------------------------------------------------------------------------
