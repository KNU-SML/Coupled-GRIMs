!
   subroutine chgr_horizontal_wavecut(maxwv,wave,maxwvc,wavec)                 
!-------------------------------------------------------------------------------
   save                                                                      
!                                                                               
   real  ::  wave(*),wavec(*)                                                
!                                                                               
   mdim=(maxwv+1)*(maxwv+2)                                                  
!                                                                               
   nmc=0                                                                     
   do m = 1,maxwvc+1                                                        
     nend=maxwvc+1                                                             
     do n = m,nend                                                            
       nmc=nmc+1                                                                
       nm=m*(maxwv+1)-(maxwv+1-n)-(m-1)*m/2                                     
       do ii = 1,2                                                              
         if(m.gt.maxwv+1.or.n.gt.maxwv+1) then
           wavec(nmc*2-2+ii)=0.
         else
           if(nm*2-2+ii.le.0) then
             print *,'logic error'
             call abort
           endif
           wavec(nmc*2-2+ii)=wave(nm*2-2+ii)
         endif
       enddo
     enddo 
   enddo
!                                                                               
   return                                                                    
   end subroutine chgr_horizontal_wavecut
!-------------------------------------------------------------------------------
