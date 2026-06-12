!
   subroutine co2_quardratic(nlv,nlp1v,nlp2v,p,pd,trns)                              
!-------------------------------------------------------------------------------
   use comco2, only : p1, p2, trnslo, ia, ja, n
!-------------------------------------------------------------------------------
   real  ::  p(nlp1v),pd(nlp2v),trns(nlp1v,nlp1v)                            
   real  ::  wt(101)                                                         
!
   n2=2*n                                                                    
   n2p=2*n+1                                                                 
!
!  *****weights are calculated                                                  
!
   wt(1)=1.                                                                  
   do i = 1,n                                                               
     wt(2*i)=4.                                                                
     wt(2*i+1)=1.                                                              
   enddo
!
   if (n.eq.1) go to 25                                                      
!
   do i = 2,n                                                               
     wt(2*i-1)=2.                                                              
   enddo
!
   25  continue                                                                  
   trnsnb=0.                                                                 
   dp=(pd(ia)-pd(ia-1))/n2                                                   
   pfix=p(ja)                                                                
!
   do kk = 1,n2p                                                             
     pvary=pd(ia-1)+(kk-1)*dp                                                  
     if (pvary.ge.pfix) p2=pvary                                               
     if (pvary.ge.pfix) p1=pfix                                                
     if (pvary.lt.pfix) p1=pvary                                               
     if (pvary.lt.pfix) p2=pfix                                                
     call co2_interp       
     trnsnb=trnsnb+trnslo*wt(kk)                                               
   enddo
!
   trns(ia,ja)=trnsnb*dp/(3.*(pd(ia)-pd(ia-1)))                              
!
   return                                                                    
   end subroutine co2_quardratic
!-------------------------------------------------------------------------------
