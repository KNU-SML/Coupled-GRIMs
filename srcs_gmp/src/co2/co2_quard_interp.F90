!
   subroutine co2_quard_interp(xm,x0,xp,fm,f0,fp,x,f)
!-------------------------------------------------------------------------------
!....                                                                           
!     double precision fm,f0,fp,f,d1,d2,b,a,del                                 
!
   d1=(fp-f0)/(xp-x0)                                                        
   d2=(fm-f0)/(xm-x0)                                                        
   b=(d1-d2)/(xp-xm)                                                         
   a=d1-b*(xp-x0)                                                            
   del=(x-x0)                                                                
   f=f0+del*(a+del*b)                                                        
   return                                                                    
   end subroutine co2_quard_interp
!-------------------------------------------------------------------------------
