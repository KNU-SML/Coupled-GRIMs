#include <define.h>
   subroutine co2_new_sigma(ci, si, del, sl, cl, rpi)                               
!-------------------------------------------------------------------------------
   use constant, only   : rd_,cp_
   use paramodel, only  : latg=>latg_,lonf=>lonf_
   use paramco2, only   : kdim
!-------------------------------------------------------------------------------
!   save                                                                      
! ------------------********************-----------------------                 
! this routine sets coordinates for levels and use phillips method              
! to get layers. ( can be stretched levels or not )                             
!  --------------------------------------------------------------               
   integer  :: kdimp,kdimm
!                                                                               
   real  ::  ci(kdim+1),si(kdim+1),del(kdim),sl(kdim),cl(kdim),rpi(kdim)
!                                                                               
   real  ::  rk,rk1,rkinv                                                         
!                                                                               
   real, dimension(100)  ::  delmdl
!                                                                               
   data ifp/0/                                                               
!
   namelist /modlsig/ delmdl
   delmdl = 0.
   read(1,modlsig)
!
   kdimp=kdim+1
   kdimm=kdim-1
!
!  in case of shocking parameter discord snafu---                               
!
   do k = 1,kdim                                                             
     del(k)=delmdl(k)                                                          
   enddo
!                                                                               
   ci(1) = 0.                                                                
   do k = 1,kdim                                                             
     ci(k+1)=ci(k)+del(k)                                                      
   enddo
   ci(kdimp)=1.                                                              
!                                                                               
   rk  = rd_/cp_                                                             
   rk1 = rk + 1.                                                             
   rkinv=1./rk                                                               
!                                                                               
   levs=kdim                                                                 
!                                                                               
   do li = 1,kdimp                                                           
     si(li) = 1. - ci(li)                                                      
   enddo
!                                                                               
   do le = 1,kdim                                                            
     dif = si(le)**rk1 - si(le+1)**rk1                                         
     dif = dif /(rk1*(si(le)-si(le+1)))                                        
     sl(le) = dif**rkinv                                                       
     cl(le) = 1. - sl(le)                                                      
   enddo
!                                                                               
!     compute pi ratios for temp. matrix.                                       
!                                                                               
   do le = 1,kdimm                                                           
     rpi(le) = (sl(le+1)/sl(le))                                               
   enddo
!
   do le = 1,kdimm                                                          
     rpi(le) = rpi(le)**rk                                                     
   enddo
!                                                                               
   do le = 1,kdimp                                                           
     print 100, le, ci(le), si(le)                                             
     100 format (1h , 'level=', i2, 2x, 'ci=', f6.3, 2x, 'si=', f6.3)              
   enddo
!                                                                               
   print 200                                                                 
   200 format (1h0)                                                              
!
   do le = 1,kdim                                                            
     print 300, le, cl(le), sl(le), del(le)                                    
     300 format (1h , 'layer=', i2, 2x, 'cl=', f6.3, 2x, 'sl=', f6.3, 2x,      &
                   'del=', f6.3)                                               
   enddo
!                                                                               
   print 400, (rpi(le), le=1,kdimm)                                          
   400 format (1h0, 'rpi=', (18(1x,f6.3)) )                                      
!                                                                               
   return                                                                    
   end subroutine co2_new_sigma
!-------------------------------------------------------------------------------
