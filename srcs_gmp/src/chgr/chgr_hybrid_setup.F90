#include <define.h>
   subroutine chgr_hybrid_setup(ak5,bk5,ci, si, del, sl, cl, rpi)     
!-------------------------------------------------------------------------------
! 
! subprogram: chgr_hybrid_setup
!
! abstract: this routine sets coordinates for levels and use phillips method
! to get layers. ( can be stretched levels or not )                             
!
!-------------------------------------------------------------------------------
   use paramodel, only : levs_
   use constant,  only : cp_,akapa_
   use paramter, only  : kdim
   use comchgr, only   : kdimp, kdimm
!-------------------------------------------------------------------------------
   save
   real                ::  ci(kdimp),si(kdimp),del(kdim),sl(kdim),cl(kdim),    &
                           rpi(kdimm)
   real                ::  ak5(kdimp),bk5(kdimp)
   real                ::  rk1,rkinv   
   real, dimension(100)::  akmdl, bkmdl
   data ifp/0/                                                               
!-------------------------------------------------------------------------------
   namelist /modlhyb/ akmdl, bkmdl
   akmdl = 0.
   bkmdl = 0.
   read(1,modlhyb)
!
   do k = 1,kdimp
     ak5(k)=akmdl(k)/1000. ! Pa -> cb
     bk5(k)=bkmdl(k)
   enddo
!
   ci(1) = 0.                                                                
   do k = 1,kdim                                                             
     del(k)=ak5(k)/100+bk5(k)-ak5(k+1)/100-bk5(k+1)
     ci(k+1)=ci(k)+del(k)                                                   
   enddo
   ci(kdimp)=1.                                                              
!                                                                               
   rk1 = akapa_ + 1.
   rkinv=1./akapa_
!                                                                               
   levs=kdim
!                                                                               
   do li = 1,kdimp
     si(li) = 1. - ci(li)                                                      
   enddo                                                                     
!                                                                               
   do le = 1,kdim
     dif = si(le)**rk1 - si(le+1)**rk1 
     dif = dif / (rk1*(si(le)-si(le+1))) 
     sl(le) = dif**rkinv
     cl(le) = 1. - sl(le) 
   enddo
!                                                                               
!  compute pi ratios for temp. matrix.                                       
!                                                                               
   do le = 1,kdimm                                                           
     rpi(le) = (sl(le+1)/sl(le))                                               
   enddo
   do le = 1,kdimm                                                          
     rpi(le) = rpi(le)**akapa_  
   enddo
!                                                                               
   do le = 1,kdimp                                                           
     print 100, le, ci(le), si(le)                                             
     100   format (1h , 'level=', i2, 2x, 'ci=', f6.3, 2x, 'si=', f6.3)
   enddo
! 
   print 200
   200 format (1h0) 
!
   do le=1,kdim    
     print 300, le, cl(le), sl(le), del(le)                                    
     300 format (1h , 'layer=', i2, 2x, 'cl=', f6.3, 2x, 'sl=', f6.3, 2x,      &
                 'del=', f6.3)   
   enddo
!                                                                               
   print 400, (rpi(le), le=1,kdimm)                                          
   400 format (1h0, 'rpi=', (18(1x,f6.3)) ) 
!                                                                               
   return
   end  
!-------------------------------------------------------------------------------
