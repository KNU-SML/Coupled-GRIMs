#include "define.h"
#ifndef HYBRID
   subroutine chgr_write_sigma(n,fhour,idate,q,te,di,ze,chgrp,sl,si,gz         &
#else
   subroutine chgr_write_sigma(n,fhour,idate,q,te,di,ze,chgrp,ak5,bk5,gz       &
#endif
!-------------------------------------------------------------------------------
#ifdef DFS
      ,tave                                                                    &
#endif
      ,ltrn                                                                    &
      ,waves,xlayers,trun,order,realform,gencode                               &
      ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,water                               &
      ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid                    &
      ,pdryini,dummy2,gases                                                    &
      ,iret)
!-------------------------------------------------------------------------------
!
! subprogram: chgr_write_sigma         
!
!-------------------------------------------------------------------------------
   use paramodel, only : igen_,levs_,jcap_,latg_,lonf_
   use paramter, only  : kdim, kdimq
   use comchgr, only   : kdimp, indxmm, b, c, lab,                             &
#ifdef DFS
                         mdim => mdimd
#else
                         mdim
#endif
   use comcon
   save                                                                      
!-------------------------------------------------------------------------------
   real                 ::  gz(mdim),q(mdim),                                  &
                            te(mdim,kdim),di(mdim,kdim),ze(mdim,kdim),         &
                            chgrp(mdim,kdimq)
   real                 ::  si(kdimp),sl(kdim)                                 
   integer              ::  idate(4)
!
   integer              ::  kdum
   integer, parameter   ::  kdum2=21,kens=2
#ifdef HYBRID
   real                 ::  ak5(kdimp),bk5(kdimp),ak5x(kdimp),bk5x(kdimp)
#endif
   real, allocatable    ::  dummy(:)
   real                 ::  dummy2(kdum2),ensemble(kens)
#ifdef DFS
   real                 ::  tave(kdim)
#endif
!                                                                             
   logical ltrn                                                              
!-------------------------------------------------------------------------------
#ifndef HYBRID
   kdum=201-levs_-1-levs_
#else
   kdum=201-levs_-2-levs_
#endif
   allocate(dummy(kdum))
!
   call chgr_make_header(lab,idate,fhour,ifin)   
!
   write(n) lab                                                              
!
!     write(n)fhour,idate,si,sl                                                 
!    &       ,dummy,waves,xlayers,trun,order,realform,gencode                   
!
   do i = 1,kdum
     dummy(i)=0.
   enddo
!
   waves=jcap_
   xlayers=levs_
   trun=1.
   order=2.
   realform=1.
   gencode=igen_
   rlond=lonf_
   rlatd=latg_
   rlonp=lonf_
   rlatp=latg_
   rlonr=lonf_
   rlatr=latg_
   gases=0.
   water=1.
   pdryini=0.
   subcen=0.
!
   do i = 1,kens
     ensemble(i)=0.
   enddo
!
   ppid=0.
   slid=0.
   vcid=0.
   vmid=0.
   vtid=0.
!
   do k = 1,kdum2
     dummy2(k)=0.
   enddo
!
#ifndef HYBRID
   write(n)fhour,idate,si,sl                                                   &
#else
   do k = 1,levs_+1
     ak5x(k)=ak5(k)*1000. ! cb -> Pa
     bk5x(k)=bk5(k)
   enddo
!
   write(n)fhour,idate,ak5x,bk5x                                               &
#endif
      ,dummy,waves,xlayers,trun,order,realform,gencode                         &
      ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,water                               &
      ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid                    &
      ,pdryini,dummy2,gases
#ifdef DFS
   write(n) gz(:)
   write(n) q(:),tave(:)
!                                                                               
   call print_maxmin_six(te,mdim,kdim,1,kdim,'temp in chgr_write_sigma')
   do k = 1,kdim                                                            
     write(n) te(:,k)
   enddo
!                                                                               
   do k = 1,kdim                                                            
     write(n) di(:,k)
     write(n) ze(:,k)
   enddo
!                                                                               
   call print_maxmin_six(chgrp,mdim,kdimq,1,kdimq,'rq in chgr_write_sigma')
   do k = 1,kdimq                                                           
     write(n) chgrp(:,k)
   enddo
#else
   if(ltrn) then                                                             
     do i = 1, mdim                                                       
       b(indxmm(i)) = gz(i) 
       c(indxmm(i)) =  q(i)
     enddo
     write(n) b           
     write(n) c          
!                       
     do k = 1,kdim    
       do i = 1, mdim
         b(indxmm(i)) = te(i,k)
       enddo
       write(n) b   
     enddo
!                                                                               
     do k = 1,kdim                                                            
       do i = 1, mdim                                                       
         b(indxmm(i)) = di(i,k)
         c(indxmm(i)) = ze(i,k) 
       enddo
       write(n) b 
       write(n) c
     enddo
!                                                                               
     do k = 1,kdimq                                                           
       do i = 1, mdim                                                       
         b(indxmm(i)) = chgrp(i,k) 
       enddo
       write(n) b  
     enddo
!                                                                               
   else 
!                                                                               
     write(n) gz                                                               
     write(n) q                                                                
!                                                                               
     do k = 1,kdim                                                            
       write(n) (te(i,k),i=1,mdim)
     enddo
!                                                                               
     do k = 1,kdim                                                            
       write(n) (di(i,k),i=1,mdim) 
       write(n) (ze(i,k),i=1,mdim)
     enddo
!                                                                               
     do k = 1,kdimq                                                           
       write(n) (chgrp(i,k),i=1,mdim)
     enddo
!                                                                               
   endif                                                                     
#endif
!                                                                               
   deallocate(dummy)
!
   return                                                                    
   end                                                                       
!-------------------------------------------------------------------------------
