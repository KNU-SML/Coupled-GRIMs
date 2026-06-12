#include "define.h"
   subroutine phys_grid_diffusion(deltim,prsl,pslap,gt0,gq0,ntotal,            &
                                    ids,ide, jds,jde, kds,kde,                 &
                                    ims,ime, jms,jme, kms,kme,                 &
                                    its,ite, jts,jte, kts,kte)
!-------------------------------------------------------------------------------
!
! program history log:
!   1995-01-01  masao kanamitsu        development
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
#ifdef DFS
!-------------------------------------------------------------------------------
   use dfsvar, only : iope
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ntotal,                                            &
                            ids,ide, jds,jde, kds,kde,                         &
                            ims,ime, jms,jme, kms,kme,                         &
                            its,ite, jts,jte, kts,kte
   real                 ::  prsl(ims:ime,kms:kme)
   real                 ::  pslap(ims:ime) 
   real                 ::  gt0(ims:ime,kms:kme)
   real                 ::  gq0(ims:ime,kms:kme*ntotal)
   real                 ::  deltim
!                                                                               
   integer              ::  i,k,n,ki,ke,kk,kd,ku,kd1,ku1
!-------------------------------------------------------------------------------
   do k = kts,kte
     kd=max(k-1,1)                                                           
     ku=min(k+1,kte)                                                       
     do i = its,ite
       gt0(i,k)=gt0(i,k)+pslap(i)*deltim*                                      &
                (gt0(i,ku)-gt0(i,kd))*prsl(i,k)/(prsl(i,ku)-prsl(i,kd))                   
     enddo                                                                   
   enddo                                                                     
!
   do n = 1,ntotal
     ki = (n-1)*kte + 1
     ke = n*kte 
     do k = kts,kte                                                              
       kk = (n-1)*kte + k
       kd=max(kk-1,ki)                                                           
       ku=min(kk+1,ke)                                                       
       kd1=kd-(n-1)*kte
       ku1=ku-(n-1)*kte
       do i = its,ite
         gq0(i,kk)=gq0(i,kk)+pslap(i)*deltim*                                  &
                  (gq0(i,ku)-gq0(i,kd))*prsl(i,k)/(prsl(i,ku1)-prsl(i,kd1))                   
       enddo                                                                   
     enddo                                                                     
   enddo                                                                     
!                                                                               
   return                                                                    
   end subroutine phys_grid_diffusion                                                                      
