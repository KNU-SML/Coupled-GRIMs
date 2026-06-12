#include <define.h>
   subroutine chgr_comp_rms
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [chgr_comp_rms] ----- [chgr_rms_sph] *
!                      |-- [chgr_rms_dfs] *
!                      |-- [chgr_rms_input] *
!                      |-- [chgr_rms_ouput] *
!
!-------------------------------------------------------------------------------
   end subroutine chgr_comp_rms
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_rms_sph(f, g, fgbar)
!-------------------------------------------------------------------------------
   use comchgr, only : mwavep, mdim => mdim2
!-------------------------------------------------------------------------------
   save
!
   real   ::  f(mdim)
   real   ::  g(mdim)
!-------------------------------------------------------------------------------
   fgbar = 0.
   l=0
   do i = 1,mwavep
     l=l+2
     fgbar = fgbar + f(l-1)*g(l-1)*0.5
     nnmax=mwavep-i+1
     do j = 2,nnmax
       l=l+2
       fgbar = fgbar+f(l-1)*g(l-1)+f(l  )*g(l  )
     enddo
   enddo
   fgbar =  sqrt  (fgbar)
!
   return
   end subroutine chgr_rms_sph
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_rms_dfs(f, g, fgbar)
!-------------------------------------------------------------------------------
#ifdef DFS
   use comchgr, only : mwavep, mdim => mdimd
#else
   use comchgr, only : mwavep, mdim => mdim2
#endif
!-------------------------------------------------------------------------------
   save
! 
   real               ::  f(mdim), g(mdim)
!-------------------------------------------------------------------------------
   fgbar = 0.
#ifdef DFS
   do l=1,mdim
     fgbar = fgbar + f(l)*g(l)
   enddo
!
#else
   l=0
   do i = 1,mwavep
     l=l+2
     fgbar = fgbar + f(l-1)*g(l-1)*0.5
     nnmax=mwavep-i+1
     do j = 2,nnmax
       l=l+2
       fgbar = fgbar+f(l-1)*g(l-1)+f(l  )*g(l  )
     enddo
   enddo
!
#endif
   fgbar =  sqrt  (fgbar)
!
   return
   end subroutine chgr_rms_dfs
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_rms_input(q,x,y,w,del,r) 
!-------------------------------------------------------------------------------
   use parmchgr, only : kdimi
   use comchgr, only : mdim
!-------------------------------------------------------------------------------
   save                                                                      
!
   real                 ::  q(mdim),x(mdim,kdimi),y(mdim,kdimi) 
   real                 ::  w(mdim,kdimi),r(mdim,kdimi) 
   real                 ::  del(kdimi) 
! 
   real                 ::  rx(kdimi),ry(kdimi),rw(kdimi),rr(kdimi)             
!-------------------------------------------------------------------------------
   vr=0.                                                                     
   vx=0.                                                                     
   vy=0.                                                                     
   vw=0.                                                                     
!
   do k = 1,kdimi                                                           
     call chgr_rms_sph(x(1,k),x(1,k),rx(k))
     call chgr_rms_sph(y(1,k),y(1,k),ry(k)) 
     call chgr_rms_sph(w(1,k),w(1,k),rw(k))  
     call chgr_rms_sph(r(1,k),r(1,k),rr(k))   
     vx=vx+rx(k)*del(k)                                                        
     vy=vy+ry(k)*del(k)                                                        
     vw=vw+rw(k)*del(k)                                                        
     vr=vr+rr(k)*del(k)                                                        
   enddo
! 
   call chgr_rms_sph(q,q,rq)   
!
      print 100,vx,vw,vy,vr,rq                                                  
100   format('chgr_rms_input div vort temp mixratio ln(ps)',5(e8.3,1x)) 
200   format(1h ,4(2x,e8.3))                                                    
!
   do k = 1,kdimi                                                           
     print 200,rx(k),rw(k),ry(k),rr(k)                                         
   enddo
!
   return                                                                    
   end subroutine chgr_rms_input   
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_rms_output(q,x,y,w,del,r)
!-------------------------------------------------------------------------------
   use paramter, only : kdim
#ifdef DFS
   use comchgr, only : mdim => mdimd
#else
   use comchgr, only : mdim
#endif
!-------------------------------------------------------------------------------
   save
!
   real                 ::  q(mdim),x(mdim,kdim),y(mdim,kdim)
   real                 ::  w(mdim,kdim),r(mdim,kdim)
   real                 ::  del(kdim)
   real                 ::  rx(kdim),ry(kdim),rw(kdim),rr(kdim)
!-------------------------------------------------------------------------------
   vr=0.
   vx=0.
   vy=0.
   vw=0.
!
   do k = 1,kdim
     call chgr_rms_dfs(x(1,k),x(1,k),rx(k))
     call chgr_rms_dfs(y(1,k),y(1,k),ry(k))
     call chgr_rms_dfs(w(1,k),w(1,k),rw(k))
     call chgr_rms_dfs(r(1,k),r(1,k),rr(k))
     vx=vx+rx(k)*del(k)
     vy=vy+ry(k)*del(k)
     vw=vw+rw(k)*del(k)
     vr=vr+rr(k)*del(k)
   enddo
   call chgr_rms_dfs(q,q,rq)
!
      print 100,vx,vw,vy,vr,rq
100   format('chgr_rms_output div vort temp mixratio ln(ps)',5(e8.3,1x))
200   format(1h ,4(2x,e8.3))
!
   do k = 1,kdim
     print 200,rx(k),rw(k),ry(k),rr(k)
   enddo
!
   return
   end subroutine chgr_rms_output
!-------------------------------------------------------------------------------
