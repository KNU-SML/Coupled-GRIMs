#include "define.h"
!-------------------------------------------------------------------------------
   subroutine dyn_comp_rms(q,x,y,w,del,r)
!-------------------------------------------------------------------------------
   use paramodel, only : levh_,levs_,ntotal_,LNT22S
#ifdef DFS
   use dfsvar, only : mt,mtg,jl,jlg,levs,levh,levsp,levhp,iope,ntotal
#endif
#ifdef MP
   use paramodel, only : lnt22_,lnt22p_, lnt2_
   use commpi
#endif
   use comio
!-------------------------------------------------------------------------------
!      
! subprogram:    dyn_comp_rms       computes root mean square.
!
! abstract: computes the root mean square in each level and of
!   the vertical integral (when appropriate) given the spectral
!   coefficients of the model variables or the tendencies.
!   the results are printed.  dyn_comp_rms is strictly diagnostic.
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call dyn_comp_rms (q, x, y, w, del, r)
!   input argument list:
!     q        - spectral coefs of ln(psfc)     or its tendency.
!     x        - spectral coefs of divergence   or its tendency.
!     y        - spectral coefs of temperature  or its tendency.
!     w        - spectral coefs of vorticity    or its tendency.
!     del      - sigma spacing at each layer.
!     r        - spectral coefs of mixing ratio or its tendency.
!
!   output files:
!     output   - print file.
!
!  ::: structure :::
!
!    [dyn_comp_rms] --- [dyn_comp_rms_sub]
!
!-------------------------------------------------------------------------------
#ifdef RMP
#undef MP
#endif
!
#ifdef DFS
#define LNT22S mt*jlg
#define LNT22_ mtg*jlg
#define LNT2_  mtg*jlg
#endif
!
#ifdef DFS
   real                 ::  q(LNT22S),x(LNT22S,levsp),y(LNT22S,levsp)
   real                 ::  w(LNT22S,levsp),r(LNT22S,levhp),del(levs_)
#else
   real                 ::  q(LNT22S),x(LNT22S,levs_),y(LNT22S,levs_)
   real                 ::  w(LNT22S,levs_),r(LNT22S,levh_),del(levs_)
#endif
!
! local
!
#ifdef MP
   real,allocatable     ::  qf(:),xf(:,:),yf(:,:),wf(:,:),rf(:,:)
#endif
   real                 ::  rx(levs_),ry(levs_),rw(levs_),rr(levh_)
!
#ifdef MP
   allocate(qf(LNT22_))
   allocate(xf(LNT22_,levs_))
   allocate(yf(LNT22_,levs_))
   allocate(wf(LNT22_,levs_))
   allocate(rf(LNT22_,levh_))
   qf=0. ; xf=0. ; yf=0. ; wf=0. ; rf=0.
#ifdef DFS
   call mpm2f (q,mt,qf,mtg,jlg)
   call mpmz2f(x,mt,jlg,levsp,xf,mtg,jlg,levs,levs)
   call mpmz2f(y,mt,jlg,levsp,yf,mtg,jlg,levs,levs)
   call mpmz2f(w,mt,jlg,levsp,wf,mtg,jlg,levs,levs)
   call mpmz2f(r,mt,jlg,levhp,rf,mtg,jlg,levs,levh)
#else
   call mpsp2f(q,lnt22p_,qf,lnt22_,1)
   call mpsp2f(x,lnt22p_,xf,lnt22_,levs_)
   call mpsp2f(y,lnt22p_,yf,lnt22_,levs_)
   call mpsp2f(w,lnt22p_,wf,lnt22_,levs_)
   call mpsp2f(r,lnt22p_,rf,lnt22_,levh_)
#endif
   if( iope ) then
!#ifdef DBG
!     print *,' In dyn_comp_rms PE0'
!     call print_maxmin_seven(qf,LNT2_,LNT22_,1,1,1,'ps')
!     call print_maxmin_seven(yf,LNT2_,LNT22_,levs_,1,levs_,'t')
!     call print_maxmin_seven(xf,LNT2_,LNT22_,levs_,1,levs_,'div')
!     call print_maxmin_seven(wf,LNT2_,LNT22_,levs_,1,levs_,'vor')
!     call print_maxmin_seven(rf,LNT2_,LNT22_,levh_,1,levh_,'q')
!#endif
#define QS qf
#define XS xf
#define YS yf
#define WS wf
#define RS rf
#else
#define QS q
#define XS x
#define YS y
#define WS w
#define RS r
#endif
     vr=0.e0
     vx=0.e0
     vy=0.e0
     vw=0.e0
     do k = 1,levs_
       call dyn_comp_rms_sub(XS(1,k),XS(1,k),rx(k))
       call dyn_comp_rms_sub(YS(1,k),YS(1,k),ry(k))
       call dyn_comp_rms_sub(WS(1,k),WS(1,k),rw(k))
       vx=vx+rx(k)*del(k)
       vy=vy+ry(k)*del(k)
       vw=vw+rw(k)*del(k)
     enddo
     do n = 1,ntotal_
       do k = 1,levs_
         nk = (n-1)*levs_+k
         call dyn_comp_rms_sub(RS(1,nk),RS(1,nk),rr(nk))
         vr=vr+rr(nk)*del(k)
       enddo
     enddo
     call dyn_comp_rms_sub(QS,QS,rq)
!
     write(6,100)vx,vw,vy,vr,rq
     do k = 1,levs_
       write(6,200)rx(k),rw(k),ry(k),rr(k)
     enddo
     if(ntotal_.ge.2) then
       write(6,300)ntotal_ - 1
       do k = 1,levs_
         ki = levs_+k
         ke = (ntotal_-1)*levs_  + k
         if(ntotal_.eq.2) write(6,402)(rr(kk),kk=ki,ke,levs_)
         if(ntotal_.eq.3) write(6,403)(rr(kk),kk=ki,ke,levs_)
         if(ntotal_.eq.4) write(6,404)(rr(kk),kk=ki,ke,levs_)
         if(ntotal_.eq.5) write(6,405)(rr(kk),kk=ki,ke,levs_)
         if(ntotal_.eq.6) write(6,406)(rr(kk),kk=ki,ke,levs_)
         if(ntotal_.eq.7) write(6,407)(rr(kk),kk=ki,ke,levs_)
         if(ntotal_.gt.7) write(6,408)(rr(kk),kk=ki,ke,levs_)
       enddo
     endif
     call flush(6)
100  format(1h0,'div vort temp mixratio ln(ps) ',5(e15.8,1x))
200  format(1h ,4(2x,e17.10))
300  format(1h0,'clouds and gases = ',i5)
402  format(1x,1(2x,e17.10))
403  format(1x,2(2x,e17.10))
404  format(1x,3(2x,e17.10))
405  format(1x,4(2x,e17.10))
406  format(1x,5(2x,e17.10))
407  format(1x,6(2x,e12.5))
408  format(1x,9(2x,e10.3))
#undef QS
#undef XS
#undef YS
#undef WS
#undef RS
#ifdef MP
   endif
!
   deallocate(qf , xf , yf , wf , rf)
#endif
!
   return
   end subroutine dyn_comp_rms
!-------------------------------------------------------------------------------
   subroutine dyn_comp_rms_sub (f, g, fgbar)      
!-------------------------------------------------------------------------------
   use paramodel, only : jcap_,lnt22_,jcap1_,jcap2_
!
#ifdef DFS
   real                 ::  f(-jcap_:jcap_,0:jcap_),g(-jcap_:jcap_,0:jcap_)
#else
   real                 ::  f(lnt22_),g(lnt22_)         
!-------------------------------------------------------------------------------
   joff(n,L)=(jcap1_)*(jcap2_)-(jcap1_-L)*(jcap2_-L)+2*(n-L)                 
!                                                                               
   l=0                                                                       
#endif
   fgbar = 0.                                                                
#ifdef DFS
   do n = 0,jcap_                                                           
     fgbar = fgbar + f(0,n)*g(0,n)*0.5
     do L = -jcap_,1
       fgbar = fgbar + f(L,n)*g(L,n)
     enddo
     do L = 1,jcap_
       fgbar = fgbar + f(L,n)*g(L,n)
     enddo
   enddo
#else
   do n = 0,jcap_                                                           
     fgbar = fgbar + f(joff(n,l)+1)*g(joff(n,l)+1)                             
   enddo                                                              
   do n = 0,jcap_                                                          
     fgbar = fgbar + f(joff(n,l)+2)*g(joff(n,l)+2)                             
   enddo                                                                  
   fgbar=fgbar*0.5                                                           
   do L = 1,jcap_                                                           
     do n = L, jcap_                                                           
       fgbar = fgbar + f(joff(n,L)+1)*g(joff(n,L)+1)                             
     enddo                                                                  
     do n = L, jcap_                                                          
       fgbar = fgbar + f(joff(n,L)+2)*g(joff(n,L)+2)                             
     enddo                                                                  
   enddo                                                                  
#endif
   fgbar = sqrt(fgbar)                                                       
!
   return                                                                    
   end subroutine dyn_comp_rms_sub
