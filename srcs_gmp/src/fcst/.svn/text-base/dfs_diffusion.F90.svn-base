#include "define.h"
!-------------------------------------------------------------------------------
!  ::: structure ::: This file contains ...
!
!        |--- [dfs_diffusion] *
!        |--- [dfs_diffusion_coef] *
!        |--- [dfs_diffusion_del_matrix] *
!        |--- [dfs_diffusion_driver] *
!        |--- [dfs_diffusion_laplacian] *
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_diffusion( DT,D1,V1,Q1,T1,PSL1,VJ,QJ,Tj,PSJ,                 &
                     D3,V3,Q3,T3,PSL3)
!-------------------------------------------------------------------------------
!
! abstract :
!  - spectral conponent of stream functions of temperature and surface pressure
!  - with leap-frog scheme
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
   use dfsvar, only : mt,jl,levs,levh,levsp,levhp,jlg,ngs,nge,trsdm,phsdm
!-------------------------------------------------------------------------------
   implicit none
!
   integer                         ::  k,mj,mx
   real                            ::  dt, dt2
   real, dimension(mt,jlg)         ::  psl1,psj,psl3
   real, dimension(mt,jlg,levsp)   ::  d1,   d3  ,&
                                       v1,vj,v3  ,&
                                       t1,tj,t3
   real, dimension(mt,jlg,levhp)   ::  q1,qj,q3
   real, dimension(2*levsp)        ::  XMAT
!
   dt2=dt*2.0
   do k = 1,levsp
     do mj = 1,JLg
       do mx = 1,MT
         V3(mx,mj,k)= V1(mx,mj,k) + DT2* VJ(mx,mj,k)
         T3(mx,mj,k)= T1(mx,mj,k) + DT2*(TJ(mx,mj,k)-trsdm(mx,mj,k))
#ifndef BACKWARD
         D3(mx,mj,k)=2.0*D3(mx,mj,k)-D1(mx,mj,k)
#endif
       enddo
     enddo
   enddo ! k=1,levsp
#ifdef BACKWARD
   call dfs_global_mean( V3, XMAT,mt,ngs,nge,  levsp, +1 )
#else
   call dfs_global_mean( D3, XMAT,mt,ngs,nge,2*levsp, +1 )
#endif
!
   do mj = 1,jlg
     do mx = 1,mt
       PSL3(mx,mj)= PSL1(mx,mj) + DT2*(PSJ(mx,mj)-phsdm(mx,mj))
     enddo
   enddo
!
#ifndef NISLQ
   do k = 1,levhp
     do mj = 1,JLg
       do mx = 1,MT
         Q3(mx,mj,k)= Q1(mx,mj,k) + DT2* QJ(mx,mj,k)
       enddo
     enddo
   enddo ! k=1,levhp
#endif
!
   RETURN
   END SUBROUTINE DFS_DIFFUSION
!
!-------------------------------------------------------------------------------
   subroutine dfs_diffusion_coef(delt,gamma,mta,opt)
!-------------------------------------------------------------------------------
   use constant, only   :   rerth_,pi_
   use dfsvar, only     :   lap_dim,iope
!-------------------------------------------------------------------------------
   implicit none
!
   real                 ::  delt
   real                 ::  gamma(3)
   integer              ::  mta
   integer              ::  opt
!
! local
!
   real                 ::  pi2,awav,fac1,fac2,fac3
!-------------------------------------------------------------------------------
   pi2=pi_*2.
!
!-----------------------------------------------   PAI  delt  -----!
!----for filtering with Real 5-diag matrix
!       gamma0= (mta*(mta+1.0D0))**LAP_DIM
!       gamma0=  delt*2/(0.1D0*PI2*gamma0)
!----for filtering with Complex 3-diag matrix
!
!awav=mta*2/3
!awav=mta*0.8
#ifdef BARO_TEST
   awav=42
   gamma(1)= ( (awav*(awav+1.0D0))/(rerth_**2) )**LAP_DIM
   gamma(1)=  delt*2/(8640.00D0*gamma(1))
   gamma(1)=  gamma(1)**(1.D0/LAP_DIM)
#ifdef DBG
   if (iope) then
     write(6,'(A,2E12.5,I3)')'gamma(1),delt,LAP_DIM=',                         &
                              gamma(1),delt,LAP_DIM
   endif
#endif
#else
   if (opt.eq.1) then          ! for vorticity
     fac1=0.75
     fac2=0.65
     fac3=0.37
   elseif (opt.eq.2) then      ! for divergence
     fac1=0.75
     fac2=0.65
     fac3=0.37
   elseif (opt.eq.3) then      ! for temperture
     fac1=0.75
     fac2=0.65
     fac3=0.37
   elseif (opt.eq.4) then      ! for specific humidity
     fac1=0.75
     fac2=0.65
     fac3=0.37
   endif
#ifdef DCMIP
   fac1=0.65
#endif
#ifdef DFS_2009
   fac1=fac1+0.05
#endif
!
   awav=mta*fac1
   gamma(1)= ( (awav*(awav+1.0D0))/(rerth_**2) )**LAP_DIM
   gamma(1)=  delt*2/(8640.00D0*gamma(1))
   gamma(1)=  gamma(1)**(1.D0/LAP_DIM)
!
   awav=mta*fac2
   gamma(2)= ( (awav*(awav+1.0D0))/(rerth_**2) )**LAP_DIM
   gamma(2)=  delt*2/(8640.00D0*gamma(2))
   gamma(2)=  gamma(2)**(1.D0/LAP_DIM)
!
   awav=mta*fac3
   gamma(3)= ( (awav*(awav+1.0D0))/(rerth_**2) )**LAP_DIM
   gamma(3)=  delt*2/(8640.00D0*gamma(3))
   gamma(3)=  gamma(3)**(1.D0/LAP_DIM)
!
#ifdef DBG
   if (iope) then
     write(6,'(A,4E12.5,I3)')'gamma(1),gamma(2),gamma(3),delt,LAP_DIM=',       &
                               gamma(1),gamma(2),gamma(3),delt,LAP_DIM
     call flush(6)
   endif
#endif
#endif
!
!----for filtering with Complex 3-diag matrix
!

!                   filtcoef(0)= dcmplx(  1.d0,0.d0)
!                   filtcoef(2)= dcmplx(gamma2,0.d0)
!                   filtcoef(4)= dcmplx(gamma4,0.d0)
!
!      CALL ZROOTS( filtcoef, LAP_DIM, roots, .TRUE. )
!
  return
  end subroutine dfs_diffusion_coef
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_diffusion_del_matrix ( zalpha, ZDMATmi_ )
!-------------------------------------------------------------------------------
   use dfsvar, only : mt,jlha,midx,                                            &
                      aMSQUAR,AMATm,DMATm
!-------------------------------------------------------------------------------
!                                                                    
!     To solve for psi; COMPLEX VARIABLES                            
!     psi + alpha* ( LAPLACian of psi ) = X_i,j ,                    
!     which arises from the high-order spectral filter               
!                                                                    
!    [CALL ZMAT200S] -> [CALL ZLAP200S] : COMPLEX-variable version   
!                                                                    
!     01 MAR 2005                                                    
!                                                                    
!-------------------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------!
   COMPLEX*16  ::  zalpha,ZDMATmi_(MT,3,0:jlha,2)
!-----local variable-------------------------------------------------!
   COMPLEX*16  ::  ZYMATD(MT,3,0:jlha),zrat
   integer     ::  js,k,i,mx,mc,mc1
   real        ::  amx2
!-------------------------------------------------------------------------------
   if(zalpha.eq.0.d0) then
     write(*,*) 'The zalpha is zero: SUB ZMAT200S_VEC'
     call exit(99)
   endif
!
!************************************************!
!BEG  MATRIECS for LAPLACIAN OPERATOR            * SEMI-IMPLICIT
!************************************************!
!
!beg-- Matrices for Y= ( D^2/Dlambda^2 + COS(LAT) D/Dphi COS(LAT) D/Dphi ) X
!
   do js = 1,2
     do k = 0,jlha
       do i = 1,3
         do mx = 1,MT
           ZDMATmi_(mx,i,K,js)= AMATm(mx,i,K,js)+zalpha*DMATm(mx,i,K,js)
         enddo
       enddo
     enddo
   enddo
!
!************************************************!
!END  MATRIECS for LAPLACIAN OPERATOR            * SEMI-IMPLICIT
!************************************************!
!
!
!
!************************************************!
!BEG  Time saving, prehandling of the MATRIX     * SEMI-IMPLICIT
!************************************************!
!
!beg-------------------------------------------------- DMAT -> ZDMAT_
!-----m0,m1,m2,(n1) and m0,m1,m2,(n2)
!
   do js = 1,2
!
     do mc = 1,jlha-1
       do mx = 1,MT
         amx2= -aMSQUAR(mx)
         ZYMATD( mx,1,mc)= ZDMATmi_(mx,1,mc,js)
         ZYMATD( mx,2,mc)= ZDMATmi_(mx,2,mc,js) - amx2*zalpha  !  mx2<0
         ZYMATD( mx,3,mc)= ZDMATmi_(mx,3,mc,js)
       enddo
     enddo
     mc=0
     do mx = 1,MT
       amx2= -aMSQUAR(mx)
       ZYMATD( mx,1,mc)= ZDMATmi_(mx,1,mc,js) - amx2*zalpha  !
       ZYMATD( mx,2,mc)= ZDMATmi_(mx,2,mc,js)
       ZYMATD( mx,3,mc)=0.
     enddo
     mc=jlha
     do mx = 1,MT
       amx2= -aMSQUAR(mx)
       ZYMATD( mx,1,mc)= ZDMATmi_(mx,1,mc,js)
       ZYMATD( mx,2,mc)= ZDMATmi_(mx,2,mc,js) - amx2*zalpha  !
       ZYMATD( mx,3,mc)=0.
     enddo
!-------------------
     do mc = 0,jlha-1  
       mc1= mc+1   
       do mx = 1,MT
!
         IF(mc.eq.0) then
           if(midx(mx).eq.0 .AND. js.eq.2) then
             zrat= ZYMATD( mx,1,mc1) / ZYMATD( mx,1,mc)
           else
             zrat= 0.d0
           endif
         ELSE
           zrat= ZYMATD( mx,1,mc1) / ZYMATD( mx,1,mc)
         ENDIF
!
         ZYMATD( mx,2,mc1)= ZYMATD( mx,2,mc1)-ZYMATD( mx,2,mc )*zrat   
         ZYMATD( mx,1,mc1)= ZYMATD( mx,2,mc1)
         ZYMATD( mx,2,mc1)= ZYMATD( mx,3,mc1)
!
         IF(mc.eq.0) then
           if(midx(mx).eq.0 .AND. js.eq.2) then
             ZDMATmi_( mx,1,mc ,js)= 1.d0/ZYMATD( mx,1,mc )
           else
             ZDMATmi_( mx,1,mc ,js)= 0.d0
           endif
         ELSE
           ZDMATmi_( mx,1,mc ,js)= 1.d0/ZYMATD( mx,1,mc )
         ENDIF
!
         ZDMATmi_( mx,2,mc ,js)= ZYMATD( mx,2,mc )       
         ZDMATmi_( mx,3,mc ,js)= zrat
       enddo ! par
     enddo ! mc

     do mx = 1,MT
       ZDMATmi_( mx,1,jlha,js)= 1.d0/ZYMATD( mx,1,jlha)
     enddo ! par
!
   enddo ! js ! js=1,2
!end-------------------------------------------------- DMAT -> DMAT_
!
!************************************************!
!END  Time saving, prehandling of the MATRIX     * SEMI-IMPLICIT
!************************************************!
!
   RETURN
   END SUBROUTINE DFS_DIFFUSION_DEL_MATRIX
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_diffusion_driver( V3, D3, gamma, LAP_DIM, sl, opt)
!-------------------------------------------------------------------------------
!
! abstract: VISCOSITY-FILTERING(SMOOTHING) BY INVERSION OF ELLIPTIC EQS.
!           with COMPLEX TRIDIAGONAL SYSTEM
!
! [1-LAP], [1+LAP^2], [1-LAP^3], [1+LAP^4]
!
! remarks : Optimization must be done by inputting two variables
!
! 01 MAR 2005
!
! program history log:
!   2000-01-01  heongbin cheong        initial development
!   2005-01-01  hoon park              grims implementation
!   2006-04-01  hoon park              mpi
!   2010-12-21  myung-seo koo          3-layer & variable seperated
!
!-------------------------------------------------------------------------------
   use constant, only : pi=>pi_
   use dfsvar, only   : mt,jla,jlha,pi2,iope                                  ,&
                        ZDMATmi_a_RE_v,ZDMATmi_a_IM_v                         ,&
                        ZDMATmi_a_RE_d,ZDMATmi_a_IM_d                         ,&
                        ZDMATmi_a_RE_t,ZDMATmi_a_IM_t                         ,&
                        ZDMATmi_a_RE_q,ZDMATmi_a_IM_q
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer          ::  LAP_DIM, opt
   real             ::  gamma(3),sl
   real             ::  V3(1:MT,0:jla),D3(1:MT,0:jla)
!
! local variables
!
   real             ::  AJX_RE(1:MT,0:jla)
   real             ::  AJX_IM(1:MT,0:jla),gamm(LAP_DIM*3)
   COMPLEX*16       ::  zalpha(LAP_DIM*3),ZZDMATmi_a(1:MT,3,0:jlha,2)
   real,dimension(:,:,:,:,:),pointer  ::  ZDMATmi_a_RE, ZDMATmi_a_IM
!
   real, pointer    ::  gam_old
   real, target     ::  gam_old_v, gam_old_d, gam_old_t, gam_old_q
   DATA gam_old_v/0.d0/, gam_old_d/0.d0/, gam_old_t/0.d0/, gam_old_q/0.d0/
   integer          ::  kl,k2,mx,mj,i,k,k3,k1,lapend,lapst
   real             ::  xx,ccc,sss
   real             ::  gam_new
!-------------------------------------------------------------------------------
#include "abort.h"
!
   IF(LAP_DIM.gt.4) then
     if (iope)write(6,*)' in dfs_diffusion_driver:                             &
                          LAP_DIM must not be larger than 4'
     call MPABORT
   endif
!
   gam_new=gamma(1)
   if (opt.eq.1) then
     gam_old=>gam_old_v
     ZDMATmi_a_RE=>ZDMATmi_a_RE_v
     ZDMATmi_a_IM=>ZDMATmi_a_IM_v
   elseif (opt.eq.2) then
     gam_old=>gam_old_d
     ZDMATmi_a_RE=>ZDMATmi_a_RE_d
     ZDMATmi_a_IM=>ZDMATmi_a_IM_d
   elseif (opt.eq.3) then
     gam_old=>gam_old_t
     ZDMATmi_a_RE=>ZDMATmi_a_RE_t
     ZDMATmi_a_IM=>ZDMATmi_a_IM_t
   elseif (opt.eq.4) then
     gam_old=>gam_old_q
     ZDMATmi_a_RE=>ZDMATmi_a_RE_q
     ZDMATmi_a_IM=>ZDMATmi_a_IM_q
   endif
!
   IF(gam_old.NE.gam_new) THEN
     gam_old=gam_new        
     gamm(1:LAP_DIM*3)=gam_old
     lapend=LAP_DIM
     if(LAP_DIM.EQ.4) then
       zalpha(1)= -dcmplx(+1.d0,+1.d0)/DSQRT(2.d0)
       zalpha(2)= -dcmplx(-1.d0,+1.d0)/DSQRT(2.d0)
       zalpha(3)= -dcmplx(-1.d0,-1.d0)/DSQRT(2.d0)
       zalpha(4)= -dcmplx(+1.d0,-1.d0)/DSQRT(2.d0)
#ifndef BARO_TEST
       lapend=LAP_DIM*3
       gamm(5:8)=gamma(2)
       zalpha(5)= -dcmplx(+1.d0,+1.d0)/DSQRT(2.d0)
       zalpha(6)= -dcmplx(-1.d0,+1.d0)/DSQRT(2.d0)
       zalpha(7)= -dcmplx(-1.d0,-1.d0)/DSQRT(2.d0)
       zalpha(8)= -dcmplx(+1.d0,-1.d0)/DSQRT(2.d0)
       gamm(9:12)=gamma(3)
       zalpha(9)= -dcmplx(+1.d0,+1.d0)/DSQRT(2.d0)
       zalpha(10)= -dcmplx(-1.d0,+1.d0)/DSQRT(2.d0)
       zalpha(11)= -dcmplx(-1.d0,-1.d0)/DSQRT(2.d0)
       zalpha(12)= -dcmplx(+1.d0,-1.d0)/DSQRT(2.d0)
#endif
     elseif(LAP_DIM.EQ.8) then
       do k = 1,8
         xx=(PI+(k-1)*PI2)/8.d0
#ifdef IBMSP
         ccc= COS(dble(xx)) 
         sss= SIN(dble(xx)) 
#else
         ccc= DCOS(dble(xx)) 
         sss= DSIN(dble(xx)) 
#endif
         zalpha(k)= -dcmplx(ccc,sss)
       enddo
     elseif(LAP_DIM.EQ.3) then
       zalpha(1)= -dcmplx(+1.d0, 0.d0       )      
       zalpha(2)= -dcmplx(-1.d0, DSQRT(3.d0))*0.5
       zalpha(3)= -dcmplx(-1.d0,-DSQRT(3.d0))*0.5
#ifndef BARO_TEST
       lapend=LAP_DIM*3
       gamm(4:6)=gamma(2)
       zalpha(4)= -dcmplx(+1.d0, 0.d0       )
       zalpha(5)= -dcmplx(-1.d0, DSQRT(3.d0))*0.5
       zalpha(6)= -dcmplx(-1.d0,-DSQRT(3.d0))*0.5
       gamm(7:9)=gamma(3)
       zalpha(7)= -dcmplx(+1.d0, 0.d0       )
       zalpha(8)= -dcmplx(-1.d0, DSQRT(3.d0))*0.5
       zalpha(9)= -dcmplx(-1.d0,-DSQRT(3.d0))*0.5
#endif
     elseif(LAP_DIM.EQ.2) then
       zalpha(1)= -dcmplx( 0.d0,+1.D0)      
       zalpha(2)= -dcmplx( 0.d0,-1.D0)
     elseif(LAP_DIM.EQ.1) then
       zalpha(1)= -dcmplx( 1.d0, 0.D0)      
     endif
!
     do i = 1,lapend 
       zalpha(i)=  zalpha(i)*gamm(i)
       CALL dfs_diffusion_del_matrix( zalpha(i), ZZDMATmi_a(1,1,0,1) ) 
       do k3 = 1,2
         do k2 = 0,jlha
           do k1 = 1,3
             do mx = 1,MT
               ZDMATmi_a_RE(mx,k1,k2,k3,i)= REal(ZZDMATmi_a(mx,k1,k2,k3))
               ZDMATmi_a_IM(mx,k1,k2,k3,i)= IMag(ZZDMATmi_a(mx,k1,k2,k3))
             enddo
           enddo
         enddo
       enddo
     enddo ! i=1,LAP_DIM
!
   ENDIF ! JUSTDOIT
!
!-beg--filtering---------
!
   lapend=LAP_DIM
   lapst=1
#ifndef BARO_TEST
   if (sl.le.0.1) then
     lapst =2*LAP_DIM+1
     lapend=3*LAP_DIM
   elseif (sl.le.0.7) then
     lapst =  LAP_DIM+1
     lapend=2*LAP_DIM
   endif
#endif
!
   do mj = 0,jla
     do mx = 1,MT
       AJX_RE(mx,mj)= V3(mx,mj)
       AJX_IM(mx,mj)= D3(mx,mj)
     enddo
   enddo
!
   do i=lapst,lapend
     CALL dfs_diffusion_laplacian( AJX_RE, AJX_IM,                             &
                        ZDMATmi_a_RE(1,1,0,1,i), ZDMATmi_a_IM(1,1,0,1,i) )
   enddo ! i=1,LAP_DIM
!
   do mj = 0,jla
     do mx = 1,MT
       V3(mx,mj)= AJX_RE(mx,mj)
       D3(mx,mj)= AJX_IM(mx,mj)
     enddo
   enddo
!
!-end--filtering---------
!
   RETURN
   END SUBROUTINE DFS_DIFFUSION_DRIVER
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_diffusion_laplacian                                          &
                            ( AJX_RE, AJX_IM, ZDMATmi_RE, ZDMATmi_IM )
!-------------------------------------------------------------------------------
   use dfsvar, only : mt,jla,jlha,jcol2js,AMATm,mls
!-------------------------------------------------------------------------------
!                                                                    
!     To get psi or Div from : COMPLEX VARIABLES                     
!                                                                    
!     psi + zalpha * ( Laplacian of psi ) = Y  : Global Ave[Y] !=0   
!     Div + zalpha * ( Laplacian of Div ) = Z  : Global Ave[Z]  =0   
!                                               because [Div]  =0    
!    [CALL ZMAT200S] -> [CALL ZLAP200S] : COMPLEX-variable version   
!                                                                    
!     01 MAR 2005                                                    
!                                                                    
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real     ::  AJX_RE(MT,0:JLA),ZDMATmi_RE(MT,3,0:JLHA,2)
   real     ::  AJX_IM(MT,0:JLA),ZDMATmi_IM(MT,3,0:JLHA,2)
   real     ::  ZYMAT_RE(MT,0:JLHA),ZYMAT_IM(MT,0:JLHA)
   real     ::  zaaar,ajx0r,tc10r,parRE,zaaai,ajx0i,tc10i,parIM,dd
   real     ::  a4,a5,a6
   integer  ::  mj,js,jcol,jja,jjb,jjc,mx,jjj
   integer  ::  mc1,mc,mco,mcop1,mzr
!-------------------------------------------------------------------------------
   if (mls.eq.0) then
     mzr=mt/2+1
     ajx0r= 0.d0
     ajx0i= 0.d0
     do mj = 0,JLA,2
       dd= DFLOAT(mj)**2-1.D-00
       ajx0r= ajx0r - AJX_RE(mzr,mj)/dd
       ajx0i= ajx0i - AJX_IM(mzr,mj)/dd
     enddo
   endif
!
!-------------- mx>0, mx<0 ------------c
!
   do js = 1,2
!
!---- 1  To get YMAT, YMAT= AMAT*ZETA
!
     do jcol = 1,JLHA-1
       jja= jcol2js(jcol-1,js)
       jjb= jcol2js(jcol  ,js)
       jjc= jcol2js(jcol+1,js)
       do mx = 1,mt
         a4= AMATm(mx,1,jcol,js)
         a5= AMATm(mx,2,jcol,js)
         a6= AMATm(mx,3,jcol,js)
         ZYMAT_RE(mx,jcol)= a4*AJX_RE(mx,jja)+a5*AJX_RE(mx,jjb)                &
                           +a6*AJX_RE(mx,jjc)
         ZYMAT_IM(mx,jcol)= a4*AJX_IM(mx,jja)+a5*AJX_IM(mx,jjb)                &
                           +a6*AJX_IM(mx,jjc)
       enddo ! par
     enddo
     jcol=0
     jjb= jcol2js(jcol  ,js)
     jjc= jcol2js(jcol+1,js)
     do mx = 1,mt
       a4= AMATm(mx,1,jcol,js)
       a5= AMATm(mx,2,jcol,js)
       ZYMAT_RE(mx,jcol)= a4*AJX_RE(mx,jjb)+a5*AJX_RE(mx,jjc)                   
       ZYMAT_IM(mx,jcol)= a4*AJX_IM(mx,jjb)+a5*AJX_IM(mx,jjc)                   
     enddo ! par
     jcol=JLHA
     jja= jcol2js(jcol-1,js)
     jjb= jcol2js(jcol  ,js)
     do mx = 1,mt
       a4= AMATm(mx,1,jcol,js)
       a5= AMATm(mx,2,jcol,js)
       ZYMAT_RE(mx,jcol)= a4*AJX_RE(mx,jja)+a5*AJX_RE(mx,jjb)
       ZYMAT_IM(mx,jcol)= a4*AJX_IM(mx,jja)+a5*AJX_IM(mx,jjb)
     enddo ! par
!
!---- 2  To get PSI, DMAT*PSI= YMAT
!
     do mc = 0,JLHA-1 
       mc1= mc+1
#ifdef LINUX_INTEL
!DEC$ NOVECTOR
#endif
       do mx = 1,mt
         ZYMAT_RE(mx,mc1)=  ZYMAT_RE(mx,mc1)                                   &
                           -ZYMAT_RE(mx,mc )*ZDMATmi_RE(mx,3,mc,js)            &
                           +ZYMAT_IM(mx,mc )*ZDMATmi_IM(mx,3,mc,js)
         ZYMAT_IM(mx,mc1)=  ZYMAT_IM(mx,mc1)                                   &
                           -ZYMAT_RE(mx,mc )*ZDMATmi_IM(mx,3,mc,js)            &
                           -ZYMAT_IM(mx,mc )*ZDMATmi_RE(mx,3,mc,js)
       enddo ! par
     enddo
!
     do mx = 1,mt
       zaaar =  ZYMAT_RE(mx,JLHA)*ZDMATmi_RE(mx,1,JLHA,js)                     &
               -ZYMAT_IM(mx,JLHA)*ZDMATmi_IM(mx,1,JLHA,js) !! (1/all)
       zaaai =  ZYMAT_RE(mx,JLHA)*ZDMATmi_IM(mx,1,JLHA,js)                     &
               +ZYMAT_IM(mx,JLHA)*ZDMATmi_RE(mx,1,JLHA,js) !! (1/all)
       ZYMAT_RE(mx,JLHA)=  zaaar
       ZYMAT_IM(mx,JLHA)=  zaaai
     enddo ! par
!
     do mco = JLHA-1,1,-1
       mcop1=mco+1
       do mx = 1,mt
         zaaar=  ZYMAT_RE(mx,mcop1)*ZDMATmi_RE(mx,2,mco,js)                    &
                -ZYMAT_IM(mx,mcop1)*ZDMATmi_IM(mx,2,mco,js)
         zaaai=  ZYMAT_RE(mx,mcop1)*ZDMATmi_IM(mx,2,mco,js)                    &
                +ZYMAT_IM(mx,mcop1)*ZDMATmi_RE(mx,2,mco,js)
         parRE=  ZYMAT_RE(mx,mco)-zaaar
         parIM=  ZYMAT_IM(mx,mco)-zaaai
         ZYMAT_RE(mx,mco)= parRE*ZDMATmi_RE(mx,1,mco,js)                       &
                          -parIM*ZDMATmi_IM(mx,1,mco,js)
         ZYMAT_IM(mx,mco)= parRE*ZDMATmi_IM(mx,1,mco,js)                       &
                          +parIM*ZDMATmi_RE(mx,1,mco,js)
       enddo ! par
     enddo
!
     do jcol = 1,JLHA
       jjj= jcol2js(jcol,js)
       do mx = 1,mt
         AJX_RE(mx,jjj)= ZYMAT_RE(mx,jcol)
         AJX_IM(mx,jjj)= ZYMAT_IM(mx,jcol)
       enddo ! par
     enddo
!
   enddo ! js=1,2
!
!----AJX->TC1      m=0, Y-even    ( mx=0, js=0 )
!
   if (mls.eq.0) then
     AJX_RE( mzr, 0 )= 0.d0
     AJX_IM( mzr, 0 )= 0.d0
!
     tc10r= 0.d0
     tc10i= 0.d0
     do mj = 0,JLA,2
       dd= DFLOAT(mj)**2-1.D-00
       tc10r= tc10r - AJX_RE(mzr,mj)/dd
       tc10i= tc10i - AJX_IM(mzr,mj)/dd
     enddo
!
     AJX_RE( mzr, 0 )= ajx0r - tc10r   ! global-mean correction
     AJX_IM( mzr, 0 )= ajx0i - tc10i   ! global-mean correction
   endif
!
   RETURN
   END SUBROUTINE dfs_diffusion_laplacian
!-------------------------------------------------------------------------------
