#include <define.h>
   subroutine dfs_advection
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [dfs_advection]
!           |
!           |-- [dfs_advection_y] *
!           |-- [dfs_dynamics_advection] *
!                   |-- [file_write_wave] *
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
   end subroutine dfs_advection
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_advection_y( SANMV, SANM2 ,mt, JLa ,lot)
!-------------------------------------------------------------------------------
   use constant, only : rrerth_
   use dfsvar, only   : mls,midx,nev,nod,mevs,mods
!-------------------------------------------------------------------------------
!                                                                    
!     Y (:SANM2) = 1/COS(LAT))*(D/DPhi) X (:SANMV)                   
!                                                                    
!     Assumes meridional coefficients are arranged in global         
!                                                                    
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  mt,jla,lot
   integer  ::  mj,mjp1,mjm1,mx,mjp2,mjm2,k,m,mx0
   real     ::  SANM2(mt,0:jla,lot),SANMV(mt,0:(jla+1),lot)
!-------------------------------------------------------------------------------
   do k = 1,lot
!
! set zero at n=0
!
     SANM2(1:mt,0,k)=0.
!
!-m=1,3,...
!
     do mj = 1,jla             ! (1/COS(LAT))*(D/DPhi) (v*COS(LAT)*Vor)
       mjp1=mj+1
       mjm1=mj-1
       do m = 1,nod
         mx=mods(m)
         SANM2(mx,mj,k)= ((mjp1)*SANMV(mx,mjm1,k)-                             &
                     (mjm1)*SANMV(mx,mjp1,k))*0.5 
       enddo ! par
     enddo
!
!-m=2,4,...
!
     do mj = 1,jla             ! (1/COS(LAT))*(D/DPhi) (v*COS(LAT)*Vor)
       mjp1=mj+1
       mjm1=mj-1
       mjp2=mj+2
       mjm2=mj-2
       do m = 1,nev
         mx=mevs(m)
         if (midx(mx).ne.0) then
           SANM2(mx,mj,k)=  ((mjp2)*SANMV(mx,mjm1,k)-                          &
                         (mjm2)*SANMV(mx,mjp1,k))*0.5 
         else             ! mx=0  
           mx0=mx
           SANM2(mx,mj,k)=  ((mjp1)*SANMV(mx,mjm1,k)-                          &
                         (mjm1)*SANMV(mx,mjp1,k))*0.5 
         endif
       enddo ! par
     enddo
!
!-m=0
!
     if (mls.eq.0) then
       SANM2(mx0,0,k )=         SANMV(mx0,0+1,k) *0.5 
       SANM2(mx0,1,k )=  ( 1+1)*SANMV(mx0,  0,k)                      
     endif
!
! Dimensionalize by Earth's parameter
! 
     do mj = 0,jla
#ifdef LINUX_PGI
!pgi$ novector
#endif
       do mx = 1,mt
         SANM2(mx,mj,k)= -SANM2(mx,mj,k)*rrerth_     ! (-) for N.P starting
       enddo
     enddo
   enddo   ! k
!
!
   RETURN
   END SUBROUTINE dfs_advection_y
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE dfs_dynamics_advection                                           &
#ifdef HYBRID
              ( V1,D1,VOR,DIO,TAI,QAI,PRS,XPSL,YPSL,APSN                      ,&
                spdmx,VJ,DJ,TJ,QJ,PSJ,ak5,bk5,kdt                             ,&
#ifndef NISLQ
                deltim)
#else
                deltim,pdot)
#endif
#else /* SIGMA */
              ( V1,D1,VOR,DIO,TAI,QAI,XPSL,YPSL,APSN,                          &
                spdmx,VJ,DJ,TJ,QJ,PSJ,kdt)
#endif /* HYBRID end */
!-------------------------------------------------------------------------------
   use constant, only      : akap=>akapa_,omega_,rd_,rrerth_                   
#ifdef HYBRID
   use constant, only      : cvap_,cp_
#endif
   use dfsvar, only        : mt,jlg,levsp,ib,jbw,levhp,jgs,jls,jle,vxc,vyc    ,&
                             jla,jlha,levs,levh,ntotal,midx,mop,sinlat,coslat ,&
#ifdef HDSZ
                             T_VIS,t_ref,z_vis,                                &
#endif
#ifdef CLM_CWF
                             cgs,                                              &
#endif
                             TAVEXY,latdef,ncldg                              ,&
#ifndef HYBRID
                             gama1,gama2,sigdot,divn,apsk,adt1,delsig         ,&
#endif
                             amatm,AMATm_,dmatm,DMATm_,amsquar,jcol2js
#ifdef DCMIP
   use dcmip_grims, only : icase,pwi,ptop
#endif
!-------------------------------------------------------------------------------
!                                                                    
!     ADVECTION OF VORTICITY, DIVERGENCE AND TEMPERATURE             
!                                                                    
! if ntotal=2,
!    nwg =4
!    ndel=5
!    nw0 =10
!    nw1 =8
!
!  ---------------------------------------------------------
!  GRID [ib,jbw,levs]    ->  WAVE [mt,jlg+1,levsp]
!  ---------------------------------------------------------
!  workg1(nw1)  AD1X         work1a(nw1)  AD1X
!               AD1Y                      AD1Y
!               AD2Y                      AD2Y
!               AD3Y-1                    AD3Y-1
!               AD3Y-2                    AD3Y-2
!               AD2X                      AD2X
!               AD3X-1                    AD3X-1
!               AD3X-2                    AD3X-2
!  ---------------------------------------------------------
!  GRID [ib,jbw,levs]    ->  WAVE [mt,jlg  ,levsp]
!  ---------------------------------------------------------
!  workg0(nwg)  AD1Z     ->  work0a(nw0)  AD1Z
!               AD2Z                      AD2Z
!               AD3Z-1                    AD3Z-1
!               AD3Z-2                    AD3Z-2
!                                         d(AD1X)  /dy
!                                         d(AD1Y)  /dy
!                                         d(AD2Y)  /dy
!                                         d(AD3Y-1)/dy
!                                         d(AD1X-2)/dy
!                                         laplacian of AD1Z
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
   integer                                        ::  kdt
#ifdef HYBRID
   real                                           ::  deltim
   real   , dimension(ib,jbw)      , intent(in)   ::  prs
   real   , dimension(levs+1)      , intent(in)   ::  ak5,bk5
#endif
   real   , dimension(mt,jlg,levsp), intent(in)   ::  v1,d1
   real   , dimension(ib,jbw,levs) , intent(in)   ::  vor,dio,tai
   real   , dimension(ib,jbw,levh)                ::  qai
   real   , dimension(ib,jbw)      , intent(in)   ::  xpsl,ypsl,apsn
   real   , dimension(levs)        , intent(out)  ::  spdmx
   real   , dimension(mt,jlg)      , intent(out)  ::  PSJ
   real   , dimension(mt,jlg,levsp), intent(out)  ::  VJ,DJ,TJ
   real   , dimension(mt,jlg,levhp), intent(out)  ::  QJ
#ifdef NISLQ
   real   , dimension(ib,jbw,levs+1),intent(out)  ::  pdot
#endif
!
! local variables
!
#ifdef CLM_CWF
   real   , dimension(mt,jlg,levsp)               ::  wcon
#endif
   real   , dimension(:,:,:), target, allocatable ::  workg1,workg0
   real   , dimension(:,:,:), target, allocatable ::  work1a,work0a
   real   , dimension(:,:,:), pointer             ::  SANMU,SANM1             ,&
                                                      SANMV,SANM2,SANM3
   integer                                        ::  k,n1,nk,i,j,mj,mx,kk    ,&
                                                      kdum,jj,jb,m,nv         ,&
#ifdef HYBRID
                                                      js                      ,&
#else
                                                      en1,enk                 ,&
#endif
                                                      lotg1,lotg0,lotw1,lotw0 ,&
                                                      nw0,nwg,ndel,nw1,ksdel
   real                                           ::  vivor,vidiv,adiv,e1     ,&
                                                      avor,wmax,UF3,UF2,UF1   ,&
#ifdef HYBRID
                                                      UF0,VF0                 ,&
#else
                                                      UF4,VF4                 ,&
#endif
                                                      WF0,VF3,VF2,VF1         ,&
                                                      ADVX,ADVY,ADVZ          ,&
                                                      atem,WF5,WF4,WF3,WF2,WF1,&
                                                      zeta,dup,dum,slat
#ifdef HYBRID
   integer                                        ::  n,kt,kt1,nvcn
   integer,save                                   ::  ifirst
   real                                           ::  cons0,cons1,cons2       ,&
                                                      cons0p5,xvcn,psi0
   real   ,save                                   ::  clog2,delta,delta1
#ifdef NISLQ
   real, dimension(ib,jbw,levs,3        )         ::  zadv
#else
   real, dimension(ib,jbw,levs,3+ntotal )         ::  zadv
#endif
   real, dimension(ib,jbw,levs+1)                 ::  pk5,dot,dotinv,si,psi2
   real, dimension(ib,jbw,levs)                   ::  dpk,rdel,rdel2          ,&
                                                      alfa,cofb,rlnp          ,&
                                                      cg,cb,db                ,&
                                                      worka,workb,workc       ,&
                                                      psi,tt,sl
   real, dimension(ib,jbw)                        ::  exprs,dqdt
   real, dimension(levs)                          ::  dbk,ck
   data ifirst/1/
#endif /* HYBRID end */
!-------------------------------------------------------------------------------
#ifdef HYBRID
! constant
!
   cons0   = 0.d0      !constant
   cons0p5 = 0.5d0     !constant
   cons1   = 1.d0      !constant
   cons2   = 2.d0      !constant
!
   if(ifirst.eq.1) then
     ifirst=0
     clog2=log(cons2)
     delta=cvap_/cp_
     delta1=delta-cons1
   endif
#endif /* HYBRID end */
!
! vertical layers
!
#ifdef NISLQ
   nwg=2
   ndel=3
   nw0=nwg+ndel+1
   nw1=ndel+1
   lotw0=nwg*levsp
   lotg0=nwg*levs
   lotw1=nw1*levsp
   lotg1=nw1*levs
#else
   nwg=ntotal+2
   ndel=ntotal+3
   nw0=nwg+ndel+1
   nw1=ndel+ntotal+1
   lotw0=nwg*levsp
   lotg0=nwg*levs
   lotw1=nw1*levsp
   lotg1=nw1*levs
#endif /* NISLQ end */
   jb=jbw/2
   allocate(workg1(ib,jbw,lotg1))
   allocate(workg0(ib,jbw,lotg0))
   allocate(work1a(mt,jlg+1,lotw1))
   allocate(work0a(mt,jlg,nw0*levsp))
!
! initialize
!
   workg1=0. ; workg0=0. ; work1a=0. ; work0a=0.
   zadv=0.   ; pk5=0.    ; dot=0.    ; dotinv=0.
   si=0.     ; psi2=0.   ; dpk=0.    ; rdel=0.
   rdel2=0.  ; alfa=0.   ; cofb=0.   ; rlnp=0.
   cg=0.     ; cb=0.     ; db=0.
   worka=0.  ; workb=0.  ; workc=0.
   psi=0.    ; tt=0.     ; sl=0.
   exprs=0.  ; dqdt=0.   ; dbk=0.    ; ck=0.
   vj=0.     ; dj=0.     ; tj=0.     ; qj=0.
   spdmx(1:levs)=0.
#ifdef NISLQ
   pdot=0.
#endif
#ifdef HYBRID
!
! exprs
!
   forall(i=1:ib,j=1:jbw) exprs(i,j)=exp(prs(i,j))
!
! [b2t] si
! [t2b] pk5
!
   do k = 1,levs+1
     do j = 1,jbw
       do i = 1,ib
         si(i,j,levs+2-k)=ak5(k)+bk5(k)*exprs(i,j)
         pk5(i,j,k)=ak5(k)+bk5(k)*exprs(i,j)
       enddo
     enddo
   enddo
!
! [b2t] sl
! [t2b] tt,dbk,ck,dpk,rdel,rdel2,cg
!
   do k = 1,levs
     dbk(k)=bk5(k+1)-bk5(k)
     ck(k)=ak5(k+1)*bk5(k)-ak5(k)*bk5(k+1)
     do j = 1,jbw
       do i = 1,ib
         sl(i,j,k)=cons0p5*(si(i,j,k)+si(i,j,k+1))
         tt(i,j,k)=tai(i,j,k)+tavexy(k)
         dpk(i,j,k)=pk5(i,j,k+1)-pk5(i,j,k)
         rdel(i,j,k)=cons1/dpk(i,j,k)
         rdel2(i,j,k)=cons0p5/dpk(i,j,k)
       enddo
     enddo
     do j = 1,jb
       jj=latdef(j+jls-jgs)
       js=jbw-j+1
       do i = 1,ib
         cg(i,j ,k)=( VXC(i,j ,levs+1-k)*xpsL(i,j )                            &
                      +VYC(i,j ,levs+1-k)*ypsL(i,j ) )*COSLAT(JJ,3) ! NH
         cg(i,js,k)=( VXC(i,js,levs+1-k)*xpsL(i,js)                            &
                      +VYC(i,js,levs+1-k)*ypsL(i,js) )*COSLAT(JJ,3) ! SH
       enddo
     enddo
   enddo
!
! [t2b] rlnp,alfa,cofb,cb,db
!
   k=1
   do j = 1,jbw
     do i = 1,ib
       rlnp(i,j,k)=99999.99
       alfa(i,j,k)=clog2
       cofb(i,j,k)=rdel(i,j,k)*(alfa(i,j,k)*dbk(k))
       db(i,j,k)=dio(i,j,levs)*dpk(i,j,k)
       cb(i,j,k)=cg(i,j,k)*dbk(k)
     enddo
   enddo
!
   do k = 2,levs
     do j = 1,jbw
       do i = 1,ib
         rlnp(i,j,k)=log( pk5(i,j,k+1)/pk5(i,j,k) )
         alfa(i,j,k)=cons1-( pk5(i,j,k)/dpk(i,j,k) )*rlnp(i,j,k)
         cofb(i,j,k)=rdel(i,j,k)*( bk5(k)*rlnp(i,j,k)+alfa(i,j,k)*dbk(k) )
         db(i,j,k)=db(i,j,k-1)+dio(i,j,levs+1-k)*dpk(i,j,k)
         cb(i,j,k)=cb(i,j,k-1)+cg(i,j,k)*dbk(k)
       enddo
     enddo
   enddo
!
! dqdt
! [t2b] dot,workb,workc
! [t2b?] psi2 (interfacial layer)
!
   k=1
   do j = 1,jbw
     do i = 1,ib
       dqdt(i,j)=-db(i,j,levs)/exprs(i,j)-cb(i,j,levs)
       !
       dot(i,j,1)        =cons0
       dot(i,j,levs+1)   =cons0
       dotinv(i,j,1)     =cons0
       dotinv(i,j,levs+1)=cons0
       !
       workb(i,j,k)=alfa(i,j,k)* ( dio(i,j,levs)*dpk(i,j,k)                    &
                                    + exprs(i,j)*cb(i,j,k)*dbk(k)  )
       workc(i,j,k)=exprs(i,j)*cg(i,j,k)*dbk(k)
       !
       psi2(i,j,levs+1)=cons0
       psi2(i,j,levs)=rd_*tt(i,j,1)*rlnp(i,j,levs)
     enddo
   enddo
!
#ifdef DCMIP
   do k = 2,levs
     do j = 1,jbw
       do i = 1,ib
         if(icase.eq.11 .or. icase.eq.12) then
           dot(i,j,k)=pwi(i,j,levs+2-k)*1.e-3
#ifndef HYBRID
           dot(i,j,k)=dot(i,j,k)/(exprs(i,j)-ptop*1.e-3)
#endif
         else if(icase.eq.13) then
           dot(i,j,k)=pwi(i,j,levs+2-k)*1.e-4
#ifndef HYBRID
           dot(i,j,k)=dot(i,j,k)/(exprs(i,j)-ptop*1.e-3)
#endif
         else
           dot(i,j,k)=-exprs(i,j)*( bk5(k)*dqdt(i,j)+cb(i,j,k-1) )-db(i,j,k-1)
         endif
       enddo
     enddo
   enddo
   if(icase.eq.13) then
     do n = 1,ntotal
       kk=(n-1)*levs
       do j = 1,jbw
         do i = 1,ib
           qai(i,j,kk+1)=max(qai(i,j,kk+1),cons0)
           qai(i,j,kk+2)=max(qai(i,j,kk+2),cons0)
           qai(i,j,kk+3)=max(qai(i,j,kk+3),cons0)
           qai(i,j,kk+4)=max(qai(i,j,kk+4),cons0)
           qai(i,j,kk+5)=max(qai(i,j,kk+5),cons0)
         enddo
       enddo
     enddo
   endif
#else
   do k = 2,levs
     do j = 1,jbw
       do i = 1,ib
         dot(i,j,k)=-exprs(i,j)*( bk5(k)*dqdt(i,j)+cb(i,j,k-1) )-db(i,j,k-1)
       enddo
     enddo
   enddo
#endif /* DCMIP end */
!
! [b2t] dotinv
!
   do k = 1,levs+1
     do j = 1,jbw
       do i = 1,ib
         dotinv(i,j,k)=dot(i,j,levs+2-k)
       enddo
     enddo
   enddo
#ifdef NISLQ
!
! [t2b] pdot: pressure velocity for nislq
!
   do k = 2,levs
     do j = 1,jbw
       do i = 1,ib
#ifdef DCMIP
         if(icase.eq.11 .or. icase.eq.12) then
           pdot(i,j,k)=pwi(i,j,levs+2-k)*1.e-3
         else if(icase.eq.13) then
           pdot(i,j,k)=pwi(i,j,levs+2-k)*1.e-4
         else
           pdot(i,j,k)=dot(i,j,k)
         endif
#else
         pdot(i,j,k)=dot(i,j,k)
#endif
       enddo
     enddo
   enddo
   pdot(1:ib,1:jbw,1)     =0.
   pdot(1:ib,1:jbw,levs+1)=0.
#endif /* NISLQ end */
!
   do k=2,levs
     do j=1,jbw
       do i=1,ib
         !
         workb(i,j,k)=rlnp(i,j,k)*( db(i,j,k-1)+exprs(i,j)*cb(i,j,k-1) )       &
                + alfa(i,j,k)* ( dio(i,j,levs+1-k)*dpk(i,j,k)                  &
                                +exprs(i,j)*cg(i,j,k)*dbk(k)   )
         workc(i,j,k)=exprs(i,j)*cg(i,j,k)                                     &
                              *( dbk(k)+ck(k)*rlnp(i,j,k)*rdel(i,j,k) )
         !
         psi2(i,j,levs+1-k)=psi2(i,j,levs+2-k)+rd_*tt(i,j,k)*rlnp(i,j,levs+1-k)
       enddo
     enddo
   enddo
!
! [t2b]  worka
! [t2b?] psi (model layer)
!
   do k = 1,levs
     do j = 1,jbw
       do i = 1,ib
         worka(i,j,k)=akap*tt(i,j,levs+1-k)                                    &
                      /(cons1+delta1*qai(i,j,levs+1-k))*rdel(i,j,k)
         psi(i,j,k)=psi2(i,j,k+1)+alfa(i,j,k)*rd_*tt(i,j,levs+1-k)
       enddo
     enddo
   enddo
!
! do vertical advection
!
!
! u,v,t
!
   do j = 1,jbw
     do i = 1,ib
       !k=1
       zadv(i,j,levs,1)=rdel2(i,j,1)*dot(i,j,2)*(vxc(i,j,levs-1)-vxc(i,j,levs))
       zadv(i,j,levs,2)=rdel2(i,j,1)*dot(i,j,2)*(vyc(i,j,levs-1)-vyc(i,j,levs))
       zadv(i,j,levs,3)=rdel2(i,j,1)*dot(i,j,2)*(tai(i,j,levs-1)-tai(i,j,levs))
       !k=levs
       zadv(i,j,1,1)=rdel2(i,j,levs)*dot(i,j,levs)*(vxc(i,j,1)-vxc(i,j,2))
       zadv(i,j,1,2)=rdel2(i,j,levs)*dot(i,j,levs)*(vyc(i,j,1)-vyc(i,j,2))
       zadv(i,j,1,3)=rdel2(i,j,levs)*dot(i,j,levs)*(tai(i,j,1)-tai(i,j,2))
     enddo
   enddo
   !
   do k = 2,levs-1
     kk=levs+1-k
     do j = 1,jbw
       do i = 1,ib
         zadv(i,j,kk,1)=                                                       &
               rdel2(i,j,k)*( dot(i,j,k+1)*( vxc(i,j,kk-1)-vxc(i,j,kk) ) +     &
                              dot(i,j,k  )*( vxc(i,j,kk  )-vxc(i,j,kk+1) ) )
         zadv(i,j,kk,2)=                                                       &
               rdel2(i,j,k)*( dot(i,j,k+1)*( vyc(i,j,kk-1)-vyc(i,j,kk) ) +     &
                              dot(i,j,k  )*( vyc(i,j,kk  )-vyc(i,j,kk+1) ) )
         zadv(i,j,kk,3)=                                                       &
               rdel2(i,j,k)*( dot(i,j,k+1)*( tai(i,j,kk-1)-tai(i,j,kk) ) +     &
                              dot(i,j,k  )*( tai(i,j,kk  )-tai(i,j,kk+1) ) )
       enddo
     enddo
   enddo
#ifndef NISLQ
!
! q
!
   do n = 1,ntotal
     kt=n*levs
     kt1=kt-levs+1
     !
     do j = 1,jbw
       do i = 1,ib
         !k=1
         zadv(i,j,levs,3+n)=rdel2(i,j,1)*dot(i,j,2)*(qai(i,j,kt-1)-qai(i,j,kt))
         !k=levs
         zadv(i,j,1   ,3+n)=rdel2(i,j,levs)*dot(i,j,levs)*(qai(i,j,kt1)-qai(i,j,kt1+1))
       enddo
     enddo
     !
     do k = 2,levs-1
       kk=levs+1-k
       do j = 1,jbw
         do i = 1,ib
           zadv(i,j,kk,3+n)= rdel2(i,j,k)*                                     &
                    ( dot(i,j,k+1)*( qai(i,j,kt  -k)-qai(i,j,kt+1-k) ) +       &
                      dot(i,j,k  )*( qai(i,j,kt+1-k)-qai(i,j,kt+2-k) ) )
         enddo
       enddo
     enddo
!
   enddo
#endif /* ~NISLQ end */
!
#ifndef DCMIP
#ifdef NISLQ
    call vcnhyb(ib*jbw,levs,3       ,deltim,si,sl,dotinv,zadv,nvcn,xvcn)
#else
    call vcnhyb(ib*jbw,levs,3+ntotal,deltim,si,sl,dotinv,zadv,nvcn,xvcn)
#endif
#endif /* ~DCMIP end */
!
#else  /* ~HYBRID */
#define ENK e1
#define EN1 e1
#endif /* HYBRID end */
!=============================  ADVECTION TERM OF VORTICITY  ===!
!=============================                and DIVERGENCE ===!
   do k = 1,levs
#ifndef HYBRID
     n1=1
     nk=1
     if(k.eq. 1  ) n1=0
     if(k.eq.levs) nk=0
!
! only for equally spaced levels
!
     e1= 0.5D0/ delsig(k)
!
! [dif_+0.5 +dif_-0.5]/2
!
     en1= 1.0d0/(delsig(k-n1)+delsig(k))
     enk= 1.0d0/(delsig(k+nk)+delsig(k))
#endif
     wmax=-10.
!
! -[ d ADVX/d lambda + d ADVY/ d the ]
!
     do j = 1,jbw
       if (j.gt.jb) then
         jj=latdef(jle-(j-jb-1))
         slat=-SINLAT(JJ)
       else
         jj=latdef(jls+j-1)
         slat= SINLAT(JJ)
       endif
       do i = 1,ib
!------------------
#ifdef DCMIP
         if(icase.eq.200) then
           ZETA =  VOR(I,J,k)
         else
           ZETA =  VOR(I,J,k)+slat*2.d0 *omega_
         endif
#else
         ZETA =  VOR(I,J,k)+slat*2.d0 *omega_
#endif /* DCMIP end */
!
#ifdef HYBRID
         UF0 = VXC(I,J,k)*ZETA
         UF1 = zadv(i,j,k,2)
         UF2 = cofb(i,j,levs+1-k)*rd_*tt(I,J,k)*exprs(i,j)*YPSL(I,J)
         ADVX= UF0+UF1+UF2
#else
         UF1 =  VXC(I,J,k)*ZETA
         UF2 = sigdot(i,j,k  )*(-VYC(I,J,k-n1)+VYC(I,J,k) )*EN1*(-1)  !(-1) up K
         UF3 = sigdot(i,j,k+1)*( VYC(I,J,k+nk)-VYC(I,J,k) )*ENK*(-1)  !(-1) up K
         UF4 = (TAI(I,J,k)*YPSL(I,J)) * rd_
!        ADVX   =   ( UF1+UF4 ) + e1*( UF2+UF3 )
         ADVX   =   ( UF1+UF4 ) +    ( UF2+UF3 )
#endif
!
#ifdef HYBRID
         VF0 = VYC(I,J,k)*ZETA
         VF1 = -zadv(i,j,k,1)
         VF2 = (-1.)*cofb(i,j,levs+1-k)*rd_*tt(I,J,k)*exprs(i,j)*XPSL(I,J)
         ADVY= VF0+VF1+VF2
#else
         VF1 =  VYC(I,J,k)*ZETA
         VF2 = sigdot(i,j,k  )*(-VXC(I,J,k-n1)+VXC(I,J,k) )*EN1*(-1)  !(-1) up K
         VF3 = sigdot(i,j,k+1)*( VXC(I,J,k+nk)-VXC(I,J,k) )*ENK*(-1)  !(-1) up K
         VF4 = -TAI(I,J,k)*XPSL(I,J)*rd_
         ADVY   =   ( VF1+VF4 ) -    ( VF2+VF3 )
#endif
!
#ifdef HYBRID
         WF0 = (VXC(I,J,k)**2+VYC(I,J,k)**2)* 0.5d0*coslat(jj,3)
         WF1 = psi(i,j,levs+1-k)
         ADVZ= (-1.)*(WF0+WF1)
#else
         ADVZ   =-(VXC(I,J,k)**2+VYC(I,J,k)**2)* 0.5d0*COSLAT(JJ,3)
#endif
!
         workg1(I,J,k      ) = ADVX * COSLAT(JJ,3)
         workg1(I,J,k+levs ) = ADVY * COSLAT(JJ,3)
         workg0(I,J,k      ) = ADVZ
#ifdef HYBRID
         wmax=max(wmax,WF0*2.0)
#else
         wmax=max(wmax,-workg0(I,J,k)*2.0)
#endif
       enddo
     enddo
     spdmx(k)=wmax
   enddo         ! k
!=============================  ADVECTION TERM OF TEMPERATURE ==!
   do k = 1,levs
#ifndef HYBRID
     n1=1
     nk=1
     if(k.eq. 1  ) n1=0
     if(k.eq.levs) nk=0
!
!  only for equally spaced levels
!
     e1=-0.5D0/ delsig(k)
!
! [dif_+0.5 +dif_-0.5]/2
!
     en1=-1.0d0/(delsig(k-n1)+delsig(k))
     enk=-1.0d0/(delsig(k+nk)+delsig(k))
!
!-------------------------------------------------------------!
!
#endif
     do j = 1,jbw
       if (j.gt.jb) then
         jj=latdef(jle-(j-jb-1))
       else
         jj=latdef(jls+j-1)
       endif
       do i = 1,ib
         ADVX   =  VXC(I,J,k)*TAI(I,J,k)
         ADVY   =  VYC(I,J,k)*TAI(I,J,k)
#ifdef HYBRID
         WF0    = TAI(I,J,k)* DIO(I,J,k)
         WF1    = -zadv(i,j,k,3)
         WF2    = worka(i,j,levs+1-k)*(workc(i,j,levs+1-k)-workb(i,j,levs+1-k))
         ADVZ   = WF0+WF1+WF2
#ifdef HDSZ
         ADVZ   = WF0 +WF1  &
                  -T_VIS(J,k)*( TAI(I,J,k)-(T_REF(I,J,k)-TAVEXY(k)) )
#endif
#else
         WF0    = TAI(I,J,k)*( DIO(I,J,k)-akap*DIVN(I,J) )
         WF1    = sigdot(I,J,k  )* (-TAI(I,J,k-n1)*gama1(k)+TAI(I,J,k) )*(-1)
         WF2    = sigdot(I,J,k+1)* ( TAI(I,J,k+nk)*gama2(k)-TAI(I,J,k) )*(-1)
         WF3    =   ADT1(I,J,k  )* (-TAVEXY(k-n1) *gama1(k)+TAVEXY(k) ) *(-1)
         WF4    =   ADT1(I,J,k+1)* ( TAVEXY(k+nk) *gama2(k)-TAVEXY(k) ) *(-1)
         WF5    = akap* (TAVEXY(k)+TAI(I,J,k))*(APSK(I,J,k)+APSN(I,J))
!        ADVZ   = (WF0+WF5) +e1*( (WF1+WF2)+(WF3+WF4) )
         ADVZ   = (WF0+WF5) +   ( (WF1+WF3)*EN1 + (WF2+WF4)*ENK )
#ifdef HDSZ
         ADVZ   = (WF0+WF5) +   ( (WF1+WF3)*EN1 + (WF2+WF4)*ENK )              &
                  -T_VIS(J,k)*( TAI(I,J,k)-(T_REF(I,J,k)-TAVEXY(k)) )
#endif
#endif /* HYBRID end */
!------------------
         workg1(I,J,k+ndel*levs) = ADVX * COSLAT(JJ,3)
         workg1(I,J,k+2*levs)    = ADVY * COSLAT(JJ,3)
         workg0(I,J,k+  levs)    = ADVZ
!------------------
       enddo
     enddo
   enddo
#ifndef NISLQ
!================================  ADVECTION TERM OF MOISTURE ==!
!
   do kk=1,ntotal
     do kdum=1,levs
       k=(kk-1)*levs+kdum
!--------------------------------------- for each level --------!
#ifndef HYBRID
       n1=1
       nk=1
       if(kdum.eq. 1  ) n1=0
       if(kdum.eq.levs) nk=0
!
! only for equally spaced levels
!
       e1=-0.5D0/ delsig(kdum)       
!
! [dif_+0.5 +dif_-0.5]/2
!
       en1=-1.0d0/(delsig(kdum-n1)+delsig(kdum))
       enk=-1.0d0/(delsig(kdum+nk)+delsig(kdum))
!
!-------------------------------------------------------------!
!
#endif
       do j = 1,jbw
         if (j.gt.jb) then
           jj=latdef(jle-(j-jb-1))
         else
           jj=latdef(jls+j-1)
         endif
         do I = 1,ib
           ADVX  =  VXC(I,J,kdum)*QAI(I,J,k)
           ADVY  =  VYC(I,J,kdum)*QAI(I,J,k)
!
#ifdef HYBRID
           WF0   =  QAI(I,J,k)* DIO(I,J,kdum)
           WF1   =  zadv(i,j,kdum,3+kk)
           ADVZ  =  WF0 - WF1
#else
           WF0   =  QAI(I,J,k)* DIO(I,J,kdum)
           WF1   = sigdot(i,j,kdum  )*(-QAI(I,J,k-n1)+QAI(I,J,k) )*EN1*(-1.)
           WF2   = sigdot(i,j,kdum+1)*( QAI(I,J,k+nk)-QAI(I,J,k) )*ENK*(-1.)
           ADVZ  = WF0 + WF1 + WF2
#endif
!------------------
           workg1(I,J,k+(NDEL+1)*levs) = ADVX * COSLAT(JJ,3)
           workg1(I,J,k+3*levs)        = ADVY * COSLAT(JJ,3)
           workg0(I,J,k+2*levs)        = ADVZ
!------------------
!
         enddo
       enddo
     enddo
   enddo         ! kk
#endif /* ~NISLQ end */
#ifndef HYBRID
#undef ENK
#undef EN1
#endif
!
   call dfs_fft_driver(1,workg1,ib,jbw,lotg1, work1a,mt,jlg+1,lotw1,jlg+1,     &
                        levsp,levs,nw1,coslat,1)
   call dfs_fft_driver(1,workg0,ib,jbw,lotg0,work0a,mt,jlg,lotw0,jlg,          &
                        levsp,levs,nwg,coslat,1)
   deallocate(workg1,workg0)
!
! Laplacian of K.E
!
   CALL dfs_divergence_laplacian1( work0a, work0a(1,1,(nw0-1)*levsp+1),        &
                     amatm,AMATm_,dmatm,DMATm_,amsquar,jcol2js,mt,jlha,levsp,0 )    
   CALL dfs_advection_y( work1a, work0a(1,1,nwg*levsp+1),MT,JLA,ndel*levsp)
! 
! compute dynamic tendency
!
!======= DIV & VOR ===================
!
   ksdel=nwg*levsp
   SANMU=>WORK1A(1:mt,1:jlg+1,1:levsp)
   SANMV=>WORK1A(1:mt,1:jlg+1,1+levsp:2*levsp)
   SANM1=>WORK0A(1:mt,1:jlg  ,1+ksdel:ksdel+levsp)
   SANM2=>WORK0A(1:mt,1:jlg  ,ksdel+levsp+1:ksdel+2*levsp)
   SANM3=>WORK0A(1:mt,1:jlg  ,1+(nw0-1)*levsp:nw0*levsp)
#ifdef BIN_DBG
   if (kdt.eq.2) then
     call file_write_wave(335,SANMU,mt*jlg,2*levsp)
     call file_write_wave(335,SANM1,mt*jlg,2*levsp)
     call file_write_wave(335,SANM3,mt*jlg,  levsp)
   endif
#endif
   do k = 1,levsp
     do mj = 1,JLg       ! VJ: VORTicity, DJ: Divergence
       do m = 1,MT
         mx=midx(m)
#ifdef HDSZ
         vivor= -Z_VIS(k)*V1(m,mj,k)
         vidiv= -Z_VIS(k)*D1(m,mj,k)
#else
         vivor= 0.
         vidiv= 0.
#endif
         avor= -(SANMU(mop(m),mj,k)*(-mx)*rrerth_+SANM2(m,mj,k))+vivor
         adiv= -(SANMV(mop(m),mj,k)*( mx)*rrerth_+SANM1(m,mj,k))+              &
                 SANM3(m,mj,k)+vidiv
         VJ(m,mj,k)= avor
         DJ(m,mj,k)= adiv
       enddo
     enddo
   enddo! k=1,levs

!======= Temperature =================
   SANMU=>WORK1A(1:mt,1:jlg+1,ndel*levsp+1:(ndel+1)*levsp)
   SANM1=>WORK0A(1:mt,1:jlg  ,1+levsp:2*levsp)
   SANM2=>WORK0A(1:mt,1:jlg  ,1+(nwg+2)*levsp:(nwg+3)*levsp)
   do k = 1,levsp
     do mj = 1,JLg       ! TJ: Temprature
       do m = 1,MT      
         mx=midx(m)
         atem= -(SANMU(mop(m),mj,k)*(-mx)*rrerth_+SANM2(m,mj,k))+SANM1(m,mj,k)
         TJ(m,mj,k)= atem
       enddo
     enddo
   enddo ! k=1,levsp
#ifndef NISLQ
!
!======= Moisture =================
!
   SANMU=>WORK1A(1:mt,1:jlg+1,1+(ndel+1)*levsp:nw1*levsp)
   SANM1=>WORK0A(1:mt,1:jlg  ,1+2*levsp:nwg*levsp)
   SANM2=>WORK0A(1:mt,1:jlg  ,1+(nwg+3)*levsp:(nwg+ndel)*levsp)
   do k = 1,levhp
     do mj = 1,JLg       ! QJ: Moisture
       do m = 1,MT      
         mx=midx(m)
         QJ(m,mj,k)= -(SANMU(mop(m),mj,k)*(-mx)*rrerth_+SANM2(m,mj,k))       &
                      +SANM1(m,mj,k)
       enddo
     enddo
   enddo
#endif /* ~NISLQ end */
#ifdef CLM_CWF
   wcon(1:mt,1:jlg,1:levsp)=0.0
   do nv = 1,ncldg+1
     do k = 1,levsp
       kk=(nv-1)*levsp+k
       do mj = 1,JLg       ! QJ: Moisture
         do m = 1,MT      
           mx=midx(m)
           wcon(m,mj,k)=wcon(m,mj,k)-                                          &
                      (SANMU(mop(m),mj,kk)*(-mx)*rrerth_+SANM2(m,mj,kk))
         enddo
       enddo
     enddo
   enddo
   call dfs_fft_driver(-1,cgs,ib,jbw,levs,wcon,mt,jlg,levsp,jlg,               &
                        levsp,levs,1,coslat,1)
#endif
!
!=========================  ADVECTION TERM OF SURFACE PRESSURE ===!
!
#ifdef HYBRID
   call dfs_fft_driver(1,dqdt,ib,jbw,1,PSJ,mt,jlg,1,jlg,1,1,0,coslat,1)
#else
   call dfs_fft_driver(1,APSN,ib,jbw,1,PSJ,mt,jlg,1,jlg,1,1,0,coslat,1)
#endif
!
   RETURN
   end subroutine dfs_dynamics_advection
!-------------------------------------------------------------------------------
!
#ifdef HYBRID
!-------------------------------------------------------------------------------
   subroutine vcnhyb(im,km,nm,dt,zint,zmid,zdot,zadv,nvcn,xvcn)
!-------------------------------------------------------------------------------
!
! subprogram:    vcnhyb      vertical advection instability filter
!
! abstract: filters vertical advection tendencies
!   in the dynamics tendency equation in order to ensure stability
!   when the vertical velocity exceeds the cfl criterion.
!   the vertical velocity in this case is sigmadot.
!   for simple second-order centered eulerian advection,
!   filtering is needed when vcn=zdot*dt/dz>1.
!   the maximum eigenvalue of the linear advection equation
!   with second-order implicit filtering on the tendencies
!   is less than one for all resolvable wavenumbers (i.e. stable)
!   if the nondimensional filter parameter is nu=(vcn**2-1)/4.
!
! program history log:
!   1997-07-30  iredell
!   2000-01-01  hann-ming henry juang  mpi
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call vcnhyb(im,km,nm,dt,zint,zmid,zdot,zadv,nvcn,xvcn)
!
!   input argument list:
!     im       - integer number of gridpoints to filter
!     km       - integer number of vertical levels
!     nm       - integer number of fields
!     dt       - real timestep in seconds
!     zint     - real (im,km+1) interface vertical coordinate values
!     zmid     - real (im,km) midlayer vertical coordinate values
!     zdot     - real (im,km+1) vertical coordinate velocity
!     zadv     - real (im,km,nm) vertical advection tendencies
!
!   output argument list:
!     zadv     - real (im,km,nm) vertical advection tendencies
!     nvcn     - integer number of points requiring filtering
!     xvcn     - real maximum vertical courant number
!
!   subprograms called:
!     tridim_hyb   - tridiagonal matrix solver
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer,intent(in):: im,km,nm
   real,intent(in):: dt,zint(im,km+1),zmid(im,km),zdot(im,km+1)
   real,intent(inout):: zadv(im,km,nm)
   integer,intent(out):: nvcn
   real,intent(out):: xvcn
   integer i,j,k,n,ivcn(im)
   logical lvcn(im)
   real zdm,zda,zdb,vcn(im,km-1)
   real rnu,cm(im,km),cu(im,km-1),cl(im,km-1)
   real rr(im,km,nm)
!-------------------------------------------------------------------------------
!
!  compute vertical courant number
!  increase by 10% for safety
!
   nvcn=0
   xvcn=0.
   lvcn=.false.
   do k = 1,km-1
     do i = 1,im
       zdm=zmid(i,k)-zmid(i,k+1)
       vcn(i,k)=abs(zdot(i,k+1)*dt/zdm)*1.1
       lvcn(i)=lvcn(i).or.vcn(i,k).gt.1
       xvcn=max(xvcn,vcn(i,k))
     enddo
   enddo
!
!  determine points requiring filtering
!
   if(xvcn.gt.1) then
!
     do i = 1,im
       if(lvcn(i)) then
         ivcn(nvcn+1)=i
         nvcn=nvcn+1
       endif
     enddo
!
!  compute tridiagonal matrim
!
     do j = 1,nvcn
       cm(j,1)=1
     enddo
     do k = 1,km-1
       do j = 1,nvcn
         i=ivcn(j)
         if(vcn(i,k).gt.1) then
          zdm=zmid(i,k)-zmid(i,k+1)
          zda=zint(i,k+1)-zint(i,k+2)
          zdb=zint(i,k)-zint(i,k+1)
           rnu=(vcn(i,k)**2-1)/4
           cu(j,k)=-rnu*zdm/zdb
           cl(j,k)=-rnu*zdm/zda
           cm(j,k)=cm(j,k)-cu(j,k)
           cm(j,k+1)=1-cl(j,k)
         else
           cu(j,k)=0
           cl(j,k)=0
           cm(j,k+1)=1
         endif
       enddo
     enddo
!
!  fill fields to be filtered
!
     do n = 1,nm
       do k = 1,km
         do j = 1,nvcn
           i=ivcn(j)
           rr(j,k,n)=zadv(i,k,n)
         enddo
       enddo
     enddo
!
!  solve tridiagonal system
!
        call tridim_hyb(nvcn,im,km,km,nm,cl,cm,cu,rr,cu,rr)
!
!  replace filtered fields
!
     do n = 1,nm
       do k = 1,km
         do j = 1,nvcn
           i=ivcn(j)
           zadv(i,k,n)=rr(j,k,n)
         enddo
       enddo
     enddo
!
   endif
!
   end subroutine vcnhyb
!
!-------------------------------------------------------------------------------
   subroutine tridim_hyb(l,lx,n,nx,m,cl,cm,cu,r,au,a)
!-------------------------------------------------------------------------------
!
! subprogram:    tridim_hyb      solves tridiagonal matrix problems.
!
! abstract: this routine solves multiple tridiagonal matrix problems
!   with multiple right-hand-side and solution vectors for every matrix.
!   the solutions are found by eliminating off-diagonal coefficients,
!   marching first foreward then backward along the matrix diagonal.
!   the computations are vectorized around the number of matrices.
!   no checks are made for zeroes on the diagonal or singularity.
!
! program history log:
!   07-30-1997    mark iredell    development
!   04-30-2007    jung-eun kim    grims implementation
!
! usage:    call tridim_hyb(l,lx,n,nx,m,cl,cm,cu,r,au,a)
!
!   input argument list:
!     l        - integer number of tridiagonal matrices
!     lx       - integer first dimension (lx>=l)
!     n        - integer order of the matrices
!     nx       - integer second dimension (nx>=n)
!     m        - integer number of vectors for every matrix
!     cl       - real (lx,2:n) lower diagonal matrix elements
!     cm       - real (lx,n) main diagonal matrix elements
!     cu       - real (lx,n-1) upper diagonal matrix elements
!                (may be equivalent to au if no longer needed)
!     r        - real (lx,nx,m) right-hand-side vector elements
!                (may be equivalent to a if no longer needed)
!
!   output argument list:
!     au       - real (lx,n-1) work array
!     a        - real (lx,nx,m) solution vector elements
!
!-------------------------------------------------------------------------------
   real  ::  cl(lx,2:n),cm(lx,n),cu(lx,n-1),r(lx,nx,m),                        &
                                 au(lx,n-1),a(lx,nx,m)
!-------------------------------------------------------------------------------
!  march up
!
   do i = 1,l
     fk=1./cm(i,1)
     au(i,1)=fk*cu(i,1)
   enddo
!
   do j = 1,m
     do i = 1,l
       fk=1./cm(i,1)
       a(i,1,j)=fk*r(i,1,j)
     enddo
   enddo
!
   do k = 2,n-1
     do i = 1,l
       fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
       au(i,k)=fk*cu(i,k)
     enddo
     do j = 1,m
       do i = 1,l
         fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
         a(i,k,j)=fk*(r(i,k,j)-cl(i,k)*a(i,k-1,j))
       enddo
     enddo
   enddo
!-------------------------------------------------------------------------------
!  march down
!
   do j = 1,m
     do i = 1,l
       fk=1./(cm(i,n)-cl(i,n)*au(i,n-1))
       a(i,n,j)=fk*(r(i,n,j)-cl(i,n)*a(i,n-1,j))
     enddo
   enddo
!
   do k = n-1,1,-1
     do j = 1,m
       do i = 1,l
         a(i,k,j)=a(i,k,j)-au(i,k)*a(i,k+1,j)
       enddo
     enddo
   enddo
!
   end subroutine tridim_hyb
!-------------------------------------------------------------------------------
#endif /* ~HYBRID end */
!
!-------------------------------------------------------------------------------
   subroutine file_write_wave(nn,var,mt,jl,lot)
!-------------------------------------------------------------------------------
   use dfsvar, only : mtg,jlg,iope
!-------------------------------------------------------------------------------
   integer,intent(in)      ::  nn,mt,jl,lot
   real,target,intent(in)  ::  var(mt*jl,lot)
!
   real*4                  ::  work4(mtg*jlg)
#ifdef MP
   real                    ::  full(mtg*jlg,lot)
!-------------------------------------------------------------------------------
   call mpsp2f(var,mt,jl,full,mtg,jlg,lot)
#define VAR full
#endif
   if (iope) then
     do k = 1,lot
       work4(:)=VAR(:,k)
       write(nn)work4(:)
     enddo
   endif
#undef VAR
!
   return
   end subroutine file_write_wave
!-------------------------------------------------------------------------------
