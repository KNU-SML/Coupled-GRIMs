#include <define.h>
   subroutine dfs_semi_implicit
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [dfs_semi_implicit]
!           |
!           |-- [dfs_semi_impl_integrate] *
!           |-- [dfs_semi_matrix] *
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
   end subroutine dfs_semi_implicit
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_semi_impl_integrate(d3,t3,psl2,dj,tj,psj,dt)
!-------------------------------------------------------------------------------
   use constant, only : rd_
   use dfsvar, only   : MT,JL,jlg,jla,JLHA,levs,levsp,nls,nle,ngs,nge,         &
                        grs,trs,phs,EIGVAL,EIGVEC,AEIGVEC,                     &
                        TAVEXY,jcol2js,                                        &
                        AMATm,DMATm,aMSQUAR,DMATmi_,amatm_,dmatm_
!-------------------------------------------------------------------------------
!                                                                    
!     A Do-loop is deployed for efficient vector-processing          
!                                                                    
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real                             ::  dt
   real, dimension(mt,jlg,levsp)     ::  d3,t3,dj,tj
   real, dimension(mt,jlg)          ::  psl2,psj
!
! local-variable
!
   real, dimension(mt,jlg,levsp)     ::  WNM3D1, WNM3D2
#ifdef MP
   real, dimension(mt,jl,levs)     ::  WNM3DL1, WNM3DL2          ,&
                                        WTJ,WT3
   real, dimension(mt,jl)           ::  psl2l,psjl
#ifdef BACKWARD
   real, dimension(mt,jl,levs)     ::  WD3
#endif
#endif
#ifndef BACKWARD
   real, dimension(mt,jl,levs)     ::  dum
#ifdef MP
   real, dimension(mt,jlg,levsp)   ::  dum2
#endif
#endif
   real                             ::  DTS2,DT2
   integer                          ::  k,i4,i3,i2,mx,mj,kk
   real                             ::  DTold,sdtterm,alpha,dtterm
   real                             ::  XMATg(levs),XMATl(levsp)
   DATA DTold/-1234.d0/   ! arbitrary negative value
!
#ifdef BACKWARD
   DT2 =   2.*DT
#ifdef MP
#define DUM WD3
#else
#define DUM D3
#endif
#else /* ~BACKWARD */
   DT2 =   DT
#endif
   DTS2 =  DT2*DT2
!
#ifdef MP
   call mpmz2mn(TJ,jlg,levsp,WTJ,jl,mt,levs,1)
#define TJ WTJ
#define WNM3D1 WNM3DL1
   call mpm2mn(psj ,jlg,psjl ,jl,mt,1)
   call mpm2mn(psl2,jlg,psl2l,jl,mt,1)
#define PSJ psjl
#define PSL2 psl2l
#endif
   do k = 1,levs
     !add hoon
     WNM3D1(:,:,k)=0.0
     do mj = 1,jl
       do mx = 1,MT
         sdtterm=TAVEXY(k)*PSJ(mx,mj)*rd_
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
         do kk = 1,levs
           sdtterm= sdtterm + GRS(k,kk)*TJ(mx,mj,kk)
         enddo
         WNM3D1(mx,mj,k)= -sdtterm*DTS2
       enddo ! mj
     enddo ! mx
   enddo
!
#ifdef MP
#undef WNM3D1
   call mpmn2mz(WNM3DL1,jl,levs,WNM3D1,jlg,mt,levsp,1)
#endif
   CALL dfs_divergence_laplacian1 ( WNM3D1, WNM3D2 ,                           &
           amatm,AMATm_,dmatm,DMATm_,amsquar,jcol2js,mt,jlha,levsp,0 )
   do k = 1,levsp
     do mj = 1,JLG
       do mx = 1,MT
         WNM3D1(mx,mj,k)= DT2*DJ(mx,mj,k)+WNM3D2(mx,mj,k)
       enddo ! mj
     enddo ! mx
   enddo ! k=1,levsp
   CALL dfs_global_mean( WNM3D1, XMATl, mt,ngs,nge ,levsp,+1 )
!
#ifdef MP
   call mpmz2mn(WNM3D1,jlg,levsp,WNM3DL1,jl,mt,levs,1)
#define WNM3D1 WNM3DL1
#define WNM3D2 WNM3DL2
#endif
   do k = 1,levs
     do mj = 1,JL
       do mx = 1,MT
         WNM3D2(mx,mj,k)=0.
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
         do kk = 1,levs
           WNM3D2(mx,mj,k)=WNM3D2(mx,mj,k)+                                    &
                        AEIGVEC(k,kk)*WNM3D1(mx,mj,kk)
         enddo
       enddo    ! mj
     enddo    ! mx
   enddo    ! k=1,levs
#ifdef MP
#undef WNM3D1
#undef WNM3D2
   call mpmn2mz(WNM3DL2,jl,levs,WNM3D2,jlg,mt,levsp,1)
#endif
   CALL dfs_divergence_laplacian2( WNM3D2, DMATmi_,mt,jla,jlha,levsp)
#ifdef MP
   call mpmz2mn(WNM3D2,jlg,levsp,WNM3DL2,jl,mt,levs,1)
#define WNM3D2 WNM3DL2
#define T3 WT3
#endif
!
   do k = 1,levs
     do mj = 1,JL
       do mx = 1,MT
         DUM(mx,mj,k)=0.
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
         do kk = 1,levs
           DUM(mx,mj,k)=DUM(mx,mj,k)+EIGVEC(k,kk)*WNM3D2(mx,mj,kk)
         enddo
       enddo ! mj
     enddo ! mx
   enddo ! k=1,levs
   CALL dfs_global_mean( DUM, XMATg, mt,nls,nle,levs,+1 )
!
   do k = 1,levs
     do mj = 1,JL
       do mx = 1,MT
         WNM3D2(mx,mj,k)=0.
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
         do kk = 1,levs
           WNM3D2(mx,mj,k)=WNM3D2(mx,mj,k)+TRS(k,kk)*DUM(mx,mj,kk)
         enddo
         T3(mx,mj,k)= T3(mx,mj,k) +DT2*(tj(mx,mj,k)-WNM3D2(mx,mj,k))
       enddo ! mj
     enddo ! mx
   enddo ! k=1,levs
!                                                          ! Level_split
   do mj = 1,JL
     do mx = 1,MT
       WNM3D2(mx,mj,1)=0.
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
       do kk = 1,levs
         WNM3D2(mx,mj,1)=WNM3D2(mx,mj,1)+PHS(kk)*DUM(mx,mj,kk)
       enddo
       PSL2(mx,mj)= PSL2(mx,mj) + DT2*(PSJ(mx,mj)-WNM3D2(mx,mj,1))
     enddo
   enddo
#ifdef MP
#undef T3
#undef PSL2
#undef PSJ
   call mpmn2m(psl2l,jl,psl2,jlg,mt,1)
   call mpmz2mn(WT3,jl,levs,T3,jlg,mt,levsp,1)
#endif
!
#ifdef BACKWARD
#undef DUM
#ifdef MP
   call mpmn2mz(WD3,jl,levs,D3,jlg,mt,levsp,1)
#endif
#else /* ~BACKWARD */
#ifdef MP
   call mpmn2mz(dum,jl,levs,dum2,jlg,mt,levsp,1)
#define DUM dum2
#else
#define DUM dum
#endif
!
   do k = 1,levsp
     do mj = 1,JLg
       do mx = 1,MT
         D3(mx,mj,k)=D3(mx,mj,k)+DUM(mx,mj,k)
       enddo ! mj
     enddo ! mx
   enddo ! k=1,levs
   CALL dfs_global_mean( D3, XMATl, mt,ngs,nge,levsp,+1 )
#undef DUM
#endif /* BACKWARD end */
!
   RETURN
   END SUBROUTINE dfs_semi_impl_integrate
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_semi_matrix(AMATm,DMATm,aMSQUAR,DMATm_,AMATm_,jcol2js,       &
                          mt,jlha)
!-------------------------------------------------------------------------------
   use constant, only : rrerth_
   use dfsvar, only   : midx,nev,nod,mevs,mods,mls
!-------------------------------------------------------------------------------
!                                                                    
!     SUBROUTINE : MATRICES and INVERSION                            
!                                                                    
!     Matrices for Laplacian of PSI                                  
!                                                                    
!    [CALL MAT200E] -> [CALL LAP200E]                                
!                                                                    
!     AMATm, DMATm, BMATm: Size=(0:JLHP,0:JLHP)                      
!                          Used=(0:JLH ,0:JLH ) for all inversions   
!                          except for 5-diagonal-matrix construction 
!                          in SUB MATS_BIH                           
!                                                                    
!     01 MAR 2005                                                    
!                                                                    
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  mt,jlha
   real     ::  AMATm (mt,3,0:jlha,2),DMATm (mt,3,0:jlha,2)
   real     ::  DMATm_(mt,3,0:jlha,2),AMATm_(mt,3,0:jlha,2)
   real     ::  aMSQUAR(mt)
   integer  ::  jcol2js(0:jlha,2)
!
! local variable
!
   integer  ::  jl
   real     ::  YMATD(mt,3,0:jlha)
   integer  ::  mt1
   integer  ::  js,k,i,mx,mc,mc1,jsmo,jcol,jcoljsm,j,m
   real     ::  rat,amx2,reis
!
   jl=jlha*2
!
   do js = 1,2         ! meridional index used in SUB LAP---
     jsmo=MOD(js,2)
     do jcol = 0,jlha                     ! mj=0,2,4 : js=even
       jcoljsm = jcol*2 - jsmo         ! mj=1,3,5 : js=odd
       if(jcoljsm.LT.0) jcoljsm= 1  !! 0th row for m1 m2 m0:n2 are unnece.
       jcol2js(jcol,js)= jcoljsm     !! and to avoid -1 array index
     enddo
   enddo
!
   reis=rrerth_*rrerth_
   do mx = 1,mt
     aMSQUAR(mx)= -midx(mx)**2 *reis
   enddo
!
   AMATm (1:mt,1:3,0:jlha,1:2)=0.
   DMATm (1:mt,1:3,0:jlha,1:2)=0.
   AMATm_(1:mt,1:3,0:jlha,1:2)=0.
   DMATm_(1:mt,1:3,0:jlha,1:2)=0.
!
!************************************************!
!BEG  MATRIECS for LAPLACIAN OPERATOR            !
!************************************************!
!beg-- Matrices for Y= ( D^2/Dlambda^2 + COS(LAT) D/Dphi COS(LAT) D/Dphi ) X
!                  js
!        DMATm(1,K, 1 ) :  m=1,3,..., n=1,3,...
!        DMATm(1,K, 2 ) :  m=1,3,..., n=2,4,...
!        DMATm(1,K, 3 ) :  m=2,4,..., n=1,3,...      ! Same for AMATm
!        DMATm(1,K, 4 ) :  m=2,4,..., n=2,4,...
!        DMATm(1,K, 5 ) :  m=0      , n=1,3,...
!        DMATm(1,K, 6 ) :  m=0      , n=0,2,4,.
!
!      Be careful of the identity that DMAT0(0,1,2)= DMAT0(1,1,2)= 0
!
   Do m = 1,nod
     mx=mods(m)
!
!----- m1,n2
!
     do j = 4,jl-2,2                   ! m=1,3,...
       K= (j+0)/2                              ! n=2,4,...
       DMATm(mx,1,K,2)=  (j-2.d0)*(j-1.d0)*0.25
       DMATm(mx,2,K,2)= -(dfloat(j)**2)*0.5 
       DMATm(mx,3,K,2)=  (j+1.d0)*(j+2.d0)*0.25
     enddo
     j=2
     K= (j+0)/2                              ! n=2,4,...
     DMATm(mx,2,K,2)= -(dfloat(j)**2)*0.5    ! 1st row of 3-diag MATRIX
     DMATm(mx,3,K,2)=  (j+1.d0)*(j+2.d0)*0.25
     DMATm(mx,1,K,2)=   0.
     j=JL
     K= (j+0)/2                              ! n=2,4,...
     DMATm(mx,1,K,2)=  (j-2.d0)*(j-1.d0)*0.25   ! last row of 3-diag MATRIX
     DMATm(mx,2,K,2)= -(dfloat(j)**2)*0.50
     DMATm(mx,3,K,2)=   0.
!
!----- m1,n1
!
     do j = 3,JL-3,2                   ! m=1,3,...
       K= (j+1)/2                              ! n=1,3,...
       DMATm(mx,1,K,1)=  (j-2.d0)*(j-1.d0)*0.25
       DMATm(mx,2,K,1)= -(dfloat(j)**2)*0.50
       DMATm(mx,3,K,1)=  (j+1.d0)*(j+2.d0)*0.25
     enddo
     j=1
     K= (j+1)/2                              ! n=1,3,...
     DMATm(mx,2,K,1)= -(dfloat(j)**2)*0.50      ! 1st row of 3-diag MATRIX
     DMATm(mx,3,K,1)=  (j+1.d0)*(j+2.d0)*0.25
     DMATm(mx,1,K,1)=   0.
     j=JL-1
     K= (j+1)/2                              ! n=1,3,...
     DMATm(mx,1,K,1)=  (j-2.d0)*(j-1.d0)*0.25   ! last row of 3-diag MATRIX
     DMATm(mx,2,K,1)= -(dfloat(j)**2)*0.50
     DMATm(mx,3,K,1)=   0.
   ENDDo ! mx=1,3,...
!
!-----m2,n2
!
   Do m = 1,nev
     mx=mevs(m)
     if(midx(mx).NE.0) then
       do j = 4,JL-2,2                   ! m=2,4,...
         K= (j+0)/2                              ! n=2,4,...
         DMATm(mx,1,K,2)=   dfloat(j)*(j-1.d0)*0.25
         DMATm(mx,2,K,2)= -(dfloat(j)**2)*0.50
         DMATm(mx,3,K,2)=   dfloat(j)*(j+1.d0)*0.25
       enddo
       j=2
       K= (j+0)/2                              ! n=2,4,...
       DMATm(mx,2,K,2)= -(dfloat(j)**2)*0.50      ! 1st row of 3-diag MATRIX
       DMATm(mx,3,K,2)=   dfloat(j)*(j+1.d0)*0.25
       DMATm(mx,1,K,2)=   0.
       j=JL
       K= (j+0)/2                              ! n=2,4,...
       DMATm(mx,1,K,2)=   dfloat(j)*(j-1.d0)*0.25 ! last row of 3-diag MATRIX
       DMATm(mx,2,K,2)= -(dfloat(j)**2)*0.50
       DMATm(mx,3,K,2)=   0.
!
!-----m2,n1
!
       do j = 3,JL-3,2                   ! m=2,4,...
         K= (j+1)/2                              ! n=1,3,...
         DMATm(mx,1,K,1)=   dfloat(j)*(j-1.d0)*0.25
         DMATm(mx,2,K,1)= -(dfloat(j)**2)*0.50
         DMATm(mx,3,K,1)=   dfloat(j)*(j+1.d0)*0.25
       enddo
       j=1
       K= (j+1)/2                              ! n=1,3,...
       DMATm(mx,2,K,1)= -(dfloat(j)**2)*0.50    ! 1st row of 3-diag MATRIX
       DMATm(mx,3,K,1)=   dfloat(j)*(j+1.d0)*0.25
       DMATm(mx,1,K,1)=   0.
       j=JL-1
       K= (j+1)/2                              ! n=1,3,...
       DMATm(mx,1,K,1)=   dfloat(j)*(j-1.d0)*0.25 ! last row of 3-diag MATRIX
       DMATm(mx,2,K,1)= -(dfloat(j)**2)*0.50
       DMATm(mx,3,K,1)=   0.
     else                  ! m=0
!
!-----m0,n1
!
       do j = 3,JL-3,2           ! n=1,3,...
         K= (j+1)/2 
         DMATm(mx,1,K,1)=  (j-2.d0)*(j-1.d0)*0.25
         DMATm(mx,2,K,1)= -(dfloat(j)**2)*0.50
         DMATm(mx,3,K,1)=  (j+1.d0)*(j+2.d0)*0.25
       enddo
       j=1
       K= (j+1)/2                              ! n=1,3,...
       DMATm(mx,2,K,1)= -(dfloat(j)**2)*0.5000    ! 1st row of 3-diag MATRIX
       DMATm(mx,3,K,1)=  (j+1.d0)*(j+2.d0)*0.2500
       DMATm(mx,1,K,1)=   0.
       j=JL-1
       K= (j+1)/2                              ! n=1,3,...
       DMATm(mx,1,K,1)=  (j-2.d0)*(j-1.d0)*0.25   ! last row of 3-diag MATRIX
       DMATm(mx,2,K,1)= -(dfloat(j)**2)*0.50
       DMATm(mx,3,K,1)=   0.
!
!-----m0,n2
!
       do j = 2,JL-2,2         ! m=0  & n=0,2,...
         K= (j+0)/2 
         DMATm(mx,1,K,2)=  (j-2.d0)*(j-1.d0)*0.25
         DMATm(mx,2,K,2)= -(dfloat(j)**2)*0.50
         DMATm(mx,3,K,2)=  (j+1.d0)*(j+2.d0)*0.25
       enddo
!
!----j=0
!
       j=0
       K= (j+0)/2                              ! n=0,2,...
       DMATm(mx,1,K,2)= -(dfloat(j)**2)*0.50    ! 1st row of 3-diag MATRIX
       DMATm(mx,2,K,2)=  (j+1.d0)*(j+2.d0)*0.25
       DMATm(mx,3,K,2)=   0.
       j=JL
       K= (j+0)/2                              ! n=0,2,...
       DMATm(mx,1,K,2)=  (j-2.d0)*(j-1.d0)*0.25   ! last row of 3-diag MATRIX
       DMATm(mx,2,K,2)= -(dfloat(j)**2)*0.50
       DMATm(mx,3,K,2)=   0.
     endif
   ENDDo ! m=2,4,...
!
!end-- Matrices for Y= ( D^2/Dlambda^2 + COS(LAT) D/Dphi COS(LAT) D/Dphi ) X
!
!
!beg-- Matrices for Z= VORTicity * COS(LAT)^2
!
!
!-----m1
!
   Do m = 1,nod
     mx=mods(m)
     do K = 2,jlha-1                            ! n=2,4,...
       AMATm(mx,1,K,2)=  -0.2500
       AMATm(mx,2,K,2)=   0.5000
       AMATm(mx,3,K,2)=  -0.2500
     enddo
     K=1                                  ! n=2,4,...
     AMATm(mx,2,K,2)=   0.5000
     AMATm(mx,3,K,2)=  -0.2500
     AMATm(mx,1,K,2)=   0.
     K=jlha                                ! n=2,4,...
     AMATm(mx,1,K,2)=  -0.2500
     AMATm(mx,2,K,2)=   0.5000
     AMATm(mx,3,K,2)=   0.
     do K = 1,jlha                              ! n=1,3,...
       do I = 1,3
         AMATm(mx,I,K,1)= AMATm(mx,I,K,2)
       enddo
     enddo
   ENDDo ! mx=1,3,...

   Do m = 1,nod
     mx=mods(m)
!    AMATm(mx,1,1,1)=   0.5000 + 0.2500  ! note
     AMATm(mx,2,1,1)=   0.5000 + 0.2500  ! note
   ENDDo
!
!-----m2
!
! Any odd m
!
   mt1=mods(1)
   Do m = 1,nev
     mx=mevs(m)
     if(midx(mx).NE.0) then
       do K = 1,jlha
         do I = 1,3                ! m1  n1n2
           AMATm(mx,I,K,1)= AMATm(mt1,I,K,1)     ! m2 : n=1,3,5...
           AMATm(mx,I,K,2)= AMATm(mt1,I,K,2)     ! m2 : n=2,4,6...
         enddo
       enddo
!
!-----m0
!
     else
       do K = 1,jlha
         do I = 1,3                ! m1    n2
           AMATm(mx,I,K,1)= AMATm(mt1,I,K,2)     ! m=0, n=1,3,...
           AMATm(mx,I,K,2)= AMATm(mt1,I,K,2)     ! m=0, n=0,2,4,.
         enddo
       enddo
!      AMATm(mx,1,1,1)=   0.5000 - 0.2500  ! note
       AMATm(mx,2,1,1)=   0.5000 - 0.2500  ! note
       AMATm(mx,1,0,2)=   0.5000                    
       AMATm(mx,2,0,2)=  -0.2500                    
       AMATm(mx,3,0,2)=   0.                          
       AMATm(mx,1,1,2)=  -0.5000                    
       AMATm(mx,2,1,2)=   0.5000                    
       AMATm(mx,3,1,2)=  -0.2500                   
     endif
   ENDDo ! mx=2,4,...
!
!end-- Matrices for Z= VORTicity * COS(LAT)^2
!
!************************************************!
!END  MATRIECS for LAPLACIAN OPERATOR            !
!************************************************!
!
! (1/(radius**2)] Dimensionalize by Earth's parameter
!
   do js = 1,2
     do k = 0,jlha
       do i = 1,3
         do mx = 1,mt
           DMATm(mx,i,k,js)= DMATm(mx,i,k,js)*reis
         enddo
       enddo
     enddo
   enddo
!
!************************************************!
!BEG  Time saving, prehandling of the MATRIX     !
!************************************************!
!
!---- AJX*COS(LAT)^2                  = AJX*AMAT = SNM   -A-
!---- COSLAT d/dphi COSLAT d/dphi TC1 = DMAT*TC1 = SNM   -B-
!---- To get AJX, AJX= Laplacian of TC1, 1) SNM from -B- 2) AJX= AMAT_*SNM
!---- To get TC1, AJX= Laplacian of TC1, 1) SNM from -A- 2) TC1= DMAT_*SNM
!
!
!beg-------------------------------------------------- DMAT -> DMAT_
!
!                12......                12......     
!                123.....                .12.....
!        DMATx   .123....       DMATx_   ..12....
!                ........                ...12...
!                .....123                ........
!                ......12                ......12
!
!-----m0,m1,m2,(n1) and m0,m1,m2,(n2)
!
   DO js = 1,2
     do mc = 1,jlha-1
       Do mx = 1,mt
         amx2= -aMSQUAR(mx)
         YMATD( mx,1,mc)= DMATm(mx,1,mc,js)
         YMATD( mx,2,mc)= DMATm(mx,2,mc,js) - amx2
         YMATD( mx,3,mc)= DMATm(mx,3,mc,js)
       ENDDo
     enddo
     mc=0
     Do mx = 1,mt
       amx2= -aMSQUAR(mx)
       YMATD( mx,1,mc)= DMATm(mx,1,mc,js) - amx2
       YMATD( mx,2,mc)= DMATm(mx,2,mc,js)
     ENDDo
     mc=jlha
     Do mx = 1,mt
       amx2= -aMSQUAR(mx)
       YMATD( mx,1,mc)= DMATm(mx,1,mc,js)
       YMATD( mx,2,mc)= DMATm(mx,2,mc,js) - amx2
     ENDDo
!
     do mc = 0,jlha-1                    
       mc1= mc+1                                 
       Do mx = 1,mt
         IF(mc.eq.0) then
           rat= 0.d0
         ELSE
           rat= YMATD( mx,1,mc1) / YMATD( mx,1,mc)
         ENDIF
!
!next-sub      YMAT(mc )= YMAT(mc )*rat   
!
         YMATD( mx,2,mc1)= YMATD( mx,2,mc1)-YMATD( mx,2,mc )*rat    
!
!next-sub      YMAT(mc1)= YMAT(mc1)-YMAT(mc )
!
         YMATD( mx,1,mc1)= YMATD( mx,2,mc1)
         YMATD( mx,2,mc1)= YMATD( mx,3,mc1)
!
         IF(mc.eq.0) then
           DMATm_( mx,1,mc ,js)= 0.d0
         ELSE
           DMATm_( mx,1,mc ,js)= 1.d0/YMATD( mx,1,mc ) !! (1/all)
         ENDIF
         DMATm_( mx,2,mc ,js)=      YMATD( mx,2,mc )                       
         DMATm_( mx,3,mc ,js)= rat
       ENDDo ! par
     enddo ! mc
!
     Do mx = 1,mt
       DMATm_( mx,1,jlha,js)= 1.d0/YMATD( mx,1,jlha) !! (1/all)
     ENDDo ! par
!
   ENDDO ! js=1,2
!
!end-------------------------------------------------- DMAT -> DMAT_
!
!beg-------------------------------------------------- AMAT -> AMAT_
!-----m0,m1,m2,(n1) and m0,m1,m2,(n2)
!
   DO js = 1,2
!
     do mc = 1,jlha-1
       Do mx = 1,mt
         YMATD( mx,1,mc)= AMATm(mx,1,mc,js)
         YMATD( mx,2,mc)= AMATm(mx,2,mc,js)
         YMATD( mx,3,mc)= AMATm(mx,3,mc,js)
       ENDDo
     enddo
!
     mc=0
     Do mx = 1,mt
       YMATD( mx,1,mc)= AMATm(mx,1,mc,js)
       YMATD( mx,2,mc)= AMATm(mx,2,mc,js)
     END Do
     mc=jlha
     Do mx = 1,mt
       YMATD( mx,1,mc)= AMATm(mx,1,mc,js)
       YMATD( mx,2,mc)= AMATm(mx,2,mc,js)
     END Do
!
     do mc = 0,jlha-1                    
       mc1= mc+1                               
       Do mx = 1,mt
         IF(mc.eq.0) then
           if(midx(mx).eq.0 .AND. js.eq.2) then
             rat= YMATD( mx,1,mc1) / YMATD( mx,1,mc)
           else
             rat= 0.d0
           endif 
         ELSE
           rat= YMATD( mx,1,mc1) / YMATD( mx,1,mc)
         ENDIF
         YMATD( mx,2,mc1)= YMATD( mx,2,mc1)-YMATD( mx,2,mc)*rat
         YMATD( mx,1,mc1)= YMATD( mx,2,mc1)
         YMATD( mx,2,mc1)= YMATD( mx,3,mc1)
         IF(mc.eq.0) then
           if(midx(mx).eq.0 .AND. js.eq.2) then
             AMATm_( mx,1,mc ,js)= 1.d0/YMATD( mx,1,mc ) !! (1/all)
           else
             AMATm_( mx,1,mc ,js)= 0.d0
           endif
         ELSE
           AMATm_( mx,1,mc ,js)= 1.d0/YMATD( mx,1,mc ) !! (1/all)
         ENDIF
!
         AMATm_( mx,2,mc ,js)=      YMATD( mx,2,mc )
         AMATm_( mx,3,mc ,js)= rat
       ENDDo
     enddo ! mc
!
     Do mx = 1,mt
       AMATm_( mx,1,jlha,js)= 1.d0/YMATD( mx,1,jlha) !! (1/all)
     ENDDo

   ENDDO ! js
!
!end-------------------------------------------------- AMAT -> AMAT_
!
!************************************************!
!END  Time saving, prehandling of the MATRIX     !
!************************************************!
!
   RETURN
   END SUBROUTINE dfs_semi_matrix
!-------------------------------------------------------------------------------
