#include "define.h"
!-------------------------------------------------------------------------------
!  ::: structure ::: This file contains ...
!
!        |--- [dfs_divergence] *
!        |--- [dfs_divergence_del_matrix] *
!        |--- [dfs_divergence_hydro] *
!        |--- [dfs_divergence_laplacian1] *
!        |--- [dfs_divergence_laplacian2] *
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!-------------------------------------------------------------------------------
#ifndef HYBRID
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_divergence                                                   &
   ( DT,SF1,D1,T1,PSL1,DJ,TJ,PSJ, EIGVAL,EIGVEC,AEIGVEC,D3)
!-------------------------------------------------------------------------------
   use constant, only : rd_
   use dfsvar, only   : MT,JL,levsp,levs,jla,jlha,JB,jls,jlg,kls,              &
                        nls,nle,ngs,nge,iope,                                  &
                        grs,trs,phs,                                           &
                        PHSDM,TRSDM,TAVEXY,jcol2js,                            &
#ifdef MP
                        SF1l,                                                  &
#endif
                        AMATm,DMATm,aMSQUAR,DMATmi_,amatm_,dmatm_,             &
                        tai,ib,jbw,coslat
!-------------------------------------------------------------------------------
!                                                                    
!     A Do-loop is deployed for efficient vector-processing          
!                                                                    
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real                              ::  dt
   real, dimension(mt,jlg)           ::  SF1,PSJ,PSL1
   real, dimension(mt,jlg,levsp)     ::  D1,DJ,D3,T1,TJ
   real                              ::  EIGVAL(levs),EIGVEC(levs,levs),       &
                                         AEIGVEC(levs,levs)
!-----local-variable-------------------------------------------------!
   real, dimension(mt,jlg,levsp)     ::  WNMJG1,WNMJG2
   real, dimension(mt,jl,levs)       ::  WNMJL1,WNMJL2
#ifdef MP
   real, dimension(mt,jlg,levsp)     ::  WNMJGA
   real, dimension(mt,jl,levs)       ::  WT1,WTJ
   real, dimension(mt,jl,levs)       ::  WD3,WTRSDM
   real, dimension(mt,jl)            ::  psl1l, phsdml, psjl
#endif
   real                              ::  DT2,DT2S
   integer                           ::  k,i4,i3,i2,mx,mj,kk,iu
   real                              ::  DTold,sdtterm,alpha,dtterm
   real                              ::  XMATl(levsp),XMATg(levs)
   DATA DTold/-1234.d0/   ! arbitrary negative value
!--------------------------------------------------------------------!
!
#ifdef BACKWARD
   DT2  = DT*2.0
#else
   DT2  = DT
#endif
   DT2S = DT2*DT2
!
   IF(DT2S.ne.DTold) THEN
     DTold= DT2S
     do k = 1,levsp
       kk=kls+k-1
       alpha= -EIGVAL(kk)*DT2S
       CALL dfs_divergence_del_matrix (alpha,mt,jlha,AMATm,DMATm,aMSQUAR,      &
                         DMATmi_(1,1,0,1,k) )
     enddo
   ENDIF
!
   WNMJG2(1:mt,1:jlg,1:levsp)=0.0
!
!---------------------------------------
!
#ifdef MP
   call mpmz2mn(T1,jlg,levsp,WT1,jl,mt,levs,1)
   call mpmz2mn(TJ,jlg,levsp,WTJ,jl,mt,levs,1)
   call mpm2mn(psl1,jlg,psl1l,jl,mt,1)
   call mpm2mn(psj ,jlg,psjl ,jl,mt,1)
#define SF1 SF1l
#define PSL1 psl1l
#define PSJ psjl
#define T1 WT1
#define TJ WTJ
#endif
   CALL dfs_divergence_hydro( SF1,T1, WNMJL2, mt,jl,levs, +1 )
!
!  Level_split
!
   do k = 1,levs
     do mj = 1,jl
       do mx = 1,mt
         sdtterm=TAVEXY(k)*PSJ(mx,mj)*rd_
         dtterm= WNMJL2(mx,mj,k) + TAVEXY(k)*PSL1(mx,mj)*rd_
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
         do kk = 1,levs
           sdtterm= sdtterm + GRS(k,kk)*TJ(mx,mj,kk)
         enddo
         WNMJL1(mx,mj,k)= -dtterm*DT2 - sdtterm*DT2S
       enddo ! mj
     enddo ! mx
   enddo ! levs
!
#ifdef MP
#undef SF1
#undef PSL1
#undef PSJ
#undef T1
#undef TJ
!
! to transpose divergent field
!
   call mpmn2mz(WNMJL1,jl,levs,WNMJGA,jlg,mt,levsp,1)
#define WNMJL1 WNMJGA
#endif
   CALL dfs_divergence_laplacian1(  WNMJL1, WNMJG2,                            &
                     amatm,AMATm_,dmatm,DMATm_,                                &
                     amsquar,jcol2js,mt,jlha,levsp,0 )
!
   do k = 1,levsp
     do mj = 1,jlg
       do mx = 1,mt
         WNMJG1(mx,mj,k)= D1(mx,mj,k) + DJ(mx,mj,k)*DT2+ WNMJG2(mx,mj,k)
       enddo ! mj
     enddo ! mx
   enddo ! k=1,levsp
!
   call dfs_global_mean( WNMJG1,XMATl,mt,ngs,nge,levsp,+1 )
!
!-------------------
!
#ifdef MP
#undef WNMJL1
   call mpmz2mn(WNMJG1,jlg,levsp,WNMJL1,jl,mt,levs,1)
#define WNMJG1 WNMJL1
#endif
   do k = 1,levs
     do mj = 1,JL
       do mx = 1,MT
         WNMJL2(mx,mj,k)=0.
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
         do kk = 1,levs
           WNMJL2(mx,mj,k)=WNMJL2(mx,mj,k)+                                    &
                       AEIGVEC(k,kk)*WNMJG1(mx,mj,kk)
         enddo
       enddo ! mj
     enddo ! mx
   enddo ! levs
!
#ifdef MP
#undef WNMJG1
   call mpmn2mz(WNMJL2,jl,levs,WNMJGA,jlg,mt,levsp,1)
#define WNMJL2 WNMJGA
#endif
   CALL dfs_divergence_laplacian2 (WNMJL2,DMATmi_,mt,jla,jlha,levsp)
#ifdef MP
#undef WNMJL2
   call mpmz2mn(WNMJGA,jlg,levsp,WNMJL2,jl,mt,levs,1)
   call mpm2mn(phsdm,jlg,phsdml,jl,mt,1)
#define D3 WD3
#define TRSDM WTRSDM 
#define PHSDM PHSDMl
#endif
!
!-------------------
!  Level_split
!
   do k = 1,levs
     do mj = 1,JL
       do mx = 1,MT
         D3(mx,mj,k)=0.
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
         do kk = 1,levs
           D3(mx,mj,k)=D3(mx,mj,k)+EIGVEC(k,kk)*WNMJL2(mx,mj,kk)
         enddo
       enddo ! mj
     enddo ! mx
   enddo ! k=1,levs
!
   CALL dfs_global_mean( D3,XMATg,mt,nls,nle,levs,+1 )
!
!-------------------
!  Level_split
!
   do k = 1,levs
     do mj = 1,JL
       do mx = 1,MT
         TRSDM(mx,mj,k)=0.
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
         do kk = 1,levs
           TRSDM(mx,mj,k)=TRSDM(mx,mj,k)+TRS(k,kk)*D3(mx,mj,kk)
         enddo
       enddo ! mj
     enddo ! mx
   enddo ! k=1,levs
!
!-------------------
!
   do mj = 1,JL
     do mx = 1,MT
       PHSDM(mx,mj)=0.
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
       do kk = 1,levs
         PHSDM(mx,mj)=PHSDM(mx,mj)+PHS(kk)*D3(mx,mj,kk)
       enddo
     enddo
   enddo
!
#ifdef MP
#undef D3
#undef TRSDM
#undef PHSDM
   call mpmn2mz(WD3   ,jl,levs,D3   ,jlg,mt,levsp,1)
   call mpmn2mz(WTRSDM,jl,levs,TRSDM,jlg,mt,levsp,1)
   call mpmn2m (phsdml,jl      ,phsdm,jlg,mt     ,1)
#endif
!
   RETURN
   END SUBROUTINE dfs_divergence
!-------------------------------------------------------------------------------
#else /* HYBRID */
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_divergence_hybrid                                            &
   ( DT,SF1,D1,T1,PSL1,DJ,TJ,PSJ, EIGVAL,EIGVEC,AEIGVEC,D2,T2,PSL2,D3)
!-------------------------------------------------------------------------------
   use constant, only : rd_
   use dfsvar, only   : MT,JL,levsp,levs,jla,jlha,JB,jls,jlg,kls,              &
                        nls,nle,ngs,nge,iope,                                  &
                        grs,trs,phs,gv,                                        &
                        PHSDM,TRSDM,TAVEXY,jcol2js,                            &
#ifdef MP
                        SF1l,                                                  &
#endif
                        AMATm,DMATm,aMSQUAR,DMATmi_,amatm_,dmatm_,             &
                        tai,ib,jbw,coslat
!-------------------------------------------------------------------------------
!                                                                    
!     A Do-loop is deployed for efficient vector-processing          
!                                                                    
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real                                        ::  dt
   real, dimension(mt,jlg)                     ::  SF1,PSJ,PSL1,PSL2
   real, dimension(mt,jlg,levsp)               ::  D1,DJ,D2,D3,T1,T2,TJ
   real, dimension(levs)                       ::  EIGVAL
   real, dimension(levs,levs)                  ::  EIGVEC,AEIGVEC
!-----local-variable-------------------------------------------------!
   real, dimension(mt,jlg,levsp)               ::  WNMJG1,WNMJG2
   real, dimension(mt,jl,levs)                 ::  WNMJL1,WNMJL2
#ifdef MP
   real, dimension(mt,jlg,levsp)               ::  WNMJGA
!
   real, dimension(mt,jlg,3*levsp)             ::  WORKG
   real, dimension(mt,jlg,2)                   ::  PWORKG
   real, dimension(mt,jl,3*levs)  , target     ::  WORKL
   real, dimension(mt,jl,2)       , target     ::  PWORKL
   real, dimension(:,:,:)         , pointer    ::  WT0,WD2,WTJ
   real, dimension(:,:)           , pointer    ::  ps0l,psjl
!
   real, dimension(mt,jl,levs)                 ::  WD3,WTRSDM
   real, dimension(mt,jl)                      ::  phsdml
#endif
   integer                                     ::  i,j,k,i4,i3,i2,mx,mj,kk,iu
   real                                        ::  DT2,DT2S
   real                                        ::  DTold,sdtterm,alpha,dtterm
   real, dimension(levsp)                      ::  XMATl
   real, dimension(levs)                       ::  XMATg
   real, dimension(mt,jlg)                     ::  PS0
   real, dimension(mt,jlg,levsp)               ::  T0
   DATA DTold/-1234.d0/   ! arbitrary negative value
!--------------------------------------------------------------------!
!
#ifdef BACKWARD
   DT2  = DT*2.0
#else
   DT2  = DT
#endif
   DT2S = DT2*DT2
!
   if ( DT2S.ne.DTold ) then
     DTold= DT2S
     do k = 1,levsp
       kk=kls+k-1
       alpha= -EIGVAL(kk)*DT2S
       CALL dfs_divergence_del_matrix (alpha,mt,jlha,AMATm,DMATm,aMSQUAR,      &
                         DMATmi_(1,1,0,1,k) )
     enddo
   endif
!
! initialize
!
   WNMJG2(1:mt,1:jlg,1:levsp)=0.0
!
!---------------------------------------
!
   forall(i=1:mt,j=1:jlg,k=1:levsp)  t0(i,j,k)=t1(i,j,k)-t2(i,j,k)
   forall(i=1:mt,j=1:jlg)           ps0(i,j)  =psl1(i,j)-psl2(i,j)
#ifdef MP
!
! MPI transpose
!
   do k = 1,levsp
     do j = 1,jlg
       do i = 1,mt
         workg(i,j,        k)=t0(i,j,k)
         workg(i,j,  levsp+k)=tj(i,j,k)
         workg(i,j,2*levsp+k)=d2(i,j,k)
       enddo
     enddo
   enddo
!
   do j = 1,jlg
     do i = 1,mt
       pworkg(i,j,1)=ps0(i,j)
       pworkg(i,j,2)=psj(i,j)
     enddo
   enddo
!
! transpose
!
   call mpmz2mn(WORKG ,jlg,levsp,WORKL ,jl,mt,levs,3)
   call mpm2mn (PWORKG,jlg,      PWORKL,jl,mt,     2)
!
! pointing
!
   wt0=>WORKL(1:mt,1:jl,       1:levs)
   wtj=>WORKL(1:mt,1:jl,  levs+1:2*levs)
   wd2=>WORKL(1:mt,1:jl,2*levs+1:3*levs)
   ps0l=>PWORKL(1:mt,1:jl,1)
   psjl=>PWORKL(1:mt,1:jl,2)
#define SF1 SF1l
#define PS0 ps0l
#define PSJ psjl
#define T0 wt0
#define TJ wtj
#define D2 wd2
#endif /* MP end */
!
!  subtract off linear dependence on divergence
!
   do k = 1,levs
     do mj = 1,jl
       do mx = 1,mt
         do kk = 1,levs
           TJ(mx,mj,k)=TJ(mx,mj,k)+TRS(k,kk)*D2(mx,mj,kk)
         enddo
       enddo
     enddo
   enddo
!
   do k = 1,levs
     do mj = 1,jl
       do mx = 1,mt
         PSJ(mx,mj)=PSJ(mx,mj)+PHS(k)*D2(mx,mj,k)
       enddo
     enddo
   enddo
!
   CALL dfs_divergence_hydro( SF1,T0, WNMJL2, mt,jl,levs, +1 )
!
!  Level_split
!
   do k = 1,levs
     do mj = 1,jl
       do mx = 1,mt
         sdtterm=GV(k)*PSJ(mx,mj)
         dtterm= WNMJL2(mx,mj,k) + GV(k)*PS0(mx,mj)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#else
!$dir unroll_and_jam(unroll_factor=100)
#endif
         do kk = 1,levs
           sdtterm= sdtterm + GRS(k,kk)*TJ(mx,mj,kk)
         enddo
         WNMJL1(mx,mj,k)= -dtterm*DT2 - sdtterm*DT2S
       enddo ! mj
     enddo ! mx
   enddo ! levs
!
#ifdef MP
#undef SF1
#undef PS0
#undef PSJ
#undef T0
#undef TJ
#undef D2
   call mpmn2mz(WNMJL1,jl,levs,WNMJGA,jlg,mt,levsp,1) ! transpose div. field
#define WNMJL1 WNMJGA
#endif
   CALL dfs_divergence_laplacian1(  WNMJL1, WNMJG2,                            &
                     amatm,AMATm_,dmatm,DMATm_,                                &
                     amsquar,jcol2js,mt,jlha,levsp,0 )
!
   forall(mx=1:mt,mj=1:jlg,k=1:levsp)
     WNMJG1(mx,mj,k)= D1(mx,mj,k) + DJ(mx,mj,k)*DT2+ WNMJG2(mx,mj,k)
   end forall
!
   call dfs_global_mean( WNMJG1,XMATl,mt,ngs,nge,levsp,+1 )
!
!-------------------
!
#ifdef MP
#undef WNMJL1
   call mpmz2mn(WNMJG1,jlg,levsp,WNMJL1,jl,mt,levs,1)
#define WNMJG1 WNMJL1
#endif /* MP end */
!
   do k = 1,levs
     do mj = 1,JL
       do mx = 1,MT
         WNMJL2(mx,mj,k)=0.
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#else
!$dir unroll_and_jam(unroll_factor=100)
#endif
         do kk = 1,levs
           WNMJL2(mx,mj,k)=WNMJL2(mx,mj,k)+                                    &
                         AEIGVEC(k,kk)*WNMJG1(mx,mj,kk)
         enddo
       enddo ! mj
     enddo ! mx
   enddo ! levs
!
#ifdef MP
#undef WNMJG1
   call mpmn2mz(WNMJL2,jl,levs,WNMJGA,jlg,mt,levsp,1)
#define WNMJL2 WNMJGA
#endif /* MP end */
   CALL dfs_divergence_laplacian2 (WNMJL2,DMATmi_,mt,jla,jlha,levsp)
#ifdef MP
#undef WNMJL2
   call mpmz2mn(WNMJGA,jlg,levsp,WNMJL2,jl,mt,levs,1)
#define D3 WD3
#define TRSDM WTRSDM 
#define PHSDM phsdml
#endif /* MP end */
!
!-------------------
!  Level_split
!
   do k = 1,levs
     do mj = 1,JL
       do mx = 1,MT
         D3(mx,mj,k)=0.
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#else
!$dir unroll_and_jam(unroll_factor=100)
#endif
         do kk = 1,levs
           D3(mx,mj,k)=D3(mx,mj,k)+EIGVEC(k,kk)*WNMJL2(mx,mj,kk)
         enddo
       enddo ! mj
     enddo ! mx
   enddo ! k=1,levs
!
   CALL dfs_global_mean( D3,XMATg,mt,nls,nle,levs,+1 )
!
!-------------------
!  Level_split
!
   do k = 1,levs
     do mj = 1,JL
       do mx = 1,MT
         TRSDM(mx,mj,k)=0.
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#else
!$dir unroll_and_jam(unroll_factor=100)
#endif
         do kk = 1,levs
           TRSDM(mx,mj,k)=TRSDM(mx,mj,k)+TRS(k,kk)*D3(mx,mj,kk)
         enddo
       enddo ! mj
     enddo ! mx
   enddo ! k=1,levs
!
!-------------------
!
   do mj = 1,JL
     do mx = 1,MT
       PHSDM(mx,mj)=0.
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#else
!$dir unroll_and_jam(unroll_factor=100)
#endif
       do kk = 1,levs
         PHSDM(mx,mj)=PHSDM(mx,mj)+PHS(kk)*D3(mx,mj,kk)
       enddo
     enddo
   enddo
#ifdef MP
#undef D3
#undef TRSDM
#undef PHSDM
!
! parallel update tj,d3,trsdm,psj,phsdm
!
   forall(i=1:mt,j=1:jlg,k=1:levsp)  workg(i,j,k)=0.0
   forall(i=1:mt,j=1:jlg,k=1:2)     pworkg(i,j,k)=0.0
!
   do k = 1,levs
     do j = 1,jl
       do i = 1,mt
         workl(i,j,       k)=wtj(i,j,k)
         workl(i,j,  levs+k)=wd3(i,j,k)
         workl(i,j,2*levs+k)=wtrsdm(i,j,k)
       enddo
     enddo
   enddo
!
   do j = 1,jl
     do i = 1,mt
       pworkl(i,j,1)=psjl(i,j)
       pworkl(i,j,2)=phsdml(i,j)
     enddo
   enddo
!
! transpose
!
   call mpmn2mz(workl ,jl,levs, workg,jlg,mt,levsp,3)
   call mpmn2m (pworkl,jl,     pworkg,jlg,mt,      2)
!
! update
!
   do k = 1,levsp
     do j = 1,jlg
       do i = 1,mt
         tj   (i,j,k)=workg(i,j,        k)
         d3   (i,j,k)=workg(i,j,  levsp+k)
         trsdm(i,j,k)=workg(i,j,2*levsp+k)
       enddo
     enddo
   enddo
!
   do j = 1,jlg
     do i = 1,mt
       psj  (i,j)=pworkg(i,j,1)
       phsdm(i,j)=pworkg(i,j,2)
     enddo
   enddo
#endif /* MP end */
!
   RETURN
   END SUBROUTINE dfs_divergence_hybrid
!-------------------------------------------------------------------------------
#endif /* ~HYBRID end */
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_divergence_del_matrix                                        &
                           (alpha,mt,jlha,AMATm,DMATm,aMSQUAR,DMATmi_)
!-------------------------------------------------------------------------------
   use dfsvar, only : midx
!-------------------------------------------------------------------------------
!                                                                    
!     To solve for psi; Real    VARIABLES                            
!     psi + alpha* ( LAPLACian of psi ) = X_i,j ,                    
!     which arises from the high-order spectral filter               
!                                                                    
!    [CALL  MAT200S] -> [CALL  LAP200S] : Real-variable version      
!                                                                    
!     01 MAR 2005                                                    
!                                                                    
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real              ::  alpha
   integer           ::  jlha,mt
   real              ::  AMATm(MT,3,0:jlha,2)
   real              ::  DMATm(MT,3,0:jlha,2)
   real              ::  aMSQUAR(MT)
   real, intent(out) ::  DMATmi_(MT,3,0:jlha,2)         ! outout
!-----local variable-------------------------------------------------!
   real              ::  YMATD(MT,3,0:jlha),YMATDi      ! add inverse
   integer           ::  js,k,i,mx,mc,mc1
   real              ::  rat,rei,amx2
!-------------------------------------------------------------------------------
!
   if(alpha.eq.0.d0) then
     write(*,*)'The alpha is zero: SUB MAT200S_VEC'
     call exit(33)
   endif
!
!************************************************!
!BEG  MATRIECS for LAPLACIAN OPERATOR            * SEMI-IMPLICIT
!************************************************!
!
   DO js = 1,2
     DO K = 0,jlha
       DO i = 1,3
         DO mx = 1,mt
           DMATmi_(mx,i,K,js)= AMATm(mx,i,K,js) + alpha*DMATm(mx,i,K,js)
         END DO
       END DO
     END DO
   END DO
!
!************************************************!
!END  MATRIECS for LAPLACIAN OPERATOR            * SEMI-IMPLICIT
!************************************************!
!
!************************************************!
!BEG  Time saving, prehandling of the MATRIX     * SEMI-IMPLICIT
!************************************************!
!
!beg-------------------------------------------------- DMAT -> ZDMAT_
!-----m0,m1,m2,(n1) and m0,m1,m2,(n2)
!
   DO js = 1,2
!
     do mc = 1,jlha-1
       Do mx = 1,mt
         amx2= -aMSQUAR(mx)
         YMATD( mx,1,mc)= DMATmi_(mx,1,mc,js)
         YMATD( mx,2,mc)= DMATmi_(mx,2,mc,js) - amx2*alpha 
         YMATD( mx,3,mc)= DMATmi_(mx,3,mc,js)
       END Do
     end do
!
     mc=0
     Do mx = 1,mt
       amx2= -aMSQUAR(mx)
       YMATD( mx,1,mc)= DMATmi_(mx,1,mc,js) - amx2*alpha 
       YMATD( mx,2,mc)= DMATmi_(mx,2,mc,js)
     END Do
!
     mc=jlha
     Do mx = 1,mt
       amx2= -aMSQUAR(mx)
       YMATD( mx,1,mc)= DMATmi_(mx,1,mc,js)
       YMATD( mx,2,mc)= DMATmi_(mx,2,mc,js) - amx2*alpha 
     END Do
!
!-------------------
!
     do mc = 0,jlha-1  
       mc1= mc+1   
       Do mx = 1,mt
         YMATDi=1./YMATD( mx,1,mc)
         IF(mc.eq.0) then
           if(midx(mx).eq.0 .AND. js.eq.2) then
             rat= YMATD( mx,1,mc1) *YMATDi
           else
             rat= 0.d0
           endif
         ELSE
           rat= YMATD( mx,1,mc1) *YMATDi
         ENDIF

!
!next-sub         YMAT(mc )= YMAT(mc )*rat
!
         YMATD( mx,2,mc1)= YMATD( mx,2,mc1)-YMATD( mx,2,mc )*rat
!
!next-sub         YMAT(mc1)= YMAT(mc1)-YMAT(mc )
!
         YMATD( mx,1,mc1)= YMATD( mx,2,mc1)
         YMATD( mx,2,mc1)= YMATD( mx,3,mc1)
         !
         IF(mc.eq.0) then
           if(midx(mx).eq.0 .AND. js.eq.2) then
             DMATmi_( mx,1,mc ,js)= YMATDi !! (1/all)
           else
             DMATmi_( mx,1,mc ,js)= 0.d0
           endif
         ELSE
           DMATmi_( mx,1,mc ,js)= YMATDi !! (1/all)
         ENDIF

         DMATmi_( mx,2,mc ,js)=      YMATD( mx,2,mc )       
         DMATmi_( mx,3,mc ,js)= rat
       END Do ! par
     end do ! mc
!
     Do mx = 1,mt
       DMATmi_( mx,1,jlha,js)= 1.d0/YMATD( mx,1,jlha) !! (1/all)
     END Do ! par
!
   END DO ! js ! js=1,2
!
!end-------------------------------------------------- DMAT -> DMAT_
!
!************************************************!
!END  Time saving, prehandling of the MATRIX     * SEMI-IMPLICIT
!************************************************!
!
   RETURN
   END SUBROUTINE dfs_divergence_del_matrix
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_divergence_hydro( SF1,T1, F1,mt,jl,levs,isign ) 
                                                ! isign=1 : temp->geop
!-------------------------------------------------------------------------------
   use constant, only : g_            ! gravity
   use dfsvar, only   : GRS,GRSINV
!-------------------------------------------------------------------------------
   implicit none
!
   integer  ::  isign,mt,jl,levs
   real     ::  SF1(MT,JL)
   real     ::  T1(MT,JL,levs),F1(MT,JL,levs)
!---  local variables -----------------------------------------------!
   integer  ::  k,mj,mx,kk
!-------------------------------------------------------------------------------
!
   IF(isign.eq.+1) THEN 
!
! Level_split
!
     do k = 1,levs
       do mj = 1,JL
         do mx = 1,MT
           F1(mx,mj,k)= 0.
           F1(mx,mj,k)= SF1(mx,mj)*g_
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
           do kk = 1,levs
             F1(mx,mj,k)= F1(mx,mj,k)+GRS(k,kk) * T1(mx,mj,kk)
           enddo
         end do ! mj
       end do ! mx
     end do ! k=1,levs
!
   ELSEIF(isign.eq.-1) THEN 
!
! Level_split
!
     do k = 1,levs
       do mj = 1,JL
         do mx = 1,MT
           T1(mx,mj,k)= -SF1(mx,mj)
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
           do kk = 1,levs
             T1(mx,mj,k)= T1(mx,mj,k)+GRSINV(k,kk)*F1(mx,mj,kk)
           enddo
         end do ! mj
       end do ! mx
     end do ! k=1,levs
!
   ENDIF
!
   RETURN
   END SUBROUTINE dfs_divergence_hydro
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_divergence_laplacian1 ( TC1, AJX ,amatm,AMATm_,dmatm,DMATm_, &
                       amsquar,jcol2js,mt,jlha,levsp,MSEL )
!-------------------------------------------------------------------------------
   use dfsvar, only : mls,nls,nle,ngs,nge
!-------------------------------------------------------------------------------
!                                                                     
!     Laplacian of X = Y , Inverse Laplacian of Y = X                 
!     on the Spherical Domain including the Poles                     
!                                                                     
!    [CALL MAT200E] -> [CALL LAP200E]                                 
!                                                                     
!     01 MAR 2005                                                     
!                                                                     
!    MSEL=0 : Laplacian TC1 -> AJX   
!         1 : Inverse                
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  mt,jlha,levsp,MSEL
   real     ::  TC1(mt,0:jlha*2,levsp), AJX(mt,0:jlha*2,levsp)
   real     ::  AMATm (mt,3,0:jlha,2),DMATm (mt,3,0:jlha,2)
   real     ::  DMATm_(mt,3,0:jlha,2),AMATm_(mt,3,0:jlha,2)
   real     ::  aMSQUAR(mt)
   integer  ::  jcol2js(0:jlha,2)
!
   real     ::  YMAT(mt,0:jlha)
   integer  ::  js,jja,jjb,jjc,jcol,mx,jjj,mc,mco,mcop1,mc1
   integer  ::  jla,k,mx0
   real     ::  a4,a5,a6,ace,aaa,XMATl(levsp)
!-------------------------------------------------------------------------------
   jla=jlha*2
   mx0=0
   if (mls.eq.0) then
     mx0=mt/2+1
   endif
!
!--------------------------------------------------------------------!
   IF(MSEL.eq.0) THEN           ! Laplacian of TC1 = AJX
!--------------------------------------------------------------------!
!
     DO k = 1,levsp
       AJX(1:mt,0,k)=0.
!
       DO js = 1,2
!
!---- 1  To get YMAT, YMAT= DMAT*TC1
!
         do jcol = 1,jlha-1
           jja= jcol2js(jcol-1,js) ! For jcol=1 js=1, jja=-1. This problem
           jjb= jcol2js(jcol  ,js) ! was cured by letting jcol2js(0,1)=+1,
           jjc= jcol2js(jcol+1,js) ! noticing [*MATm(1,1,js)=0] except for mx=0 js=2
           Do mx = 1,mt
             a4= DMATm(mx,1,jcol,js)
             a5= DMATm(mx,2,jcol,js)
             a6= DMATm(mx,3,jcol,js) 
             ace=(a5+amsquar(mx))*TC1(mx,jjb,k)
             YMAT(mx,jcol)= a4*TC1(mx,jja,k) + ace + a6*TC1(mx,jjc,k)
           END Do ! par
         end do
!
         jcol=0
         jjb= jcol2js(jcol  ,js)
         jjc= jcol2js(jcol+1,js)
         Do mx = 1,mt
           a4= DMATm(mx,1,jcol,js)
           a5= DMATm(mx,2,jcol,js)
           YMAT(mx,jcol)=(a4+amsquar(mx))*TC1(mx,jjb,k)+                       &
                    a5*TC1(mx,jjc,k)
         END Do ! par
!
         jcol=jlha
         jja= jcol2js(jcol-1,js)
         jjb= jcol2js(jcol  ,js)
         Do mx = 1,mt
           a4= DMATm(mx,1,jcol,js)
           a5= DMATm(mx,2,jcol,js)
           YMAT(mx,jcol)= a4*TC1(mx,jja,k)+                                    &
                   (a5+amsquar(mx))*TC1(mx,jjb,k)
         END Do ! par
!
!---- 2  To get ZETA, AMAT*ZETA= YMAT
!
         do mc = 0,jlha-1
           mc1= mc+1
           Do mx = 1,mt
             YMAT(mx,mc1)= YMAT(mx,mc1)-YMAT(mx,mc )*AMATm_(mx,3,mc,js)
           END Do ! par
         end do
!
         Do mx = 1,mt
           YMAT(mx,jlha)=  YMAT(mx,jlha) * AMATm_(mx,1,jlha,js) !! (1/all)
         END Do ! par
!
         do mco = jlha-1,1,-1
           mcop1=mco+1
           Do mx = 1,mt
             aaa =  YMAT(mx,mcop1)*AMATm_(mx,2,mco,js)
             YMAT(mx,mco)= (YMAT(mx,mco)-aaa) * AMATm_(mx,1,mco,js)
           END Do ! par
         end do
!
         do jcol = 1,jlha
           jjj= jcol2js(jcol,js)
           Do mx = 1,mt
             AJX(mx,jjj,k)= YMAT(mx,jcol)
           END Do ! par
         end do
!
       END DO ! js=1,2
!
     enddo      ! K
!
     CALL dfs_global_mean( AJX, XMATl, mt,ngs,nge,levsp,+1 ) ! MAK=+1: glo-mean->0
!
!--------------------------------------------------------------------!
   ELSEIF(MSEL.eq.1) THEN           ! Inverse Laplacian of AJX= TC1
!--------------------------------------------------------------------!
!
     CALL dfs_global_mean( AJX, XMATl,mt,ngs,nge, levsp,+1 ) ! MAK=+1: glo-mean->0
!
     DO k = 1,levsp                              ! To end
!
       tc1(1:mt,0,k)=0.
!
       DO js = 1,2
!
!---- 1  To get YMAT, YMAT= AMAT*ZETA
!
         do jcol = 1,jlha-1
           jja= jcol2js(jcol-1,js)
           jjb= jcol2js(jcol  ,js)
           jjc= jcol2js(jcol+1,js)
           Do mx = 1,mt
             a4= AMATm(mx,1,jcol,js)
             a5= AMATm(mx,2,jcol,js)
             a6= AMATm(mx,3,jcol,js)
             YMAT(mx,jcol)= a4*AJX(mx,jja,k)+a5*AJX(mx,jjb,k)+                 &
               a6*AJX(mx,jjc,k)
           END Do ! par
         end do
!
         jcol=0
         jjb= jcol2js(jcol  ,js)
         jjc= jcol2js(jcol+1,js)
         Do mx = 1,mt
           a4= AMATm(mx,1,jcol,js)
           a5= AMATm(mx,2,jcol,js)
           YMAT(mx,jcol)= a4*AJX(mx,jjb,k)+a5*AJX(mx,jjc,k)    
         END Do ! par
!
         jcol=jlha
         jja= jcol2js(jcol-1,js)
         jjb= jcol2js(jcol  ,js)
         Do mx = 1,mt
           a4= AMATm(mx,1,jcol,js)
           a5= AMATm(mx,2,jcol,js)
           YMAT(mx,jcol)= a4*AJX(mx,jja,k)+a5*AJX(mx,jjb,k)
         END Do ! par
!
!---- 2  To get PSI, DMAT*PSI= YMAT
!
         do mc = 0,jlha-1
           mc1= mc+1
           Do mx = 1,mt
             YMAT(mx,mc1)= YMAT(mx,mc1)                                        &
                     -YMAT(mx,mc )*DMATm_(mx,3,mc,js)
           END Do ! par
         end do
!
         Do mx = 1,mt
           YMAT(mx,jlha)=  YMAT(mx,jlha) * DMATm_(mx,1,jlha,js) !! (1/all)
         END Do ! par
!
         do mco = jlha-1,1,-1
           mcop1=mco+1
           Do mx = 1,mt
             aaa =  YMAT(mx,mcop1)*DMATm_(mx,2,mco,js)
             YMAT(mx,mco)= (YMAT(mx,mco)-aaa) * DMATm_(mx,1,mco,js)
           END Do ! par
         end do
!
         do jcol = 1,jlha
           jjj= jcol2js(jcol,js)
           Do mx = 1,mt
             TC1(mx,jjj,k)= YMAT(mx,jcol)
           END Do ! par
         end do
!
       END DO ! js=1,2
!
       if (mx0.ne.0) TC1( mx0, 0 ,k)= 0.             
!
     enddo      ! k=1,levsp
!
!--------------------------------------------------------------------!
!
   ENDIF
!
   RETURN
   END SUBROUTINE dfs_divergence_laplacian1
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_divergence_laplacian2 ( AJX, DMATmi_ ,mt,jla,jlha,levsp)
!-------------------------------------------------------------------------------
   use dfsvar, only : mls,jcol2js,AMATm
!-------------------------------------------------------------------------------
!                                                                    
!     To get psi or Div from : Real    VARIABLES                     
!                                                                    
!     psi + alpha * ( Laplacian of psi ) = Y  : Global Ave[Y] !=0    
!     Div + alpha * ( Laplacian of Div ) = Z  : Global Ave[Z]  =0    
!                                               because [Div]  =0    
!    [CALL  MAT200S] -> [CALL  LAP200S] : Real-variable version      
!                                                                    
!     01 MAR 2005                                                    
!
!     The meridional wave should be global
!                                                                    
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  mt,jlha,jla,levsp
   real     ::  AJX(MT,0:jla,levsp),DMATmi_(MT,3,0:jlha,2,levsp)
!
   real     ::  ymat(mt,0:jlha)
   real     ::  suma,dd,ajx0,a4,a5,a6,tc10,aaa
   integer  ::  mj,js,jcol,mx,jja,jjb,jjc,mc,mcop1,jjj
   integer  ::  mc1,mco,mzr,k
!-------------------------------------------------------------------------------
!
   do k = 1,levsp
     if (mls.eq.0) then
       mzr=mt/2+1
       suma= 0.                      ! global-mean
       do mj = 0,jla,2
         dd= DFLOAT(mj)**2-1.D-00
         suma= suma - AJX(mzr,mj,k)/dd
       end do
       ajx0= suma
     endif
!
!-------------- mx>0, mx<0 ------------c
!
     DO js = 1,2
!
!---- 1  To get YMAT, YMAT= AMAT*ZETA
!
       do jcol = 1,jlha-1
         jja= jcol2js(jcol-1,js)
         jjb= jcol2js(jcol  ,js)
         jjc= jcol2js(jcol+1,js)
         Do mx = 1,mt
           a4= AMATm(mx,1,jcol,js)
           a5= AMATm(mx,2,jcol,js)
           a6= AMATm(mx,3,jcol,js)
           YMAT(mx,jcol)= a4*AJX(mx,jja,k)+a5*AJX(mx,jjb,k)                    &
                         +a6*AJX(mx,jjc,k)
         END Do ! par
       end do
!
       jcol=0
       jjb= jcol2js(jcol  ,js)
       jjc= jcol2js(jcol+1,js)
       Do mx = 1,mt
         a4= AMATm(mx,1,jcol,js)
         a5= AMATm(mx,2,jcol,js)
         YMAT(mx,jcol)= a4*AJX(mx,jjb,k)+a5*AJX(mx,jjc,k)                   
       END Do ! par
!
       jcol=jlha
       jja= jcol2js(jcol-1,js)
       jjb= jcol2js(jcol  ,js)
       Do mx = 1,mt
         a4= AMATm(mx,1,jcol,js)
         a5= AMATm(mx,2,jcol,js)
         YMAT(mx,jcol)= a4*AJX(mx,jja,k)+a5*AJX(mx,jjb,k)
       END Do ! par
!
!---- 2  To get PSI, DMAT*PSI= YMAT
!
       do mc = 0,jlha-1 
         mc1= mc+1
         Do mx = 1,mt
           YMAT(mx,mc1)=  YMAT(mx,mc1)                                         &
                         -YMAT(mx,mc )*DMATmi_(mx,3,mc,js,k)
         END Do ! par
       end do
!
       Do mx = 1,mt
         YMAT(mx,jlha)=  YMAT(mx,jlha) * DMATmi_(mx,1,jlha,js,k) !! (1/all)
       END Do ! par
!
       do mco = jlha-1,1,-1
         mcop1=mco+1
         Do mx = 1,mt
           aaa =  YMAT(mx,mcop1)*DMATmi_(mx,2,mco,js,k)
           YMAT(mx,mco)= (YMAT(mx,mco)-aaa)*DMATmi_(mx,1,mco,js,k) !! (1/all)
         END Do ! par
       end do
!
       do jcol = 1,jlha
         jjj= jcol2js(jcol,js)
         Do mx = 1,mt
           AJX(mx,jjj,k)=  YMAT(mx,jcol)
         END Do ! par
       end do
!
     END DO ! js=1,2
!
     if (mls.eq.0) then
       AJX( mzr, 0,k )= 0.
       suma= 0.                      ! global-mean
       do mj = 0,jla,2
         dd= DFLOAT(mj)**2-1.D-00
         suma= suma - AJX(mzr,mj,k)/dd
       end do
       tc10= suma
       AJX( mzr, 0 ,k)= ajx0 - tc10   ! global-mean correction
     endif
   enddo         ! levsp
!
   RETURN
   END SUBROUTINE dfs_divergence_laplacian2
!-------------------------------------------------------------------------------
