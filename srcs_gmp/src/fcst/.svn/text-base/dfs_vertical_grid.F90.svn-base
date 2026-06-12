#include <define.h>
   subroutine dfs_vertical_grid
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [dfs_vertical_grid]
!           |
!           |-- [dfs_vertical_descr] *
!           |-- [dfs_eigen_vector] *
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
   end subroutine dfs_vertical_grid
!-------------------------------------------------------------------------------

#ifndef HYBRID
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_vertical_descr( VERMAT, EIGVAL, EIGVEC, AEIGVEC ,levs)
!-------------------------------------------------------------------------------
   use constant, only : akap=>akapa_,rd_
   use dfsvar, only   : trs,grs,phs,tavexy,GRSINV,                             &
                        delsig,VERTCORD,sigmafull,sigma,iope,                  &
                        gama1,gama2,tsigma,valpha,beta,FAVEXY
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                          ::  levs
   real   , dimension(levs)        ::  EIGVAL
   real   , dimension(levs,levs)  ::  EIGVEC, AEIGVEC, VERMAT
!-----local-variable-------------------------------------------------!
   real   , dimension(levs,levs)  ::  WORKIN, GRS_CHARP
   integer, dimension(levs)        ::  INDX
!--------------------------------------------------------------------!
   integer                          ::  ir,is,n1,nk,k,kk,j,mm,ITER
   real                             ::  en1,enk,xtu,xtl,as2,&
                                        pa2,pa1,aglevs,e1
   real                             ::  ErrorP
!-------------------------------------------------------------------------------
   if (iope) write(6,*) ' VERTCORD=', VERTCORD
!
   do ir = 1,levs           ! reset
     do is = 1,levs
       TRS(ir,is)= 0.
       GRS(ir,is)= 0.
       GRS_CHARP(ir,is)= 0.d0
     end do
   end do
!
!----------------------------- row  PHS matrix -----!
   do is = 1,levs
     PHS(is)= delsig(is)
   end do
!==========================================================c
#define EN1 e1
#define ENK e1
!---------------------------------- TRS matrix -----!
   DO ir = 1,levs
!-----------------------------------------  ir -----c
     n1=1
     nk=1
     if(ir.eq. 1  ) n1=0
     if(ir.eq.levs) nk=0
!
     e1= 0.5d0/ delsig(ir)     ! only for equally spaced levels
     en1= 1.0d0/(delsig(ir-n1)+delsig(ir))  ! [dif_+0.5 +dif_-0.5]/2
     enk= 1.0d0/(delsig(ir+nk)+delsig(ir))
!
!    xtu= e1*( -tavexy(ir-n1)*gama1(ir)+tavexy(ir) )
!    xtl= e1*(  tavexy(ir+nk)*gama2(ir)-tavexy(ir) )
     xtu=    ( -tavexy(ir-n1)*gama1(ir)+tavexy(ir) )*EN1*(-1) ! (-1) up K
     xtl=    (  tavexy(ir+nk)*gama2(ir)-tavexy(ir) )*ENK*(-1) !      up K
     pa1= (xtu*tsigma(ir-n1)*n1+ xtl*tsigma(ir))        *(-1) !      up K
     pa2= akap*TAVEXY(ir)
!------------------------------------------------  is -----c
     do is = 1,levs
       TRS(ir,is)= delsig(is)* ( pa1 + pa2 )
     end do
!------------------------------------------------  is -----c
     do is = 1,ir-n1
       TRS(ir,is)= TRS(ir,is) - xtu*delsig(is)*n1 *(-1) ! (-1) up K
     end do
!------------------------------------------------  is -----c
     do is = 1,ir
       TRS(ir,is)= TRS(ir,is) - xtl*delsig(is)    *(-1) ! (-1) up K
     end do
!------------------------------------------------  is -----c
   END DO ! ir=1,levs

!==========================================================c

   IF(VERTCORD.EQ.'LORENZ') THEN
!
!---------------------------------- GRS matrix --------!
!     do 501 ir=1,levs
!     do 501 is=1,levs
! 501        VERMAT(ir,is) = 0
!--------------------
!     do 500 ir=1,levs-1
!            VERMAT(ir  ,ir  ) =valpha(ir)
! 500        VERMAT(ir  ,ir+1) =  beta(ir)
!                    ir  =  levs
!             GRS(ir,ir) =  beta(levs)+valpha(levs)
!     do 530 ir=levs-1,1,-1
!     do 530 is=     1,levs
!------------------------=----------------------------c
! 530         GRS(ir,is) = GRS(ir+1,is) + VERMAT(ir,is)
!------------------------=----------------------------c
!---------------------------------- GRS matrix --------!
!
!----conserving  angular momentum-- GRS matrix --------!
!
!            VERMAT : HERE, work storage 
!-------------------------------
     do 5500 ir=1,levs
     do 5500 is=1,levs
 5500        VERMAT(ir,is) = 0.
!-------------------------------
     do ir =   1,levs-1
       do is = ir+1,levs-1
         VERMAT(ir,is  ) = valpha(is) + beta(is  )
       enddo
       VERMAT(ir,ir  ) = valpha(ir)
       VERMAT(ir,levs) =              beta(levs)
     enddo
!
     do is = 1,levs
       aglevs = 0.
       do ir = 1,levs-1
         aglevs = aglevs + VERMAT(ir,is)*delsig(ir)
       enddo
       GRS(levs,is) = -aglevs + delsig(is)
     enddo
!
     do ir = levs-1,1,-1
        do is = 1,levs
           GRS(ir,is  ) = GRS(ir+1,is  )
        enddo
        GRS(ir,ir  ) = GRS(ir  ,ir  ) + valpha(ir  )
        GRS(ir,ir+1) = GRS(ir  ,ir+1) +   beta(ir+1)
     enddo
!
! conserving  angular momentum-- GRS matrix
!
     do is = 1,levs
       do ir = 1,levs
!
!       *Rgas Dimensionalize by Earth's parameter
!
         GRS(ir,is)= GRS(ir,is) *rd_
       end do
     end do

!========================================================c

   ELSEIF(VERTCORD.EQ.'CHARNEY_PHILLIPS') THEN
!
!-beg-Charney-Phillips grids : sigdot,Phi at full level--!
!
!          GRS_CHARP(ir,is): matrix for full level GEOPOTENTIAL
!
     do ir = 1,levs
       do is = ir,levs
         as2= sigmafull(is+1)**akap - sigmafull(is)**akap
         GRS_CHARP(ir,is)= as2 * sigma(is)**(-akap) / akap
       end do
     end do
!
!          GRS(ir,is): matrix for half level GEOPOTENTIAL
!
     do ir = 1,levs-1
       do is = ir,levs
         GRS(ir,is)= GRS_CHARP(ir,is) *0.5
       end do
       do is = ir+1,levs
         GRS(ir,is)= GRS(ir,is) + GRS_CHARP(ir+1,is) *0.5
       end do
     end do
     do ir = levs,levs
       do is=ir,levs
         GRS(ir,is)= GRS_CHARP(ir,is) *0.5
       end do
     end do
!
!-end-Charney-Phillips grids : sigdot,Phi at same level--!
!
   ENDIF ! IF(VERTCORD)
!
!========================================================c

!
!-b---horizontal mean geopotential without topography
!                                                          ! Level_split
   do k=1,levs
     FAVEXY(k)= 0.
!$dir unroll_and_jam(unroll_factor=100)
#ifdef LINUX_PGI
!pgi$ unroll=n:100
#endif
     do kk=1,levs
       FAVEXY(k)= FAVEXY(k) +GRS(k,kk) * TAVEXY(kk)
     enddo
!             a1= GRS(k, 1) * TAVEXY( 1) + GRS(k, 2) * TAVEXY( 2)
!     &         + GRS(k, 3) * TAVEXY( 3) + GRS(k, 4) * TAVEXY( 4)
!             a2= GRS(k, 5) * TAVEXY( 5) + GRS(k, 6) * TAVEXY( 6)
!     &         + GRS(k, 7) * TAVEXY( 7) + GRS(k, 8) * TAVEXY( 8)
!             a3= GRS(k, 9) * TAVEXY( 9) + GRS(k,10) * TAVEXY(10)
!     &         + GRS(k,11) * TAVEXY(11) + GRS(k,12) * TAVEXY(12)
!             a4= GRS(k,13) * TAVEXY(13) + GRS(k,14) * TAVEXY(14)
!     &         + GRS(k,15) * TAVEXY(15) + GRS(k,16) * TAVEXY(16)
!             a5= GRS(k,17) * TAVEXY(17) + GRS(k,18) * TAVEXY(18)
!     &         + GRS(k,19) * TAVEXY(19) + GRS(k,20) * TAVEXY(20)
!             a6= GRS(k,21) * TAVEXY(21) + GRS(k,22) * TAVEXY(22)
!     &         + GRS(k,23) * TAVEXY(23) + GRS(k,24) * TAVEXY(24)
!             o7= GRS(k,25) * TAVEXY(25) + GRS(k,26) * TAVEXY(26)
!     &         + GRS(k,27) * TAVEXY(27) + GRS(k,28) * TAVEXY(28)
!            FAVEXY(k)= a1 + a2 + a3 + a4 + a5 + a6 + a7
   end do ! k=1,levs
!
!-e---horizontal mean geopotential
! 
   do j=1,levs
     do k=1,levs
       WORKIN(j,k)= GRS(j,k)
     end do
   end do
!
   CALL dfs_inverse_matrix( WORKIN, levs, levs, INDX, AEIGVEC )
!
   do j=1,levs
     do k=1,levs
       GRSINV(j,k)= AEIGVEC(j,k) ! AEIGVEC(j,k)= temporary storage
     end do
   end do
!
!----- VERMAT for semi-implicit method -----------
!
   do ir=1,levs
     do is=1,levs
       VERMAT(ir,is) = 0
       do mm=1,levs
         VERMAT(ir,is) = VERMAT(ir,is) + GRS(ir,mm)*TRS(mm,is)
       end do
       VERMAT(ir,is) = VERMAT(ir,is) + TAVEXY(ir)*PHS(is)*rd_
                           ! *Rgas Dimensionalize by Earth's parameter
     end do
   end do
!
!----- EigenValues and -Vectors of VERMAT --------
!
   ErrorP= 1.0D-12
   ITER= 100    
   CALL DFS_EIGEN_VECTOR(VERMAT, levs, EIGVEC, AEIGVEC, ErrorP, ITER)  
!
   do ir=1,levs
     EIGVAL(ir)= VERMAT(ir,ir)
   end do
!
!===========
!
   RETURN
   END SUBROUTINE dfs_vertical_descr
!-------------------------------------------------------------------------------
#else /* HYBRID */
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_vertical_descr_hybrid(VERMAT,EIGVAL,EIGVEC,AEIGVEC,levs,ak5,bk5)
!-------------------------------------------------------------------------------
   use constant, only : akap=>akapa_,rd_
   use dfsvar, only   : trs,grs,phs,gv,tavexy,GRSINV,VERTCORD,iope
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer                          ::  levs
   real   , dimension(levs)        ::  EIGVAL
   real   , dimension(levs,levs)  ::  EIGVEC, AEIGVEC, VERMAT
!
!-----local-variable-------------------------------------------------!
!
   real   , dimension(levs,levs)  ::  WORKIN
   integer, dimension(levs)        ::  INDX(levs)
!--------------------------------------------------------------------!
   integer                          ::  ir,is,n1,nk,k,kk,j,mm,ITER
   real                             ::  ErrorP
!
   integer                          ::  icol, irow, icolbeg
   real                             ::  psref, factor
   real, dimension(levs+1)         ::  ak5,bk5,pk5ref
   real, dimension(levs)           ::  dpkref                                 ,&
                                        alfaref                               ,&
                                        tref                                  ,&
                                        vecm
   real, dimension(levs,levs)     ::  yecm                                    ,&
                                        tecm
!-------------------------------------------------------------------------------
   if (iope) write(6,*) ' VERTCORD=', VERTCORD
!
   do ir = 1,levs           ! reset
     do is = 1,levs
       TRS(ir,is)= 0.
       GRS(ir,is)= 0.
     end do
   end do
!
! psref,tref
!
   psref=80.
   do k = 1,levs
     tref(k)=300.
   enddo
!
! pk5ref
!
   do k = 1,levs+1
     pk5ref(k)=ak5(k)+bk5(k)*psref
   enddo
   do k = 1,levs
     dpkref(k)=pk5ref(k+1)-pk5ref(k)
   enddo
!
! alfaref
!
   alfaref(1)=log(2.)
   do k = 2,levs
     alfaref(k)=1.-(pk5ref(k)/dpkref(k))*log(pk5ref(k+1)/pk5ref(k))
   enddo
!
! vecm
!
   do k = 1,levs
     vecm(k)=dpkref(k)/psref
   enddo
!
! tecm
!
   tecm=0.
   do irow = 1,levs
     tecm(irow,irow)=akap*tref(irow)*alfaref(irow)
     do icol = 1,irow-1
       factor=(akap*tref(irow)/dpkref(irow))*            &
              log(pk5ref(irow+1)/pk5ref(irow))
       tecm(irow,icol)=factor*dpkref(icol)
     enddo
   enddo
!
! yecm
!
   yecm=0.
   do irow = 1,levs
     yecm(irow,irow)=alfaref(irow)*rd_
     icolbeg=irow+1
     if(icolbeg.le.levs)then
       do icol = icolbeg,levs
         yecm(irow,icol)=rd_*log( pk5ref(icol+1)/pk5ref(icol) )
       enddo
     endif
   enddo
!
! PHS, TRS
!
   do j = 1,levs
   PHS(j)= vecm(levs+1-j)
   GV(j) = rd_*tref(j)
     do k = 1,levs
       GRS(k,j)=yecm(levs+1-k,levs+1-j)
       TRS(k,j)=tecm(levs+1-k,levs+1-j)
     enddo
   enddo
!
! inverse matrix
!
   do j = 1,levs
     do k = 1,levs
       WORKIN(j,k)= GRS(j,k)
     end do
   end do
!
   CALL dfs_inverse_matrix( WORKIN, levs, levs, INDX, AEIGVEC )
!
   do j = 1,levs
     do k = 1,levs
       GRSINV(j,k)= AEIGVEC(j,k) ! AEIGVEC(j,k)= temporary storage
     end do
   end do
!
!----- VERMAT for semi-implicit method -----------
!
   do ir = 1,levs
     do is = 1,levs
       VERMAT(ir,is) = 0.
       do mm = 1,levs
         VERMAT(ir,is) = VERMAT(ir,is) + GRS(ir,mm)*TRS(mm,is)
       end do
       VERMAT(ir,is) = VERMAT(ir,is) + GV(ir)*PHS(is)
                           ! *Rgas Dimensionalize by Earth's parameter
     end do
   end do
!
!----- EigenValues and -Vectors of VERMAT --------
!
   ErrorP= 1.0D-12
   ITER= 100
   CALL DFS_EIGEN_VECTOR(VERMAT, levs, EIGVEC, AEIGVEC, ErrorP, ITER)
   do ir = 1,levs
     EIGVAL(ir)= VERMAT(ir,ir)
   end do
!
!===========
!
   RETURN
   END SUBROUTINE dfs_vertical_descr_hybrid
!-------------------------------------------------------------------------------
#endif /* ~HYBRID end */

!-------------------------------------------------------------------------------
   SUBROUTINE DFS_EIGEN_VECTOR(A,N,TR,TL,EP ,ITER)
!-------------------------------------------------------------------------------
   use dfsvar, only : iope
!-------------------------------------------------------------------------------
!
!  - calculate eigrnvectior TR and inversion TL of given matrix A(N,N).
!  - input argument 
!      EP     : 1.0D-12 or..
!      ITER   : max iteration number to get EP
!
!-------------------------------------------------------------------------------
   integer  ::  n
   real     ::  A(N,N)          ! on entry input matrix
                                ! on output : diagnal eigvalue
   real     ::  TL(N,N),TR(N,N) ! eigvec,inv_eigvec
!-------------------------------------------------------------------------------
   DO J = 1,N                                                               
     DO I = 1,N                                                               
       TL(I,J)=0.0D0
       TR(I,J)=0.0D0
     enddo
   enddo
!
   DO I = 1,N
     TL(I,I)=1.0D0
     TR(I,I)=1.0D0
   enddo
!
   ITERX=0                                                                   
   DO IT = 1,ITER                                                          
!                                                                               
!------ CONVERGENCE CHECK-----------------------------------------------        
!                                                                               
     DO 200 I = 1,N-1                                                            
     DO 200 J = I+1,N                                                            
       IF(ABS(A(I,J)+A(J,I)).GT.EP)GOTO 101
       IF(ABS(A(I,J)-A(J,I)).LE.EP)GOTO 200
       IF(ABS(A(I,I)+A(J,J)).GT.EP)GOTO 101
!
     200 CONTINUE                                                                  
!
     GOTO 999                                                                  
!
     101 CONTINUE                                                                  
!
!
     ITERX=ITERX+1                                                             
     DO K = 1,N-1                                                            
       DO M = K+1,N                                                            
         H=0.;G=0.;HJ=0.;YH=0.
         DO 110 I = 1,N
           TE=A(I,K)*A(I,K);TEE=A(I,M)*A(I,M)
           YH=YH+TE-TEE
!
           IF(I.EQ.K.OR.I.EQ.M)GOTO 110
!
           H  =H  +A(K,I)*A(M,I)-A(I,K)*A(I,M)
           TEP=TE +A(M,I)**2
           TEM=TEE+A(K,I)**2
           G  =G  +TEP+TEM
           HJ =HJ -TEP+TEM
!
         110 CONTINUE
!                                                                               
         H=H*2.
         D=A(K,K)-A(M,M)
         C=A(K,M)+A(M,K)
         E=A(K,M)-A(M,K)
!
         IF(ABS(C).EQ.0.)THEN
           CX=1.0D0;sx=0.0D0
         ELSE
           COT2X=D/C
           SIG=1.0D0
           IF(COT2X.LT.0.0D0)SIG=-1.0
           COTX=COT2X+SIG*SQRT(1.0+COT2X**2)
           sx  =SIG/SQRT(1.0+COTX**2)
           CX  =sx*COTX
         END IF
!                                                                               
         IF(YH.LT.0.)THEN
           TEM=CX
           CX =sx
           sx =-TEM
         END IF
!
         COS2X=CX*CX-sx*sx
         SIN2X=2.0*sx*CX
         D=D*COS2X+C *SIN2X
         H=H*COS2X-HJ*SIN2X
         TANHY=(E*D-H/2.)/(G+2.0*(E*E+D*D))
         CHY=1./SQRT(1.-TANHY**2)
         SHY=CHY*TANHY
         C1=CHY*CX-SHY*sx
         C2=CHY*CX+SHY*sx
         S1=CHY*sx+SHY*CX
         S2=-CHY*sx+SHY*CX
!
         DO I=1,N
           AKI=A(K,I)
           AMI=A(M,I)
           A(K,I)=C1*AKI+S1*AMI
           A(M,I)=S2*AKI+C2*AMI
           TKI=TL(K,I)
           TMI=TL(M,I)
           TL(K,I)=C1*TKI+S1*TMI
           TL(M,I)=S2*TKI+C2*TMI
         enddo
!
         DO I=1,N
           AKI=A(I,K)
           AMI=A(I,M)
           A(I,K)=C2*AKI-S2*AMI
           A(I,M)=-S1*AKI+C1*AMI
           TKI=TR(I,K)
           TMI=TR(I,M)
           TR(I,K)=C2*TKI-S2*TMI
           TR(I,M)=-S1*TKI+C1*TMI
         enddo
       enddo
     enddo
   enddo ! IT loop
#ifdef DBG
   WRITE(6,*)'ITERATIONS ARE OVER ',ITER
#endif
   RETURN
 !
 999  continue
 !
#ifdef DBG
   if (iope) WRITE(6,210)ITERX
 210 FORMAT(1H ,'ITERATION=',I10)
#endif
   RETURN
   END SUBROUTINE DFS_EIGEN_VECTOR
!-------------------------------------------------------------------------------
