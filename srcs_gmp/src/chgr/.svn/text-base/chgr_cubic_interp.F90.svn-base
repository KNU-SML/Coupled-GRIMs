#include <define.h>
   subroutine chgr_cubic_interp
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [chgr_cubic_interp]
!      |
!      |--- [chgr_cubic_coeff] *
!      |--- [chgr_cubic_compute] *
!      |--- [chgr_cubic_spline] *
!
!-------------------------------------------------------------------------------
   end subroutine chgr_cubic_interp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_cubic_coeff(n,p,q,y)
!-------------------------------------------------------------------------------
!
! subprogram: chgr_cubic_coeff
! 
! abstract: calculate weights for cubic spline interpolation
!
! usage:   call chgr_cubic_coeff(n,p,q,y)
!   input argument list:            
!     n        - number ofc     n        - number of observations             
!     y        - function y(p) to be interpolated                            
!     p        - coor (independent variables)                               
!     q,ovh,sh - contain weights                                                
!
!-------------------------------------------------------------------------------
   save
   common/spl/ ovh(100),sh(100),iflag,jflag,theta
   real  ::  c(100),b(100),p(n),q(n),d(100),y(n),                              &
             tchs(100),tb(100)   
!-------------------------------------------------------------------------------
   nu=n-1                                                                    
!
   if (iflag.ne.0) go to 99                                                  
!
   do i = 2,n                                                                
     im=i-1                                                                    
     him=p(i)-p(im)                                                            
     if(him .eq. 0.) him=0.001                                                 
     ovh(im)=1/him                                                             
     th=theta*him                                                              
     sh(im)=1/sinh(th)                                                        
     tchs(im)=(theta*cosh(th))*sh(im)-ovh(im)                                  
   enddo
!
   do i = 2,nu                                                               
     im=i-1                                                                    
     c(im)=ovh(im)-theta*sh(im)                                                
     b(i)=tchs(im)+tchs(i)                                                     
   enddo
!
   c(n-1)=ovh(n-1)-theta*sh(n-1)                                             
   b(1)=tchs(1)                                                              
   b(n)=tchs(n-1)                                                            
   tb(1)=c(1)/b(1)                                                           
   do i = 2,nu                                                               
     tb(i)=c(i)/(b(i)-c(i-1)*tb(i-1))                                          
   enddo
   d(1)=0.                                                                   
   d(n)=0.                                                                   
!
!   iflag=1                                                                   
!
 99   yp=(y(2)-y(1))*ovh(1)                                                     
!
   do i = 2,nu                                                               
     ynow=(y(i+1)-y(i))*ovh(i)                                                 
     d(i)=ynow-yp                                                              
     yp=ynow                                                                   
   enddo
   q(1)=d(1)/b(1)                                                            
   do i = 2,n                                                                
     im=i-1                                                                    
     q(i)=(d(i)-c(im)*q(im))/(b(i)-c(im)*tb(im))                              
   enddo
   do i = 1,nu                                                               
     ii=n-i                                                                    
     q(ii)=q(ii)-tb(ii)*q(ii+1)                                                
   enddo
   return                                                                    
   end subroutine chgr_cubic_coeff 
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_cubic_compute(q,vst,vinc,p,y,m,n,yn)        
!-------------------------------------------------------------------------------
!
! abstract: calculates interpolated values of function using
!   weights calcuated by chgr_cubic_coeff
! 
! usage:    call chgr_cubic_compute(q,vst,vinc,p,y,m,n,yn)
!   input argument list:            
!     m        - number of grid points desired                             
!     yn       - interpolated value of y(p) at grid points                
!     dyn      - analytic value of derivative at grid point               
!     vst      - starting location of interval                           
!     vinc     - increment of regular grid                              
!     n        - number of points input from chgr_cubic_coeff      
!     q        - weight function from chgr_cubic_coeff     
!     ovh,sh   - weights from chgr_cubic_coeff         
!
!-------------------------------------------------------------------------------
   save
   common/spl/ ovh(100),sh(100),iflag,jflag,theta
   real  ::  q(n),p(n),y(n),yn(m),vinc(m)                                    
!-------------------------------------------------------------------------------
   vs=vst                                                                    
   j=1                                                                       
   i=1                                                                       
!
1  if(vs .gt. p(1)) go to 9                                                  
!    
   yn(i)=99999.9                                                             
   i=i+1                                                                     
   vs =vs +vinc(i-1)                                                         
!
   go to 1                                                                   
!                                                                               
9  j=j+1                                                                     
!
   if (j.gt.n) go to 90                                                      
8  if (vs .gt.p(j)) go to 9                                                  
!
   jm=j-1                                                                    
   yn(i)=q(jm)*sinh(theta*(p(j)-vs))*sh(jm)+(y(jm)-q(jm))*(p(j)-vs)            &
         *ovh(jm)+q(j)*sinh(theta*(vs-p(jm)))*sh(jm)+(y(j)-q(j))*(vs-          &
         p(jm))*ovh(jm)                                                           
   i=i+1                                                                     
   vs =vs +vinc(i-1)                                                         
!
   if(m.ge.i) go to 8                                                        
!
   return                                                                    
!
90 yn(i)=99999.9                                                             
!
   i=i+1                                                                     
!
   if (m.ge.i) go to 90                                                      
!
   return                                                                    
   end subroutine chgr_cubic_compute
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE chgr_cubic_spline(ZI,PILN,KMAX,Z,PLN,LMAX,IDX)                   
!-------------------------------------------------------------------------------
!
! subprogram: spline
!
! abstract: INTERPOLATION USING A CUBIC NONPERIODIC SPLINE
!   
!     ZI(I)........INTERPOLATED VALUE AT PILN(I)                       
!     PILN(I)......LN(P) COORDINATE IN DECREASING ORDER                
!     KMAX.........NUMBER OF INTERPOLATION POINTS                      
!     Z(I).........DATA AT PLN(I)                                      
!     PLN(I).......LN(P) COORDINATE  IN DECREASING ORDER               
!     LMAX.........NUMBER OF DATA POINTS                               
!     SM(I)........SECOND DERIVATIVES AT DATA POINTS                   
!     *******   1   INTERPOLATION FOR HEIGHT                           
!     * IDX *   2   INTERPOLATION FOR WIND                             
!     *******   3   INTERPOLATION FOR TEMPERATURE ( CUBIC SPLINE FOR   
!                   HEIGHT IS DIFFERENCIATED WITH RESPECT TO LOG P     
!                   --- HYDROSTATIC RELATION )                         
!
!-------------------------------------------------------------------------------
   REAL   ::   ZI(KMAX),PILN(KMAX),Z(LMAX),PLN(LMAX)                  
   REAL   ::   SM(LMAX),H(LMAX),AL(LMAX),AM(LMAX),AP(LMAX),C(LMAX)    
!-------------------------------------------------------------------------------
!
   LM1=LMAX-1                                                       
   G=9.8
   R=287.04                                                         
   GR=-G/R                                                          
   DO I = 2,LMAX                                                  
     H(I)=PLN(I)-PLN(I-1)                                             
   ENDDO
   DO I = 2,LM1                                                   
     AL(I)=0.5*H(I+1)/(H(I)+H(I+1))                                   
     AM(I)=0.5-AL(I)                                                  
   ENDDO
!
   IF( IDX.EQ.2 ) GO TO 80                                          
!
!     END CONDITION FOR HEIGHT AND TEMPERATURE (LAPSE RATE IS CONSTANT)
!     SM(1)=SM(2)     ;     SM(LMAX-1)=SM(LMAX)                        
!
     AL(1)=-1.                                                        
!ORG AM(LMAX)=-1.                                                     
     AM(LMAX)=-0.5                                                    
!
   GO TO 90                                                         
!
80 CONTINUE                                                         
!
!     END CONDITION FOR WIND ( WIND SHEAR  IS CONSTANT )               
!     SM(1)=0.     ;     SM(LMAX)=0.                                  
!
   AL(1)=0.                                                          
   AM(LMAX)=0.                                                      
90 CONTINUE                                                     
   AL(LMAX)=0.0                                                   
   DO I = 2,LMAX                                               
     AP(I)=1.0/(1.0-AL(I-1)*AM(I))                                
     AL(I)=AL(I)*AP(I)                                        
   ENDDO
!                                                             
   C(1)=0.                                                   
   C(LMAX)=0.                                               
   DO I = 2,LM1                                          
     C(I)=3.0*((Z(I+1)-Z(I))/H(I+1)-(Z(I)-Z(I-1))/H(I))                        &
     /(H(I)+H(I+1))                                       
   ENDDO
!
!     FORWARD SUBSTITUTION                              
!
   DO I = 2,LMAX                                     
     C(I)=(C(I)-C(I-1)*AM(I))*AP(I)                  
   ENDDO
   SM(LMAX)=C(LMAX)                                  
!
!     BACKWARD SUBSTUTUTION                         
!
   DO K = 1,LM1                                  
     I=LMAX-K                                       
     SM(I)=C(I)-AL(I)*SM(I+1)                   
   ENDDO
!
!     INTERPOLATION                             
!
   IB=2                                        
   DO 500 L = 1,KMAX                            
     X=PILN(L)                                 
     DO I=IB,LMAX                         
!
       IF(X.GE.PLN(I)) GOTO 310                
!
     ENDDO
     I=LMAX                                
310  IB=I                              
!
     IF(IDX.EQ.3) GO TO 400              
     ZI(L)=(PLN(I)-X)/H(I)*                                                    &
            (Z(I-1)-SM(I-1)/6.*(X-PLN(I-1))*(H(I)+PLN(I)-X))                   &
            +(X-PLN(I-1))/H(I)*                                                &
            (Z(I)-SM(I)/6.*(PLN(I)-X)*(H(I)+X-PLN(I-1)))                 
!
     GO TO 500
400  CONTINUE                                                
!
!     DIFFERENTIAL CALCULUS OF CUBIC SPLINE FOR HEIGHT       
!     GR=-G/R ...  COEFFICIENT IN HYDROSTATIC EQUATION      
!
     ZI(L)=SM(I-1)*(-(PLN(I)-X)**2/(2.*H(I))+H(I)/6.)                          &
           +SM(I)*((X-PLN(I-1))**2/(2.*H(I))-H(I)/6.)                          &
           +(Z(I)-Z(I-1))/H(I)                              
     ZI(L)=ZI(L)*GR
!
500  CONTINUE                                         
!
   RETURN                                             
   END SUBROUTINE chgr_cubic_spline
!-------------------------------------------------------------------------------
