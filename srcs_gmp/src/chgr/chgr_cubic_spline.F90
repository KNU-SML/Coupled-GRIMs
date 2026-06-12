!
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
     AL(1)=-1.                                                        
!ORG AM(LMAX)=-1.                                                     
     AM(LMAX)=-0.5                                                    
!
   GO TO 90                                                         
!
80 CONTINUE                                                         
!     END CONDITION FOR WIND ( WIND SHEAR  IS CONSTANT )               
!     SM(1)=0.     ;     SM(LMAX)=0.                                  
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
     DO I = IB,LMAX                         
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
      ZI(L)=SM(I-1)*(-(PLN(I)-X)**2/(2.*H(I))+H(I)/6.)                         &
           +SM(I)*((X-PLN(I-1))**2/(2.*H(I))-H(I)/6.)                          &
           +(Z(I)-Z(I-1))/H(I)                              
      ZI(L)=ZI(L)*GR                                       
!
500 CONTINUE                                         
!
   RETURN                                             
   END SUBROUTINE chgr_cubic_spline
!-------------------------------------------------------------------------------
