!
   SUBROUTINE chgr_linear_interp(ZINT,PINTLG,KMAX,ZDAT,PDATLG,LMAX)             
!-------------------------------------------------------------------------------
!
! subprogram: chgr_linear_interp         
!
! abstract: 1-DIMENSIONAL LINEAR INTERPOLATION                             
! ***  PDATLG MUST BE IN DECREASING ORDER
!
!-------------------------------------------------------------------------------
   REAL      ::  ZINT(KMAX),PINTLG(KMAX)                                
   REAL      ::  ZDAT(LMAX),PDATLG(LMAX)                               
   INTEGER   ::   LLL(KMAX)
!-------------------------------------------------------------------------------
!                                                                 
   DO K = 1,KMAX
     LLL(K)=LMAX+1
   ENDDO
!                                                                
   DO L = LMAX,1,-1
     DO K=1,KMAX
       IF(PINTLG(K).GE.PDATLG(L)) LLL(K)=L
     ENDDO
   ENDDO
!                                                              
   DO 10 K = 1,KMAX
     LL = LLL(K)
     IF(LL.EQ.LMAX+1) THEN
       ZINT(K)=ZDAT(LMAX)
       GOTO 10       
     END IF
     IF(LL.EQ.1) THEN                                           
       ZINT(K)=ZDAT(1)                                       
       GOTO 10                                              
     END IF                                                  
!                                                         
     LL1=LL-1                                              
     C1=PDATLG(LL)-PDATLG(LL1)                            
     C2=PINTLG(K)-PDATLG(LL1)                            
     C3=PDATLG(LL)-PINTLG(K)                            
     ZINT(K)=(ZDAT(LL)*C2+ZDAT(LL1)*C3)/C1             
!                                                  
   10 CONTINUE                                    
!                                                
   RETURN                                       
   END SUBROUTINE chgr_linear_interp                                        
!-------------------------------------------------------------------------------
