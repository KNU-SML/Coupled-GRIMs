!
   subroutine co2_tau_interp(rat,ir)
!-------------------------------------------------------------------------------
! **********************************************************************        
!                                                                               
!                                                                               
!            the transmission function between p1 and p2 is assumed to          
!       the  functional form                                                    
!                     tau(p1,p2)= 1.0-sqrt(c*log(1.0+x*path)),                  
!               where                                                           
!                     path(p1,p2)=((p1-p2)**2)*(p1+p2+core)/                    
!                                 (eta*(p1+p2+core)+(p1-p2))                    
!                                                                               
!                                                                               
!        the parameters c and x are functions of p2, and are to be deter        
!        while core is a prespecified number.eta is a function of the th        
!        product (cx);it is obtaited iteratively. the derivation of all         
!        values will be explained in a forthcoming paper.                       
!            subroutine co2_tau_interp determines c(i) and x(i) by using the act  
!        values of tau(p(i-2),p(i)) and tau(p(i-1),p(i)) and the previou        
!        iteration value of eta.                                                
!             define:                                                           
!                patha=path(p(i),p(i-2),core,eta)                               
!                pathb=path(p(i),p(i-1),core,eta);                              
!        then                                                                   
!                r=(1-tau(p(i),p(i-2)))/(1-tau(p(i),p(i-1)))                    
!                 = sqrt(log(1+x*patha)/log(1+x*pathb)),                        
!        so that                                                                
!                r**2= log(1+x*patha)/log(1+x*pathb).                           
!        this equation can be solved by newton s method for x and then t        
!        result used to find c. this is repeated for each value of i gre        
!        than 2 to give the arrays x(i) and c(i).                               
!             newton s method for solving the equation                          
!                 f(x)=0                                                        
!        makes use of the loop xnew= xold-f(xold)/fprime(xold).                 
!        this is iterated 20 times, which is probably excessive.                
!        the first guess for eta is 3.2e-4*exp(-p(i)/1000),which has            
!        been found to be fairly realistic by experiment; we iterate 5 t        
!        (again,probably excessively) to obtain the values for c,x,eta t        
!        used for interpolation.                                                
!           there are several possible pitfalls:                                
!              1) in the course of iteration, x may reach a value which         
!                 1+x*patha negative; in this case the iteration is stop        
!                 and an error message is printed out.                          
!              2) even if (1) does not occur, it is still possible that         
!                 be negative and large enough to make 1+x*path(p(i),0,c        
!                 negative. this is checked for in a final loop, and if         
!                 a warning is printed out.                                     
!                                                                               
!  *********************************************************************        
!....                                                                           
!     implicit double precision (a-h,o-z)                                       
!
!-------------------------------------------------------------------------------
   use comco2,  only : pa
   use comco2,  only : transa
   use comco2,  only : xa, ca, eta, sexpv, core, uexp, sexp
!-------------------------------------------------------------------------------
!  real  ::  rat                                                            
!  real  ::  pa,core,transa,path,uexp,sexp,eta,sexpv                              
!  real  ::  pa2                                                                  
!  real  ::  path0(109),etap(109),xap(109),cap(109)                          
!  real  ::  sinv(4)                                                         
   real  ::  rat                                                            
   real  ::  path
   real  ::  pa2                                                                  
   real  ::  path0(109),etap(109),xap(109),cap(109)                          
   real  ::  sinv(4)                                                         
   data sinv/2.74992,2.12731,4.38111,0.0832926/                              
!
!nov89   dimension sinv(3)                                                      
!nov89   data sinv/2.74992,2.12731,4.38111/                                     
!o222  old code used 2.7528 rather than 2.74992 ---k.a.c. october 1988          
!o222   when 2.7528 was used,we exactly reproduced the mrf co2 arrays           
!
   core=5.000                                                                
   uexp=0.90                                                                 
   p0=0.7                                                                    
   do i = 1,109                                                            
     pa2=pa(i)*pa(i)                                                           
     sexpv(i)=.505+2.0e-5*pa(i)+.035*(pa2-.25)/(pa2+.25)                       
   enddo
   do i = 1,109                                                            
     eta(i)=3.2e-4*exp(-pa(i)/500.)                                            
     etap(i)=eta(i)                                                            
   enddo
   do np = 1,10                                                           
     do i = 3,109                                                           
       sexp=sexpv(i)                                                             
       r=(1.0-transa(i,i-2))/(1.0-transa(i,i-1))                                 
       rexp=r**(uexp/sexp)                                                       
       patha=(path(pa(i),pa(i-2),core,eta(i)))**uexp                             
       pathb=(path(pa(i),pa(i-1),core,eta(i)))**uexp                             
       xx=2.0*(pathb*rexp-patha)/(pathb*pathb*rexp-patha*patha)                  
       do ll = 1,20                                                           
         f1=log(1.0+xx*patha)                                                      
         f2=log(1.0+xx*pathb)                                                      
         f=f1/f2-rexp                                                              
         fprime=(f2*patha/(1.0+xx*patha)-f1*pathb/(1.0+xx*pathb))/             &
                   (f2*f2)                                                               
         xx=xx-f/fprime                                                            
         check=1.0+xx*patha                                                        
         if (check) 1020,1020,1025                                                 
         1020 continue                                                                  
         write (6,360) i,ll,check                                                  
         360   format (' error,i=',i3,'ll=',i3,'check=',f20.10)                          
         stop                                                                      
         1025  continue                                                                  
       enddo
       ca(i)=(1.0-transa(i,i-2))**(uexp/sexp)/                                 &
               (log(1.0+xx*patha)+1.0e-20)                                              
       xa(i)=xx                                                                  
     enddo
     xa(2)=xa(3)                                                               
     xa(1)=xa(3)                                                               
     ca(2)=ca(3)                                                               
     ca(1)=ca(3)                                                               
     do i = 3,109                                                           
       path0(i)=(path(pa(i),0.,core,eta(i)))**uexp                               
       path0(i)=1.0+xa(i)*path0(i)                                               
       if (path0(i).lt.0.) write (6,361) i,path0(i),xa(i)                        
     enddo
     do i = 1,109                                                           
       sexp=sexpv(i)                                                             
       etap(i)=eta(i)                                                            
       eta(i)=(sinv(ir)/rat)**(1./sexp)*                                       &
                (ca(i)*xa(i))**(1./uexp)                                                
     enddo
!                                                                               
!     the eta formulation is detailed in schwarzkopf and fels(1985).            
!        the quantity sinv=(g*deltanu)/(rco2*d*s)                               
!      in cgs units,with d,the diffusicity factor=2, and                        
!      s,the sum of co2 line strengths over the 15um co2 band                   
!       also,the denominator is multiplied by                                   
!      1000 to permit use of mb units for pressure.                             
!        s is actually weighted by b(250) at 10 cm-1 wide intervals,in          
!      order to be consistent with the methods used to obtain the lbl           
!      1-band consolidated trancmission functions.                              
!      for the 490-850 interval (deltanu=360,ir=1) sinv=2.74992.                
!      (slightly different from 2.7528 used in earlier versions)                
!      for the 490-670 interval (ir=2) sinv=2.12731                             
!      for the 670-850 interval (ir=3) sinv=4.38111                             
!      for the 2270-2380 interval (ir=4) sinv=0.0832926                         
!      sinv has been obtained using the 1982 afgl catalog for co2               
!        rat is the actual co2 mixing ratio in units of 330 ppmv,               
!      letting use of this formulation for any co2 concentration.               
!                                                                               
!     write (6,366) (np,i,ca(i),xa(i),eta(i),sexpv(i),i=1,109)                  
!366   format (2i4,4e20.12)                                                     
!
   enddo
!
361  format (' **warning:** 1+xa*path(pa(i),0) is negative,i= ',i3,/           &
                20x,'path0(i)=',f16.6,' xa(i)=',f16.6)                                   
!
   return                                                                    
   end subroutine co2_tau_interp
!-------------------------------------------------------------------------------
