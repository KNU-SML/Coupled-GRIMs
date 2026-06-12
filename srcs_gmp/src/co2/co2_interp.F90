!
   subroutine co2_interp                                                         
!-------------------------------------------------------------------------------
   use comco2,  only : p1, p2, trnslo, ia, ja, n
   use comco2,  only : pa
   use comco2,  only : transa
   use comco2,  only : xa, ca, eta, sexpv, core, uexp, sexp
!-------------------------------------------------------------------------------
!....                                                                           
!     implicit double precision (a-h,o-z)                                       
!
   do 70 l = 1,109                                                             
     ip1=l                                                                     
     if (p2-pa(l)) 65,65,70                                                    
   70 continue                                                                  
!
   65 i=ip1-1                                                                   
   if (ip1.eq.1) ip1=2                                                       
   if (i.eq.0) i=1                                                           
!
   do 80 l = 1,109                                                             
     jp1=l                                                                     
     if (p1-pa(l)) 75,75,80                                                    
   80 continue                                                                  
!
   75 j=jp1-1                                                                   
   if (jp1.eq.1) jp1=2                                                       
   if (j.eq.0) j=1                                                           
   jjj=j                                                                     
   iii=i                                                                     
   j=jjj                                                                     
   jp1=j+1                                                                   
   i=iii                                                                     
   ip1=i+1                                                                   
!
!  determine etap,the value of eta to use by linear interpolation               
!    for peta(=0.5*(p1+p2))                                                     
!
   peta=p2                                                                   
   do 90 l = 1,109                                                             
     ietap1=l                                                                  
     if (peta-pa(l)) 85,85,90                                                  
   90    continue                                                                  
!
   85    ieta=ietap1-1                                                             
   if (ietap1.eq.1) ietap1=2                                                 
   if (ieta.eq.0) ieta=1                                                     
   etap=eta(ieta)+(peta-pa(ieta))*(eta(ietap1)-eta(ieta))/                     &
        (pa(ietap1)-pa(ieta))                                                    
   sexp=sexpv(ieta)+(peta-pa(ieta))*(sexpv(ietap1)-                            &
        sexpv(ieta))/ (pa(ietap1)-pa(ieta))                                      
   pipmpi=pa(ip1)-pa(i)                                                      
   up2p1=(path(p2,p1,core,etap))**uexp                                       
   if (i-j) 126,126,127                                                      
   126 continue                                                                  
   trip=(ca(ip1)*log(1.0+xa(ip1)*up2p1))**(sexp/uexp)                        
   tri=(ca(i)*log(1.0+xa(i)*up2p1))**(sexp/uexp)                             
   trnslo=1.0-((pa(ip1)-p2)*tri+(p2-pa(i))*trip)/pipmpi                      
   go to 128                                                                 
   127 tij=transa(i,j)                                                           
   tipj=transa(i+1,j)                                                        
   tijp=transa(i,j+1)                                                        
   tipjp=transa(i+1,j+1)                                                     
   uij=(path(pa(i),pa(j),core,etap))**uexp                                   
   uipj=(path(pa(i+1),pa(j),core,etap))**uexp                                
   uijp=(path(pa(i),pa(j+1),core,etap))**uexp                                
   uipjp=(path(pa(i+1),pa(j+1),core,etap))**uexp                             
   prodi=ca(i)*xa(i)                                                         
   prodip=ca(i+1)*xa(i+1)                                                    
   prod=((pa(i+1)-p2)*prodi+(p2-pa(i))*prodip)/pipmpi                        
   xint=((pa(i+1)-p2)*xa(i)+(p2-pa(i))*xa(i+1))/pipmpi                       
   cint=prod/xint                                                            
   aij=(cint*log(1.0+xint*uij))**(sexp/uexp)                                 
   aijp=(cint*log(1.0+xint*uijp))**(sexp/uexp)                               
   aipj=(cint*log(1.0+xint*uipj))**(sexp/uexp)                               
   aipjp=(cint*log(1.0+xint*uipjp))**(sexp/uexp)                             
   eij=tij+aij                                                               
   eipj=tipj+aipj                                                            
   eijp=tijp+aijp                                                            
   eipjp=tipjp+aipjp                                                         
   dtdj=(eijp-eij)/(pa(j+1)-pa(j))                                           
   dtdpj=(eipjp-eipj)/(pa(j+1)-pa(j))                                        
   epip1=eij+dtdj*(p1-pa(j))                                                 
   epipp1=eipj+dtdpj*(p1-pa(j))                                              
   epp2p1=((pa(i+1)-p2)*epip1+(p2-pa(i))*epipp1)/pipmpi                      
   trnslo=epp2p1-(cint*log(1.0+xint*up2p1))**(sexp/uexp)                     
   if (i.ge.108.or.j.ge.108) go to 350                                       
   if (i-j-2) 350,350,355                                                    
   355   continue                                                                  
   tip2j=transa(i+2,j)                                                       
   tip2jp=transa(i+2,j+1)                                                    
   ti2j2=transa(i+2,j+2)                                                     
   tijp2=transa(i,j+2)                                                       
   tipjp2=transa(i+1,j+2)                                                    
   uip2j=(path(pa(i+2),pa(j),core,etap))**uexp                 
   uijp2=(path(pa(i),pa(j+2),core,etap))**uexp                  
   uipjp2=(path(pa(i+1),pa(j+2),core,etap))**uexp                
   ui2j2=(path(pa(i+2),pa(j+2),core,etap))**uexp                  
   uip2jp=(path(pa(i+2),pa(j+1),core,etap))**uexp                  
   aijp2=(cint*log(1.0+xint*uijp2))**(sexp/uexp)                             
   aipjp2=(cint*log(1.0+xint*uipjp2))**(sexp/uexp)                           
   aip2j=(cint*log(1.0+xint*uip2j))**(sexp/uexp)                             
   aip2jp=(cint*log(1.0+xint*uip2jp))**(sexp/uexp)                           
   ai2j2=(cint*log(1.0+xint*ui2j2))**(sexp/uexp)                             
   eip2j=tip2j+aip2j                                                         
   eip2jp=tip2jp+aip2jp                                                      
   eijp2=tijp2+aijp2                                                         
   eipjp2=tipjp2+aipjp2                                                      
   ei2j2=ti2j2+ai2j2                                                         
   call co2_quard_interp(pa(j),pa(j+1),pa(j+2),eij,eijp,eijp2,p1,ei)   
   call co2_quard_interp(pa(j),pa(j+1),pa(j+2),eipj,eipjp,eipjp2,p1,ep) 
   call co2_quard_interp(pa(j),pa(j+1),pa(j+2),eip2j,eip2jp,ei2j2,p1,ep2)
   call co2_quard_interp(pa(i),pa(i+1),pa(i+2),ei,ep,ep2,p2,epsil)        
   trnslo=epsil-(cint*log(1.0+xint*up2p1))**(sexp/uexp)                      
   350 continue                                                                  
   128 continue                                                                  
   205 continue                                                                  
!
   return                                                                    
   end subroutine co2_interp
!-------------------------------------------------------------------------------
