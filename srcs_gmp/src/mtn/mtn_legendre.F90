!
   subroutine mtn_legendre(ipd,cosclt,maxwv,iromb,ex,px,pnm)                     
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:  mtn_legendre   compute legendre polynomials.                         
!                                                                               
! abstract: evaluates the normalized associated legendre polynomials            
!   at a given cosine of colatitude.                                            
!                                                                               
! program history log:                                                          
!   1991-01-01  iredell                initial mrf
!   2000-01-01  hann-ming henry juang  mpi
!   2000-01-01  song-you hong          usgs topo, kagwd options
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call mtn_legendre(ipd,cosclt,maxwv,iromb,ex,px,pnm)                     
!                                                                               
!   input argument list:                                                        
!     ipd      - ipd = 0 to evaluate polynomial itself,                         
!                ipd = 1 to evaluate polynomial derivative.                     
!     cosclt   - cosine of colatitude at which to evalutate.                    
!     maxwv    - spectral truncation                                            
!     iromb    - iromb = 0 for triangular truncation                            
!                iromb = 1 for rhomboidal truncation                            
!                                                                               
!   work argument list:                                                         
!     ex       - real (maxwv+2) work area                                       
!     px       - real (maxwv+3) work area                                       
!                                                                               
!   output argument list:                                                       
!     pnm      - real (kmax) polynomial values,                                 
!                where kmax=(maxwv+1)*(iromb+1)*maxwv+2)/2.                     
!
!-------------------------------------------------------------------------------
   save                                                                      
   real                 ::  pnm((maxwv+1)*((iromb+1)*maxwv+2)/2)                            
   real                 ::  ex(0:maxwv+1),px(-1:maxwv+1)                                    
!
   sinclt=sqrt(1.-cosclt**2)                                                 
   k=0                                                                       
   ex(0)=0.                                                                  
   px(-1)=0.                                                                 
   px(0)=sqrt(0.5)                                                           
!
   do m = 0,maxwv                                                              
     if(m.gt.0) px(0)=px(0)*sinclt/(ex(1)*sqrt(float(2*m)))                  
     do n = m+1,maxwv+iromb*m+1                                                
       ex(n-m)=sqrt(float(n**2-m**2)/float(4*n**2-1))                        
     enddo                                                                   
     do n = m+1,maxwv+iromb*m+1                                                
       px(n-m)=(cosclt*px(n-m-1)-ex(n-m-1)*px(n-m-2))/ex(n-m)                
     enddo                                                                   
     if(ipd.eq.0) then                                                       
       do n = m,maxwv+iromb*m                                                  
         k=k+1                                                               
         pnm(k)=px(n-m)                                                      
       enddo                                                                 
     else                                                                    
       do n = m,maxwv+iromb*m                                                  
         k=k+1                                                               
         pnm(k)=n*ex(n-m+1)*px(n-m+1)-(n+1)*ex(n-m)*px(n-m-1)                
       enddo                                                                 
     endif                                                                   
   enddo                                                                     
!
   return                                                                    
   end subroutine mtn_legendre
!-------------------------------------------------------------------------------
