!
   subroutine post_grad_scalar(m,enn1,elonn1,eon,eontop,q,qdx,qdy,qdytop) 
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:    post_grad_scalar       compute gradient in spectral space                 
!                                                                               
! abstract: computes the horizontal vector gradient of a scalar field           
!           in spectral space. 
!           subprogram post_spectral_field should be called already.        
!           if l is the zonal wavenumber, n is the total wavenumber,            
!           eps(l,n)=sqrt((n**2-l**2)/(4*n**2-1)) and a is earth radius,        
!           then the zonal gradient of q(l,n) is simply i*l/a*q(l,n)            
!           while the meridional gradient of q(l,n) is computed as              
!           eps(l,n+1)*(n+2)/a*q(l,n+1)-eps(l,n+1)*(n-1)/a*q(l,n-1).            
!           extra terms are computed over top of the spectral triangle.         
!           advantage is taken of the fact that eps(l,l)=0                      
!           in order to vectorize over the entire spectral triangle.            
!                                                                               
! program history log:                                                          
!    1992-10-31  iredell                development
!    2000-03-09  songyou hong           cvs verion setup, prognostic clouds
!    2009-10-01  jung-eun kim           f90 format with standard physics modules
!    2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call post_grad_scalar(m,enn1,elonn1,eon,eontop,q,          
!    &                 qdx,qdy,qdytop)                                          
!                                                                               
!   input argument list:                                                        
!     m        - integer spectral truncation                                    
!     enn1     - real ((m+1)*(m+2)/2) n*(n+1)/a**2                              
!     elonn1   - real ((m+1)*(m+2)/2) l/(n*(n+1))*a                             
!     eon      - real ((m+1)*(m+2)/2) epsilon/n*a                               
!     eontop   - real (m+1) epsilon/n*a over top                                
!     q        - real ((m+1)*(m+2)) scalar field                                
!                                                                               
!   output argument list:                                                       
!     qdx      - real ((m+1)*(m+2)) zonal gradient (times coslat)               
!     qdy      - real ((m+1)*(m+2)) merid gradient (times coslat)               
!     qdytop   - real (2*(m+1)) merid gradient (times coslat) over top          
!                                                                               
!-------------------------------------------------------------------------------
   real  ::  enn1((m+1)*(m+2)/2),elonn1((m+1)*(m+2)/2)                            
   real  ::  eon((m+1)*(m+2)/2),eontop(m+1)                                       
   real  ::  q((m+1)*(m+2))                                                       
   real  ::  qdx((m+1)*(m+2)),qdy((m+1)*(m+2)),qdytop(2*(m+1))                    
!-------------------------------------------------------------------------------
!
!  take zonal and meridional gradients                                          
!
   i=1                                                                       
   qdx(2*i-1)=0.                                                             
   qdx(2*i)=0.                                                               
   qdy(2*i-1)=eon(i+1)*enn1(i+1)*q(2*i+1)                                    
   qdy(2*i)=eon(i+1)*enn1(i+1)*q(2*i+2)                                      
!
   do i = 2,(m+1)*(m+2)/2-1                                                    
     qdx(2*i-1)=-elonn1(i)*enn1(i)*q(2*i)                                    
     qdx(2*i)=elonn1(i)*enn1(i)*q(2*i-1)                                     
     qdy(2*i-1)=eon(i+1)*enn1(i+1)*q(2*i+1)-eon(i)*enn1(i-1)*q(2*i-3)        
     qdy(2*i)=eon(i+1)*enn1(i+1)*q(2*i+2)-eon(i)*enn1(i-1)*q(2*i-2)          
   enddo                                                                     
!
   i=(m+1)*(m+2)/2                                                           
   qdx(2*i-1)=-elonn1(i)*enn1(i)*q(2*i)                                      
   qdx(2*i)=elonn1(i)*enn1(i)*q(2*i-1)                                       
   qdy(2*i-1)=-eon(i)*enn1(i-1)*q(2*i-3)                                     
   qdy(2*i)=-eon(i)*enn1(i-1)*q(2*i-2)                                       
! 
!  take meridional gradient over top                                            
!
   do l = 0,m                                                                  
     i=l*(2*m+1-l)/2+m+1                                                     
     qdytop(2*l+1)=-eontop(l+1)*enn1(i)*q(2*i-1)                             
     qdytop(2*l+2)=-eontop(l+1)*enn1(i)*q(2*i)                               
   enddo                                                                     
!
   return                                                                    
   end                                                                       
