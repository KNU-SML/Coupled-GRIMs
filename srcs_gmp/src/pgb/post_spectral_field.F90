!
   subroutine post_spectral_field(m,eps,epstop,enn1,elonn1,eon,eontop)   
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:    post_spectral_field        compute utility spectral fields                    
!                                                                               
! abstract: computes constant fields indexed in the spectral triangle           
!           in "ibm order" (zonal wavenumber is the slower index).              
!           if l is the zonal wavenumber and n is the total wavenumber          
!           and a is the earth radius, then the fields returned are:            
!           (1) normalizing factor epsilon=sqrt((n**2-l**2)/(4*n**2-1))         
!           (2) laplacian factor n*(n+1)/a**2                                   
!           (3) zonal derivative/laplacian factor l/(n*(n+1))*a                 
!           (4) meridional derivative/laplacian factor epsilon/n*a              
!                                                                               
! program history log:                                                          
!    1992-10-31  iredell                development
!    2000-03-09  songyou hong           cvs verion setup
!    2009-10-01  jung-eun kim           f90 format with standard physics modules
!    2010-07-01  myung-seo koo          dimension allocatable with namelist input


!                                                                               
! usage:    call post_spectral_field(m,eps,epstop,enn1,elonn1,eon,eontop)  
!                                                                               
!   input argument list:                                                        
!     m        - integer spectral truncation                                    
!                                                                               
!   output argument list:                                                       
!     eps      - real ((m+1)*(m+2)/2) sqrt((n**2-l**2)/(4*n**2-1))              
!     epstop   - real (m+1) sqrt((n**2-l**2)/(4*n**2-1)) over top               
!     enn1     - real ((m+1)*(m+2)/2) n*(n+1)/a**2                              
!     elonn1   - real ((m+1)*(m+2)/2) l/(n*(n+1))*a                             
!     eon      - real ((m+1)*(m+2)/2) epsilon/n*a                               
!     eontop   - real (m+1) epsilon/n*a over top                                
!                                                                               
!-------------------------------------------------------------------------------
   real             ::  eps((m+1)*(m+2)/2),epstop(m+1)                                       
   real             ::  enn1((m+1)*(m+2)/2),elonn1((m+1)*(m+2)/2)                            
   real             ::  eon((m+1)*(m+2)/2),eontop(m+1)                                       
   real, parameter  ::  rerth=6.3712e6,ra2=1./rerth**2                                 
!-------------------------------------------------------------------------------
!
   do l = 0,m                                                                  
     ill=l*(2*m+3-l)/2+1                                                     
     eps(ill)=0.                                                             
     enn1(ill)=ra2*l*(l+1)                                                   
     elonn1(ill)=rerth/(l+1)                                                 
     eon(ill)=0.                                                             
   enddo  
!                                                                   
   do l = 0,m                                                                  
     is=l*(2*m+1-l)                                                          
     ip=is/2+1                                                               
     do n = l+1,m                                                              
       eps(ip+n)=sqrt(float(n**2-l**2)/float(4*n**2-1))                      
       enn1(ip+n)=ra2*n*(n+1)                                                
       elonn1(ip+n)=rerth*l/(n*(n+1))                                        
       eon(ip+n)=rerth/n*eps(ip+n)                                           
     enddo                                                                   
   enddo    
!                                                                 
   do l = 0,m                                                                  
     epstop(l+1)=sqrt(float((m+1)**2-l**2)/float(4*(m+1)**2-1))              
     eontop(l+1)=rerth/(m+1)*epstop(l+1)                                     
   enddo                                                                     
! 
   return                                                                    
   end                                                                       
