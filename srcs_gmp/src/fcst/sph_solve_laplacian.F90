#include <define.h>
   subroutine sph_solve_laplacian(q,qlap,llstr,llens,lwvdef)
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:  sph_solve_laplacian compute modified laplacian in spectral space
!                                                                               
! abstract: computes the modified horizontal laplacian in spectral space        
!           for horizontal diffusion calculation of temp and hum on             
!           grid space.                                                         
!           if n is the total wavenumber, modified laplacian is expressed by:   
!             -(n-n0)*(n-n0+1)/a**2*spc(n,m), n-n0>0                            
!           where n0=0.55*jcap                                                  
!                                                                               
! program history log:
!   1997-03-10  masao kanamitsu        development                                           
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call sph_solve_laplacian(q,rq,llstr,llens,lwvdef)  
!                                                                               
!   input argument list:                                                        
!     q        - real (lnt22) ln(ps) spectral coeffs                           
!                                                                               
!   output argument list:                                                       
!     qlap     - real (lnt22) qusi-laplacian of ln(ps) spectral coeffs         
! 
!-------------------------------------------------------------------------------
   use constant,  only : rerth_
   use paramodel, only : LLN22S,jcap_,jcap1_
!-------------------------------------------------------------------------------
   real                 ::  q(LLN22S)
   real                 ::  qlap(LLN22S)
   integer              ::  lwvdef(jcap1_)
!                                                                               
   integer,parameter    ::  lefres=80,jdel=2
   real,parameter       ::  difcof=3.e15
!-------------------------------------------------------------------------------
   np=jcap_                                                                  
   n0=0.55*jcap_                                                             
   rtnp=difcof/(rerth_**4)*float(lefres*(lefres+1))**2                       
   jdelh=jdel/2                                                              
   npd=max(np-n0,0)                                                          
   dn1=2.*rtnp/float(npd*(npd+1))**jdelh                                     
!                    
   i=0                                                                       
   do ll = 1,llens
      nm=lwvdef(llstr+ll)
      do mm = 0,jcap_-nm                                                        
         nd=max(nm+mm-n0,0)                                                    
         dn=dn1*float(nd*(nd+1))**jdelh                                        
         qlap(i+1)=q(i+1)*dn                                                   
         qlap(i+2)=q(i+2)*dn                                                   
         i=i+2                                                                 
      enddo                                                                   
   enddo                                                                     
!                                                                               
   return                                                                    
   end subroutine sph_solve_laplacian
