#include <define.h>
   subroutine sph_comp_coeff
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [phys_main_solver]
!        |
!        |--- [sph_comp_coeff1] *
!        |--- [sph_comp_coeff2] *
!        |--- [sph_sum_coeff] *
!        |--- [sph_sum_coeff_top] *
!        |--- [sph_sum_wind_coeff] *
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   end subroutine sph_comp_coeff
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_comp_coeff1(fp,fm,fln,qln,llstr,llens,lwvdef,lota)
!-------------------------------------------------------------------------------
!
! abstract:
!     compute the coefficients of the expansion in spherical harmonics          
!     of the field at each level                                                
! 
!                                                                               
!     statement function                                                        
!     ------------------                                                        
!                                                                               
!     offset(n,l) is the offset in words                                        
!     to the (n,l)-element of a lower                                           
!     triangular matrix of complex numbers                                      
!     in an array containing the matrix                                         
!     packed in column-major order,                                             
!     where l and n range from 0 to jcap,                                       
!     inclusive                                                                 
!                                                                               
!          lower triangular matrix of complex numbers:                          
!                                                                               
!                     l -->                                                     
!                                                                               
!                   x                                                           
!               n   x x                                                         
!                   x x x                                                       
!               |   x x x x                                                     
!               v   x x x x x                                                   
!                   x x x x x x                                                 
!                                                                               
!          order of the matrix elements in memory:                              
!                                                                               
!          (0,0), (1,0), (2,0), ..., (jcap,0), (1,1), (2,1), (3,1), ...         
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only   :  lnt_,jcap_,jcap1_
#ifndef RMP
   use paramodel, only   :  JCAPS       , LLN2S       , LLN22S
#else
   use paramodel, only   :  JCAPS=>jcap_, LLN2S=>lnt2_, LLN22S=>lnt22_
#endif
! 
   integer,intent(in)   ::  llstr,llens,lota
   real                 ::  fp(2,0:JCAPS,lota), fm(2,0:JCAPS,lota)
   real                 ::  qln(LLN2S), fln(LLN22S,lota)
   integer              ::  lwvdef(jcap1_)
!                                
!     local scalars               
!                                                                               
   integer              ::  n, l, k   
   integer              ::  offset                                      
!
   offset(n,l) = (jcap_+1)*(jcap_+2) - (jcap_+1-l)*(jcap_+2-l) + 2*(n-l) 
!                                                                               
   lntx=0
!
   do lx = 1,llens
     l=lwvdef(llstr+lx)
     lnt0=offset(l,l)-lntx
     lntx=lntx+offset(l+1,l+1)-offset(l,l)
     ll=lx-1
!
#define DEFAULT
#ifdef SGERX1
#undef DEFAULT
     ls=l*((2*jcap_+3)-l)                                                    
     lls=ls-lnt0
#endif
!                                                                               
!        compute the even (n-l) expansion coefficients for each level           
!        ------------------------------------------------------------           
!        real part                                                              
!                                                                               
#ifdef SGERX1
     call sgerx1((jcap_+2-l)/2,lota,1.,qln(lls+1),4,                           &
                  fp(1,ll,1),(jcap_+1)*2,fln(lls+1,1),4,LLN22S)
#endif
#ifdef DEFAULT
     do n = l,jcap_,2                                                      
       do k = 1,lota                                                      
         nl=offset(n,l)-lnt0
         fln(nl+1,k) = fln(nl+1,k)+ fp(1,ll,k)*qln(nl+1)
       end do                                                              
     end do                                                                 
#endif
!                                                                               
!        imaginary part                                                         
!                                                                               
#ifdef SGERX1
     call sgerx1((jcap_+2-l)/2,lota,1.,qln(lls+2),4,                           &
                  fp(2,ll,1),(JCAPS+1)*2,fln(lls+2,1),4,LLN22S) 
#endif
#ifdef DEFAULT
     do n = l,jcap_,2                                                      
       do k = 1,lota                                                      
         nl=offset(n,l)-lnt0
         fln(nl+2,k) = fln(nl+2,k)+ fp(2,ll,k)*qln(nl+2)
       end do                                                              
     end do                                                                 
#endif
!                                                                               
!        compute the odd (n-l) expansion coefficients for each level            
!        -----------------------------------------------------------            
!
#ifdef SGERX1
     if(l.lt.jcap_) then                                                     
#endif
       !                                                                               
       ! real part                                                              
       !                                                                               
#ifdef SGERX1
       call sgerx1((jcap_+1-l)/2,lota,1.,qln(lls+3),4,                         &
                    fm(1,ll,1),(JCAPS+1)*2,fln(lls+3,1),4,LLN22S)               
#endif
#ifdef DEFAULT
!
       do n = l+1,jcap_,2                                                  
         do k = 1,lota                                                    
           nl=offset(n,l)-lnt0
           fln(nl+1,k) = fln(nl+1,k)+ fm(1,ll,k)*qln(nl+1)          
         end do                                                            
       end do                                                               
#endif
       !                                                                               
       ! imaginary part                                                         
       !                                                                               
#ifdef SGERX1
       call sgerx1((jcap_+1-l)/2,lota,1.,qln(lls+4),4,                         &
                    fm(2,ll,1),(JCAPS+1)*2,fln(lls+4,1),4,LLN22S)               
#endif
#ifdef DEFAULT
!
       do n = l+1,jcap_,2                                                  
         do k = 1,lota                                                    
           nl=offset(n,l)-lnt0
           fln(nl+2,k) = fln(nl+2,k)+ fm(2,ll,k)*qln(nl+2)
         end do                                                            
       end do                                                               
#endif
!                                                                               
#ifdef SGERX1
     endif                                                                  
#endif
!
   enddo
!
   return                                                                    
   end subroutine sph_comp_coeff1
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_comp_coeff2(fp,fm,fln,qln,llstr,llens,lwvdef,lota)
!-------------------------------------------------------------------------------
#ifndef RMP
!
! abstract:
!     compute the coefficients of the expansion in spherical harmonics          
!     of the field at each level                                                
!                                                                               
!
!     statement function                                                        
!     ------------------                                                        
!                                                                               
!     offset(n,l) is the offset in words                                        
!     to the (n,l)-element of a lower                                           
!     triangular matrix of complex numbers                                      
!     in an array containing the matrix                                         
!     packed in column-major order,                                             
!     where l and n range from 0 to jcap,                                       
!     inclusive                                                                 
!                                                                               
!          lower triangular matrix of complex numbers:                          
!                                                                               
!                     l -->                                                     
!                                                                               
!                   x                                                           
!               n   x x                                                         
!                   x x x                                                       
!               |   x x x x                                                     
!               v   x x x x x                                                   
!                   x x x x x x                                                 
!                                                                               
!          order of the matrix elements in memory:                              
!                                                                               
!          (0,0), (1,0), (2,0), ..., (jcap,0), (1,1), (2,1), (3,1), ...         
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only : JCAPS,LLN2S,LLN22S,lnt_,jcap_,jcap1_
!-------------------------------------------------------------------------------
!                                                                               
! parallelize at k level.                                                       
!                                                                               
   integer,intent(in)   ::  llstr,llens,lota
   real                 ::  fp(2,0:JCAPS,lota),fm(2,0:JCAPS,lota)
   real                 ::  qln(LLN2S), fln(LLN22S,lota)
   integer              ::  lwvdef(jcap1_)
   
!                                                                               
! local scalars                                                             
! -------------                                                             
   integer n, l, k                                                           
   integer              ::  offset                                                            
!
   offset(n,l)=(jcap_+1)*(jcap_+2)-(jcap_+1-l)*(jcap_+2-l)+2*(n-l)
!
#ifdef ORIGIN_THREAD
!$doacross share(fp,fm,fln,qln,llstr,llens,lwvdef,jcap_,lota)
!$& local(k,l,n,ll,nl,lntx,lnt0,lx)                                                                
#endif
#ifdef OPENMP
!$omp parallel do private(k,l,n,ll,nl,lntx,lnt0,lx)
#endif
!                                                                               
   do k = 1,lota  ! mj 5/8/1998                                                
!
     lntx=0
!
     do lx = 1,llens
       l=lwvdef(llstr+lx)
       lnt0=offset(l,l)-lntx
       lntx=lntx+offset(l+1,l+1)-offset(l,l)
       ll=lx-1
!                                                                               
!      compute the even (n-l) expansion coefficients for each level           
!      ------------------------------------------------------------           
!                                                                               
!      real part                                                              
!                                                                               
       do n = l,jcap_,2
         nl=offset(n,l)-lnt0
         fln(nl+1,k) = fln(nl+1,k)+ fp(1,ll,k)*qln(nl+1)
       end do                                                                 
!                                                                               
!      imaginary part                                                         
!                                                                               
       do n = l,jcap_,2 
         nl=offset(n,l)-lnt0
         fln(nl+2,k) = fln(nl+2,k)+ fp(2,ll,k)*qln(nl+2)
       end do                                                                 
!                                                                               
!      compute the odd (n-l) expansion coefficients for each level            
!      -----------------------------------------------------------            
!                                                                               
!      real part                                                              
!                                                                               
       do n = l+1,jcap_,2                                                  
         nl=offset(n,l)-lnt0
         fln(nl+1,k) = fln(nl+1,k)+ fm(1,ll,k)*qln(nl+1)
       end do                                                               
!                                                                               
!      imaginary part                                                         
!                                                                               
       do n = l+1,jcap_,2                                                  
         nl=offset(n,l)-lnt0
         fln(nl+2,k) = fln(nl+2,k)+ fm(2,ll,k)*qln(nl+2)
       end do                                                               
!                                                                               
     enddo
!
   enddo     !mj 5/8/1998                                                    
#endif
!                                                                               
   return                                                                    
   end subroutine sph_comp_coeff2
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_sum_coeff(fln,ap,qln,llstr,llens,lwvdef,levs)
!-------------------------------------------------------------------------------
   use paramodel, only   :  jcap_,jcap1_
#ifndef RMP
   use paramodel, only   :  LLN2S       ,LLN22S        ,len0=>LCAPS
#else
   use paramodel, only   :  LLN2S=>lnt2_,LLN22S=>lnt22_,len0=>lonf_
#endif
!-------------------------------------------------------------------------------
!
! if RMP, run se
!
   integer              ::  lenh
!
   integer,intent(in)   ::  llstr,llens,levs
   real                 ::  ap(2,0:len0,levs), qln(LLN2S), fln(LLN22S,levs)
   integer              ::  lwvdef(jcap1_)
!
!     local scalars  
!
   integer              ::  i, n, l, k
   real                 ::  evenr, eveni    
!                                                                               
!     statement functions                                                       
!     -------------------                                                       
!                                                                               
!     offset(n,l) is the offset in words                                        
!     to the (n,l)-element of a lower                                           
!     triangular matrix of complex numbers                                      
!     in an array containing the matrix                                         
!     packed in column-major order,                                             
!     where l and n range from 0 to jcap,                                       
!     inclusive                                                                 
!                                                                               
!          lower triangular matrix of complex numbers:                          
!                                                                               
!                     l -->                                                     
!                                                                               
!                   x                                                           
!               n   x x                                                         
!                   x x x                                                       
!               |   x x x x                                                     
!               v   x x x x x                                                   
!                   x x x x x x                                                 
!                                                                               
!          order of the matrix elements in memory:                              
!                                                                               
!          (0,0), (1,0), (2,0), ..., (jcap,0), (1,1), (2,1), (3,1), ...         
!                                                                               
   integer              ::  offset
!-------------------------------------------------------------------------------
!
   offset(n,l)=(jcap_+1)*(jcap_+2)-(jcap_-l+1)*(jcap_-l+2)+2*(n-l)
!                                                                               
!     term(1,n,l,k) and term(2,n,l,k) are                                       
!     the real and imaginary part, resp.,                                       
!     of exp((0,1)*l*phi) times the (n,l) term                                  
!     in the expansion in spherical                                             
!     harmonics of the field at level k,                                        
!     where phi is the azimuthal angle                                          
!                                                                               
   term(i,n,l,k) = qln(offset(n,l)-lnt0+i)*fln(offset(n,l)-lnt0+i,k)
!
   lenh=len0/2
!                                                                               
!     zero the accumulators                                                     
!     ---------------------                                                     
!                                                                               
!     lens=l1e-l1s+l2e-l2s+1
!
   do k = 1,levs                                                            
     do l = 0,len0
       ap(1,l,k) = 0.                                                      
       ap(2,l,k) = 0.                                                      
     end do                                                                 
   end do                                                                    
!                                                                               
!     compute the even and odd (n-l) components                                 
!     of the fourier coefficients                                               
!     ---------------------------------------------------------                 
!                                                                               
   lntx=0
   DO lx = 1,llens
     l=lwvdef(llstr+lx)
     lnt0=offset(l,l)-lntx
     lntx=lntx+offset(l+1,l+1)-offset(l,l)
     ll=lx-1
!
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
     ls=l*((2*jcap_+3)-l)                                                    
     lls=ls-lnt0
#endif
!                                                                               
!        compute the sum of the even (n-l) terms for each level                 
!        ------------------------------------------------------                 
!                                                                               
!        real part                                                              
!                                                                               
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
     call sgemvx1(levs,(jcap_+2-l)/2,1.,fln(lls+1,1),LLN22S,4,                 &
                  qln(lls+1),4,1.,ap(1,ll,1),(len0+1)*2)
#endif
#ifdef DEFAULT
!
     do n = l,jcap_,2                                                      
       do k = 1,levs
         ap(1,ll,k) = ap(1,ll,k) + term(1,n,l,k)
       end do                                                              
     end do                                                                 
#endif
!                                                                               
!        imaginary part                                                         
!                                                                               
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
     call sgemvx1(levs,(jcap_+2-l)/2,1.,fln(lls+2,1),LLN22S,4,                 &
                  qln(lls+2),4,1.,ap(2,ll,1),(len0+1)*2)
#endif
#ifdef DEFAULT
!
     do n = l,jcap_,2                                                      
       do k = 1,levs                                                      
         ap(2,ll,k) = ap(2,ll,k) + term(2,n,l,k)
       end do                                                              
     end do                                                                 
#endif
!                                                                               
!        compute the sum of the odd (n-l) terms for each level                  
!        -----------------------------------------------------                  
!                                                                               
!          real part                                                            
!                                                                               
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
     if(l.lt.jcap_) then                                                     
       call sgemvx1(levs,(jcap_+1-l)/2,1.,fln(lls+3,1),LLN22S,4,               &
                    qln(lls+3),4,1.,ap(1,lenh+ll,1),(len0+1)*2)
#endif
#ifdef DEFAULT
!
       do n = l+1,jcap_,2                                                  
         do k = 1,levs                                                    
           ap(1,lenh+ll,k) = ap(1,lenh+ll,k) + term(1,n,l,k)
         end do                                                            
       end do                                                               
#endif
!                                                                               
!          imaginary part                                                       
!                                                                               
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
       call sgemvx1(levs,(jcap_+1-l)/2,1.,fln(lls+4,1),LLN22S,4,               &
                    qln(lls+4),4,1.,ap(2,lenh+ll,1),(len0+1)*2)
#endif
#ifdef DEFAULT
!
       do n = l+1,jcap_,2                                                  
         do k = 1,levs                                                    
           ap(2,lenh+ll,k) = ap(2,lenh+ll,k) + term(2,n,l,k)
         end do                                                            
       end do                                                               
#endif
!                                                                               
#define DEFAULT
#ifdef SGEMVX1
#undef DEFAULT
     endif                                                                  
#endif
!
!     compute the fourier coefficients for each level                           
!     -----------------------------------------------                           
!                                                                               
     do k = 1,levs                                                            
       evenr = ap(1,ll,k)                                                   
       eveni = ap(2,ll,k)                                                   
       ap(1,ll,k) = ap(1,ll,k) + ap(1,lenh+ll,k)                              
       ap(2,ll,k) = ap(2,ll,k) + ap(2,lenh+ll,k)                              
       ap(1,lenh+ll,k) = evenr - ap(1,lenh+ll,k)                             
       ap(2,lenh+ll,k) = eveni - ap(2,lenh+ll,k)                             
     end do                                                                    
!
   ENDDO
#undef DEFAULT
!                                                                               
   return                                                                    
   end subroutine sph_sum_coeff
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_sum_coeff_top(ap,top,qvv,llstr,llens,lwvdef,kall) 
!-------------------------------------------------------------------------------
   use paramodel, only   : jcap1_
#ifndef RMP
   use paramodel, only   : JCAPS       , LEVSS       , len0=>LCAPS
#else
   use paramodel, only   : JCAPS=>jcap_, LEVSS=>levs_, len0=>lonf_
#endif
!-------------------------------------------------------------------------------
   integer              ::  lenh
! 
   real                 ::  ap(2,0:len0,kall)            
   real                 ::  top(2,0:JCAPS,kall)                           
   real                 ::  qvv(2,0:JCAPS)                                 
   integer              ::  lwvdef(jcap1_)                             
!
! local array
!
   real                 ::  qtop(2,0:JCAPS)
!-------------------------------------------------------------------------------
   lenh=len0/2
!
   do lx = 1,llens
     l=lwvdef(llstr+lx)
     ll=lx-1
     do k = 1,kall                                                           
       qtop(1,ll) = top(1,ll,k) * qvv(1,ll)
       qtop(2,ll) = top(2,ll,k) * qvv(2,ll)
       ap(1,ll,k) = ap(1,ll,k) + qtop(1,ll)
       ap(2,ll,k) = ap(2,ll,k) + qtop(2,ll)
       if( mod(l,2).eq.0 ) then
         !
         ! odd wavenumber  for l=0,2,4,.....
         !
         ap(1,ll+lenh,k) = ap(1,ll+lenh,k) - qtop(1,ll)
         ap(2,ll+lenh,k) = ap(2,ll+lenh,k) - qtop(2,ll)
       else
         !
         ! even wavenumber for l=1,3,5,......
         !
         ap(1,ll+lenh,k) = ap(1,ll+lenh,k) + qtop(1,ll)
         ap(2,ll+lenh,k) = ap(2,ll+lenh,k) + qtop(2,ll)
       endif
     enddo
   enddo
!
   return                                                                    
   end subroutine sph_sum_coeff_top
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_sum_wind_coeff(fpu,fmu,fpv,fmv,topuln,topvln,qvv,wgt,        &
                     llstr,llens,lwvdef,levs)             
!-------------------------------------------------------------------------------
!                                                                               
!            this sr assumes jcap is even                                       
!            this sr assumes jcap is even                                       
!            this sr assumes jcap is even                                       
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only   : JCAPS,jcap1_
!-------------------------------------------------------------------------------
   integer              ::  llstr,llens,levs
   real                 ::  fpu(2,0:JCAPS,levs)
   real                 ::  fpv(2,0:JCAPS,levs)
   real                 ::  fmu(2,0:JCAPS,levs)
   real                 ::  fmv(2,0:JCAPS,levs)
   real                 ::  topuln(2,0:JCAPS,levs)
   real                 ::  topvln(2,0:JCAPS,levs)
   real                 ::  toppln(0:JCAPS)  
   real                 ::  qvv(2,0:JCAPS)  
   integer              ::  lwvdef(jcap1_)
!-------------------------------------------------------------------------------
!                                                                               
!     ----------------------------------------------------------------          
!     compute expansion coeffs. for top rows of u and v                         
!     ----------------------------------------------------------------          
!
   do l = 0,JCAPS                                                           
     toppln(l) = qvv(1,l)*wgt  
   enddo
!    
   do lx = 1,llens
     l=lwvdef(llstr+lx)
     ll=lx-1
!
     if( mod(l,2).eq.1 ) then
!                                                                               
!  even wavenumber for l=1,3,5,7,.........
!
!         compute the even (n-l) expansion coefficients for each level    
!
!         ------------------------------------------------------------    
!                                                                               
!         real part                                                     
!                                                                               
       do k = 1,levs
         topuln(1,ll,k) = topuln(1,ll,k)+fpu(1,ll,k)*toppln(ll)
         topvln(1,ll,k) = topvln(1,ll,k)+fpv(1,ll,k)*toppln(ll)
       enddo
!                                                                               
!         imaginary part                            
!                                                                               
       do k = 1,levs
         topuln(2,ll,k) = topuln(2,ll,k)+fpu(2,ll,k)*toppln(ll)
         topvln(2,ll,k) = topvln(2,ll,k)+fpv(2,ll,k)*toppln(ll)
       enddo
!                                                                               
!  odd wavenumber for l=0,2,4,6,.........
!
     else
!
!         compute the odd (n-l) expansion coefficients for each level      
!         -----------------------------------------------------------   
!
!         real part                                          
!
       do k = 1,levs
            topuln(1,ll,k) = topuln(1,ll,k)+fmu(1,ll,k)*toppln(ll)
            topvln(1,ll,k) = topvln(1,ll,k)+fmv(1,ll,k)*toppln(ll)
       enddo
!                                                                               
!         imaginary part                                                   
!                                                                               
       do k = 1,levs
         topuln(2,ll,k) = topuln(2,ll,k)+fmu(2,ll,k)*toppln(ll)
         topvln(2,ll,k) = topvln(2,ll,k)+fmv(2,ll,k)*toppln(ll)
       enddo
!
!
     endif
!
   enddo                                                                    
!
   return                                                                    
   end subroutine sph_sum_wind_coeff
!-------------------------------------------------------------------------------
