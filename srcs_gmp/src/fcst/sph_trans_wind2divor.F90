#include <define.h>
   subroutine sph_trans_wind2divor
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!    [sph_trans_wind2divor]
!      |
!      |--- [sph_wind2divor] *
!      |--- [sph_divor2wind] *
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   end subroutine sph_trans_wind2divor
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_wind2divor(uln,vln,di,ze,topuln,topvln,llstr,llens,lwvdef)
!-------------------------------------------------------------------------------
#ifndef RMP
   use paramodel, only   :  JCAPS,LLN2S,LLN22S,LEVSS,jcap_,jcap1_,lnt22_
   use constant, only    :  rerth_
!-------------------------------------------------------------------------------
   integer,intent(in )  ::  llstr,llens
   real,intent(out)     ::  di(LLN22S,LEVSS)
   real,intent(out)     ::  ze(LLN22S,LEVSS)
   real,intent(in )     ::  uln(LLN22S,LEVSS)
   real,intent(in )     ::  vln(LLN22S,LEVSS)
   real,intent(in )     ::  topuln(2,0:JCAPS,LEVSS)
   real,intent(in )     ::  topvln(2,0:JCAPS,LEVSS)
   integer,intent(in )  ::  lwvdef(jcap1_)        
!
!  local
!
   real,save,allocatable  ::  topeps(:)        
   real,save,allocatable  ::  eps(:)      
   integer,save           ::  ifirst
!                                                                               
!  local scalars                                                             
!                                                                               
   integer i, n, l, k                                                        
!-------------------------------------------------------------------------------
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
!               ?   x x x x                                                     
!               v   x x x x x                                                   
!                   x x x x x x                                                 
!                                                                               
!          order of the matrix elements in memory:                              
!                                                                               
!          (0,0), (1,0), (2,0), ..., (jcap,0), (1,1), (2,1), (3,1), ...         
!                                                                               
!-------------------------------------------------------------------------------
   integer              ::  offset    
!-------------------------------------------------------------------------------
   offset(n,l) = (jcap_+1)*(jcap_+2) - (jcap_-l+1)*(jcap_-l+2) + 2*(n-l) 
!
   if(.not.allocated(topeps)) allocate(topeps(0:jcap_))
   if(.not.allocated(eps))    allocate(eps(lnt22_))
!                                                                               
!     ---                                                                       
!                                                                               
!     term(1,n,l,k) and term(2,n,l,k) are                                       
!     the real and imaginary part, resp.,                                       
!     of exp((0,1)*l*phi) times the (n,l) term                                  
!     in the expansion in spherical                                             
!     harmonics of the field at level k,                                        
!     where phi is the azimuthal angle                                          
!                                                                               
!     term(i,n,l,k) = di(offset(n,l)+i,k)                                       
!                                                                               
!                                                                               
   data ifirst/0/                                                            
   if(ifirst.eq.0) then
     do l = 0,jcap_  
       do n = l,jcap_
         temp=((n*n-l*l)/(4.*n*n-1.))                                        
         if(n.eq.0) temp=0.                                                  
         temp=sqrt(temp)                                                     
         eps(offset(n,l)+1)=temp                                             
         eps(offset(n,l)+2)=temp                                             
       end do                                                                 
       n=jcap_+1                                                            
       temp=((n*n-l*l)/(4.*n*n-1.))                                        
       topeps(l)=sqrt(temp)                                                
     end do                                                                    
     ifirst=1                                                                  
   endif
!
! --------------------------------------------------------------
!
   do k = 1,LEVSS                                        
     lntx=0
     do lx = 1,llens
       l=lwvdef(llstr+lx)
       lnt0=offset(l,l)-lntx
       lntx=lntx+offset(l+1,l+1)-offset(l,l)
       ll=lx-1
       !
       !        the case n=l      
       !-----do l = 1,jcap_-1   
       ! 
       if( l.ge.1 .and. l.le.jcap_-1 ) then
         rl= l
         n = l                                                               
         rn= n                 
         nl=offset(n,l)-lnt0
         npl=offset(n+1,l)-lnt0
         ze(nl+1,k)=-rl*vln(nl+2,k)-rn*eps(offset(n+1,l)+1)*uln(npl+1,k) 
         ze(nl+2,k)= rl*vln(nl+1,k)-rn*eps(offset(n+1,l)+2)*uln(npl+2,k)
         di(nl+1,k)=-rl*uln(nl+2,k)+rn*eps(offset(n+1,l)+1)*vln(npl+1,k)
         di(nl+2,k)= rl*uln(nl+1,k)+rn*eps(offset(n+1,l)+2)*vln(npl+2,k)
       endif
       ! 
       !---  do l=0
       ! 
       if( l.eq.0 ) then
         ze(1,k)=0.                                                    
         ze(2,k)=0.                                                   
         di(1,k)=0.                                                  
         di(2,k)=0.                                                 
       endif
       !
       !---  do l = 0,jcap_ 
       !  
       rl=l
       do n = l+1,jcap_-1    
         rn=n               
         nl=offset(n,l)-lnt0
         npl=offset(n+1,l)-lnt0
         nml=offset(n-1,l)-lnt0
         ze(nl+1,k)=-rl*vln(nl+2,k)                                            &
            -rn*eps(offset(n+1,l)+1)*uln(npl+1,k)                              &
            +(rn+1.)*eps(offset(n  ,l)+1)*uln(nml+1,k)                 
         ze(nl+2,k)= rl*vln(nl+1,k)                                            &
            -rn*eps(offset(n+1,l)+2)*uln(npl+2,k)                              &
            +(rn+1.)*eps(offset(n  ,l)+2)*uln(nml+2,k)                 
         di(nl+1,k)=-rl*uln(nl+2,k)                                            &
            +rn*eps(offset(n+1,l)+1)*vln(npl+1,k)                              &
            -(rn+1.)*eps(offset(n  ,l)+1)*vln(nml+1,k)                 
         di(nl+2,k)= rl*uln(nl+1,k)                                            &
            +rn*eps(offset(n+1,l)+2)*vln(npl+2,k)                              &
            -(rn+1.)*eps(offset(n  ,l)+2)*vln(nml+2,k)                 
       end do                                                                 
       !   
       ! do top row involving u,v at n=jcap_+1
       !
       n =  jcap_                                                             
       rn=n                                                                  
       !
       !-----do l = 0, jcap_
       !
       rl=l                
       nl=offset(n,l)-lnt0
       nml=offset(n-1,l)-lnt0
       !     
       ze(nl+1,k)=-rl*vln(nl+2,k)+(rn+1.)*eps(offset(n,l)+1)*uln(nml+1,k)  
       ze(nl+2,k)= rl*vln(nl+1,k)+(rn+1.)*eps(offset(n,l)+2)*uln(nml+2,k) 
       !
       di(nl+1,k)=-rl*uln(nl+2,k)-(rn+1.)*eps(offset(n,l)+1)*vln(nml+1,k)
       di(nl+2,k)= rl*uln(nl+1,k)-(rn+1.)*eps(offset(n,l)+2)*vln(nml+2,k) 
       !                                                                 
       !-----do l = 0, jcap_                                            
       !  
       !         nl=offset(n,l)-lnt0
       !
       ze(nl+1,k)=ze(nl+1,k)-rn*topeps(l)*topuln(1,ll,k)
       ze(nl+2,k)=ze(nl+2,k)-rn*topeps(l)*topuln(2,ll,k)
       !                                                               
       di(nl+1,k)=di(nl+1,k)+rn*topeps(l)*topvln(1,ll,k)
       di(nl+2,k)=di(nl+2,k)+rn*topeps(l)*topvln(2,ll,k)
     enddo                                                                    
   enddo
!
   do k = 1,LEVSS
     do j = 1,LLN2S
       di(j,k)=di(j,k)/rerth_                                            
       ze(j,k)=ze(j,k)/rerth_                                            
     end do                                                             
   enddo
! 
#endif 
   return                                                                    
   end subroutine sph_wind2divor                                                                      
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_divor2wind(d,z,u,v,utop,vtop,llstr,llens,lwvdef)
!-------------------------------------------------------------------------------
   use constant, only  : rerth_
   use paramodel, only : jcap_,jcap1_
#ifndef RMP
   use paramodel, only : JCAP1S,LEVSS,LLN22S
#else
   use paramodel, only : JCAP1S=>jcap1_,LEVSS=>levs_,LLN22S=>lnt22_
#endif
!-------------------------------------------------------------------------------
   integer              ::  lnep
   integer,intent(in)   ::  llstr,llens
   real,intent(in)      ::  d(LLN22S,LEVSS),z(LLN22S,LEVSS)
   real,intent(out)     ::  u(LLN22S,LEVSS),v(LLN22S,LEVSS)
   real,intent(out)     ::  utop(2,JCAP1S,LEVSS)
   real,intent(out)     ::  vtop(2,JCAP1S,LEVSS)
   integer,intent(in)   ::  lwvdef(jcap1_)
!
   real,save,allocatable::  e(:)
   integer,save         ::  ifirst
!-------------------------------------------------------------------------------
!                     
!   array e =eps/n   
!   array e =eps/n  
!            eps/n=0. for n=l        
!   array e =eps/n                  
!   array e =eps/n                 
!                                 
   je(n,l) =((jcap_+2)*(jcap_+3)-(jcap_+2-l)*(jcap_+3-l))/2+n-l       
!                                                                 
   jc(n,l) = (jcap_+1)*(jcap_+2)-(jcap_+1-l)*(jcap_+2-l)+2*(n-l)    
!                                                               
   data ifirst/1/                                                            
!
! ------------------ initial setting --------------------      
!
   lnep=(jcap_+2)*(jcap_+3)/2
   if(.not.allocated(e)) allocate(e(lnep))
!
   if(ifirst.eq.1) then
     do l = 0,jcap_                                                           
       n=l                                                                
       e(je(n,l)+1)=0.
     enddo
     do l= 0,jcap_                                                         
       do n = l+1,jcap_+1                                        
         rn=n                                                
         rl=l                                               
         a=(rn*rn-rl*rl)/(4.*rn*rn-1.)                     
         e(je(n,l)+1)=sqrt(a) / rn                        
       enddo
     enddo
     ifirst=0                                                                  
   endif
!
! kang initialize u and v 
!
   do k = 1,LEVSS
     do l = 1,LLN22S
       u(l,k) = 0.
       v(l,k) = 0.
     enddo
   enddo
!
! -------------------- end of initial setting ---------------
!
   lntx=0
   do ll = 1,llens
     l=lwvdef(llstr+ll)
     lnt0=jc(l,l)-lntx
     lntx=lntx+jc(l+1,l+1)-jc(l,l)
!
#ifdef ORIGIN_THREAD
!$doacross share(d,z,u,v,utop,vtop,e,lnt0,l,ll),
!$&        local(k,n,rl,rn,j,jc0,je0)
#endif
#ifdef OPENMP
!$omp parallel do private(k,n,rl,rn,j,jc0,je0)
#endif
     do k = 1,LEVSS                                                        
       !
       if(l.eq.0) then
         do n = 0,jcap_
           jc0=jc(n,l)
           u(jc0+1,k)=0.0               
           u(jc0+2,k)=0.0              
           v(jc0+1,k)=0.0 
           v(jc0+2,k)=0.0
         enddo
       endif
       !                                
       !----- l=1,jcap_
       ! 
       if( l.ge.1 ) then
         rl=l                         
         do n = l,jcap_                
           rn=n                    
           jc0=jc(n,l)-lnt0
           u(jc0+2,k)=-rl*d(jc0+1,k)/(rn*(rn+1.))      
           u(jc0+1,k)= rl*d(jc0+2,k)/(rn*(rn+1.))     
           v(jc0+2,k)=-rl*z(jc0+1,k)/(rn*(rn+1.))    
           v(jc0+1,k)= rl*z(jc0+2,k)/(rn*(rn+1.))   
         enddo
       endif
       !                                             
       !----- l=  0,jcap_-1
       ! 
       if( l.le.jcap_-1 ) then
         !
         do n = l+1,jcap_                            
           jc0=jc(n,l)-lnt0
           je0=je(n,l)
           !
           u(jc0+1,k)=u(jc0+1,k)-e(je0+1)*z(jc0-1,k)                 
           u(jc0+2,k)=u(jc0+2,k)-e(je0+1)*z(jc0  ,k)                 
           !                                                                               
           v(jc0+1,k)=v(jc0+1,k)+e(je0+1)*d(jc0-1,k)                 
           v(jc0+2,k)=v(jc0+2,k)+e(je0+1)*d(jc0  ,k)                 
         enddo
         !
         do n = l,jcap_-1    
           jc0=jc(n,l)-lnt0
           je0=je(n+1,l)
           !
           u(jc0+1,k)=u(jc0+1,k)+e(je0+1)*z(jc0+3,k)               
           u(jc0+2,k)=u(jc0+2,k)+e(je0+1)*z(jc0+4,k)               
           !                                                                               
           v(jc0+1,k)=v(jc0+1,k)-e(je0+1)*d(jc0+3,k)               
           v(jc0+2,k)=v(jc0+2,k)-e(je0+1)*d(jc0+4,k)               
         enddo
         ! 
       endif
       !
       n=jcap_+1            
       !
       !----- l=0,jcap_
       !     
       jc0=jc(n,l)-lnt0
       je0=je(n,l)
       !
       utop(1,ll,k)=-e(je0+1)*z(jc0-1,k)                                
       utop(2,ll,k)=-e(je0+1)*z(jc0  ,k)                                
       !                   
       vtop(1,ll,k)= e(je0+1)*d(jc0-1,k)                                
       vtop(2,ll,k)= e(je0+1)*d(jc0  ,k)                                
       !                  
       !----- l=0,jcap_
       !  
       utop(1,ll,k)=utop(1,ll,k)*rerth_      
       utop(2,ll,k)=utop(2,ll,k)*rerth_     
       vtop(1,ll,k)=vtop(1,ll,k)*rerth_    
       vtop(2,ll,k)=vtop(2,ll,k)*rerth_   
       !
     enddo
     !
   enddo
   !
   do k = 1,LEVSS
     do j = 1,LLN22S
       u(j,k)=u(j,k)*rerth_ 
       v(j,k)=v(j,k)*rerth_
     enddo
   enddo
!                                                                               
   return                                                                    
   end subroutine sph_divor2wind
