#include <define.h>
   subroutine chgr_sph_gaussian(lgghaf,colrad,wgt,wgtcs,rcs2)     
!-------------------------------------------------------------------------------
   use constant, only : pi_
   use module_sph_legendre, only : sph_legendre_polinomial 
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   save                                     
!
   integer         iter,k,k1,l2,lgghaf,iprint
   real*8 drad,dradz,eps,p1,p2,phi,pi,rad,rc
   real*8 rl2,scale,si,sn,w,x
   real*8 colrad(lgghaf)
   real*8 wgt(lgghaf)
   real*8 wgtcs(lgghaf)
   real*8 rcs2(lgghaf)
   real*8 cons0,cons0p25,cons1          !constant
   real*8 cons2,cons4,cons180,cons360   !constant
   cons0    =   0.d0    !constant
   cons0p25 =   0.25d0  !constant
   cons1    =   1.d0    !constant
   cons2    =   2.d0    !constant
   cons4    =   4.d0    !constant
   cons180  = 180.d0    !constant
   cons360  = 360.d0    !constant
!
   eps=1.d-12           !constant
!-------------------------------------------------------------------------------
!                                        
   si = cons1                  !constant
   l2=2*lgghaf
   rl2=l2
   scale = cons2/(rl2*rl2)     !constant
   k1=l2-1
   pi = atan(si)*cons4         !constant
   dradz = pi / cons360 / 10.0 !constant
   rad = cons0                 !constant
!
   do k = 1,lgghaf
     iter=0
     drad=dradz
1    call sph_legendre_polinomial(l2,rad,p2)
2    p1 =p2
     iter=iter+1
     rad=rad+drad
     call sph_legendre_polinomial(l2,rad,p2)
     if(sign(si,p1).eq.sign(si,p2)) go to 2
     if(drad.lt.eps)go to 3
     rad=rad-drad
     drad = drad * cons0p25           !constant
     go to 1
3    continue
     colrad(k)=rad
     phi = rad * cons180 / pi         !constant
     call sph_legendre_polinomial(k1,rad,p1)
     x = cos(rad)
     w = scale * (cons1 - x*x)/ (p1*p1)     !constant
     wgt(k) = w
     sn = sin(rad)
     w=w/(sn*sn)
     wgtcs(k) = w
     rc=cons1/(sn*sn)              !constant
     rcs2(k) = rc
     call sph_legendre_polinomial(l2,rad,p1)
!                                        
     print 102,k,phi,colrad(k),wgt(k),wgtcs(k),iter,p1             
102  format(1h ,i2,2x,f6.2,2x,f10.7,2x,e13.7,2x,e13.7,2x,i4,2x,d13.7)       
!                                        
   enddo
!                                        
   return                                   
   end subroutine chgr_sph_gaussian 
!-------------------------------------------------------------------------------
