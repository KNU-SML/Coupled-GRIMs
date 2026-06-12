#include <define.h>
   subroutine sph_derivative_init(epsi)
!-------------------------------------------------------------------------------
!   
! subroutine:    sph_derivative_init  
!
! abstract:
! - sets common for subroutine sph_derivative.
! - initializes the constant variables and arrays
!   of a common for subroutine sph_derivative.
!
! program history log:
!   1991-03-14  joseph sela
!
! usage:    call sph_derivative_init ( eps )
!
! remarks: call subroutine once before calls to sph_derivative.
!          refer to sph_derivative for additional documentation.
!
!  ::: structure :::
!
!    [sph_derivative_init] *
!        |
!        |--- [sph_matrix_trans] *
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only : jcap1_,jcap2_,jcap_,lnt2_,twoj1_
   use comfcst, only : dxa, dxb
!-------------------------------------------------------------------------------
   real  ::  epsi(jcap2_,jcap1_)
!
   real  ::  dxint(lnt2_)
   real  ::  dx(twoj1_,jcap2_)
   real  ::  deps(twoj1_,jcap2_)
!-------------------------------------------------------------------------------
!
   do ll = 1,twoj1_
     dxint(2*ll-1) = ll
     dxint(2*ll  ) = ll
   enddo
!
   lp = 0
   do i = 1,jcap2_
     do ll = 1,twoj1_
       dx(ll,i) = dxint(ll+lp)
     enddo
     lp = lp + 2
   enddo
!
   do i = 1,jcap2_
     do ll = 1,jcap1_
       deps(2*ll-1,i) = epsi(i,ll)
       deps(2*ll  ,i) = epsi(i,ll)
     enddo
   enddo
!
   do ll = 1,twoj1_
      dxa(ll) = 0.0
   enddo
!
   lp1 = twoj1_
   len = twoj1_ - 2
   do i = 1,jcap_
     do ll = 1,len
       dxa(ll+lp1) =  dx(ll,i+1) * deps(ll,i+1)
       dxb(ll+lp1) = -dx(ll,i  ) * deps(ll,i+2)
     enddo
     lp1 = lp1 + len
     len = len - 2
   enddo
!
   do i = 1,jcap2_
     do ll = 1,twoj1_
       dx(ll,i) = dx(ll,i) - 1.e0
     enddo
   enddo
!
   do ll = 1,twoj1_
      dxb(ll) = -dx(ll,1) * deps(ll,2)
   enddo
!
!  transpose scalar arrays dxa, dxb, dxc, dxd
!  from cray order to ibm order.
!
   call sph_matrix_trans (dxa,1)
   call sph_matrix_trans (dxb,1)
!
   return
   end subroutine sph_derivative_init
!
!-------------------------------------------------------------------------------
   subroutine sph_matrix_trans(a,kmax)   
!-------------------------------------------------------------------------------
   use paramodel, only : mdim=>lnt2_
   use comfcst, only   : indxnn, indxmm
!-------------------------------------------------------------------------------
   implicit none
!                                                                               
   integer, intent(in)   ::  kmax
   real   , intent(out)  ::  a(mdim,kmax)
!                                                                               
   real                  ::  b(mdim)
   integer               ::  k,m
!-------------------------------------------------------------------------------
   do k = 1,kmax
     do m = 1,mdim                                                            
       b(indxmm(m))=a(m,k)                                                       
     enddo
     do m = 1,mdim                                                            
       a(m,k)=b(m)
     enddo
   enddo
!                                                                               
   return                                                                    
   end subroutine sph_matrix_trans
!-------------------------------------------------------------------------------
