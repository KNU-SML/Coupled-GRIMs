#include <define.h>
   module module_sph_semi_implicit
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [sph_semi_implicit]
!      |
!      |--- [sph_semi_gwave] *
!      |--- [sph_semi_hydro] *
!      |--- [sph_semi_hydro_hybrid] *
!      |--- [sph_semi_impl_integrate] *
!      |--- [sph_semi_inverse_matrix] *
!
!-------------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------------
!
#ifndef DFS
!-------------------------------------------------------------------------------
   subroutine sph_semi_gwave(deltim,am,bm,gv,sv,cm)
!-------------------------------------------------------------------------------
!
! subprogram:    sph_semi_gwave   setup for semi-implicit time integration.
!
! abstract: computes matrix inverse of rhs divergence autodependence
!           in the semi-implicit treatment of the gravity wave modes.
!
! program history log:
!   1993-03-15  mark iredell           initial mrf
!   2000-01-01  hann-ming henry juang  mpi
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call sph_semi_gwave(deltim,am,bm,gv,sv,cm)
!   input argument list:
!     deltim   - timestep
!     am       - div dependence on temp (hydrostatic)
!     bm       - temp dependence on div (energy conversion)
!     gv       - div dependence on lnps (pressure gradient)
!     sv       - lnps dependence on div (continuity)
!     cm       - div autodependence (cm=sv*gv+am*bm)
!
!-------------------------------------------------------------------------------
   use paramodel, only : jcap_,km=>levs_
   use constant, only  : rd_,rerth_
   use comfcst, only   : dt, gvdt, svdt, amdt, bmdt, dm
!-------------------------------------------------------------------------------
   real                 ::  am(km,km),bm(km,km),sv(km),gv(km),cm(km,km)   
   real,parameter       ::  raa=rd_/(rerth_**2)
   real,parameter       ::  tol=1.e-12
#define DEFAULT
#ifdef MINV
#undef DEFAULT
   real                 ::  work(2*km)
#endif
#ifdef DEFAULT
   integer              ::  iwork(2*km)
#endif
!
   dt=deltim
   do k = 1,km
     gvdt(k)=dt*gv(k)
     svdt(k)=dt*sv(k)
   enddo
   do kj = 1,km
     do ki = 1,km
       amdt(ki,kj)=dt*am(ki,kj)
       bmdt(ki,kj)=dt*bm(ki,kj)
       dm(ki,kj,0)=0.
     enddo
   enddo
   do kj = 1,km
     dm(kj,kj,0)=1.
   enddo
#ifdef ORIGIN_THREAD
!$doacross share(dt,dm,cm,km,tol),
!$&        local(n,dt2nn1,kj,work,iwork,det)
#endif
#ifdef OPENMP
!$omp parallel do private(n,dt2nn1,kj,work,iwork,det)
#endif
   do n = 1,jcap_
     dt2nn1=dt**2*(n*(n+1))                                                  
     do kj = 1,km
       do ki = 1,km
#ifndef HYBRID
         dm(ki,kj,n)=dm(ki,kj,0)-dt2nn1*cm(ki,kj)
#else
         dm(ki,kj,n)=dm(ki,kj,0)+dt2nn1*cm(ki,kj)
#endif
       enddo
     enddo
#define DEFAULT
#ifdef MINV
#undef DEFAULT
     call minv(dm(1,n),km,km,work,det,tol,0,1)
#endif
#ifdef DEFAULT
     call sph_semi_inverse_matrix(dm(1,1,n),km,det,iwork(1),iwork(km+1))
#endif
   enddo
#undef DEFAULT
!
   return
   end subroutine sph_semi_gwave
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_semi_hydro(d,t,q,x,y,z,u,v,lnt2,lnoffset)
!-------------------------------------------------------------------------------
!
! subprogram:    sph_semi_hydro      semi-implicit time integration.
!
! abstract: integrates divergence, temperature and log surface pressure
!           semi-implicitly in time.
!
! program history log:
!   89-03-15  joseph sela
!   93-03-15  mark iredell   linear matrices passed in common
!   99-03-15  hann-ming henry juang  pass offset and length
!
! usage:    call sph_semi_hydro(d,t,q,x,y,z,u,v)
!   input argument list:
!     d        - divergence at time t-dt
!     t        - temperature at time t-dt
!     q        - ln(psfc) at time t-dt
!     x        - divergence nonlinear tendency at time t
!     y        - temperature nonlinear tendency at time t
!     z        - ln(psfc) nonlinear tendency at time t
!
!   output argument list:
!     x        - divergence at time t+dt
!     y        - temperature at time t+dt
!     z        - ln(psfc) at time t+dt
!     u        - work array
!     v        - work array
!
!-------------------------------------------------------------------------------
   use paramodel, only : km=>levs_,jcap=>jcap_,lnt22=>LNT22S
   use comio
   use comfcst, only : dt, gvdt, svdt, amdt, bmdt, dm
!-------------------------------------------------------------------------------
   integer              ::  lnt2,lnoffset
   real                 ::  d(lnt22,km),t(lnt22,km),q(lnt22)
   real                 ::  x(lnt22,km),y(lnt22,km),z(lnt22)
   real                 ::  u(lnt22,km),v(lnt22,km)
   integer              ::  i,j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  explicitly integrate lnps and temperature halfway in time.
!
   do i = 1,lnt2
     z(i)=q(i)+dt*z(i)
   enddo
   do k = 1,km
     do i = 1,lnt2
       y(i,k)=t(i,k)+dt*y(i,k)
     enddo
   enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  compute linear dependence of divergence on lnps and temperature.
!  explicitly integrate divergence halfway including linear terms.
!
   do k = 1,km
     do i = 1,lnt2
       v(i,k)=0.
     enddo
     do j = 1,km
       do i = 1,lnt2
         v(i,k)=v(i,k)+amdt(k,j)*y(i,j)
       enddo
     enddo
     do i = 1,lnt2
       ii=lnoffset+i
       u(i,k)=d(i,k)+dt*x(i,k)+snnp1(ii)*(v(i,k)+gvdt(k)*z(i))
     enddo
   enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  solve helmholz equation for semi-implicit divergence.
!
   do k = 1,km
     do i = 1,lnt2
       v(i,k)=0.
     enddo
!
     do j = 1,km
       do i = 1,lnt2,2
         ii=lnoffset+i
         n=ndex(ii)
         v(i,k)=v(i,k)+dm(k,j,n)*u(i,j)
         v(i+1,k)=v(i+1,k)+dm(k,j,n)*u(i+1,j)
       enddo
     enddo
   enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  back solve for lnps.
!
   do j = 1,km
     do i = 1,lnt2
       z(i)=z(i)+svdt(j)*v(i,j)
     enddo
   enddo
   do i = 1,lnt2
     z(i)=2*z(i)-q(i)
   enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  back solve for temperature and divergence.
!
   do k = 1,km
     do j = 1,km
       do i = 1,lnt2
         y(i,k)=y(i,k)+bmdt(k,j)*v(i,j)
       enddo
     enddo
     do i = 1,lnt2
       y(i,k)=2*y(i,k)-t(i,k)
       x(i,k)=2*v(i,k)-d(i,k)
     enddo
   enddo
!
   return
   end subroutine sph_semi_hydro
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_semi_hydro_hybrid(d,t,q,x,y,z,di,te,qe,u,v,lnt2,lnoffset)
!-------------------------------------------------------------------------------
!
! subprogram:    sph_semi_hydro_hybrid      semi-implicit time integration.
!
! abstract: integrates divergence, temperature and log surface pressure
!           semi-implicitly in time.
!
! program history log:
!   89-03-15  joseph sela
!   93-03-15  mark iredell  linear matrices passed in common
!   99-03-15  hann-ming henry juang  pass offset and lengt
!
! usage:    call sph_semi_hydro_hybrid(d,t,q,x,y,z,u,v)
!   input argument list:
!     d        - divergence at time t-dt
!     t        - temperature at time t-dt
!     q        - ln(psfc) at time t-dt
!     x        - divergence nonlinear tendency at time t
!     y        - temperature nonlinear tendency at time t
!     z        - ln(psfc) nonlinear tendency at time t
!
!   output argument list:
!     x        - divergence at time t+dt
!     y        - temperature at time t+dt
!     z        - ln(psfc) at time t+dt
!     u        - work array
!     v        - work array
!
!-------------------------------------------------------------------------------
   use paramodel, only : lnt22=>LNT22S,jcap=>jcap_,km=>levs_
   use comio
   use comfcst, only : dt, gvdt, svdt, amdt, bmdt, dm
!-------------------------------------------------------------------------------
   real d(lnt22,km),t(lnt22,km),q(lnt22)
   real di(lnt22,km),te(lnt22,km),qe(lnt22)
   real x(lnt22,km),y(lnt22,km),z(lnt22),u(lnt22,km),v(lnt22,km)
!-------------------------------------------------------------------------------
!
!  explicitly integrate lnps and temperature halfway in time.
!kei because full term is included in x (d(div)/dt), q and T time n terms
!kei should be removed before. they will be added back.
!
   do i = 1,lnt2
     z(i)=q(i)+dt*z(i)-qe(i)
   enddo
   do k = 1,km
     do i = 1,lnt2
       y(i,k)=t(i,k)+dt*y(i,k)-te(i,k)
     enddo
   enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  compute linear dependence of divergence on lnps and temperature.
!  explicitly integrate divergence halfway including linear terms.
!
   do k = 1,km
     do i = 1,lnt2
       v(i,k)=0.
     enddo
     do j = 1,km
       do i = 1,lnt2
         v(i,k)=v(i,k)+amdt(k,j)*y(i,j)
       enddo
     enddo
     do i = 1,lnt2
       ii=lnoffset+i
       u(i,k)=d(i,k)+dt*x(i,k)+snnp1(ii)*(v(i,k)+gvdt(k)*z(i))
     enddo
   enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  solve helmholz equation for semi-implicit divergence.
!
   do k = 1,km
     do i = 1,lnt2
       v(i,k)=0.
     enddo
     do j = 1,km
       do i = 1,lnt2,2
         ii=lnoffset+i
         n=ndex(ii)
         v(i,k)=v(i,k)+dm(k,j,n)*u(i,j)
         v(i+1,k)=v(i+1,k)+dm(k,j,n)*u(i+1,j)
       enddo
     enddo
   enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  back solve for lnps.
!
   do j = 1,km
     do i = 1,lnt2
       z(i)=z(i)-svdt(j)*v(i,j)
     enddo
   enddo
   do i = 1,lnt2
     z(i)=z(i)+qe(i)
   enddo
!
   do i = 1,lnt2
     z(i)=2*z(i)-q(i)
   enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  back solve for temperature and divergence.
!
   do k = 1,km
     do j = 1,km
       do i = 1,lnt2
         y(i,k)=y(i,k)-bmdt(k,j)*v(i,j)
       enddo
     enddo
     do i = 1,lnt2
       y(i,k)=y(i,k)+te(i,k)
     enddo
     do i = 1,lnt2
       y(i,k)=2*y(i,k)-t(i,k)
       x(i,k)=2*v(i,k)-d(i,k)
     enddo
   enddo
!
   return                                                                    
   end subroutine sph_semi_hydro_hybrid
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_semi_impl_integrate(d,t,q,x,y,z,u,v,lnt2,lnoffset)
!-------------------------------------------------------------------------------
!
! subprogram:    sph_semi_impl_integrate
!   - implicit adjustment of physics tendencies.
!
! abstract: extends the semi-implicit time integration to include
!           the phyical forcing terms computed in phys_main_driver.
!
! program history log:
!   91-03-15  mark iredell
!   93-03-15  mark iredell   change argument list
!   93-06-22  hann-ming henry juang  pass offset and length for mpp
!
! usage:    call sph_semi_impl_integrate(d,t,q,x,y,z,u,v)
!   input argument list:
!     d        - divergence before adjustment
!     t        - temperature before adjustment
!     q        - ln(psfc) before adjustment
!     x        - divergence tendency adjustment
!     y        - temperature tendency adjustment
!     z        - ln(psfc) tendency adjustment
!
!   output argument list:
!     d        - divergence adjusted
!     t        - temperature adjusted
!     q        - ln(psfc) adjusted
!     u        - work array
!     v        - work array
!
!-------------------------------------------------------------------------------
   use paramodel, only : km=>levs_,lnt22=>LNT22S
   use comio
   use comfcst, only : dt, gvdt, svdt, amdt, bmdt, dm
!-------------------------------------------------------------------------------
   real                 ::  d(lnt22,km),t(lnt22,km),q(lnt22)
   real                 ::  x(lnt22,km),y(lnt22,km),z(lnt22)
   real                 ::  u(lnt22,km),v(lnt22,km)
!
!  compute linear dependence of divergence on lnps and temperature.
!
   do k = 1,km
     do i = 1,lnt2
       v(i,k)=0.
     enddo
#ifdef CRAY_THREAD
!fpp$ unroll l
#endif
     do j = 1,km
       do i = 1,lnt2
         v(i,k)=v(i,k)+amdt(k,j)*y(i,j)
       enddo
     enddo
     do i = 1,lnt2
       ii=lnoffset+i
       u(i,k)=x(i,k)+snnp1(ii)*(v(i,k)+gvdt(k)*z(i))
     enddo
   enddo
!
!  solve helmholz equation for semi-implicit divergence.
!
#ifdef CRAY_THREAD
!mic$ do all autoscope
#endif
   do k = 1,km
     do i = 1,lnt2
       v(i,k)=0.
     enddo
#ifdef CRAY_THREAD
!fpp$ unroll l
#endif
     do j = 1,km
       do i = 1,lnt2,2
         ii=lnoffset+i
         n=ndex(ii)
         v(i,k)=v(i,k)+dm(k,j,n)*u(i,j)
         v(i+1,k)=v(i+1,k)+dm(k,j,n)*u(i+1,j)
       enddo
     enddo
   enddo
!
!  back solve for lnps.
!
#ifdef CRAY_THREAD
!fpp$ unroll l
#endif
   do j = 1,km
     do i = 1,lnt2
#ifndef HYBRID
       q(i)=q(i)+svdt(j)*v(i,j)
#else
       q(i)=q(i)-svdt(j)*v(i,j)
#endif
     enddo
   enddo
   do i = 1,lnt2
     q(i)=q(i)+z(i)
   enddo
!
!  back solve for temperature and divergence.
!
#ifdef CRAY_THREAD
!mic$ do all autoscope
#endif
   do k = 1,km
#ifdef CRAY_THREAD
!fpp$ unroll l
#endif
     do j = 1,km
       do i = 1,lnt2
#ifndef HYBRID
         t(i,k)=t(i,k)+bmdt(k,j)*v(i,j)
#else
         t(i,k)=t(i,k)-bmdt(k,j)*v(i,j)
#endif
       enddo
     enddo
     do i = 1,lnt2
       t(i,k)=t(i,k)+y(i,k)
       d(i,k)=d(i,k)+v(i,k)
     enddo
   enddo
!
   return
   end subroutine sph_semi_impl_integrate
#endif  /* not DFS */
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_semi_inverse_matrix (a,n,d,l,m)
!-------------------------------------------------------------------------------
!
! subprogram:    sph_semi_inverse_matrix  computes inverse of matrix in place.
!
! abstract: computes inverse of matrix in place using
!   gauss-jordan reduction with max pivot.
!
! program history log:
!   88-04-21  ibm scientific subroutine package.
!
! usage:    call sph_semi_inverse_matrix (a, n, d, l, m)
!   input argument list:
!     a        - square matrix which will be inverted.
!                matrix a will be destroyed and replaced by inverse.
!     n        - order of matrix a.
!
!   output argument list:
!     a        - inverse of input matrix a.
!     d        - determinant of a inverse.
!                if d=0.0, matrix a is singular.
!     l        - work vector of length n.
!     m        - work vector of length n.
!
!-------------------------------------------------------------------------------
!
!     ..................................................................
!
!        ................
!
!        purpose
!           invert a matrix
!                                                                               
!        usage
!           call sph_semi_inverse_matrix (a,n,d,l,m)
!
!        description of parameters
!           a - input matrix, destroyed in computation and replaced by
!               resultant inverse.
!           n - order of matrix a
!           d - resultant determinant
!           l - work vector of length n
!           m - work vector of length n
!
!        remarks
!           matrix a must be a general matrix
!
!        .............................................
!           none
!
!        method
!           the standard gauss-jordan method is used. the determinant
!           is also calculated. a determinant of zero indicates that
!           the matrix is singular.
!
!     ..................................................................
!
!-------------------------------------------------------------------------------
   integer              ::  n
   real                 ::  a(n*n),d
   integer              ::  l(n),m(n)
!
!        ...............................................................
!
!        if a double precision version of this routine is desired, the
!        c in column 1 should be removed from the double precision
!        statement which follows.
!
!     double precision a, d, biga, hold
!
!        the c must also be removed from double precision statements
!        appearing in other routines used in conjunction with this
!        routine.
!
!        the double precision version of this sr........ must also
!        contain double precision fortran functions.  abs in statement
!        10 must be changed to dabs.
!
!        ...............................................................
!
!        search for largest element
!
   d=1.0e0
   nk=-n
   do k = 1,n
     nk=nk+n
     l(k)=k
     m(k)=k
     kk=nk+k
     biga=a(kk)
     do j = k,n
       iz=n*(j-1)
       do 20 i = k,n
         ij=iz+i
!
!  10 if (dabs(biga)-dabs(a(ij))) 15,20,20
!
   10    if(abs(biga)-abs(a(ij))) 15,20,20
   15    biga=a(ij)
         l(k)=i
         m(k)=j
   20  continue
     enddo
!
!        interchange rows
!
     j=l(k)
     if(j-k) 35,35,25
  25 ki=k-n
     do i = 1,n
       ki=ki+n
       hold=-a(ki)
       ji=ki-k+j
       a(ki)=a(ji)
       a(ji) =hold
     enddo
!
!        interchange columns
!
 35  i=m(k)
     if(i-k) 45,45,38
 38  jp=n*(i-1)
     do j = 1,n
       jk=nk+j
       ji=jp+j
       hold=-a(jk)
       a(jk)=a(ji)
       a(ji) =hold
     enddo
!
!        divide column by minus pivot (value of pivot element is
!        contained in biga)
!
 45  if(biga) 48,46,48
 46  d=0.0e0
     return
!
 48  do 55 i = 1,n
       if(i-k) 50,55,50
 50    ik=nk+i
       a(ik)=a(ik)/(-biga)
 55  continue
!
!        reduce matrix
!
     do 65 i = 1,n
       ik=nk+i
       ij=i-n
       do 65 j = 1,n
         ij=ij+n
         if(i-k) 60,65,60
  60     if(j-k) 62,65,62
  62     kj=ij-i+k
         a(ij)=a(ik)*a(kj)+a(ij)
  65 continue
!
!        divide row by pivot
!
     kj=k-n
     do 75 j = 1,n
       kj=kj+n
       if(j-k) 70,75,70
70     a(kj)=a(kj)/biga
75   continue
!
!        product of pivots
!
     d=d*biga
!
!        replace pivot by reciprocal
!
     a(kk)=1.0e0/biga
   enddo
!
!        final row and column interchange
!
   k=n
100 k=(k-1)
   if(k) 150,150,105
105 i=l(k)
   if(i-k) 120,120,108
108 jq=n*(k-1)
   jr=n*(i-1)
   do j = 1,n
     jk=jq+j
     hold=a(jk)
     ji=jr+j
     a(jk)=-a(ji)
     a(ji) =hold
   enddo
!
120 j=m(k)
   if(j-k) 100,100,125
125 ki=k-n
   do i = 1,n
     ki=ki+n
     hold=a(ki)
     ji=ki-k+j
     a(ki)=-a(ji)
     a(ji) =hold
   enddo
!
   go to 100
!
150 return
   end subroutine sph_semi_inverse_matrix
!-------------------------------------------------------------------------------
   end module module_sph_semi_implicit
