#include <define.h>
   module module_sph_legendre
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [module_sph_legendre]
!      |
!      |--- [sph_gaussian_lat] *
!               |--- [sph_legendre_polinomial] *
!      |--- [sph_poly_funct_init] *
!      |--- [sph_poly_funct] *
!      |--- [sph_poly_epsilon1] *
!      |--- [sph_poly_epsilon2] *
!
!-------------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_gaussian_lat(lgghaf,colrad,wgt,wgtcs,rcs2)
!-------------------------------------------------------------------------------
!
! subroutine: sph_gaussian_lat
!
! abstract:
!   computes the location of the gaussian latitudes for the
!   input lgghaf.  the latitudes are determined by finding
!   the zeros of the legendre polynomials.
!
! program history log:
!   1988-04-05  joseph sela
!   2000-01-01  henry juang            mpi
!   2000-01-01  song-you hong          cvs version
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call sph_gaussian_lat (lgghaf, colrad, wgt, wgtcs, rcs2)
!   input argument list:
!     lgghaf   - number of gaussian latitudes in a hemisphere.
!
!   output argument list:
!     colrad   - array of colatitude of gaussian latitudes
!                in northern hemisphere.
!     wgt      - array of weights at each gaussian latitude
!                required for gaussian quadrature.
!     wgtcs    - array of gaussian weight/sin of colatitude squared.
!     rcs2     - array of reciprocal  of  sin of colatitude squared.
!
!   output files:
!     output   - printout file.
!
!-------------------------------------------------------------------------------
   use paramodel, only  :  latg_,latg2_
   use constant, only   :  pi_
!-------------------------------------------------------------------------------
   real  ::  colrad( latg2_ ),wgt( latg2_ ),wgtcs( latg2_ )                  
   real  ::  rcs2( latg2_ )                                                  
!
!hoon eps=1.e-12
   eps=1.e-14
   si = 1.0
   l2=2*lgghaf
   rl2=l2
   scale = 2.0/(rl2*rl2)
   k1=l2-1
!hoon dradz = pi_ / 360.
   dradz = pi_ / (4.*lgghaf)
   rad = 0.0
!
   do k = 1,lgghaf
     iter=0
     drad=dradz
!
1    call sph_legendre_polinomial(l2,rad,p2)
!
2    p1 =p2
     iter=iter+1
     rad=rad+drad
!
     call sph_legendre_polinomial(l2,rad,p2)
!
     if(sign(si,p1).eq.sign(si,p2)) go to 2
     if(drad.lt.eps)go to 3
     rad=rad-drad
!hoon drad = drad * 0.25
     drad = drad * 0.1
     go to 1
3    continue
     colrad(k)=rad
     phi = rad * 180 / pi_
!
     call sph_legendre_polinomial(k1,rad,p1)
!
     x = cos(rad)
     w = scale * (1.0 - x*x)/ (p1*p1)
     wgt(k) = w
     sn = sin(rad)
     w=w/(sn*sn)
     wgtcs(k) = w
     rc=1./(sn*sn)
     rcs2(k) = rc
!
     call sph_legendre_polinomial(l2,rad,p1)
!
   enddo
!
   return
   end subroutine sph_gaussian_lat
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_legendre_polinomial(n,rad,p)
!-------------------------------------------------------------------------------
!
! abstract: evaluates the unnormalized legendre polynomial
!   of specified degree at a given colatitude using a standard
!   recursion formula.  real arithmetic is used.
!
! program history log:
!   88-04-01  joseph sela
!
! usage:    call sph_legendre_polinomial (n, rad, p)
!   input argument list:
!     n        - degree of legendre polynomial.
!     rad      - real colatitude in radians.
!
!   output argument list:
!     p        - real value of legendre polynomial.
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer, intent(in)   ::  n
   real,    intent(in)   ::  rad
   real,    intent(out)  ::  p
   real     ::  x,y1,y2,y3,g
   integer  ::  i
!-------------------------------------------------------------------------------
   x = cos(rad)
   y1 = 1.0
   y2=x
   do i = 2,n
     g=x*y2
     y3=g-y1+g-(g-y1)/float(i)
     y1=y2
     y2=y3
   enddo
   p=y3
!
   return
   end subroutine sph_legendre_polinomial
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_poly_funct_init
!-------------------------------------------------------------------------------
!
! subprogram:    gpln2t      sets common for subroutine pln2t.
!
! abstract: initializes the constant variables and arrays
!   of a common for subroutine pln2t.
!
! program history log:
!   91-03-14  joseph sela
!
! usage:    call gpln2t
!
! remarks: call subroutine once before calls to pln2t.
!          refer to pln2t for additional documentation.
!
!-------------------------------------------------------------------------------
   use paramodel, only : jcap1_,jcap2_,jcap_,lnut2_,lnut_,twoj1_
   use comfcst, only : deps, rdeps, dx, y, indxmv
!-------------------------------------------------------------------------------
   real     ::  x(jcap1_)
!
!   data ifir /0/
!
!   if  (ifir .eq. 1)  go to 500
!        ifir = 1
!
   do ll = 1,jcap1_
     rdeps(ll) = 0.0
   enddo
   lplus = jcap1_
   len   = jcap1_
   do inde = 2,jcap2_
     do ll = 1,len
       l = ll - 1
       n = l + inde - 1
       rdeps(ll+lplus) = (n*n - l*l) / (4.0 * n*n - 1.0)
     enddo
     lplus = lplus + len
     len = len - 1
   enddo
   do i = jcap2_,lnut_
     rdeps(i) = sqrt(rdeps(i))
   enddo
   do i = 1,lnut_
     deps(2*i-1) = rdeps(i)
     deps(2*i  ) = rdeps(i)
   enddo
!
   ibegin = twoj1_ + 1
   do i = ibegin,lnut2_
     rdeps(i) = 1.0/deps(i)
   enddo
   do ll = 1,jcap1_
     x(ll) = ll*2+1
   enddo
   do ll = 1,jcap1_
     y(ll) = x(ll)/(x(ll)-1.)
   enddo
   do ll = 1,jcap1_
     x(ll) = sqrt(x(ll))
   enddo
   do ll = 1,jcap1_
     dx(2*ll-1) = x(ll)
     dx(2*ll  ) = x(ll)
   enddo
!
!    set index array for transposing vector array
!    from cray order to ibm order.
!
   l=0
   do nn = 1,jcap2_
     lln=min0(jcap2_-nn+1,jcap1_)
     do ll = 1,lln
       indx=((jcap_+3)*(ll-1)-(ll-1)*ll/2+nn)*2
       l=l+2
       indxmv(l-1)=indx-1
       indxmv(l  )=indx
     enddo
   enddo
!
   return
   end subroutine sph_poly_funct_init
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_poly_funct(qlnt,qlnv,colrad,lat)
!-------------------------------------------------------------------------------
!
! subprogram:    pln2t       evaluates associated legendre functions.
!
! abstract: evaluates the required values of the normalized
!   associated legendre function at a prescribed colatitude.
!   a standard recursion relation is used with real arithmetic.
!
! program history log:
!   88-10-25  joseph sela
!
! usage:    call pln2t (qlnt, qlnv, colrad, lat)
!   input argument list:
!     colrad   - half precision colatitudes in radians for which
!                the associated legendre functions are to be
!                computed.
!     lat      - index which indicates the current latitude.
!
!   output argument list:
!     qlnt     - doubled scalar triangle of
!                half precision associated legendre functions.
!     qlnv     - doubled vector triangle of
!                half precision associated legendre functions.
!
!-------------------------------------------------------------------------------
   use paramodel, only : jcap1_,jcap2_,lnut2_,latg2_,twoj1_
   use comfcst, only  : deps, rdeps, dx, y, indxmv
!-------------------------------------------------------------------------------
   real                 ::  qlnv(lnut2_)
   real                 ::  qlnt(lnut2_)
   real                 ::  colrad(latg2_)
   real                 ::  x(jcap1_)
   real                 ::  dpln(lnut2_)
!
!   data ifir /0/
!         part between guards made into sr sph_poly_funct_init.
!         7 dec 1990      m. rozwodoski
!
   colr   = colrad(lat)
   sinlat = cos(colr)
   cos2   = 1.0 - sinlat * sinlat
   prod   = 1.0
   do ll = 1,jcap1_
     x(ll) = 0.5*prod
!    if (prod .lt. 1.0e-75)  prod=0.0
     prod = prod*cos2*y(ll)
   enddo
   do ll = 1,jcap1_
     x(ll) = sqrt(x(ll))
   enddo
   do ll = 1,jcap1_
     dpln(2*ll-1) = x(ll)
     dpln(2*ll  ) = x(ll)
   enddo
   lplus = twoj1_
   do ll = 1,twoj1_
     dpln(ll+lplus) = dx(ll) * sinlat * dpln(ll)
   enddo
   lp2 = 0
   lp1 =     twoj1_
   lp0 = 2 * twoj1_
   len =     twoj1_ - 2
   do n = 3,jcap2_
     do ll = 1,len
       dpln(ll+lp0) = (sinlat * dpln(ll+lp1)                                   &
                  - deps(ll+lp1) * dpln(ll+lp2)) * rdeps(ll+lp0)
     enddo
     lp2 = lp1
     lp1 = lp0
     lp0 = lp0 + len
     len = len - 2
   enddo
! 
!     transpose vector dpln array from cray order to ibm order.
!
   do i = 1,lnut2_
     qlnv(indxmv(i)) = dpln(i)
   enddo
! 
   lpv = 0
   lpt = 0
   len = twoj1_
   do n = 1,jcap1_
     do ll = 1,len
       qlnt(ll+lpt) = qlnv(ll+lpv)
     enddo
     lpv = lpv + len + 2
     lpt = lpt + len
     len = len - 2
   enddo
!
   return
   end subroutine sph_poly_funct
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_poly_epsilon1(eps,jcap)                                               
!-------------------------------------------------------------------------------
!
! subprogram:    sph_poly_epsilon1      computes eps, a function of wave number.
!
! abstract: computes eps, a function of wave number.
!   eps is used in calculating legendre polys. and their derivatives.
!   eps is also used in computing winds from divergence and vorticity.
!
! program history log:
!   88-04-01  joseph sela
!
! usage:    call sph_poly_epsilon1 (eps, jcap)
!   input argument list:
!     jcap     - index indicating the spectral truncation used.
!
!   output argument list:
!     eps      - array computed from sqrt((n**2-l**2)/(4*n**2-1)).
!
!-------------------------------------------------------------------------------
   use paramodel, only : jcap1_,jcap2_
!-------------------------------------------------------------------------------
   real     ::  eps(jcap1_,jcap2_)
   integer  ::  jcap
!
   jcap1 = jcap + 1
   jcap2 = jcap + 2
   do ll = 1,jcap1
     l = ll - 1
     do inde = 2,jcap2
       n = l + inde - 1
       a = (n*n - l*l) / (4.0 * n*n - 1.0)
       eps(ll,inde)=sqrt(a)
     enddo
   enddo
   do ll = 1,jcap1_
     eps(ll,1) = 0.0e0
   enddo
!
   return
   end subroutine sph_poly_epsilon1
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_poly_epsilon2(eps,jcap)
!-------------------------------------------------------------------------------
!
! subprogram:    sph_poly_epsilon2  computes eps, a function of wave number.
!
! abstract: computes eps, a function of wave number.
!   eps is used in calculating legendre polys. and their derivatives.
!   eps is also used in computing winds from divergence and vorticity.
!
! program history log:
!   88-04-01  joseph sela
!
! usage:    call sph_poly_epsilon2 (eps, jcap)
!   input argument list:
!     jcap     - index indicating the spectral truncation used.
!
!   output argument list:
!     eps      - array computed from sqrt((n**2-l**2)/(4*n**2-1)).
!
!-------------------------------------------------------------------------------
   use paramodel, only : jcap1_,jcap2_
!-------------------------------------------------------------------------------
   integer, intent(in)   ::  jcap
   real,    intent(out)  ::  eps(jcap2_,jcap1_)
!
   jcap1 = jcap + 1
   jcap2 = jcap + 2
   do ll = 1,jcap1
     l = ll - 1
     do inde = 2,jcap2
       n = l + inde - 1
       a = (n*n - l*l) / (4.0 * n*n - 1.0)
       eps(inde,ll)=sqrt(a)
     enddo
   enddo
!
   do ll = 1,jcap1_
     eps(1,ll) = 0.0e0
   enddo
!
   return
   end subroutine sph_poly_epsilon2
!-------------------------------------------------------------------------------
   end module module_sph_legendre
