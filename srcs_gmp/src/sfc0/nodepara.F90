!
   subroutine nodepara(msl, nsl, dph, dphnd,                                   &
                      smx, bub, expt, smxnd,                                   &
                      expnd, bubnd, alpha, beta,                               &
                      gamma)
!-------------------------------------------------------------------------------
!
! subprogram: nodepara
!
! abstract: This subroutine sets the thermal node soil parameters to
!   constant values based on those defined for the current grid
!   cells soil type. Thermal node propertiers for the energy
!   balance solution are also set (these constants are used to
!   reduce the solution time required within each iteration).
!
!-------------------------------------------------------------------------------
!
! program history log:
!   2005-02-06  ji chen                development
!   2008-03-09  kyeong hee seol        cvs version setup
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! input variables 
!
   integer     ::  msl          ! number of soil moisture layers
   integer     ::  nsl          ! number of soil thermal nodes
   real        ::  dph(msl)     ! soil moisture layer thickness (m)
   real        ::  dphnd(nsl)   ! thermal node thicknes (m)
   real        ::  smx(msl)     ! maximum layer moisture content (mm)
   real        ::  bub(msl)     ! soil layer bubbling pressure (cm)
   real        ::  expt(msl)    ! soil moisture exponential (N/A)
!
! return variables 
!
   real        ::  smxnd(nsl)   ! maximum moisture thermal node (m3/m3)
   real        ::  expnd(nsl)   ! exponential at thermal node (N/A)
   real        ::  bubnd(nsl)   ! bubbling pressure at thermal node (cm)
   real        ::  alpha(nsl)   ! first thermal eqn term
   real        ::  beta(nsl)    ! second thermal eqn term
   real        ::  gamma(nsl)   ! third thermal eqn term
!
! local variables
!
   integer     ::  nidx, lidx    ! loop indexes
   real        ::  Lsum, Zsum
   logical     ::  Past_Bottom   ! check whether pass soil column bottom
!-------------------------------------------------------------------------------
!
! initialize variables
!
   Past_Bottom = .FALSE.
   lidx = 1
   Lsum = 0.
   Zsum = 0.
!
! set node parameters
!
   do nidx = 1,nsl
     Zsum = Zsum + dphnd(nidx)
     if(Zsum.gt.Lsum.and..not.Past_Bottom) then
       Lsum = Lsum + dph(lidx)
       lidx = lidx + 1
       if( lidx .ge. msl ) then
         Past_Bottom = .TRUE.
         lidx = msl
       endif
       do while(Zsum.gt.Lsum)
         Lsum = Lsum + dph(lidx)
         lidx = lidx + 1
         if( lidx .ge. msl ) then
           Past_Bottom = .TRUE.
           lidx = msl
         endif
       enddo
     endif
!
! node on layer boundary
!
     if(Zsum.eq.Lsum.and.nidx.ne.1.and.lidx.ne.1) then
       smxnd(nidx) = (smx(lidx-1)/dph(lidx-1) +                                &
                      smx(lidx)/dph(lidx))/2.0/1000.0
       expnd(nidx) = (expt(lidx-1)+expt(lidx))/2.0
       bubnd(nidx) = (bub(lidx-1) + bub(lidx))/2.0
!
! node completely in layer
!
     else
       smxnd(nidx) = smx(lidx)/dph(lidx)/1000.0
       expnd(nidx) = expt(lidx)
       bubnd(nidx) = bub(lidx)
     endif
   enddo
!
! compute constant parameters for thermal calculations
!
   do nidx = 1,nsl
     alpha(nidx) = 0.0
     beta(nidx)  = 0.0
     gamma(nidx) = 0.0
   enddo
   do nidx = 1,nsl-1
     alpha(nidx)= dphnd(nidx)+dphnd(nidx+1)
     beta(nidx) = dphnd(nidx)*dphnd(nidx)+                                     &
                  dphnd(nidx+1)*dphnd(nidx+1)
     gamma(nidx)= dphnd(nidx) - dphnd(nidx+1)
   enddo
!
! no flux bottom boundary activated
!
   alpha(nsl) = 2 * dphnd(nsl)
   beta(nsl)  = 2 * dphnd(nsl)*dphnd(nsl)
   gamma(nsl) = 0.0
!
   return
   end
!-------------------------------------------------------------------------------
