#include <define.h>
   subroutine exmoist(stdprs,qin,t,qou,ijdim,kdimi,kdimqi,xlat)
!-------------------------------------------------------------------------------
!
! subprogram: exmoist
!
! abstract: 
!  this code extrapolates moisture up into stratosphere
!  using exponential decrease up to a specified q at 80 mb,
!  then use lims profile above(see middle atmos program handbook,
!  vol 22,sep 86,ed. j.m.russell ---> j alpert has copy)....
!
!  input - slin  - sigma level 
!        - pstar - surface pressure (cb)
!        - qin   - specific humidity
!        - t     - virtual temperature
!        - xlat  - latitude in radians
!        - ijdim - number of grid points
!        - kdimqi- number of input moisture levels
!        - kdimi - number of output sigma levels
!
!  output - qou - moisture (use lims data abuv 8 cb)
!         - t   - thermodynamic temperature back into t
!
!-------------------------------------------------------------------------------
   use paramodel, only : pi_,psat_,rd_,rv_,t0c_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer ijdim
   integer kdimi,kdimqi
   real stdprs(kdimi)
   real qin(ijdim,kdimqi)
   real t(ijdim,kdimi)
   real xlat(ijdim)
   real qou(ijdim,kdimi)
!
   real qout(kdimi)
   real qsat(kdimi),prs(kdimi)
!
   integer i,j,k,ij
   integer jdx,ldry,kdy
   real pmin,pmm,xlnpm
   real expon,es,dlat,rh,dx,q8,q1,xlnqm,dy
   real qqout,xlnpq,xlnqq,xlamb
!
!  qq stores lims (satellite-1986-see j.alpert) observed h2o above
!  80 mb...(k,j)-k=1,8 for 80,70,..,10 mb
!  j=1,19 for every 10 deg lat (starting at n.pole)
!
   real qq(8,19)
   save qq
   data qq /32*5.  ,  4.8,4.55,4.5,4.6,4.7,4.75,4.75,4.7,                      &
    4.,3.85,3.9,4.,4.25,4.5,4.5,4.6, 3.,2.75,3.,3.3,3.7,4.2,4.3,4.5,           &
              2.1,2.,2.,2.5,3.,3.9,4.,4. ,                                     &
    2.,2.,2.,2.3,2.8,3.5,3.7,3.75 , 2.,2.,2.,2.3,2.7,3.25,3.45,3.75,           &
    2.,2.,2.,2.3,2.8,3.5,3.7,3.75 , 2.1,2.,2.,2.5,3.,3.9,4.,4. ,               &
              3.,2.75,3.,3.3,3.7,4.2,4.3,4.5,                                  &
         4.,3.85,3.9,4.,4.25,4.5,4.5,4.6 ,                                     &
              4.8,4.55,4.5,4.6,4.7,4.75,4.75,4.7 , 32*5.  /
!-------------------------------------------------------------------------------
!
!  first convert mdl virtual temp to thermodynmic temp
!
   pmin = 8.0
   pmm  = 1.0
   xlnpm =  log(pmin)
!
!  get correct units of lims moisture
!
   do j = 1,19
     do k = 1,8
       qq(k,j) = qq(k,j) * 1.e-6
     enddo
   enddo
!
!  get lyr pressure and then saturated moisture
!
   do k = 1,kdimi
     prs(k) = stdprs(k)*0.1      ! kpa
   enddo

   do i = 1,ijdim
!
!  get latitude in degrees..
!
     dlat = xlat(i)* 180.0/pi_
!
!  model temp is virtual (since hydrostatic from hgts initially),
!
!  do k=1,kdimqi
!    if (qin(i,k).le.0.0) go to 2
!    t(i,k)=t(i,k)*(1.0+qin(i,k))/(1.0+rv_/rd_*qin(i,k))
!  2 continue
!  enddo
!
!  compute saturation specific humidity(dimensionless) from
!  temperature t (deg k) and pressure (cb)
!  conversion to specific humidity follows from definition
!
     do k = 1,kdimi
       expon = 7.50*(t(i,k)-t0c_)/((t(i,k)-t0c_)+237.30)
       es = psat_*1.e-2 * 10.0**expon
       qsat(k)=rd_/rv_*es/(prs(k)*10.0-(1.0-rd_/rv_)*es)
     enddo
!
!  limit moisture in lowest levh layers--rh le 1. but ge .15
!          the latter avoids the negative q problem...
!
     do k = 1,kdimqi
       qout(k) = qin(i,k)
       rh = qin(i,k) / qsat(k)
     enddo
!
!  obtain 8.cb (k=1) and 1.cb (k=8) valu by horiz interpolation
!
     jdx = (90.-dlat)/10.+1.
     dx = (90.-dlat)/10.+1.-jdx
     q8  = qq(1,jdx)*(1.-dx)+qq(1,jdx+1)*dx
     q1  = qq(8,jdx)*(1.-dx)+qq(8,jdx+1)*dx
     ldry = kdimqi + 1
! 
     if( ldry .le. kdimi ) then
!
       xlnqm =  log(q8)
!
!  extrapolate moisture to min valu(q8) at pressure pmin(80 mb)
!   use exponential decrease from layer levh----
!   i.e.  q=q(levh)*(p/p(levh)) ** xlamb
!   where xlamb is computed to fit q8 at pmin and q,p at levh.
!
       qqout = qin(i,kdimqi)
       rh = qqout / qsat(kdimqi)
       if (rh.le..15) qqout = .15 * qsat(kdimqi)
       if (rh.gt. 1.0) qqout = qsat(kdimqi)
       xlnpq =  log(prs(kdimqi))
       xlnqq =  log(qqout)
       xlamb = (xlnqm-xlnqq) / (xlnpm-xlnpq)
       do k = ldry,kdimi
         if (prs(k).ge.pmin) go to 21
         if (prs(k).ge.pmm ) go to 22
!
!  above pmm(10 mb) use constant value from table
!
         qout(k) = q1
         go to 13
   21    continue
         qout(k) = qout(kdimqi)*(prs(k)/prs(kdimqi)) ** xlamb
         go to 13
   22    continue
!
!  above 8 cb so complete linear interp from table..
!
         kdy = 9. - prs(k)
         dy = 9. - prs(k) - kdy
         qout(k) = qq(kdy,jdx  )*(1.-dy)*(1.-dx)+qq(kdy+1,jdx+1)*dy*dx         &
                  +qq(kdy,jdx+1)*(1.-dy)*dx+qq(kdy+1,jdx)*dy*(1.-dx)
   13    rh = qout(k) / qsat(k)
         if(rh.gt.1.) qout(k)=qsat(k)
       enddo
     endif
!
!  store extrapolated moisture
!
     do k = 1,kdimi
       qou(i,k) = qout(k)
     enddo
!
!  temp. is restored to virtual temp.
!
!  do k=1,kdimqi
!    if (qin(i,k).gt.0.) then
!      t(i,k) = t(i,k) * (1. + rv_/rd_ * qin(i,k))/(1. + qin(i,k))
!    endif
!  enddo
!
   enddo
!
   return
   end
