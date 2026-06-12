#include <define.h>
   subroutine post_grid_interp
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!  [phys_river_main_driver] or [post_pgb_2d]
!      |
!      |--- (post_grid_interp) ----- [post_gauss2latlon] *
!                                     |-- [post_gauss2latlon_river] *
!                                     |-- [post_rmp2latlon_river] *
!  [post_pgb_2d]
!      |--- (post_grid_interp) ----- [post_sigma2p_avg] *
!                                     |-- [post_sigma2pressure] *
!                                     |-- [post_sigma2tropopause] *
!                                     |-- [post_wave2grid_fcst] *
!
! program history log:
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   end subroutine post_grid_interp
!
!-------------------------------------------------------------------------------
   subroutine post_gauss2latlon(gauin,imxin,jmxin,                             &
                                xlonw,xlatn,dxlon,dxlat,regout,imxout,jmxout)
!-------------------------------------------------------------------------------
!
!  interpolation from lat/lon grid to other lat/lon grid
!
!-------------------------------------------------------------------------------
#ifdef RIVER
   use comfcst, only   : iindx1, iindx2, jindx1, jindx2, ddx, ddy
   real                 ::  gauin (imxin,jmxin)
   real                 ::  regout(imxout,jmxout)
!
   real                 ::  gaul(jmxin),regl(jmxout)
!
   save ifp
   data ifp/0/
!
   if(ifp.ne.0) go to 111
   ifp=1
!
#ifdef DBG
   write(6,*) 'imxin=',imxin,' jmxin=',jmxin
   write(6,*) 'xlatn=',xlatn,' xlonw=',xlonw,                                  &
              'dxlat=',dxlat,' dxlon=',dxlon
   write(6,*) 'imxout=',imxout,' jmxout=',jmxout
#endif
!
   call gaulat(gaul,jmxin)
   do j = 1,jmxin
     gaul(j)=90.-gaul(j)
   enddo
!
   do j = 1,jmxout
     regl(j)=xlatn-float(j-1)*dxlat
   enddo
!
   dxin =360./float(imxin)
!
   do i = 1,imxout
     alamd=xlonw+float(i-1)*dxlon
     if(alamd.lt.0.) alamd=360.+alamd
     i1=alamd/dxin+1.001
     if(i1.gt.imxin) i1=1
     iindx1(i)=i1
     i2=i1+1
     if(i2.gt.imxin) i2=1
     iindx2(i)=i2
     ddx(i)=(alamd-float(i1-1)*dxin)/dxin
   enddo
!
   do 40 j = 1,jmxout
   aphi=regl(j)
   do 50 jj = 1,jmxin
   if(aphi.lt.gaul(jj)) go to 50
   j2=jj
   go to 42
   50 continue
   j2=jmxin
   42 continue
   if(j2.gt.2) go to 43
   j1=1
   j2=2
   go to 44
   43 continue
   if(j2.le.jmxin) go to 45
   j1=jmxin-1
   j2=jmxin
   go to 44
   45 continue
   j1=j2-1
   44 continue
   jindx1(j)=j1
   jindx2(j)=j2
   ddy(j)=(aphi-gaul(j1))/(gaul(j2)-gaul(j1))
   40 continue
!
#ifdef DBG
   write(6,*) 'iindx1'
   write(6,*) (iindx1(n),n=1,imxout)
   write(6,*) 'iindx2'
   write(6,*) (iindx2(n),n=1,imxout)
   write(6,*) 'jindx1'
   write(6,*) (jindx1(n),n=1,jmxout)
   write(6,*) 'jindx2'
   write(6,*) (jindx2(n),n=1,jmxout)
   write(6,*) 'ddy'
   write(6,*) (ddy(n),n=1,jmxout)
   write(6,*) 'ddx'
   write(6,*) (ddx(n),n=1,jmxout)
#endif
  111 continue
!
   do j = 1,jmxout
     y=ddy(j)
     j1=jindx1(j)
     j2=jindx2(j)
     do i = 1,imxout
       x=ddx(i)
       i1=iindx1(i)
       i2=iindx2(i)
       regout(i,j)=(1.-x)*(1.-y)*gauin(i1,j1)+                                 &
                   (1.-y)*x*gauin(i2,j1)+                                      &
                   (1.-x)*y*gauin(i1,j2)+x*y*gauin(i2,j2)
     enddo
   enddo
!
   sum1=0.
   sum2=0.
   do i = 1,imxin
     sum1=sum1+gauin(i,1)
     sum2=sum2+gauin(i,jmxin)
   enddo
   sum1=sum1/float(imxin)
   sum2=sum2/float(imxin)
!
   do i = 1,imxout
     if(abs(regl(1)).eq.90.) then
       regout(i,     1)=sum1
     endif
     if(abs(regl(jmxout)).eq.90.) then
       regout(i,jmxout)=sum2
     endif
   enddo
!
#endif  /* RIVER */
   return
   end subroutine post_gauss2latlon
!
!-------------------------------------------------------------------------------
      subroutine post_gauss2latlon_river(gauin,imxin,jmxin,                    &
                        xlonw,xlatn,dxlon,dxlat,regout,imxout,jmxout)
!-------------------------------------------------------------------------------
!
!  interpolation from lat/lon grid to other lat/lon grid
!
!-------------------------------------------------------------------------------
#ifdef RIVER
   use paramodel, only : io2_, jo2_
!  parameter(io2_=io_,jo2_=jo_)
!  parameter(io2_=360,jo2_=181) !! for 1-deg but not T126
!  parameter(io2_=720,jo2_=361) !! for 0.5-deg but not T248
   real ::  gauin (imxin,jmxin)
   real ::  regout(imxout,jmxout)
!
   real ::  gaul(jmxin),regl(jmxout)
!
   real,save,allocatable :: iindx1(:),iindx2(:),                               &
                            jindx1(:),jindx2(:),                               &
                               ddx(:),   ddy(:)
   save ifp
   data ifp/0/
!
   if(.not.allocated(iindx1)) allocate(iindx1(io2_))
   if(.not.allocated(iindx2)) allocate(iindx2(io2_))
   if(.not.allocated(jindx1)) allocate(jindx1(jo2_-1))
   if(.not.allocated(jindx2)) allocate(jindx2(jo2_-1))
   if(.not.allocated(ddx)   ) allocate(ddx(io2_))
   if(.not.allocated(ddy)   ) allocate(ddy(jo2_-1))
!
   if(ifp.ne.0) go to 111
   ifp=1
!
!      write(6,*) 'imxin=',imxin,' jmxin=',jmxin
!      write(6,*) 'xlatn=',xlatn,' xlonw=',xlonw,
!     1           ' dxlat=',dxlat,' dxlon=',dxlon
!      write(6,*) 'imxout=',imxout,' jmxout=',jmxout
!
   call gaulat(gaul,jmxin)
   do j = 1,jmxin
     gaul(j)=90.-gaul(j)
   enddo
!  print *,'gaul=',gaul
!
   do j = 1,jmxout
     regl(j)=xlatn-float(j-1)*dxlat
   enddo
!
   do 40 j = 1,jmxout
     aphi=regl(j)
     do 50 jj = 1,jmxin
       if(aphi.lt.gaul(jj)) go to 50
       j2=jj
       go to 42
50   continue
     j2=jmxin
42   continue
     if(j2.gt.2) go to 43
     j1=1
     j2=2
     go to 44
43   continue
     if(j2.le.jmxin) go to 45
     j1=jmxin-1
     j2=jmxin
     go to 44
45   continue
     j1=j2-1
44   continue
     jindx1(j)=j1
     jindx2(j)=j2
     ddy(j)=(aphi-gaul(j1))/(gaul(j2)-gaul(j1))
40 continue
!
   dxin =360./float(imxin)
!
   do i = 1,imxout
     alamd=xlonw+float(i-1)*dxlon
     if(alamd.lt.0.) alamd=360.+alamd
     i1=alamd/dxin+1.001
     if(i1.gt.imxin) i1=1
     iindx1(i)=i1
     i2=i1+1
     if(i2.gt.imxin) i2=1
     iindx2(i)=i2
     ddx(i)=(alamd-float(i1-1)*dxin)/dxin
!    print*,'ddx',alamd,ddx(i)
   enddo
!
!
#ifdef DBG
   write(6,*) 'iindx1'
   write(6,*) (iindx1(n),n=1,imxout)
   write(6,*) 'iindx2'
   write(6,*) (iindx2(n),n=1,imxout)
   write(6,*) 'jindx1'
   write(6,*) (jindx1(n),n=1,jmxout)
   write(6,*) 'jindx2'
   write(6,*) (jindx2(n),n=1,jmxout)
   write(6,*) 'ddy'
   write(6,*) (ddy(n),n=1,jmxout)
   write(6,*) 'ddx'
   write(6,*) (ddx(n),n=1,jmxout)
#endif
111 continue
!
   do j = 1,jmxout
     y=ddy(j)
     j1=jindx1(j)
     j2=jindx2(j)
     do i = 1,imxout
       x=ddx(i)
       i1=iindx1(i)
       i2=iindx2(i)
       regout(i,j)=(1.-x)*(1.-y)*gauin(i1,j1)+                                 &
                   (1.-y)*x*gauin(i2,j1)+                                      &
                   (1.-x)*y*gauin(i1,j2)+x*y*gauin(i2,j2)
     enddo
   enddo
!
#endif
   return
   end subroutine post_gauss2latlon_river
!
!-------------------------------------------------------------------------------
   subroutine post_rmp2latlon_river(regin,imxin,jmxin,llout,imxout,jmxout,     &
                        delxo,delyo)
!-------------------------------------------------------------------------------
!
!  interpolation from rmp grid to lat/lon grid for river model (TRIP)
!
!-------------------------------------------------------------------------------
#ifdef RMP
#ifdef RIVER
!
   use rscomloc
   use module_trans, only     :  rmp_trans2output_grid
   use constant, only : pi_
!
   integer            :: imxin,jmxin,imxout,jmxout
   integer            :: ijmxin
   real               :: delxo, delyo
   real               :: regin(imxin,jmxin), llout(imxout,jmxout)
   integer            :: i, j
   real, parameter    :: rdg = 180.0 / pi_
   real, parameter    :: undef = 0.
!
   real, allocatable  :: rlat(:,:),rlon(:,:)
   real, allocatable  :: glon(:),glat(:)
   real               :: glato, glono
   real               :: d00, d01, d10, d11
!
   real               :: rlonmin, rlatmin, rlonmax, rlatmax
   integer            :: is, js, ie, je
!
   integer            :: xlon, xlat, lon, lat
   real               :: x, y
   real               :: x00, y00
   real               :: rlon1i, rlat1i, glon1o,glat1o
!
   real               :: wsum, wsumiv
!
   ijmxin=imxin*jmxin
!
#ifdef DBG
   print*, "post_rmp2latlon_river start"
#endif
!
   allocate (rlat(imxin,jmxin),rlon(imxin,jmxin))
   allocate (glon(imxout),glat(jmxout))
!
!  Initialize
!
   do j = 1,jmxout
     do i = 1,imxout
       llout(i,j) = undef
     enddo
   enddo
!
!  compute latitude and longitude of global
!
   do j = 1,jmxout-1
     glat(j) = (delyo*j) - 90.
#ifdef DBG
     if ( j == jmxout-1 ) print*, "j, delyo, glat", j, delyo, glat(j)
#endif
   enddo
   glat1o=glat(1)
!
   do i = 1,imxout
     glon(i) = delxo * i
#ifdef DBG
     if ( i == imxout ) print*, "i, delxo, glon", i, delxo, glon(i)
#endif
   enddo
   glon1o=glon(1)
!
!  compute latitude and longitude of rmp
!
   call rmp_setup_rsm_grid(rlat,rlon,rdelx,rdely,dlamda0)
   call rmp_trans2output_grid(rlat,1)
   call rmp_trans2output_grid(rlon,1)
!
   rlonmin =  999.
   rlonmax = -999.
   rlatmin =  999.
   rlatmax = -999.
!
! radian to degree
!
   do j = 1,jmxin
     do i = 1,imxin
       rlat(i,j) = rlat(i,j) * 180. / pi_
       rlon(i,j) = rlon(i,j) * 180. / pi_
     enddo
   enddo
!
   do j = 1,jmxin
     do i = 1,imxin
       if ( rlonmin .gt. rlon(i,j) ) rlonmin = rlon(i,j)
       if ( rlonmax .lt. rlon(i,j) ) rlonmax = rlon(i,j)
       if ( rlatmin .gt. rlat(i,j) ) rlatmin = rlat(i,j)
       if ( rlatmax .lt. rlat(i,j) ) rlatmax = rlat(i,j)
       if ( i.eq.1 .and. j.eq.1 ) then
         rlon1i=rlon(i,j)
         rlat1i=rlat(i,j)
       endif
     enddo
   enddo
!
!  Find start and end grid of do loop on lat-lon coordinate;
!  this eliminates needless do loop calculations at out of rmp domain.
! 
   is = (rlonmin-mod(rlonmin,delxo))/delxo
   ie = (rlonmax-mod(rlonmax,delxo)+delxo)/delxo
!
   js = ((rlatmin-mod(rlatmin,delyo))+90.)/delyo
   je = ((rlatmax-mod(rlatmax,delyo))+delyo+90.)/delyo-2
!
#ifdef DBG
   print*, "is, ie, js, je = ", is, ie, js, je
   print*, "glon1, glon2, glat1, glon2 =", glon(is), glon(ie), glat(js), glat(je)
#endif
!
! calculate weighting factor
! 
   call sfc_interp_ll2xyr(1,rproj,rorient,rtruth,rcotru,rlat1i,rlon1i,x00,y00)
#ifdef DBG
   print *, ' input grid rlat1 rlon1 x00 y00 ',rlat1i,rlon1i,x00,y00
#endif
!
   do j = js, je
     glato=glat(j)
     do i = is, ie
       glono=glon(i)
       if(glono.lt.0.) glono = glono + 360.
!
! find nearest rmp point
!
       call sfc_interp_ll2xyr(1,rproj,rorient,rtruth,rcotru,glato,glono,x,y)
       xlon=(x-x00)/rdelx+1
       xlat=(y-y00)/rdely+1
#ifdef DBG
!       print*, "x, y = ", x, y
!       print*, "ii, jj = ", ii,jj
#endif
       if ( xlon.le.1 .or. xlon.gt.imxin-1 .or.                                &
            xlat.le.1 .or. xlat.gt.jmxin        ) then
!
! This point lies out of rmp domain.
!
#ifdef DBG
!         print*, "(",ii,jj,"point is not proper"
!         print*, "glato,glono,rlat,rlon:",glato,glono,rlat(ii,jj),rlon(ii,jj)
#endif
         llout(i,j)=undef
       else
!
! This point is included in rmp domain.
!
#ifdef DBG
!         print*, "(",ii,jj,"point is OK"
!         print*, "xlat,xlon,glato,glono,rlat,rlon:",                           &
!                  xlat,xlon,glato,glono,rlat(xlon,xlat),rlon(xlon,xlat)
#endif
!
! calculate weighting coefficient
!
         lon=max(xlon,1)
         lat=max(xlat,1)
         lon=min(lon,imxin)
         lat=min(lat,jmxin)
         d00=(1.-(xlon-lon)) * (1.-(xlat-lat))
         d10=(xlon-lon) * (1.-(xlat-lat))
         d11=(xlon-lon) * (xlat-lat)
         d01=(1.-(xlon-lon)) * (xlat-lat)
         wsum = d00 + d10 + d11 + d01
         wsumiv = 1./wsum
!
! interpolates to output grid
!
         llout(i,j) = (  d00*regin(xlon,xlat-1)+d10*regin(xlon+1,xlat-1)       &
                       + d11*regin(xlon+1,xlat)+d01*regin(xlon,xlat) ) *wsumiv
       endif
     enddo
   enddo
!
! adjust to TRIP grid format
!
   call ydir_reverse_2d(imxout,jmxout,llout)
!
#ifdef DBG
   print*, "post_rmp2latlon_river end"
#endif
#endif
#endif
   return
   end subroutine post_rmp2latlon_river
!
!-------------------------------------------------------------------------------
   subroutine ydir_reverse_2d(io,jo,dat)
!-------------------------------------------------------------------------------
#ifdef RMP
#ifdef RIVER
!
   integer, intent(in)    :: io, jo
   real, dimension(io,jo) :: dat, rdat
   integer                :: i, j, jrev
!
   do j = 1,jo
     jrev=jo-j+1
     do i = 1,io
       rdat(i,j)=dat(i,jrev)
     enddo
   enddo
   do j = 1,jo
     do i = 1,io
       dat(i,j)=rdat(i,j)
     enddo
   enddo
!
#endif
#endif
   return
   end subroutine ydir_reverse_2d
!
!-------------------------------------------------------------------------------
   subroutine post_sigma2p_avg(im,ix,km,si,sl,                                 &
                     ps,us,vs,ts,qs,qss,                                       &
                     kt,pt,upt,vpt,tpt,qpt,rpt)
!-------------------------------------------------------------------------------
!
! subprogram:    post_sigma2p_avg      sigma to pressure thickness
!
! abstract: interpolates winds, temperature and humidity
!   from the sigma coordinate system to constant pressure thicknesses
!   above the ground.
!
! program history log:
!   1994-07-08  iredell       development
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call post_sigma2p_avg(im,ix,km,si,sl,
!    &                  ps,us,vs,ts,qs,qss,
!    &                  kt,pt,upt,vpt,tpt,qpt,rpt)
!
!   input argument list:
!     im       - integer number of points
!     ix       - integer first dimension of upper air data
!     km       - integer number of sigma levels
!     si       - real (km+1) sigma interface values
!     sl       - real (km) sigma values
!     ps       - real (im) surface pressure in kpa
!     us       - real (ix,km) zonal wind in m/s
!     vs       - real (ix,km) merid wind in m/s
!     ts       - real (ix,km) temperature in k
!     qs       - real (ix,km) specific humidity in kg/kg
!     qss      - real (im,km) saturated specific humidity in kg/kg
!     kt       - integer number of pressure thickness layers
!     pt       - real pressure thickness in kpa
!
!   output argument list:
!     upt      - real (im,kt) zonal wind in m/s
!     vpt      - real (im,kt) merid wind in m/s
!     tpt      - real (im,kt) temperature in k
!     qpt      - real (im,kt) specific humidity in kg/kg
!     rpt      - real (im,kt) relative humidity in percent
!
!-------------------------------------------------------------------------------
   integer              ::  km,im,ix,kt
#ifdef HYBRID
   real                 ::  si(im,km+1),sl(im,km),ps(im)
#else
   real                 ::  si(km+1),sl(km),ps(im)
#endif /* HYBRID */
   real                 ::  us(ix,km),vs(ix,km),ts(ix,km),qs(ix,km),qss(im,km)
   real                 ::  upt(im,kt),vpt(im,kt),tpt(im,kt),qpt(im,kt),rpt(im,kt)
!-------------------------------------------------------------------------------
   rmin=1.e+33
   do i = 1,im
     rmin=min(rmin,ps(i))
   enddo
   rmin=rmin/pt
   rmax=-1.e+33
   do i = 1,im
     rmax=max(rmax,ps(i))
   enddo
   rmax=rmax/pt
!
   do k = 1,kt
     do i = 1,im
       upt(i,k)=0.
       vpt(i,k)=0.
       tpt(i,k)=0.
       qpt(i,k)=0.
       rpt(i,k)=0.
     enddo
     do ks = 1,km
       do i = 1,im
#ifdef HYBRID
         if(k-(1-si(i,ks))*rmin.gt.0..and.k-(1-si(i,ks+1))*rmax.lt.1.) then
           r=ps(i)/pt
           fks=min(k-(1-si(i,ks))*r,1.)-max(k-(1-si(i,ks+1))*r,0.)
           if(fks.gt.0.) then
             upt(i,k)=upt(i,k)+fks*us(i,ks)
             vpt(i,k)=vpt(i,k)+fks*vs(i,ks)
             tpt(i,k)=tpt(i,k)+fks*ts(i,ks)
             qpt(i,k)=qpt(i,k)+fks*qs(i,ks)
             rpt(i,k)=rpt(i,k)+fks*qss(i,ks)
           endif
         endif
#else
         if(k-(1-si(ks))*rmin.gt.0..and.k-(1-si(ks+1))*rmax.lt.1.) then
           r=ps(i)/pt
           fks=min(k-(1-si(ks))*r,1.)-max(k-(1-si(ks+1))*r,0.)
           if(fks.gt.0.) then
             upt(i,k)=upt(i,k)+fks*us(i,ks)
             vpt(i,k)=vpt(i,k)+fks*vs(i,ks)
             tpt(i,k)=tpt(i,k)+fks*ts(i,ks)
             qpt(i,k)=qpt(i,k)+fks*qs(i,ks)
             rpt(i,k)=rpt(i,k)+fks*qss(i,ks)
           endif
         endif
#endif /* HYBRID */
       enddo
     enddo
     do i = 1,im
       rpt(i,k)=min(max(qpt(i,k)/rpt(i,k),0.),1.)*100.
     enddo
   enddo
!
   return
   end subroutine post_sigma2p_avg
!
!-------------------------------------------------------------------------------
   subroutine post_sigma2pressure(im,ix,km,si,sl,ps,us,vs,os,                  &
                    zs,zi,ts,rs,qs,as,cls,prs,o3s,                             &
                    ko,po,lpcl,lppr,lpo3,up,vp,op,zp,tp,rp,                    &
                    qp,ap,clp,prp,o3p)
!-------------------------------------------------------------------------------
!
! abstract
!  - sigma to pressure interpolation
!  - interpolates winds, omega, height, temperature and humidity
!   from the sigma coordinate system to the mandatory pressure levels.
!   assumes that relative humidity, temperature, geopotential heights,
!   wind components and vertical velocity vary linearly in the vertical
!   with the log of pressure.  underground heights are obtained using
!   the shuell method and underground temperatures are obtained using
!   a constant moist adiabatic lapse rate.  heights above the top sigma
!   level are integrated hydrostatically.  otherwise fields are held
!   constant outside the sigma structure and no extrapolation is done.
!
! program history log:
!   92-10-31  sela,newell,gerrity,ballish,deaven,iredell
!
! usage:    call post_sigma2pressure(im,ix,km,si,sl,
!    &                 ps,us,vs,os,zs,zi,ts,rs,qs,
!    &                 ko,po,up,vp,op,zp,tp,rp)
!
!   input argument list:
!     im       - integer number of points
!     ix       - integer first dimension of upper air data
!     km       - integer number of sigma levels
!     si       - real (km+1) sigma interface values
!     sl       - real (km) sigma values
!     ps       - real (im) surface pressure in kpa
!     us       - real (ix,km) zonal wind in m/s
!     vs       - real (ix,km) merid wind in m/s
!     os       - real (im,km) vertical velocity in pa/s
!     zs       - real (im,km) heights on the full levels in m
!     zi       - real (im,km) heights on the interfaces in m
!     ts       - real (ix,km) temperature in k
!     rs       - real (im,km) relative humidity in percent
!     qs       - real (ix,km) specific humidity in kg/kg
!     ko       - integer number of pressure levels
!     po       - real (ko) mandatory pressures in kpa
!
!   output argument list:
!     up       - real (im,ko) zonal wind in m/s
!     vp       - real (im,ko) merid wind in m/s
!     op       - real (im,ko) vertical velocity in pa/s
!     zp       - real (im,ko) heights in m
!     tp       - real (im,ko) temperature in k
!     rp       - real (im,ko) relative humidity in percent
!
! subprograms called:
!   isrchflt - find first value in an array less than target value
!
!-------------------------------------------------------------------------------
   use paramodel, only : io_,jo_,ko_,levs_
   use constant, only  : g_,rd_,rv_
!-------------------------------------------------------------------------------
   real, parameter  ::  g= g_ ,rd= rd_ ,rv= rv_
   real, parameter  ::  rog=rd/g,fvirt=rv/rd-1.
   real, parameter  ::  gammam=-6.5e-3,zshul=75.,tvshul=290.66
!
   real             ::  si(im,km+1),sl(im,km),ps(im)
   real             ::  us(ix,km),vs(ix,km),os(im,km)
   real             ::  zs(im,km),zi(im,km),ts(ix,km),rs(im,km),qs(ix,km)
   real             ::  as(ix,km),o3s(ix,km),cls(ix,km),prs(ix,km)
   real             ::  po(ko)
!
   real             ::  up(im,ko),vp(im,ko),op(im,ko)
   real             ::  zp(im,ko),tp(im,ko),rp(im,ko)
   real             ::  qp(im,ko),ap(im,ko),o3p(im,ko),clp(im,ko),prp(im,ko)
   real             ::  asi(levs_),asl(levs_),apo(ko_),aps(im)
   logical lpo3, lpcl, lppr
!-------------------------------------------------------------------------------
!
!  compute log pressures for interpolation
!
   do i = 1,im
     do kk = 2,km
       asi(kk)=log(si(im,kk))
     enddo
     do kk = 1,km
       asl(kk)=log(sl(im,kk))
     enddo
   enddo
   do kk = 1,ko
     apo(kk)=log(po(kk))
   enddo
   do ii = 1,im
     aps(ii)=log(ps(ii))
   enddo
   apsmin=aps(1)
   apsmax=aps(1)
   do ii = 1,im
     apsmin=min(apsmin,aps(ii))
     apsmax=max(apsmax,aps(ii))
   enddo
! 
!  determine sigma layers bracketing pressure layer.
!  within sigma structure, interpolate fields linearly in log pressure
!  between bracketing full sigma layers except heights are interpolated
!  between the nearest full sigma layer and the nearest sigma interface
!
   kd1=1
   do k = 1,ko
     kd1=kd1+isrchflt(km-kd1-1,asl(kd1+1),1,apo(k)-apsmin)-1
     kd2=kd1+isrchflt(km-kd1-1,asl(kd1+1),1,apo(k)-apsmax)-1
     do kd = kd1,kd2
       ku=kd+1
       do i = 1,im
         ask=apo(k)-aps(i)
         if(ask.le.asl(kd).and.ask.gt.asl(ku)) then
           wu=(asl(kd)-ask)/(asl(kd)-asl(kd+1))
           wd=1.-wu
           up(i,k)=wu*us(i,ku)+wd*us(i,kd)
           vp(i,k)=wu*vs(i,ku)+wd*vs(i,kd)
           op(i,k)=wu*os(i,ku)+wd*os(i,kd)
           ki=kd+1
           di=asi(ki)-ask
           kl=nint(ki-0.5+sign(0.5,di))
           wl=di/(asi(ki)-asl(kl))
           wi=1.-wl
           zp(i,k)=wi*zi(i,ki)+wl*zs(i,kl)
           tp(i,k)=wu*ts(i,ku)+wd*ts(i,kd)
           rp(i,k)=wu*rs(i,ku)+wd*rs(i,kd)
           qp(i,k)=wu*qs(i,ku)+wd*qs(i,kd)
           ap(i,k)=wu*as(i,ku)+wd*as(i,kd)
           if(lpcl) clp(i,k)=wu*cls(i,ku)+wd*cls(i,kd)
           if(lppr) prp(i,k)=wu*prs(i,ku)+wd*prs(i,kd)
           if(lpo3) o3p(i,k)=wu*o3s(i,ku)+wd*o3s(i,kd)
         endif
       enddo
     enddo
! 
!  interpolate sigma to pressure outside the model domain
!
     do i = 1,im
       ask=apo(k)-aps(i)
!
!  below ground use shuell method to obtain height, constant lapse rate
!  to obtain temperature, and hold other fields constant
!
       if(ask.gt.0.) then
         up(i,k)=us(i,1)
         vp(i,k)=vs(i,1)
         op(i,k)=os(i,1)
         tvsf=ts(i,1)*(1.+fvirt*qs(i,1))-gammam*(zs(i,1)-zi(i,1))
         if(zi(i,1).gt.zshul) then
           tvsl=tvsf-gammam*zi(i,1)
           if(tvsl.gt.tvshul) then
             if(tvsf.gt.tvshul) then
               tvsl=tvshul-5.e-3*(tvsf-tvshul)**2
             else
               tvsl=tvshul
             endif
           endif
           gammas=(tvsf-tvsl)/zi(i,1)
         else
           gammas=0.
         endif
         part=rog*ask
         zp(i,k)=zi(i,1)-tvsf*part/(1.+0.5*gammas*part)
         tp(i,k)=ts(i,1)+gammam*(zp(i,k)-zs(i,1))
         rp(i,k)=rs(i,1)
         qp(i,k)=qs(i,1)
         ap(i,k)=as(i,1)
         if(lpcl) clp(i,k)=cls(i,1)
         if(lppr) prp(i,k)=prs(i,1)
         if(lpo3) o3p(i,k)=o3s(i,1)
!
!  between bottom sigma and ground interpolate height,
!  extrapolate temperature and hold other fields constant
!
       elseif(ask.ge.asl(1)) then
         up(i,k)=us(i,1)
         vp(i,k)=vs(i,1)
         op(i,k)=os(i,1)
         wl=ask/asl(1)
         wi=1.-wl
         zp(i,k)=wi*zi(i,1)+wl*zs(i,1)
         wu=(asl(1)-ask)/(asl(1)-asl(2))
         wd=1.-wu
         tp(i,k)=wu*ts(i,2)+wd*ts(i,1)
         rp(i,k)=rs(i,1)
         qp(i,k)=qs(i,1)
         ap(i,k)=as(i,1)
         if(lpcl) clp(i,k)=cls(i,1)
         if(lppr) prp(i,k)=prs(i,1)
         if(lpo3) o3p(i,k)=o3s(i,1)
!
!  above top sigma integrate height hydrostatically
!  and hold other fields constant
!
       elseif(ask.le.asl(km)) then
         up(i,k)=us(i,km)
         vp(i,k)=vs(i,km)
         op(i,k)=os(i,km)
         tvkm=ts(i,km)*(1.+fvirt*qs(i,km))
         zp(i,k)=zs(i,km)+rog*tvkm*(asl(km)-ask)
         tp(i,k)=ts(i,km)
         rp(i,k)=rs(i,km)
         qp(i,k)=qs(i,km)
         ap(i,k)=as(i,km)
         if(lpcl) clp(i,k)=cls(i,km)
         if(lppr) prp(i,k)=prs(i,km)
         if(lpo3) o3p(i,k)=o3s(i,km)
       endif
     enddo
   enddo
! 
   return
   end subroutine post_sigma2pressure
!
!-------------------------------------------------------------------------------
   subroutine post_sigma2tropopause(im,ix,km,clat,sl,ps,                       &
                                    u,v,t,ptp,utp,vtp,ttp,shtp)
!-------------------------------------------------------------------------------
!
! abstract
!   - sigma to tropopause interpolation
!   - locates the tropopause pressure level and interpolates
!   the winds and temperature and wind shear to the tropopause.
!   the tropopause is identified by the lowest level above 450 mb
!   where the temperature lapse rate -dt/dz becomes less than 2 k/km.
!   the tropopause is not allowed higher than 85 mb.
!   interpolations are done linearly in log of pressure.
!
! program history log:
!   92-10-31  mccalla,iredell
!
! usage:    call post_sigma2tropopause(im,ix,km,clat,sl,
!    &                  ps,u,v,t,
!    &                  ptp,utp,vtp,ttp,shtp)
!
!   input argument list:
!     im       - integer number of points
!     ix       - integer first dimension of upper air data
!     km       - integer number of sigma levels
!     clat     - real (im) latitude in radiance
!     sl       - real (km) sigma values
!     ps       - real (im) surface pressure in kpa
!     u        - real (ix,km) zonal wind in m/s
!     v        - real (ix,km) merid wind in m/s
!     t        - real (ix,km) temperature in k
!
!   output argument list:
!     ptp      - real (im) tropopause pressure in kpa
!     utp      - real (im) tropopause zonal wind in m/s
!     vtp      - real (im) tropopause merid wind in m/s
!     ttp      - real (im) tropopause temperature in k
!     shtp     - real (im) tropopause wind speed shear in (m/s)/m
!
!-------------------------------------------------------------------------------
   use paramodel, only : levs_
   use constant, only  : g_,rd_,rv_
!-------------------------------------------------------------------------------
   real, parameter  ::  g= g_ ,rd= rd_ ,rv= rv_
   real, parameter  ::  rog=rd/g
   real, parameter  ::  ptbotp=450.e-1,ptbote=350.e-1,pttop=85.e-1,gamt=2.e-3
   real             ::  sl(im,km),ps(im),clat(im)
   real             ::  u(ix,km),v(ix,km),t(ix,km)
   real             ::  ptp(im),utp(im),vtp(im),ttp(im),shtp(im)
   real             ::  asl(ix,levs_)
   real             ::  stp(im)
!-------------------------------------------------------------------------------
   fgamma(k)=(t(i,k-1)-t(i,k+1))/(rog*t(i,k)*(asl(i,k-1)-asl(i,k+1)))
!
!  identify tropopause as first layer above ptbot but below pttop
!  where the temperature lapse rate drops below gamt
!   stp is real interpolated sigma layer number of tropopause
!
   do k = 1,km
     do i = 1,im
       asl(i,k)=log(sl(i,k))
     enddo
   enddo
   do i = 1,im
     k=3
     pu=ps(i)*sl(i,k)
     ptbot=ptbote+(ptbotp-ptbote)*clat(i)/acos(0.)
     do while(k.lt.km-1.and.pu.gt.ptbot)
       k=k+1
       pu=ps(i)*sl(i,k)
     enddo
     gamd=fgamma(k-1)
     gamd=max(gamd,gamt)
     gamu=fgamma(k)
     do while(k.lt.km-1.and.pu.gt.pttop.and.gamu.gt.gamt)
       k=k+1
       pu=ps(i)*sl(i,k)
       gamd=gamu
       gamu=fgamma(k)
     enddo
     gamu=min(gamu,gamt)
     stp(i)=k-(gamt-gamu)/(gamd-gamu)
   enddo
!
!  interpolate tropopause pressure, temperature, winds and wind shear
!  tropopause pressure is constrained to be between ptbot and pttop
!
   do i = 1,im
     kd=stp(i)
     ku=kd+1
     wu=stp(i)-kd
     dlp=asl(i,ku)-asl(i,kd)
     ptp(i)=ps(i)*sl(i,kd)*exp(wu*dlp)
     if(ptp(i).gt.ptbot) then
       wu=wu+log(ptbot/ptp(i))/dlp
       ptp(i)=ptbot
     elseif(ptp(i).lt.pttop) then
       wu=wu+log(pttop/ptp(i))/dlp
       ptp(i)=pttop
     endif
     ttp(i)=t(i,kd)+wu*(t(i,ku)-t(i,kd))
     utp(i)=u(i,kd)+wu*(u(i,ku)-u(i,kd))
     vtp(i)=v(i,kd)+wu*(v(i,ku)-v(i,kd))
     spdd=sqrt(u(i,kd)**2+v(i,kd)**2)
     spdu=sqrt(u(i,ku)**2+v(i,ku)**2)
     shtp(i)=(spdu-spdd)/(rog*0.5*(t(i,ku)+t(i,kd))*dlp)
   enddo
!
   return
   end subroutine post_sigma2tropopause
!
!-------------------------------------------------------------------------------
   subroutine post_wave2grid_fcst(fgz,fq,fte,fdi,fvt,frq,                      &
                      ggz,ggzx,ggzy,gq,gqx,gqy,                                &
                      gte,guu,gvv,gdiv,gvot,grq)
!-------------------------------------------------------------------------------
#ifndef RMP
#ifndef DFS
!
! abstract: computes grid point values for post processing
!
! usage:    call post_wave2grid_fcst
!
! variable info.:
!    syn(1, 0*levs_+1, lan)  ze
!    syn(1, 1*levs_+1, lan)  di
!    syn(1, 2*levs_+1, lan)  te
!    syn(1, 3*levs_+1, lan)  rq
!    syn(1, 4*levs_+1, lan)  uln
!    syn(1, 5*levs_+1, lan)  vln
!    syn(1, 6*levs_+1, lan)  dpdphi
!    syn(1, 6*levs_+2, lan)  dpdlam
!    syn(1, 6*levs_+3, lan)  q
!    syn(1, 6*levs_+4, lan)  dgzdphi
!    syn(1, 6*levs_+5, lan)  dgzdlam
!    syn(1, 6*levs_+8, lan)  gz
!
!-------------------------------------------------------------------------------
#ifdef DFS
   use dfsvar, only : iope
#endif
   use paramodel, only : JCAP1S,LNT2S,LNT22S,LONF2S,LATG2S,                    &
                         LEVSS,LEVHS,LCAPS,LCAP22S,levs_,levh_,lnt22_,lnt2_,   &
                         jcap_,jcap1_,latg2_,ncpus_,ntotal_,                   &
                         lonf_,lonf2_
   use constant, only  : g_, rerth_
   use comio
   use comgpln
   use comfgrid
   use comfver
#ifdef MP
   use paramodel, only : lln22p_,lonf22_,latg2p_,lonf22p_
   use commpi
#endif
#ifdef REDUCE_GRID
   use comreduce
#endif
!-------------------------------------------------------------------------------
   integer    ::  lots   
   integer    ::  lotst   
   integer    ::  ksz     
   integer    ::  ksd     
   integer    ::  kst     
   integer    ::  ksr     
   integer    ::  ksu     
   integer    ::  kstb    
   integer    ::  ksv     
   integer    ::  kspphi  
   integer    ::  ksgzphi 
   integer    ::  ksplam  
   integer    ::  ksp     
   integer    ::  ksgzlam 
   integer    ::  ksgz    
!
   integer    ::  lotss   
   integer    ::  lotsts  
   integer    ::  kszs    
   integer    ::  ksds    
   integer    ::  ksts    
   integer    ::  ksrs    
   integer    ::  ksus    
   integer    ::  kstbs   
   integer    ::  ksvs    
   integer    ::  kspphis 
   integer    ::  ksgzphis
   integer    ::  ksplams 
   integer    ::  ksps    
   integer    ::  ksgzlams
   integer    ::  ksgzs   
!
#ifdef MP
   real, allocatable    ::  ffa(:,:)
   real, allocatable    ::  syf(:,:,:)
   real, allocatable    ::  grs(:,:,:)
#endif
   real, allocatable    ::  ff (:,:)
!
#ifdef MP
#define NCPUSS latg2_
#else
#define NCPUSS ncpus_
#endif
   real, allocatable    ::  syn(:,:,:)
   real, allocatable    ::  syntop(:,:,:)
#undef NCPUSS
!
   real                 ::  fgz(LNT22S)
   real                 ::  fq (LNT22S)
   real                 ::  fte(LNT22S,levs_)
   real                 ::  fdi(LNT22S,levs_)
   real                 ::  fvt(LNT22S,levs_)
   real                 ::  frq(LNT22S,levs_)
!
   real                 ::  ggz (LONF2S,LATG2S)
   real                 ::  ggzx(LONF2S,LATG2S)
   real                 ::  ggzy(LONF2S,LATG2S)
   real                 ::  gq  (LONF2S,LATG2S)
   real                 ::  gqx (LONF2S,LATG2S)
   real                 ::  gqy (LONF2S,LATG2S)
   real                 ::  gte (LONF2S,LATG2S,levs_)
   real                 ::  guu (LONF2S,LATG2S,levs_)
   real                 ::  gvv (LONF2S,LATG2S,levs_)
   real                 ::  gdiv(LONF2S,LATG2S,levs_)
   real                 ::  gvot(LONF2S,LATG2S,levs_)
   real                 ::  grq (LONF2S,LATG2S,levh_)
!
   real, allocatable    ::  spec(:),specp(:)
!-------------------------------------------------------------------------------
   lots    =5*levs_+levh_+6
   lotst   =2*levs_+2
   ksz     =1
   ksd     =1*levs_+1
   kst     =2*levs_+1
   ksr     =3*levs_+1
   ksu     =3*levs_+levh_+1
   kstb    =3*levs_+levh_+1
   ksv     =4*levs_+levh_+1
   kspphi  =5*levs_+levh_+1
   ksgzphi =5*levs_+levh_+2
   ksplam  =5*levs_+levh_+3
   ksp     =5*levs_+levh_+4
   ksgzlam =5*levs_+levh_+5
   ksgz    =5*levs_+levh_+6
!
   lotss   =5*LEVSS+LEVHS+6
   lotsts  =2*LEVSS+2
   kszs    =1
   ksds    =1*LEVSS+1
   ksts    =2*LEVSS+1
   ksrs    =3*LEVSS+1
   ksus    =3*LEVSS+LEVHS+1
   kstbs   =3*LEVSS+LEVHS+1
   ksvs    =4*LEVSS+LEVHS+1
   kspphis =5*LEVSS+LEVHS+1
   ksgzphis=5*LEVSS+LEVHS+2
   ksplams =5*LEVSS+LEVHS+3
   ksps    =5*LEVSS+LEVHS+4
   ksgzlams=5*LEVSS+LEVHS+5
   ksgzs   =5*LEVSS+LEVHS+6
!
#ifdef MP
   allocate(ffa(lln22p_,lotss))
   allocate(syf(lonf22_ ,lotss,latg2p_))
   allocate(grs(lonf22p_,lots ,latg2p_))
#endif
   allocate(ff (LNT22S,lots))
!
#ifdef MP
#define NCPUSS latg2_
#else
#define NCPUSS ncpus_
#endif
   allocate(syn(LCAP22S,lotss,NCPUSS))
   allocate(syntop(2,JCAP1S,lotsts))
#undef NCPUSS
   allocate (spec(lnt22_))
!
   ga2=g_/(rerth_*rerth_)
#ifdef MP
   call mpsp2f(fgz,lnt22p_,spec,lnt22_,1)
   if(iope) then
   call spcshfli(spec,lnt22_,1,jcap_,lwvdef)
#else
   do n = 1,lnt2_
     spec(n)=fgz(n)
   enddo
#endif
   spec(1)=z00
   spec(2)=0.
   do n = 3,lnt2_
     spec(n)=spec(n)/snnp1(n)/ga2
   enddo
   do i = lnt2_+1,lnt22_
     spec(i)=0.
   enddo
#ifdef MP
   call spcshflo(spec,lnt22_,1,jcap_,lwvdef)
   endif
   allocate (specp(lnt22p_))
   call mpsf2p(spec,lnt22_,specp,lnt22p_,1)
#define SPEC specp
#else
#define SPEC spec
#endif
   do n = 1,LNT22S
     ff(n,ksgz)=SPEC(n)
   enddo
   deallocate(spec)
#ifdef MP
   deallocate(specp)
#endif
#undef SPEC
!
   do n=1,LNT22S
     ff(n,ksp )=fq (n)
   enddo
   do k = 1,levs_
     do n = 1,LNT22S
       ff(n,kst+k-1)=fte(n,k)
       ff(n,ksd+k-1)=fdi(n,k)
       ff(n,ksz+k-1)=fvt(n,k)
     enddo
   enddo
   do k = 1,levh_
     do n = 1,LNT22S
       ff(n,ksr+k-1)=frq(n,k)
     enddo
   enddo
!
#ifdef MP
   llstr=lwvstr(mype)
   llens=lwvlen(mype)
   jstr=latstr(mype)
   jend=latstr(mype)+latlen(mype)-1
   lons2=lonlen(mype)*2
   lats2=latlen(mype)
#else
   llstr=0
   llens=jcap1_
   lons2=lonf2_
   lats2=latg2_
#endif
!
#ifdef MP
   call mpnn2nk(ff (1,ksz  ),lnt22p_,levs_,                                    &
                ffa(1,kszs ),lln22p_,levsp_,3+ntotal_)
   call mpnn2n (ff (1,ksp  ),lnt22p_,                                          &
                ffa(1,ksps ),lln22p_,1)
   call mpnn2n (ff (1,ksgz ),lnt22p_,                                          &
                ffa(1,ksgzs),lln22p_,1)
#define QS ffa(1,ksps)
#define DPDPHIS ffa(1,kspphis)
#define DPDLAMS ffa(1,ksplams)
#define GZS ffa(1,ksgzs)
#define DGZDPHIS ffa(1,ksgzphis)
#define DGZDLAMS ffa(1,ksgzlams)
#define DIS ffa(1,ksds)
#define ZES ffa(1,kszs)
#define ULNS ffa(1,ksus)
#define VLNS ffa(1,ksvs)
#else
#define QS ff(1,ksp)
#define DPDPHIS ff(1,kspphi)
#define DPDLAMS ff(1,ksplam)
#define GZS ff(1,ksgz)
#define DGZDPHIS ff(1,ksgzphi)
#define DGZDLAMS ff(1,ksgzlam)
#define DIS ff(1,ksd)
#define ZES ff(1,ksz)
#define ULNS ff(1,ksu)
#define VLNS ff(1,ksv)
#endif
   call sph_del_log_ps(QS ,DPDPHIS ,                                           &
               syntop(1,1,2*LEVSS+1),DPDLAMS ,                                 &
               llstr,llens,lwvdef)
   call sph_del_log_ps(GZS,DGZDPHIS,                                           &
               syntop(1,1,2*LEVSS+2),DGZDLAMS,                                 &
               llstr,llens,lwvdef)
   call sph_divor2wind(DIS,ZES,ULNS,VLNS,                                      &
               syntop(1,1,1),                                                  &
               syntop(1,1,LEVSS+1),                                            &
               llstr,llens,lwvdef)
#undef QS
#undef DPDPHIS
#undef DPDLAMS
#undef QLAPS
#undef DIS
#undef ZES
#undef ULNS
#undef VLNS
!
#ifndef MP
! compute latitude band limits
!
   last=mod(latg2_,ncpus_)
   nggs=(latg2_-last)/ncpus_
   if(last.ne.0)nggs=nggs+1
   inclat=ncpus_
   lat1=1-ncpus_
   lat2=0
   latdon=0
!
   do ngg = 1,nggs
     if((ngg.eq.nggs).and.(last.ne.0)) inclat=last
     lat1=lat1+ncpus_
     lat2=lat2+inclat
#endif
!
#ifdef MP
     lat1=1
     lat2=latg2_
     latdon=0
#define ZES ffa(1,kszs)
#define TES ffa(1,ksts)
#else
#define ZES ff(1,ksz)
#define TES ff(1,kst)
#endif
!
! first lat loop
!
#ifdef ORIGIN_THREAD
!$doacross share(syntop,syn,qtt,qvv,lat1,lat2,latdon,
!$&        colrad,ze,zea,llstr,llens,lwvdef,lcapdp,mype,lcapd)
!$&        local(lat,lan,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,llensd)
#endif
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
       llensd=lcapdp(lat,mype)
#else
       llensd=lcapd(lat)
#endif
#else
       llensd=llens
#endif
       call sph_sum_coeff(ZES,syn(1,1,lan),qtt(1,lat),                         &
                 llstr,llensd,lwvdef,lotss)
       call sph_sum_coeff_top(syn(1,kstbs,lan),syntop,qvv(1,lat),              &
                 llstr,llensd,lwvdef,lotsts)
     enddo 
!
#ifdef MP
     call mpnl2ny(syn,lcap22p_,latg2_,                                         &
                  syf,lonf22_,latg2p_,lotss,1,lotss)
!
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#define SYNS syf
#else
#define SYNS syn
#endif
!
#ifdef ORIGIN_THREAD
!$doacross share(syf,syn,lat1,lat2,latdon,latdef,lcapd,lonfd),
!$&        local(lat,lan,k,j,lcapf,lonff)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$& shared(syf,syn)
!mic$& shared(lat1,lat2,latdon,latdef,lcapd,lonfd)
!mic$& private(lat,lan,k,j,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,k,j,lcapf,lonff)
#endif
     do lat = lat1,lat2
#ifdef REDUCE_GRID
       lcapf=lcapd(latdef(lat))
       lonff=lonfd(latdef(lat))
#else
       lcapf=jcap1_
       lonff=lonf_
#endif
       lan=lat-latdon
       call sph_wave2grid(SYNS(1,1,lan),SYNS(1,1,lan),lotss*2,                 &
                    lcapf,lonff,latdef(lat),1)
     enddo
!
#undef SYNS
!
#ifdef MP
     call mpnk2nx(syf,lonf22_,lotss,                                           &
                 grs,lonf22p_,lots,latg2p_,levsp_,levs_,1,1,                   &
                 5+ntotal_)
     call mpx2nx (syf,lonf22_,lotss,                                           &
                 grs,lonf22p_,lots,latg2p_,kspphis,kspphi,6)
!
     lat1=jstr
     lat2=jend
     latdon=jstr-1
#define SYNS grs
#define LATX lan
#else
#define SYNS syn
#define LATX lat
#endif
!
#ifdef ORIGIN_THREAD
!$doacross share(grs,syn,lat1,lat2,latdon,spdlat,
!$&              rbs2,del,rdel2,ci,p1,p2,h1,h2,tov,lons2,lonfdp,mype),
!$&              local(lat,lan,j,k,lonsd2)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(syn,grs,lat1,lat2,latdon,spdlat)
!mic$1 shared(rbs2,del,rdel2,ci,p1,p2,h1,h2,tov,lons2,lonfdp,mype)
!mic$1 private(lat,lan,j,k,lonsd2)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,j,k,lonsd2)
#endif
!
     do lat = lat1,lat2
       lan=lat-latdon
#ifdef MP
       rcs=sqrt(rbs2(lan))
#else
       rcs=sqrt(rbs2(lat))
#endif
#ifdef REDUCE_GRID
#ifdef MP
       lonsd2=lonfdp(lan,mype)*2
#else
       lonsd2=lonfd(latdef(lat))*2
#endif
#else
       lonsd2=lons2
#endif
!
#ifdef MP
       if( lonsd2.gt.0 ) then
#endif
#ifndef MP
#define LAN lat
#else
#define LAN lan
#endif
         do i = 1,lonsd2
           ggz (i,LAN)=SYNS(i,ksgz   ,lan)
           ggzx(i,LAN)=SYNS(i,ksgzlam,lan)*rcs
           ggzy(i,LAN)=SYNS(i,ksgzphi,lan)*rcs
           gq  (i,LAN)=SYNS(i,ksp    ,lan)
           gqx (i,LAN)=SYNS(i,ksplam ,lan)*rcs
           gqy (i,LAN)=SYNS(i,kspphi ,lan)*rcs
         enddo
         do i = lonsd2+1,LONF2S
           ggz (i,LAN)=ggz (lonsd2,lan)
           ggzx(i,LAN)=ggzx(lonsd2,lan)
           ggzy(i,LAN)=ggzy(lonsd2,lan)
           gq  (i,LAN)=gq  (lonsd2,lan)
           gqx (i,LAN)=gqx (lonsd2,lan)
           gqy (i,LAN)=gqy (lonsd2,lan)
         enddo
         do k = 1,levs_
           do i = 1,lonsd2
             gte (i,LAN,k)=SYNS(i,kst-1+k,lan)
             guu (i,LAN,k)=SYNS(i,ksu-1+k,lan)*rcs
             gvv (i,LAN,k)=SYNS(i,ksv-1+k,lan)*rcs
             gdiv(i,LAN,k)=SYNS(i,ksd-1+k,lan)
             gvot(i,LAN,k)=SYNS(i,ksz-1+k,lan)
           enddo
           do i = lonsd2+1,LONF2S
             gte (i,LAN,k)=gte (lonsd2,lan,k)
             guu (i,LAN,k)=guu (lonsd2,lan,k)
             gvv (i,LAN,k)=gvv (lonsd2,lan,k)
             gdiv(i,LAN,k)=gdiv(lonsd2,lan,k)
             gvot(i,LAN,k)=gvot(lonsd2,lan,k)
           enddo
         enddo
         do k = 1,levh_
           do i = 1,lonsd2
             grq (i,LAN,k)=SYNS(i,ksr-1+k,lan)
           enddo
           do i = lonsd2+1,LONF2S
             grq (i,LAN,k)=grq (lonsd2,lan,k)
           enddo
         enddo
!
#ifdef MP
       endif
#endif
     enddo
#ifndef MP
     latdon=latdon+(lat2-lat1+1)
   enddo
#endif
!
#ifdef MP
   deallocate(ffa, syf, grs)
#endif
   deallocate(ff, syn, syntop)
#undef SYNS
#undef LAN
#undef LATX
#undef ZES
#undef TES
#undef NCPUSS
#undef QS
#undef DPDPHIS
#undef DPDLAMS
#undef GZS
#undef DGZDPHIS
#undef DGZDLAMS
#undef DIS
#undef ZES
#undef ULNS
#undef VLNS
!
#endif  /* not RMP */
#endif  /* not DFS */
   return
   end subroutine post_wave2grid_fcst
!
!------------------------------------------------------------------------------
