#include <define.h>
   module diag_point_module
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [diag_point_module]
!      |
!      |--- [diag_point_prepare] *
!      |         |--- [diag_point_ij] *
!      |--- [diag_point_arrange] *
!
!-------------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine diag_point_prepare(con,colrad,idim,jdim2,n1)
!-------------------------------------------------------------------------------
#ifdef DGP
   use paramodel, only           : lonf_,latg_,latg2_
   use constant, only            : pi_
   use module_trans, only        : dyn_trans2model_grid
   use module_sph_legendre, only : sph_gaussian_lat
   use comgpd
!-------------------------------------------------------------------------------
!
!  this routine computes the data grid point indices igrd,jgrd for
!   fcst grid, for the
!   npoint points,given the longitude and latitude of same (con).
!   colrad is the colatitude of the fcst grid (dimensioned jdim2,which
!   is half of the total latitudinal grid points),
!
!-------------------------------------------------------------------------------
   real                ::  colrad(jdim2),con(1700)
   real                ::  alat(nptken),alon(nptken)
   real                ::  slimsk(lonf_,latg_),blatf(latg_)
   real                ::  colrar(latg2_),wgr(latg2_),wgrcs(latg2_),rrs2(latg2_)
   integer             ::  kpoi(nptken)
!
!....   begin here..
!
   do k = 1,nstken
     do j = 1,nptken
       do i = 1,nvrken
         svdata(i,j,k) = 0.
       enddo
     enddo
   enddo
!
!--- get slmsk because we wish to get nearest point
!       of same sfc type...this array not available til after step1
!       if we were not trying to cover ourselves for out-board radi8
!       ,then this code could be called from step1...
!
   open (unit=n1,file='sfci',format='unformatted',err=999)
!
   go to 998
!
999 continue
   print *,'error opening sfci in diag_point_prepare'
#ifdef MP
#ifdef RMP
   call rmpabort
#else
   call mpabort
#endif
#else
   call abort
#endif
998 continue
   rewind n1
   read(n1)
   read(n1) ghour,id1,id2,id3,id4
#ifndef NOPRINT
99 format(1h ,'fhour, idate=',f6.2,2x,4(1x,i4))
   print *,'in diag_point_prepare read slmsk from unit=',n1
   print 99,ghour, id1,id2,id3,id4
#endif
   read(n1)
   read(n1)
   read(n1)
   read(n1)
   read(n1)
   read(n1)
!
!.....  skip cv, cvb, cvt, albedo
!
   read(n1)
   read(n1)
   read(n1)
   read(n1)
   read(n1) slimsk
   rewind n1
!
!  call dyn_trans2model_grid(slimsk,1)
   call sph_gaussian_lat (latg2_, colrar, wgr, wgrcs, rrs2)
   dxf = 360. / lonf_
   dxr = 360. / lonf_
   ilonf = lonf_
   jlatg2 = latg2_
   jlatg = latg_
   jfp1 = jlatg2 + 1
!----    get latitude of gaussian grids
   do j = 1,jlatg2
     blatf(j) =(pi_ /2. - colrad(j)) * 180. / pi_
   enddo
   blatf(jfp1) = -blatf(jlatg2)
!
!...    put lat/lon into useable arrays (max=200),where
!         npoint gt 0 implies npoint lat/lon s in con and
!           if abs(lat) between   0, 90 look for nearest point
!                       between 100,190 look for nearest land point
!                       between 200,290 look for nearest sea point
!         npoint lt 0 implies lat/lon of center of region ,only..
!                    lat,lon=con(1301),con(1501)
!           let xy=abs(npoint) and always be 2 digits
!                              and do not differentiate land/sea,
!             then x between 1,9 means create array of every x points
!                   (i.e. x=1 means every point,x=3 means every 3rd,..
!              and y between 0,9 means create (y+1,y+1) array..
!           thus xy can have values 10-99
!....
!
   npute = -1
   if (npoint.lt.0) then
     xy = abs(npoint)
     if (xy.lt.10..or.xy.gt.99.) then
       npute = 0
#ifndef NOPRINT
       print 98,npoint
98     format(1h ,' num(1300)=',i6,'out of -range, so set=1')
#endif
       npoint = 1
     else
       npoint = 1
       iskp = xy/10
       iy   = xy - iskp*10 + 1
#ifndef NOPRINT
       print 97,iy,iy,iskp
97     format(1h ,' prepare regional (',i2,',',i2,') array - every',           &
              i2,' points')
#endif
       npute = iy * iy
     endif
   endif
!
   do 5 k = 1,npoint
     ils = -1
     ylat = abs(con(k+1300))
     if (ylat.ge.100.and.ylat.le.190.) then
!
!...   land point is desired...
!
       ils = 1
       sgn = con(k+1300) / ylat
       con(k+1300) = ylat-100.
       if (sgn.lt.0.) con(k+1300) = - (ylat-100.)
     endif
     if (ylat.ge.200.and.ylat.le.290.) then
       ils = 0
       sgn = con(k+1300) / ylat
       con(k+1300) = ylat-200.
       if (sgn.lt.0.) con(k+1300) = - (ylat-200.)
     endif
     xlat = con(k+1300)
     xlon = con(k+1500)
#ifndef NOPRINT
     if (npute.lt.0.and.ils.eq.-1) print 197,k,xlat,xlon
     if (npute.lt.0.and.ils.eq.0) print 198,k,xlat,xlon
     if (npute.lt.0.and.ils.eq.1) print 199,k,xlat,xlon
197  format(1h ,' ==== station ',i4,' at latlon=',2f8.2,                       &
               ' desired as nearest point')
198  format(1h ,' ==== station ',i4,' at latlon=',2f8.2,                       &
               ' desired as ocean pt')
199  format(1h ,' ==== station ',i4,' at latlon=',2f8.2,                       &
               ' desired as land pt')
#endif
     if (npute.gt.0.and.k.gt.1) go to 195
     alat(k) = con(k+1300)
     alon(k) = con(k+1500)
195  continue
     if (xlon.lt.0) xlon = 360. + con(k+1500)
     if (npute.lt.0) then
       ils = -1
       call diag_point_ij (xlat,xlon,slimsk,blatf,dxf,                         &
                   ils,ilonf,jlatg,ki,kj)
     else
       ils = -1
       call diag_point_ij (xlat,xlon,slimsk,blatf,dxf,                         &
                   ils,ilonf,jlatg,ki,kj)
     endif
     igrd(k) = ki
     jgrd(k) = kj
     if(npute.gt.0) go to 5
     if(xlat.lt.0.) then
       igrd(k) = ki + ilonf
       jgrd(k) = jlatg + 1 - kj
     endif
5  continue
!
!....    regional block , i,j still in single latitude structure..
!
   if (npute.gt.0) then
     iback = iy/2
!
!....    if iy = 1 the all we want is 1 point
!
     if (iback.le.0) then
       npoint = 1
!
       go to 59
!
     endif
     istarf = igrd(1) - iback*iskp
     jstarf = jgrd(1) - iback*iskp
     npoint = 0
     do kyj = 1,iy
       do kxi = 1,iy
         npoint = npoint + 1
         igrd(npoint) = istarf + (kxi-1)*iskp
         jgrd(npoint) = jstarf + (kyj-1)*iskp
       enddo
     enddo
     do 32 n = 1,npoint
       kpoi(n) = 0
       if (jgrd(n).gt.jlatg.or.jgrd(n).lt.1) go to 32
       if (igrd(n).gt.ilonf) igrd(n) = igrd(n) - ilonf
       if (igrd(n).lt.1) igrd(n) = igrd(n) + ilonf
       kpoi(n) = n
32   continue
!
!...    squeeze out the out of bounds points(kpoi=0)
!
     npp = 0
     do 33 n = 1,npoint
       if (kpoi(n).le.0) go to 33
       npp = npp + 1
       igrd(npp) = igrd(kpoi(n))
       jgrd(npp) = jgrd(kpoi(n))
       if (jgrd(npp).gt.jlatg2) then
         igrd(npp) = igrd(npp) + ilonf
         jgrd(npp) = jlatg+1-jgrd(npp)
       endif
33   continue
     npoint = npp
   endif
!
!...................  debug print
!
59 continue
   do k = 1,npoint
     ig=igrd(k)
     jg=jgrd(k)
     iclnd=ig
     jclnd=jg
     if(igrd(k).le.ilonf) then
       blat=90.-colrad(jgrd(k))*180./pi_
       blon=(igrd(k)-1)*360./ilonf
       if(blon.gt.180.) blon=blon-360.
     else
       blat=colrad(jgrd(k))*180./pi_-90.
       blon=(igrd(k)-1-ilonf)*360./ilonf
       if(blon.gt.180.) blon=blon-360.
       iclnd=ig-ilonf
       jclnd=jlatg+1-jg
     endif
#ifndef NOPRINT
     write(6,61) k,alat(k),alon(k),blat,blon
     write(6,62) jgrd(k),igrd(k),slimsk(iclnd,jclnd)
   enddo
61 format(' diag_point_prepare: k,orig lat-lon,compt lat-lon=',i4,4f8.2)
62 format('          ....jgrd,igrd,slmsk=',2i6,f6.1)
#endif
#endif  /* DGP */
!
   return
   end subroutine diag_point_prepare
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine diag_point_ij (xlat,xlon,slmsk,blat,dx,                          &
                     ils,idm,jdm,ki,kj)                                      
!-------------------------------------------------------------------------------
#ifdef DGP
   use paramodel, only : LATGS,lonfp_,lonf_
#ifdef MP
#define LONFS lonfp_
#else
#define LONFS lonf_
#endif
!-------------------------------------------------------------------------------
   real                 ::  blat(jdm),dist(4),keni(4),kenj(4)
   real                 ::  slmsk(LONFS,LATGS)
   integer              ::  ipsort(4),ipskp(4)
!
   jdm2 = jdm / 2
   jdmp1 = jdm2 + 1
   if (abs(xlat).gt.blat(1)) go to 70
!
!----    get upper left gaussian point (ia,ja) on gridbox
!          surrounding the input lat/lon point........
   ia = xlon/dx + 1
   ib = ia + 1
   if (ia.ge.idm) ib = 1
   xi = xlon/dx + 1. - ia
   do 10 jak=2,jdmp1
     jb = jak - 1
     if(abs(xlat).gt.blat(jak)) go to 15
10 continue
15 continue
!
   ja = jb
!
!----   normalized distance from upper lat to gaussian lat                      
!
   xj = (blat(ja) - abs(xlat)) / (blat(ja)-blat(ja+1))
!       xout(i,lat)   = (1-xi)* xj   *xin(ia,ja+1) +
!    1                 xi *  xj  *xin(ib,ja+1) +
!    2              (1-xi)*(1-xj)*xin(ia  ,ja  ) +
!    3               xi   *(1-xj)*xin(ib,ja  )
!----    southern hemisphere
!
   if (xlat.lt.0.) then
     ja = jdm - ja
     xj = 1. - xj
   endif
!
!       xout(i,jout+1-lat)=(1-xi)* xj   *xin(ia,ja+1) +
!    1                 xi *  xj  *xin(ib,ja+1) +
!    2              (1-xi)*(1-xj)*xin(ia  ,ja  ) +
!    3               xi   *(1-xj)*xin(ib,ja  )
!...     upper left point
!
   dist(1)= sqrt(xi**2+xj**2)
   keni(1)= ia
   kenj(1)= ja
!
!...     upper right
!
   dist(2)= sqrt((1-xi)**2+xj**2)
   keni(2)= ib
   kenj(2)= ja
!
!...     lower right
!
   dist(3)= sqrt((1-xi)**2+(1-xj)**2)
   keni(3)= ib
   kenj(3)= ja+1
!
!...     lower left
!
   dist(4)= sqrt(xi**2+(1-xj)**2)
   keni(4)= ia
   kenj(4)= ja+1
!
!---     now sort the distances (by index, shortest first)
!
   npt = 4
   npt1 = npt + 1
   do kd = 1,npt
     ipskp(kd) = 0
   enddo
!
   do kd = 1,npt
     dd=100.
!
!---     find shortest dist of the remaining unsorted data...
!
     do 25 kk = 1,npt
       if(ipskp (kk).gt.0) go to 25
       if(dist(kk).lt.dd) then
         dd = dist(kk)
         jx = kk
       endif
25   continue
!
!---   store sorted index
!
     ipsort(kd) = jx
     ipskp (jx) = 1
   enddo
#ifndef NOPRINT
   print 102,(dist(kk),kk=1,npt),(ipsort(kk),kk=1,npt)
102 format(1h ,' distances=',4f8.4,' sorted indices=',4i4)
#endif
!
!---   end of distance sort
!
   xils = ils
   if (ils.lt.0) then
     ki = keni(ipsort(1))
     kj = kenj(ipsort(1))
     return
   endif
!
   if (ils.eq.0) then
!
!....     find nearest sea point                                                
!
     do kd = 1,npt
       ii = keni(ipsort(kd))
       jj = kenj(ipsort(kd))
       if (slmsk(ii,jj).le.xils) go to 46
     enddo
!
!....     no sea points so default to nearest point
!
#ifndef NOPRINT
     print 49,xlat,xlon                                                       
49   format(1h ,' asked for sea point but can t find one for',                 &
             ' lat lon=',2f9.2,'..so default to nearest')
#endif
     ki = keni(ipsort(1))
     kj = kenj(ipsort(1))
     return
46   continue
     ki = ii
     kj = jj
     return
   endif
!
   if (ils.eq.1) then
!
!....     find nearest land/ice point
!
     do kd = 1,npt
       ii = keni(ipsort(kd))
       jj = kenj(ipsort(kd))
       if (slmsk(ii,jj).ge.xils) go to 56
     enddo
!
!....     no land points so default to nearest point
!
#ifndef NOPRINT
     print 59,xlat,xlon
59   format(1h ,' asked for land sea point but can t find one for',            &
               ' lat lon=',2f9.2,'..so default to nearest')
#endif
     ki = keni(ipsort(1))
     kj = kenj(ipsort(1))
     return
56   continue
     ki = ii
     kj = jj
     return
   endif
!
!...   outside limit of gaussian polar rows so just take nearest
!        point without regard to land and sea
!
70 continue
   ia = xlon/dx + 1
   ib = ia + 1
   if (ia.ge.idm) ib = 1
   xi = xlon/dx + 1. - ia
   ja = 1
   if (xlat.lt.0.) ja = jdm
   if (xi.gt.0.5) then
     ki = ib
     kj = ja
   else
     ki = ia
     kj = ja
   endif
#endif
!
#undef LONFS
!
   return
   end subroutine diag_point_ij
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine diag_point_arrange(lat,rcl,                                      &
                     slmsk,psexp,tg3,snoweq,radsl,dlwsf1,                      &
                     tsea,qss,plantr,gflx,z0cm,cd,cdq,                         &
                     rnet,hflx,stsoil,                                         &
                     canopy,drain,smsoil,runof,cld1d,                          &
                     u10,v10,t2,q2,                                            &
                     hpbl,gamt,gamq,                                           &
                     dqsfc1,dtsfc1,dusfc1,dvsfc1,                              &
                     dusfcg,dvsfcg,                                            &
                     rainc,rainl,                                              &
                     u,v,t,q,ntotal,hsw,hlw,vvel,                              &
#ifdef EXPLICIT_CLOUDINESS
                     qcicps,qrscps,                                            &
#endif
                     snowmt,snowev,snowfl)                                   
!-------------------------------------------------------------------------------
   use paramodel, only : ix=>LONF2S, im=>LONF2S,km=>levs_,ngases_,slvark_
   use constant, only  : cal_
   use comgpd
!-------------------------------------------------------------------------------
   real                 ::  slmsk(im),psexp(im),tg3(im),snoweq(im),radsl(im),  &
                            dlwsf1(im),tsea(im),qss(im),plantr(im),gflx(im),   &
                            z0cm(im),cd(im),cdq(im),                           &
                            rnet(im),hflx(im),stsoil(im,2),                    &
                            canopy(im),drain(im),smsoil(im,2),runof(im),       &
                            cld1d(im),                                         &
                            u10(im),v10(im),t2(im),q2(im),                     &
                            hpbl(im),gamt(im),gamq(im),                        &
                            dqsfc1(im),dtsfc1(im),dusfc1(im),dvsfc1(im),       &
                            dusfcg(im),dvsfcg(im),                             &
                            rainc(im),rainl(im),                               &
                            u(ix,km),v(ix,km),t(ix,km),q(ix,km,ntotal),        &
                            hsw(im,km),hlw(im,km),vvel(im,km),                 &
#ifdef EXPLICIT_CLOUDINESS
                            qcicps(im,km),qrscps(im,km),                       &
#endif
                            snowmt(im),snowev(im),snowfl(im)
   real,parameter       ::  cnwatt=-cal_*1.e4/60.
!-------------------------------------------------------------------------------
#ifdef DGP
   do igpt = 1,npoint
     if(lat.eq.jgrd(igpt)) then
       svdata(11,igpt,itnum)=svdata(11,igpt,itnum)+rainc(igrd(igpt))
       svdata(12,igpt,itnum)=svdata(12,igpt,itnum)+rainl(igrd(igpt))
       if(isave.ne.0) then
         svdata(  1,igpt,itnum)= igrd(igpt)
         svdata(  2,igpt,itnum)= jgrd(igpt)
         svdata(  3,igpt,itnum)= slmsk (igrd(igpt))
         svdata(  4,igpt,itnum)= psexp (igrd(igpt)) *10.
         svdata(  8,igpt,itnum)= tg3 (igrd(igpt))
         svdata( 10,igpt,itnum)= snoweq(igrd(igpt))
         svdata( 13,igpt,itnum)= radsl(igrd(igpt))*cnwatt
         svdata( 14,igpt,itnum)= dlwsf1(igrd(igpt))
         svdata(  5,igpt,itnum)= tsea (igrd(igpt))
         svdata( 15,igpt,itnum)= qss  (igrd(igpt))
         svdata( 16,igpt,itnum)= plantr(igrd(igpt))
         svdata( 19,igpt,itnum)= gflx(igrd(igpt))
         svdata( 22,igpt,itnum)= z0cm (igrd(igpt))
         svdata( 23,igpt,itnum)= cd  (igrd(igpt))
         svdata( 24,igpt,itnum)= cdq (igrd(igpt))
         svdata( 62,igpt,itnum)= rnet (igrd(igpt))
#endif
!        svdata( 63,igpt,itnum)= evap (igrd(igpt))
#ifdef DGP
         svdata( 64,igpt,itnum)= hflx (igrd(igpt))
         svdata(  6,igpt,itnum)= stsoil (igrd(igpt),1)
         svdata(  7,igpt,itnum)= stsoil (igrd(igpt),2)
         svdata( 34,igpt,itnum)= u10 (igrd(igpt))
         svdata( 35,igpt,itnum)= v10 (igrd(igpt))
         svdata( 30,igpt,itnum)= t2  (igrd(igpt))
         svdata( 31,igpt,itnum)= q2  (igrd(igpt))
         svdata( 32,igpt,itnum)= canopy(igrd(igpt))
         svdata( 33,igpt,itnum)= drain(igrd(igpt))
         svdata( 17,igpt,itnum)= dqsfc1(igrd(igpt))
         svdata( 18,igpt,itnum)= dtsfc1(igrd(igpt))
         svdata( 20,igpt,itnum)= dusfc1(igrd(igpt))
         svdata( 21,igpt,itnum)= dvsfc1(igrd(igpt))
         svdata( 28,igpt,itnum)= dusfcg(igrd(igpt))
         svdata( 29,igpt,itnum)= dvsfcg(igrd(igpt))
         svdata(  9,igpt,itnum)= smsoil (igrd(igpt),1)
         svdata( 61,igpt,itnum)= smsoil (igrd(igpt),2)
#endif
!        svdata( 65,igpt,itnum)= runof (igrd(igpt))
#ifdef DGP
         svdata( 63,igpt,itnum)= hpbl(igrd(igpt))
         svdata( 64,igpt,itnum)= gamt(igrd(igpt))
         svdata( 65,igpt,itnum)= gamq(igrd(igpt))
         svdata( 66,igpt,itnum)= cld1d(igrd(igpt))
         svdata( 67,igpt,itnum)= snowev(igrd(igpt))
         svdata( 68,igpt,itnum)= snowmt(igrd(igpt))
         svdata( 69,igpt,itnum)= snowfl(igrd(igpt))
         if(ilshrt.lt.2) then
           r=sqrt(rcl)
           do k = 1,levs_
             svdata(k+slvark_+0*levs_,igpt,itnum)=u(igrd(igpt),k)*r
             svdata(k+slvark_+1*levs_,igpt,itnum)=v(igrd(igpt),k)*r
             svdata(k+slvark_+2*levs_,igpt,itnum)=t(igrd(igpt),k)
             svdata(k+slvark_+3*levs_,igpt,itnum)=q(igrd(igpt),k)
             if(ilshrt.lt.1) then
               svdata(k+slvark_+4*levs_,igpt,itnum)=hsw(igrd(igpt),k)
               svdata(k+slvark_+5*levs_,igpt,itnum)=hlw(igrd(igpt),k)
               svdata(k+slvark_+6*levs_,igpt,itnum)=vvel(igrd(igpt),k)
#ifdef EXPLICIT_CLOUDINESS
               svdata(k+slvark_+9*levs_,igpt,itnum)=qcicps(igrd(igpt),k)
               svdata(k+slvark_+10*levs_,igpt,itnum)=qrscps(igrd(igpt),k)
#endif
               if(nwater_.gt.1) then
                 do icloud = 2,nwater_
                   ic = icloud
                   ivar = 10 + icloud
                   svdata(k+slvark_+ivar*levs_,igpt,itnum)=q(igrd(igpt),k,ic)
                 enddo
               endif
             endif
           enddo
#endif
#ifdef DGP
         endif
       endif
     endif
   enddo
#endif
!
   return
   end subroutine diag_point_arrange
!-------------------------------------------------------------------------------
   end module diag_point_module
