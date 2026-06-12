#include <define.h>
   subroutine sfc_interp_grid
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [sfc_interp_grid]
!           |
!           |-- [sfc_interp_file_gmp] *
!           |-- [sfc_interp_ll2gmp] *
!           |-- [sfc_interp_ll2rmp] *
!           |-- [sfc_interp_ll2mask] *
!           |-- [sfc_interp_ll2xyr] *
!           |-- [sfc_interp_xy2llr] *
!
!-------------------------------------------------------------------------------
   end subroutine sfc_interp_grid
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_interp_file_gmp(regin,imxin,jmxin,dloin,dlain,rlon,rlat,     &
                    gauout,imxout,jmxout,lmask,rslmsk,slmask,igau,             &
                    inttyp)
#include "abort.h"
!-------------------------------------------------------------------------------
! 
! subprogram: sfc_interp_file_gmp
!
! abstract: interpolation from lat/lon or gaussian grid to other gaussian grid
!
! igau=0  lat/lon input
! igau=1  gaussian input
!
! inttyp=0  regular interpolation
! inttyp=1  take the closest point value
! inttyp=2  take the dominant type
!
!-------------------------------------------------------------------------------
   save
!-------------------------------------------------------------------------------
   integer,parameter    ::  numtype=20
   integer              ::  ntype(numtype)
!
   real                 ::  regin (imxin ,jmxin )
   real                 ::  gauout(imxout,jmxout)
!
   real                 ::  rslmsk(imxin,jmxin)
   real                 ::  slmask(imxout,jmxout)
!
   real                 ::  gaulo (jmxout)
   real                 ::  gauli (jmxin)
   real                 ::  rinlat(jmxin),outlat(jmxout)
   real                 ::  rinlon(imxin)
!
   integer              ::  iindx1(imxout)
   integer              ::  iindx2(imxout)
   integer              ::  jindx1(jmxout)
   integer              ::  jindx2(jmxout)
!
   real                 ::  ddx(imxout)
   real                 ::  ddy(jmxout)
!
   logical              ::  lmask
!
!
   print *,'inttyp,io,jo=',inttyp,imxin,jmxin,imxout,jmxout
!
   if(imxin.eq.1.or.jmxin.eq.1) then
     do j = 1,jmxout
       do i = 1,imxout
          gauout(i,j)=0.
       enddo
     enddo
     return
   endif
!
   if(dloin.eq.0..or.dlain.eq.0.) then
     print *,'dloin or dlain is zero .... check data cards'
     call MPABORT
   endif
!
   ifpi=imxin
   jfpi=jmxin
   ifpo=imxout
   jfpo=jmxout
   rfp1=rlon
   rfp2=rlat
!
   if(igau.eq.0) then
     do j = 1,jmxin
       if(rlat.gt.0.) then
         rinlat(j)=rlat-float(j-1)*dlain
       else
         rinlat(j)=rlat+float(j-1)*dlain
       endif
     enddo
   else
     call gaulat(gauli,jmxin)
     do j = 1,jmxin
       rinlat(j)=90.-gauli(j)
     enddo
     rlat=90.
   endif
!
!    compute gaussian latitude for output grid
!
#ifndef DFS
   call gaulat(gaulo,jmxout)
!
   do j = 1,jmxout
     outlat(j)=90.-gaulo(j)
   enddo
#else
   dlat = 180./jmxout
   do j = 1,jmxout
     outlat(j)=90.-(j-0.5)*dlat
   enddo
#endif
!
   do i = 1,imxin
     rinlon(i)=rlon+float(i-1)*dloin
   enddo
!
!  find i-index for interplation
!
   do i = 1,imxout
     alamd=float(i-1)*360./float(imxout)
     if(rlon.lt.0.) then
       if(alamd.gt.180.) alamd=alamd-360.
     endif
     do ii = 1,imxin
       if(alamd.gt.rinlon(ii)) cycle
       ix=ii
       go to 32
     enddo
     i1=nint(360./dloin)
     i2=1
     go to 34
     !
     32 continue
     !
     if(ix.ge.2) go to 33
     i1=nint(360./dloin)
     i2=1
     go to 34
     !
     33 continue
     !
     i2=ix
     i1=i2-1
     !
     34 continue
     !
     iindx1(i)=i1
     iindx2(i)=i2
     denom=rinlon(i2)-rinlon(i1)
     if(denom.lt.0.) denom=denom+360.
     rnume=alamd-rinlon(i1)
     if(rnume.lt.0.) rnume=rnume+360.
     ddx(i)=rnume/denom
   enddo
!
!  find j-index for interpolation
!
   jq=1
   do j = 1,jmxout
     aphi=outlat(j)
     do jj = 1,jmxin
       jx=jj
       if(rlat.lt.0.) jx=jmxin-jj+1
       if(aphi.lt.rinlat(jx)) cycle
       jq=jx
       go to 42
     enddo
     if(rlat.gt.0.) then
       j1=jmxin
       j2=jmxin
     else
       j1=1
       j2=1
     endif
     go to 44
     !
     42 continue
     !
     if(rlat.gt.0.) then
       if(jq.ge.2) go to 43
       j1=1
       j2=1
     else
       if(jq.lt.jmxin) go to 43
       j1=jmxin
       j2=jmxin
     endif
     go to 44
     !
     43 continue
     !
     if(rlat.gt.0.) then
       j2=jq
       j1=jq-1
     else
       j1=jq
       j2=jq+1
     endif
     !
     44 continue
     !
     jindx1(j)=j1
     jindx2(j)=j2
     if(j2.ne.j1) then
       ddy(j)=(aphi-rinlat(j1))/(rinlat(j2)-rinlat(j1))
     else
       if(j1.eq.1.and.rlat.gt.0..or.j1.eq.jmxin.and.rlat.lt.0.) then
         if(abs(90.-rinlat(j1)).gt.0.001) then
           ddy(j)=(aphi-rinlat(j1))/(90.-rinlat(j1))
         else
           ddy(j)=0.0
         endif
       endif
       if(j1.eq.1.and.rlat.lt.0..or.j1.eq.jmxin.and.rlat.gt.0.) then
         if(abs(-90.-rinlat(j1)).gt.0.001) then
           ddy(j)=(aphi-rinlat(j1))/(-90.-rinlat(j1))
         else
         ddy(j)=0.0
         endif
       endif
     endif
   enddo
!
#ifdef DBG
   write(6,*)'in sfc_interp_file_gmp'
   write(6,*)'iindx1=',iindx1
   write(6,*)'iindx2=',iindx2
   write(6,*)'jindx1=',jindx1
   write(6,*)'jindx2=',jindx2
#endif
!
! increments of index for dominant type search
!
   if(inttyp.eq.2) then
     ipm=imxin/imxout/2
     jpm=jmxin/jmxout/2
     print *,'ipm=',ipm,' jpm=',jpm
   endif
!
   sum1=0.
   sum2=0.
   wei1=0.
   wei2=0.
   do i = 1,imxin
     sum1=sum1+regin(i,1) * rslmsk(i,1)
     sum2=sum2+regin(i,jmxin) * rslmsk(i,jmxin)
     wei1=wei1+rslmsk(i,1)
     wei2=wei2+rslmsk(i,jmxin)
   enddo
!
   if(rlat.gt.0.) then
     if(wei1.gt.0.) then
       sumn=sum1/wei1
     else
       sumn=0.
     endif
     if(wei2.gt.0.) then
       sums=sum2/wei2
     else
       sums=0.
     endif
   else
     if(wei1.gt.0.) then
       sums=sum1/wei1
     else
       sums=0.
     endif
     if(wei2.gt.0.) then
       sumn=sum2/wei2
     else
       sumn=0.
     endif
   endif
!
!  quasi-bilinear interpolation or closest point or dominant type
!
   ifill=0
   do j = 1,jmxout
     y=ddy(j)
     j1=jindx1(j)
     j2=jindx2(j)
     jcl=j1
     if(y.gt.0.5) jcl=j2
     do i = 1,imxout
       x=ddx(i)
       i1=iindx1(i)
       i2=iindx2(i)
       icl=i1
       if(x.gt.0.5) icl=i2
!
       if(inttyp.le.1) then
         if(lmask) then
           if(slmask(i,j).eq.rslmsk(i1,j1).and.                                &
              slmask(i,j).eq.rslmsk(i2,j1).and.                                &
              slmask(i,j).eq.rslmsk(i1,j2).and.                                &
              slmask(i,j).eq.rslmsk(i2,j2)) then
             wi1j1=(1.-x)*(1.-y)
             wi2j1=    x *(1.-y)
             wi1j2=(1.-x)*      y
             wi2j2=    x *      y
           elseif(slmask(i,j).eq.1.) then
             wi1j1=(1.-x)*(1.-y)  *rslmsk(i1,j1)
             wi2j1=    x *(1.-y)  *rslmsk(i2,j1)
             wi1j2=(1.-x)*      y *rslmsk(i1,j2)
             wi2j2=    x *      y *rslmsk(i2,j2)
           elseif(slmask(i,j).eq.0.) then
             wi1j1=(1.-x)*(1.-y)  *(1.-rslmsk(i1,j1))
             wi2j1=    x *(1.-y)  *(1.-rslmsk(i2,j1))
             wi1j2=(1.-x)*      y *(1.-rslmsk(i1,j2))
             wi2j2=    x *      y *(1.-rslmsk(i2,j2))
           endif
         else
           wi1j1=(1.-x)*(1.-y)
           wi2j1=    x *(1.-y)
           wi1j2=(1.-x)*      y
           wi2j2=    x *      y
         endif
!
         if(inttyp.eq.1) then
           if(icl.eq.i1.and.jcl.eq.j1) then
             wi1j1=1.
             wi2j1=0.
             wi1j2=0.
             wi2j2=0.
           elseif(icl.eq.i2.and.jcl.eq.j1) then
             wi1j1=0.
             wi2j1=1.
             wi1j2=0.
             wi2j2=0.
           elseif(icl.eq.i1.and.jcl.eq.j2) then
             wi1j1=0.
             wi2j1=0.
             wi1j2=1.
             wi2j2=0.
           elseif(icl.eq.i2.and.jcl.eq.j2) then
             wi1j1=0.
             wi2j1=0.
             wi1j2=0.
             wi2j2=1.
           endif
         endif
!
         wsum=wi1j1+wi2j1+wi1j2+wi2j2
         if(wsum.ne.0.) then
           wsumiv = 1./wsum
         !
           if(j1.ne.j2) then
             gauout(i,j)=(wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+               &
             wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2))*wsumiv
           else
             if(j1.eq.1.and.rlat.gt.0..or.j1.eq.jmxin.and.rlat.lt.0.) then
               gauout(i,j)=(wi1j1*sumn        +wi2j1*sumn        +             &
                      wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2))*wsumiv
             endif
             if(j1.eq.1.and.rlat.lt.0..or.j1.eq.jmxin.and.rlat.gt.0.) then
               gauout(i,j)=(wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+             &
                                wi1j2*sums        +wi2j2*sums        )*wsumiv
             endif
           endif
         else
           if(.not.lmask) then
             write(6,*) ' sfc_interp_file_gmp called with lmask=.true.',       &
                            ' but rslmsk or slmask bad.'
             call MPABORT
           endif
           ifill=ifill+1
!
           do jx = j1,jmxin
             do ix = i1,imxin
               if((slmask(i,j).eq.1..and.                                      &
                   slmask(i,j).eq.rslmsk(ix,jx)).or.                           &
                  (slmask(i,j).eq.0..and.                                      &
                   slmask(i,j).eq.rslmsk(ix,jx))) then
                 gauout(i,j)=regin(ix,jx)
                 go to 71
               endif
             enddo
             do ix = i1,1,-1
               if((slmask(i,j).eq.1..and.                                      &
                   slmask(i,j).eq.rslmsk(ix,jx)).or.                           &
                  (slmask(i,j).eq.0..and.                                      &
                   slmask(i,j).eq.rslmsk(ix,jx))) then
                 gauout(i,j)=regin(ix,jx)
                 go to 71
               endif
             enddo
           enddo
           do jx = j1,1,-1
             do ix = i1,imxin
               if((slmask(i,j).eq.1..and.                                      &
                   slmask(i,j).eq.rslmsk(ix,jx)).or.                           &
                  (slmask(i,j).eq.0..and.                                      &
                   slmask(i,j).eq.rslmsk(ix,jx))) then
                 gauout(i,j)=regin(ix,jx)
                 go to 71
               endif
             enddo
             do ix = i1,1,-1
               if((slmask(i,j).eq.1..and.                                      &
                   slmask(i,j).eq.rslmsk(ix,jx)).or.                           &
                  (slmask(i,j).eq.0..and.                                      &
                   slmask(i,j).eq.rslmsk(ix,jx))) then
                 gauout(i,j)=regin(ix,jx)
                 go to 71
               endif
             enddo
           enddo
           write(6,*) 'error! no filling value found in sfc_interp_file_gmp'
           call MPABORT
         endif
       !
       !  dominant type
       !
       elseif(inttyp.eq.2) then
         isrs=icl-ipm
         isrf=icl+ipm
         jsrs=max(jcl-jpm,1)
         jsrf=min(jcl+jpm,jmxin)
         do n = 1,numtype
           ntype(n)=0
         enddo
         do jx = jsrs,jsrf
           do ixx = isrs,isrf
             ix=ixx
             if(ix.lt.1) ix=imxin-ixx
             if(ix.gt.imxin) ix=ixx-imxin
             n=regin(ix,jx)+1
             if(n.gt.numtype) then
               print *,'type .gt.numtype.',                                    &
                       ' change numtype in sfc_interp_file_gmp'
               call MPABORT
             endif
             if(n.lt.0) then
               print *,'type cannot be less than zero'
               call MPABORT
             endif
             ntype(n)=ntype(n)+1
           enddo
         enddo
         maxcount=0
         maxtyp=-1
         do n = 1,numtype
           if(ntype(n).gt.maxcount) then
             maxcount=ntype(n)
             maxtyp=n
           endif
         enddo
         gauout(i,j)=maxtyp-1
       endif
       !
       71  continue
       !
     enddo
   enddo
!
   if(ifill.gt.1) then
     write(6,*) ' unable to interpolate.  filled with nearest',                &
                ' point value at ',ifill,' points'
   endif
!
#ifdef DBG
   write(6,*) ' maxmin of model grid arry:'
   call print_maxmin_six(gauout,imxout*jmxout,1,1,1,'gauout')
!
#endif
   return
   end subroutine sfc_interp_file_gmp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_interp_ll2gmp(io,proj,glat,glon,ii,jj,delx,dely,rlat1,rlon1, &
                    rlat,rlon,xlon,ylat)
!-------------------------------------------------------------------------------
   real                 ::  glat(jj),glon(ii)
   real                 ::  rlat(io),rlon(io),xlon(io),ylat(io)
!
   integer, allocatable :: lon(:)
   logical, allocatable :: flag(:)
!
   allocate(lon (io))
   allocate(flag(io))
!
   ibeg = 1
   glon0 = glon(1)
   do ix = 1,ii
     if(glon(ix).le.glon0 ) then
       glon0 = glon(ix)
       ibeg = ix
     endif
   enddo
   iend = ibeg + ii - 1
!
   do i = 1,io
     flag(i)=.true.
   enddo
!
   do ix = ibeg,iend
     if(ix.le.ii) then
       iglon = ix
     else
       iglon = ix - ii
     endif
     gloni = glon(iglon)
     do i = 1,io
       if(flag(i).and.rlon(i).lt.gloni) then
         xlon(i) = iglon + (rlon(i)-gloni)/delx*1000.
         lon(i) = xlon(i)
         flag(i) = .false.
       endif
     enddo
   enddo
!
   do ix = ibeg,iend
     if(iend.le.ii) then
       do i = 1,io
         if((.not.flag(i)).and.lon(i).lt.1) then
           lon(i) = iend
         endif
       enddo
     else
       do i = 1,io
         if((.not.flag(i)).and.lon(i).lt.1) then
           lon(i) = iend - ii
         endif
       enddo
     endif
     do i = 1,io
       if((.not.flag(i)).and.lon(i).lt.1) then
         xlon(i) = iend + xlon(i)
       endif
     enddo
   enddo
!
   do ix = ibeg,iend
     if(iend.le.ii) then
       do i = 1,io
         if(flag(i)) then
           lon(i) = iend
         endif
       enddo
     else
       do i = 1,io
         if(flag(i)) then
           lon(i) = iend-ii
         endif
       enddo
     endif
   enddo
!
   do ix = ibeg,iend
     do i = 1,io
       if(flag(i)) then
         xlon(i) = lon(i) + (rlon(i)-glon(lon(i)))/delx*1000.
       endif
     enddo
   enddo
!
   do i = 1,io
     flag(i)=.true.
   enddo
!
   if(proj.gt.0.) then
     do j = jj,2,-1
       do i = 1,io
         if(flag(i).and.glat(j).le.rlat(i)) then
           ylat(i) = j + (rlat(i)-glat(j))/(glat(j+1)-glat(j))
           flag(i) = .false.
         endif
         if(rlat(i).le.glat(1)) then
           ylat(i) = 1
           flag(i) = .false.
         endif
       enddo
     enddo
   else
     do j = 2,jj
       do i = 1,io
         if(flag(i).and.glat(j).le.rlat(i)) then
           ylat(i) = j-1 + (glat(j-1)-rlat(i))/(glat(j-1)-glat(j))
           flag(i) = .false.
         endif
         if(rlat(i).ge.glat(1)) then
           ylat(i) = 1
           flag(i) = .false.
         endif
       enddo
     enddo
   endif
!
   deallocate (lon,flag)
!
   return
   end subroutine sfc_interp_ll2gmp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_interp_ll2rmp(rdat,imin,jmin,inttyp,rout,io,jo,              &
                     lmask,rslmsk,slmask,                                      &
                     proji,orienti,truthi,cotrui,                              &
                     delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i)
!-------------------------------------------------------------------------------
!
!  interpolation from lat/lon or gaussian grid to other lat/lon grid
!
!  inttyp ... type of interpolation.  0 = bilinear
!                                     1 = take the closest point value
!                                     2 = take the predominant type
!
!-------------------------------------------------------------------------------
   use paramodel, only    :  ILOTS => LONF2F, JLOTS => LATG2F
   use constant, only     :  pi_
#ifdef RMP
   use module_trans, only :  rmp_trans2output_grid
   use rscomloc
#endif
!------------------------------------------------------------------------------- 
   implicit none
!-------------------------------------------------------------------------------
#undef VECTOR
#ifdef NEC
#define VECTOR
#endif
#ifdef ES
#define VECTOR
#endif
#define VECTOR
#include "abort.h"
   integer                 ::  imin,jmin,io,jo,inttyp
   integer                 ::  idim,jdim
   real                    ::  proji,orienti,truthi,cotrui,delxi,delyi
   real                    ::  rlat1i,rlat2i,rlon1i,rlon2i
   real                    ::  rdat(imin*jmin)
   real                    ::  rout(io*jo)
   real, allocatable       ::  glon(:),glat(:)
   real                    ::  gaul(jmin)
#ifdef RMP
   integer                 ::  iminjmin,ijo,j,i,ipm,jpm,j1,j2,ifill,ij
   integer                 ::  isrs,isrf,jsrs,jsrf,n,jx,ix,ixjx,maxcount,maxtyp
   integer                 ::  iy
   real                    ::  dlamda0,rlat1o,rlat2o,rlon1o,rlon2o,cotruo
   real                    ::  x00,y00,x00o,y00o,dummy,sum1,sum2,wei1,wei2
   real                    ::  sumn,sums,a00,a10,a11,a01,wi1j1,wi2j1,wi1j2
   real                    ::  wi2j2,wsum,wsumiv,weimax
#endif
   integer, parameter      ::  numtype=100
   integer                 ::  ntype(numtype)
   real, allocatable       ::  rlat(:,:),rlon(:,:)
   real                    ::  rslmsk(imin*jmin)
   real                    ::  slmask(io*jo)
   logical lmask
   real,    allocatable    ::  d00(:),d10(:),d11(:),d01(:)
   integer, allocatable    ::  n00(:),n10(:),n11(:),n01(:)
   integer, allocatable    ::  ij1(:),icl(:),jcl(:)
!
   idim=ILOTS/2
   jdim=JLOTS*2
!
   allocate (rlat(idim,jdim),rlon(idim,jdim))
   allocate (d00(io),d10(io),d11(io),d01(io))
   allocate (n00(io),n10(io),n11(io),n01(io))
   allocate (ij1(io),icl(io),jcl(io))
!
#ifdef DBG
   print *,'inttyp=',inttyp
#endif
!
#ifdef RMP
   iminjmin = imin*jmin
   ijo = io*jo
   if(imin.eq.1.or.jmin.eq.1) then
     print*,' wrong imin jmin in sfc_interp_ll2rmp',imin,jmin
     call MPABORT
   endif
!
   if(delxi.eq.0..or.delyi.eq.0.) then
     print *,'delxi or delyi is zero in sfc_interp_ll2rmp'
     call MPABORT
   endif
!
!  compute latitude and longitude
!
   call rmp_setup_rsm_grid(rlat,rlon,rdelx,rdely,dlamda0)
   call rmp_trans2output_grid(rlat,1)
   call rmp_trans2output_grid(rlon,1)
!
   do j = 1,jo
     do i = 1,io
       rlat(i,j) = rlat(i,j) * 180. / pi_
       rlon(i,j) = rlon(i,j) * 180. / pi_
     enddo
   enddo
!
   if(inttyp.eq.2) then
     ipm=nint((rdelx/2.)/(delxi*111.))
     jpm=nint((rdely/2.)/(delyi*111.))
#ifdef DBG
     print *,'ipm,jpm=',ipm,jpm
#endif
   endif
!
!  find i and j index and weighting factor
!
   rlat1o = rlat(1,1)
   rlat2o = rlat(io,jo)
   rlon1o = rlon(1,1)
   rlon2o = rlon(io,jo)
!
#ifdef DBG
   print *, ' prepare grid to grid interpolation.'
   print *, ' input grid : '
   print *, ' proji=',proji,' orienti=',orienti,                               &
         ' truthi=',truthi,' cotrui=',cotrui,                                  &
         ' delxi=',delxi,' delyi=',delyi,                                      &
         ' rlat1i=',rlat1i,' rlon1i=',rlon1i,                                  &
         ' rlat2i=',rlat2i,' rlon2i=',rlon2i
   print *, ' output grid : '
   print *, ' projo=',rproj,' oriento=',rorient,                               &
         ' trutho=',rtruth,' cotruo=',rcotru,                                  &
         ' delxo=',rdelx,' delyo=',rdely,                                      &
         ' rlat1o=',rlat1o,' rlon1o=',rlon1o,                                  &
         ' rlat2o=',rlat2o,' rlon2o=',rlon2o
#endif
!
   if(abs(rproj).eq.4.) then
     print *,'no regional gaussian grid allowed in sfc_interp_ll2rmp'
     call MPABORT
   endif
!
   allocate (glon(imin),glat(jmin))
!
   call sfc_interp_init(proji,orienti,truthi,cotrui,                           &
               delxi,delyi,rlat1i,rlon1i,imin,jmin,                            &
               rproj,rorient,rtruth,rcotru,                                    &
               x00,y00,x00o,y00o,                                              &
               rlat1o,rlon1o,jo,dummy,glon,glat,gaul)
!
   sum1=0.
   sum2=0.
   wei1=0.
   wei2=0.
   j1 = 0
   j2 = imin*(jmin-1)
   do i = 1,imin
     sum1=sum1+rdat(i+j1) * rslmsk(i+j1)
     sum2=sum2+rdat(i+j2) * rslmsk(i+j2)
     wei1=wei1+rslmsk(i+j1)
     wei2=wei2+rslmsk(i+j2)
   enddo
#ifdef DBG
!
   print *,' step wei12 ',wei1,wei2
#endif
!
   if(rlat1i.gt.0.) then
     if(wei1.gt.0.) then
       sumn=sum1/wei1
     else
       sumn=0.
     endif
     if(wei2.gt.0.) then
       sums=sum2/wei2
     else
       sums=0.
     endif
   else
     if(wei1.gt.0.) then
       sums=sum1/wei1
     else
       sums=0.
     endif
     if(wei2.gt.0.) then
       sumn=sum2/wei2
     else
       sumn=0.
     endif
   endif
!
   ifill=0
   do j = 1,jo
     call sfc_interp_weight(io,j,proji,orienti,truthi,cotrui,                  &
               delxi,delyi,rlat1i,rlon1i,imin,jmin,                            &
               rproj,rorient,rtruth,rcotru,                                    &
               rdelx,rdely,rlat1o,rlon1o,dummy,                                &
               x00,y00,x00o,y00o,                                              &
               n00,n10,n11,n01,d00,d10,d11,d01,                                &
               glon,glat,gaul)
!
     if(inttyp.eq.1.or.inttyp.eq.2) then
       do i = 1,io
         weimax = 0.
         if(weimax.le.d00(i)) then
           weimax = d00(i)
           ij1(i) = n00(i)
         elseif(weimax.le.d01(i)) then
           weimax = d01(i)
           ij1(i) = n01(i)
         elseif(weimax.le.d10(i)) then
           weimax = d10(i)
           ij1(i) = n10(i)
         elseif(weimax.le.d11(i)) then
           weimax = d11(i)
           ij1(i) = n11(i)
         endif
         jcl(i)=(ij1(i)-1)/imin+1
         icl(i)=ij1(i)-(jcl(i)-1)*imin
       enddo
     endif
!
!  inttyp=1  take the closest point value
!
     if(inttyp.eq.1) then
       do i = 1,io
         ij=(j-1)*io+i
         rout(ij)=rdat(ij1(i))
       enddo
!
!  inttyp=2 take the predominant type
!
     elseif(inttyp.eq.2) then
       do i = 1,io
         ij=(j-1)*io+i
         jcl(i)=(ij1(i)-1)/imin+1
         icl(i)=ij1(i)-(jcl(i)-1)*imin
         isrs=max(icl(i)-ipm,1)
         isrf=min(icl(i)+ipm,imin)
         jsrs=max(jcl(i)-jpm,1)
         jsrf=min(jcl(i)+jpm,jmin)
         do n = 1,numtype
           ntype(n)=0
         enddo
         do jx = jsrs,jsrf
           do ix = isrs,isrf
             ixjx=(jx-1)*imin+ix
             n=rdat(ixjx)+1
#ifndef VECTOR
             if(n.gt.numtype) then
               print *,'type .gt.numtype.',                                    &
                       ' change numtype in sfc_interp_file_gmp'
               call MPABORT
             endif
             if(n.lt.0) then
               print *,'type cannot be less than zero'
               call MPABORT
             endif
#endif
             ntype(n)=ntype(n)+1
           enddo
         enddo
         maxcount=0
         maxtyp=-1
         do n = 1,numtype
           if(ntype(n).gt.maxcount) then
             maxcount=ntype(n)
             maxtyp=n
           endif
         enddo
         rout(ij)=maxtyp-1
       enddo
!
!  inttyp=0  bilinear interpolation without mask
!
     elseif(inttyp.eq.0.and..not.lmask) then
       do i = 1,io
         ij=(j-1)*io+i
         a00 = rdat(n00(i))
         a10 = rdat(n10(i))
         a11 = rdat(n11(i))
         a01 = rdat(n01(i))
         wi1j1 = d00(i)
         wi2j1 = d10(i)
         wi1j2 = d01(i)
         wi2j2 = d11(i)
         wsum  = wi1j1+wi2j1+wi1j2+wi2j2
         wsumiv = 1./wsum
#ifndef VECTOR
         if(n00(i).ne.n01(i)) then
#endif
           rout(ij)=(wi1j1*a00+wi2j1*a10+                                      &
                         wi1j2*a01+wi2j2*a11)*wsumiv
#ifndef VECTOR
         else
           if(n00(i).le.imin.and.rlat1i.gt.0..or.n00(i).gt.j2                  &
              .and.rlat1i.lt.0.) then
             rout(ij)=(wi1j1*sumn + wi2j1*sumn +                               &
             wi1j2*a01+wi2j2*a11)*wsumiv
           endif
           if(n00(i).le.imin.and.rlat1i.lt.0..or.n00(i).gt.j2                  &
              .and.rlat1i.gt.0.) then
             rout(ij)=(wi1j1*a00 + wi2j1*a10 +                                 &
                            wi1j2*sums + wi2j2*sums)*wsumiv
           endif
         endif
#endif
       enddo
!
!  inttyp=0  bilinear interpolation with mask
!
     elseif(inttyp.eq.0.and.lmask) then
       do i = 1,io
         ij=(j-1)*io+i
         a00 = rdat(n00(i))
         a10 = rdat(n10(i))
         a11 = rdat(n11(i))
         a01 = rdat(n01(i))
         if(slmask(ij).eq.rslmsk(n00(i)).and.                                  &
           slmask(ij).eq.rslmsk(n01(i)).and.                                   &
           slmask(ij).eq.rslmsk(n10(i)).and.                                   &
           slmask(ij).eq.rslmsk(n11(i)))then
           wi1j1 = d00(i)
           wi2j1 = d10(i)
           wi1j2 = d01(i)
           wi2j2 = d11(i)
         elseif(slmask(ij).eq.1.) then
           wi1j1 = d00(i) * rslmsk(n00(i))
           wi2j1 = d10(i) * rslmsk(n10(i))
           wi1j2 = d01(i) * rslmsk(n01(i))
           wi2j2 = d11(i) * rslmsk(n11(i))
         elseif(slmask(ij).eq.0.) then
           wi1j1 = d00(i) * (1.-rslmsk(n00(i)))
           wi2j1 = d10(i) * (1.-rslmsk(n10(i)))
           wi1j2 = d01(i) * (1.-rslmsk(n01(i)))
           wi2j2 = d11(i) * (1.-rslmsk(n11(i)))
         endif
         wsum  = wi1j1+wi2j1+wi1j2+wi2j2
         if(wsum.ne.0.) then
           wsumiv = 1./wsum
           if(n00(i).ne.n01(i)) then
             rout(ij)=(wi1j1*a00+wi2j1*a10+                                    &
                            wi1j2*a01+wi2j2*a11)*wsumiv
           else
             if(n00(i).le.imin.and.rlat1i.gt.0..or.n00(i).gt.j2                &
                .and.rlat1i.lt.0.) then 
               rout(ij)=(wi1j1*sumn + wi2j1*sumn +                             &
                               wi1j2*a01+wi2j2*a11)*wsumiv
             endif
             if(n00(i).le.imin.and.rlat1i.lt.0..or.n00(i).gt.j2                &
                .and.rlat1i.gt.0.) then
               rout(ij)=(wi1j1*a00 + wi2j1*a10 +                               &
                               wi1j2*sums + wi2j2*sums)*wsumiv
             endif
           endif
         else
           ifill=ifill+1
#ifdef DBG
           if(ifill.le.2) then
             write(6,*) 'n00 n01 n10 n11 = ',                                  &
                         n00(i),n01(i),n10(i),n11(i)
             write(6,*) 'rslmsk=',rslmsk(n00(i)),rslmsk(n01(i)),               &
                         rslmsk(n10(i)),rslmsk(n11(i))
             write(6,*) 'ij = ',ij,' slmask(ij)=',slmask(ij)
           endif
#endif
           ix = n11(i)
           iy = n11(i)
103        continue
           if(ix.gt.iminjmin) go to 101
           if((slmask(ij).eq.1..and.                                           &
               slmask(ij).eq.rslmsk(ix)).or.                                   &
              (slmask(ij).eq.0..and.                                           &
               slmask(ij).eq.rslmsk(ix))) then
             rout(ij) = rdat(ix)
             !
             go to 100
             !
           endif
           ix=ix+1
101        continue
           if(iy.lt.1) go to 102
           if((slmask(ij).eq.1..and.                                           &
               slmask(ij).eq.rslmsk(iy)).or.                                   &
              (slmask(ij).eq.0..and.                                           &
               slmask(ij).eq.rslmsk(iy))) then
             rout(ij) = rdat(iy)
             !
             go to 100
             !
           endif
           iy=iy-1
102        continue
           if(ix.le.iminjmin.or.iy.ge.1) go to 103
           write(6,*) ' error!!! no filling value found in sfc_interp_ll2rmp'
           print *,"debug print rslmsk"
           call sfc_quick_print(rslmsk,imin,jmin)
           print *,"debug print slmask"
           call sfc_quick_print(slmask,io,jo)
           call MPABORT
         endif
100      continue
!
       enddo
     endif
   enddo
!
   if(ifill.gt.1) then
     write(6,*) ' unable to interpolate.  filled with nearest',                &
                ' point value at ',ifill,' points'
   endif
!
   call print_maxmin_six(rout,io*jo,1,1,1,'sfc_interp_ll2rmp:rout')
!
#endif
   deallocate (rlat,rlon)
   deallocate (d00,d10,d11,d01,n00,n10,n11,n01,ij1,icl,jcl)
   deallocate (glon,glat)
!
   return
   end subroutine sfc_interp_ll2rmp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
#include <abort.h>
   subroutine sfc_interp_ll2mask(lugb,fnmskg,                                  &
                     inttyp,rout,io,jo,                                        &
                     projo,oriento,trutho,cotruo,                              &
                     delxo,delyo,rlat1o,rlat2o,rlon1o,rlon2o)
!-------------------------------------------------------------------------------
#undef VECTOR
#ifdef NEC
#define VECTOR
#endif
#ifdef ES
#define VECTOR
#endif
!fpp$ noconcur r
!-------------------------------------------------------------------------------
!
!  interpolation from lat/lon or gaussian grid to other lat/lon grid
!
!  inttyp ... type of interpolation.  0 = bilinear
!                                     1 = take the closest point value
!
!-------------------------------------------------------------------------------
   character(len=128)         ::  fnmskg
   integer, parameter         ::  imxgrmsk_max=2280,jmxgrmsk_max=1141
   real, dimension(io*jo)     ::  rout
   integer                    ::  kpds_grmsk(25)
   data kpds_grmsk/4*-1, 81, -1,  -1,18*-1/
!
   real, allocatable          :: grmsk(:)
   real, allocatable          :: glon(:),glat(:),gaul(:)
   real, allocatable          :: glao(:)
!
   real,    allocatable       ::  d00(:),d10(:),d11(:),d01(:)
   integer, allocatable       ::  n00(:),n10(:),n11(:),n01(:)
   integer, allocatable       ::  ij1(:)
!
   allocate (d00(io),d10(io),d11(io),d01(io))
   allocate (n00(io),n10(io),n11(io),n01(io))
   allocate (ij1(io))
!
!  read in high resolution global lat/lon land sea mask
!
   allocate (grmsk(imxgrmsk_max*jmxgrmsk_max))
   print *,'----------- start reading hires lat/lon mask -----'
   call sfc_read_grib(lugb,fnmskg,kpds_grmsk,                                  &
              0,0,0,0,0.,                                                      &
              imxgrmsk_max,jmxgrmsk_max,                                       &
              imxgrmsk,jmxgrmsk,                                               &
              kkpds5,kkgds1,                                                   &
              proji,orienti,truthi,cotrui,                                     &
              delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i,                         &
              dlon,dlat,wlon,rnlat,                                            &
              grmsk,iret)
!
   if(iret.gt.0) then
     print *,'high res global lat/lon land sea mask read failed.'
     call MPABORT
   endif
!
   if(imxgrmsk_max*jmxgrmsk_max.lt.imxgrmsk*jmxgrmsk) then
     print *,'dimension size overflow in sfc_interp_ll2mask'
     call MPABORT
   endif
!
   print *,'----------- finish reading hires lat/lon mask -----'
   close(lugb)
!
#ifdef DBG
   print *,' --- into sfc_interp_ll2mask --- '
   print *,'sfc_interp_ll2mask:grmsk'
   call sfc_quick_print(grmsk,imxgrmsk,jmxgrmsk)
#endif
!
   ijo = io*jo
   if(imxgrmsk.eq.1.or.jmxgrmsk.eq.1) then
     print*,' wrong imxgrmsk or jmxgrmsk in sfc_interp_ll2mask',               &
               imxgrmsk,jmxgrmsk
     call MPABORT
   endif
!
   if(delxi.eq.0..or.delyi.eq.0.) then
     print *,'delxi or delyi is zero in sfc_interp_ll2mask.'
     call MPABORT
   endif
!
#ifdef DBG
   print *,'delxi=',delxi
   print *,'delyi=',delyi
   print *,'rlon1i=',rlon1i
   print *,'rlat1i=',rlat1i
   print *,' prepare grid to grid interpolation.'
   print *,' input grid : '
   print *,' proj=',proji,' orient=',orienti,                                  &
           ' truth=',truthi,' cotru=',cotrui,                                  &
           ' delx=',delxi,' dely=',delyi,                                      &
           ' rlat1=',rlat1i,' rlon1=',rlon1i
   print *,' output grid : '
   print *,' projo=',projo,' oriento=',oriento,                                &
           ' trutho=',trutho,' cotruo=',cotruo,                                &
           ' delxo=',delxo,' delyo=',delyo,                                    &
           ' rlat1o=',rlat1o,' rlon1o=',rlon1o
#endif
!
   allocate (glon(imxgrmsk),glat(jmxgrmsk),gaul(jmxgrmsk))
   allocate (glao(jo))
!
   call sfc_interp_init(proji,orienti,truthi,cotrui,                           &
               delxi,delyi,rlat1i,rlon1i,imxgrmsk,jmxgrmsk,                    &
               projo,oriento,trutho,cotruo,                                    &
               x00,y00,x00o,y00o,                                              &
               rlat1o,rlon1o,jo,glao,glon,glat,gaul)
!
!  inttyp=1  take the closest point value
!
#ifdef DBG
   print *,'in sfc_interp_ll2mask.  inttyp=',inttyp
#endif
!
   do j = 1,jo
     call sfc_interp_weight(io,j,proji,orienti,truthi,cotrui,                  &
              delxi,delyi,rlat1i,rlon1i,imxgrmsk,jmxgrmsk,                     &
              projo,oriento,trutho,cotruo,                                     &
              delxo,delyo,rlat1o,rlon1o,glao(j),                               &
              x00,y00,x00o,y00o,                                               &
              n00,n10,n11,n01,d00,d10,d11,d01,                                 &
              glon,glat,gaul)
!
     if(inttyp.eq.1) then
       do i = 1,io
         weimax = 0.
         if(weimax.le.d00(i)) then
           weimax = d00(i)
           ij1(i) = n00(i)
         elseif(weimax.le.d01(i)) then
           weimax = d01(i)
           ij1(i) = n01(i)
         elseif(weimax.le.d10(i)) then
           weimax = d10(i)
           ij1(i) = n10(i)
         elseif(weimax.le.d11(i)) then
           weimax = d11(i)
           ij1(i) = n11(i)
         endif
       enddo
       do i = 1,io
         ij=(j-1)*io+i
         rout(ij)=grmsk(ij1(i))
       enddo
     else
       do i = 1,io
         ij=(j-1)*io+i
         a00 = grmsk(n00(i))
         a10 = grmsk(n10(i))
         a11 = grmsk(n11(i))
         a01 = grmsk(n01(i))
         wi1j1 = d00(i)
         wi2j1 = d10(i)
         wi1j2 = d01(i)
         wi2j2 = d11(i)
         wsum  = wi1j1+wi2j1+wi1j2+wi2j2
         wsumiv = 1./wsum
         rout(ij)=(wi1j1*a00+wi2j1*a10+                                        &
                   wi1j2*a01+wi2j2*a11)*wsumiv
       enddo
     endif
   enddo
!
#ifdef DBG
   call print_maxmin_six(rout,io*jo,1,1,1,'mask for input grib from grmsk')
#endif
   deallocate (grmsk,glon,glat,gaul,glao)
   deallocate (d00,d10,d11,d01,n00,n10,n11,n01,ij1)
#undef VECTOR
!
   return
   end subroutine sfc_interp_ll2mask
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_interp_ll2xyr(io,cproj,corient,ctruth,ccotru,clat,clon,cx,cy)
!-------------------------------------------------------------------------------
!
! input are all degree and output to global x y, not domain x y
!
! if proj= 0  do mercater projection
! if proj= 1  do north polar projection
! if proj=-1  do south polar projection
! if proj= 2  do north lambert projection
! if proj=-2  do south lambert projection
!
!-------------------------------------------------------------------------------
   use constant, only : pi=>pi_,rerth=>rerth_
!-------------------------------------------------------------------------------
   real,parameter       ::  twopi=2.0*pi,hfpi=0.5*pi,qtpi=0.5*hfpi
   real,parameter       ::  rad=pi/180.
   real                 ::  clat(io),clon(io),cx(io),cy(io)
!
! polar projection
!
   nproj = cproj
#ifdef DBG
   print*, 'nproj in ll2', nproj
#endif
   if( nproj.eq.1 .or. nproj.eq.-1 ) then
     truth  = ctruth * rad
     truth  = nproj * truth
     orient  = corient * rad
     cenlon = mod(orient,twopi)
     if(cenlon.lt.0.e0) cenlon = twopi + cenlon
     dlamda0 = cenlon + hfpi
     a2 =  rerth * ( 1.0 + sin(truth) )
     radlat = 90. * rad
     radlon = 0.0 * rad - dlamda0
     radlat = nproj * radlat
     radlon = nproj * radlon
     do i = 1,io
       blat = clat(i) * rad
       blon = clon(i) * rad
       rsoa2 = tan( (hfpi-blat*nproj)*0.5 )
       x2py2 = ( rsoa2 * a2 ) ** 2.0
       blon = mod(blon,twopi)
       if(blon.lt.0.e0) blon = twopi + blon
       rlon = nproj * (blon - dlamda0)
       rlon = amod(rlon,twopi)
       if( rlon.lt.0. ) rlon=twopi+rlon
       yox = tan(rlon)
       x = sqrt( x2py2/(1.+yox*yox) )
       y = sqrt( x2py2 - x*x )
       if( rlon.gt.hfpi .and. rlon.lt. pi+hfpi ) x = -x
       if( rlon.gt.pi .and. rlon.lt. twopi ) y = -y
       cx(i) = x
       cy(i) = y
     enddo
!
!  mercator
!
   else if ( nproj.eq.0 ) then
     truth  = ctruth * rad
     cenlon = corient * rad
     cenlon = mod(cenlon,twopi)
     if(cenlon.lt.0.e0) cenlon = twopi + cenlon
     a2 =  rerth * cos( truth )
     dlamda0 = 0.0
     do i = 1,io
       blat = clat(i) * rad
       blon = clon(i) * rad
       blon = mod(blon,twopi)
       if(blon.lt.0.e0) blon = twopi + blon
       x=a2*(blon-cenlon)
       y=a2*log(tan(blat/2.0+qtpi))
       cx(i) = x
       cy(i) = y
     enddo
!
!  lambert
!
   else if( nproj.eq.2 .or. nproj.eq.-2 ) then
     is = 1
     if( nproj.lt.0 ) is = -1
     truth  = ctruth * rad
     cotru  = ccotru * rad
     cenlon = corient * rad
     cenlon = mod(cenlon,twopi)
     if(cenlon.lt.0.e0) cenlon = twopi + cenlon
     if( truth.eq.cotru ) then
       cone= cos (hfpi-is*truth)
     else
       cone=(log(cos(truth))-log(cos(cotru)))/                                 &
              (log(tan(qtpi-is*truth/2))-log(tan(qtpi-is*cotru/2)))
     endif
     r00 = rerth/cone*cos(truth)/(tan(qtpi-is*truth/2))**cone
     do i = 1,io
       blat = clat(i) * rad
       blon = clon(i) * rad
       blon = mod(blon,twopi)
       if(blon.lt.0.e0) blon = twopi + blon
       r =  r00*(tan(qtpi-is*blat/2))**cone
       x =    r*sin(cone*(blon-cenlon))
       y = -is*r*cos(cone*(blon-cenlon))
       cx(i) = x
       cy(i) = y
     enddo
   endif
!
#ifdef DBG
   print*, 'cx in ll2', (cx(i),i=10,20)
#endif
   return
   end subroutine sfc_interp_ll2xyr
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_interp_xy2llr(io,cproj,corient,ctruth,ccotru,cx,cy,clat,clon)
!-------------------------------------------------------------------------------
   use constant, only : pi=>pi_,rerth=>rerth_
!-------------------------------------------------------------------------------
   real, parameter      ::  twopi=2.0*pi,hfpi=0.5*pi,qtpi=0.5*hfpi
   real,parameter       ::  rad=pi/180.
   real                 ::  cx(io),cy(io)
   real                 ::  clat(io),clon(io)
!
! input are all degree and output to global x y, not regional x y
!
! if proj=0  do mercator projection
! if proj=1  do north polar projection
! if proj=-1 do south polar projection
! if proj=2  do north lambert projection
! if proj=-2 do south lambert projection
!
! polar projection
!
   nproj = cproj
   if( nproj.eq.1 .or. nproj.eq.-1 ) then
     truth  = ctruth * rad
     truth  = nproj * truth
     orient  = corient * rad
     cenlon = mod(orient,twopi)
     if(cenlon.lt.0.e0) cenlon = twopi + cenlon
     dlamda0 = cenlon + hfpi
     a2 =  rerth * ( 1.0 + sin(truth) )
     radlat = 90. * rad
     radlon = 0.0 * rad - dlamda0
     radlat = nproj * radlat
     radlon = nproj * radlon
     do i = 1,io
       x = cx(i)
       y = cy(i)
       if( x.gt.0.0 ) then
         blon = atan(y/x)
       else if( x.lt.0.0 ) then
         blon = pi + atan(y/x)
       else
         blon = hfpi
         if( y.lt.0.0 ) blon = blon * 3.0
       endif
       blon = blon + dlamda0
       blon = mod(blon,twopi)
       blon = nproj * blon
       rsoa2 = sqrt( x*x + y*y )/a2
       blat = hfpi - 2. * atan(rsoa2)
       clat(i) = blat / rad
       clon(i) = blon / rad
     enddo
!
!  mercator
!
   else if ( nproj.eq.0 ) then
     truth  = ctruth * rad
     cenlon = corient * rad
     cenlon = mod(cenlon,twopi)
     if(cenlon.lt.0.e0) cenlon = twopi + cenlon
     a2 =  rerth * cos( truth )
     dlamda0 = 0.0
     do i = 1,io
       x = cx(i)
       y = cy(i)
       blon = x / a2 + cenlon
       blon = mod(blon,twopi)
       if(blon.lt.0.e0) blon = twopi + blon
       blat = 2.*(atan(exp(y/a2))-qtpi)
       clat(i) = blat / rad
       clon(i) = blon / rad
     enddo
!
!  lambert
!
   else if( nproj.eq.2 .or. nproj.eq.-2 ) then
     is = 1
     if( nproj.lt.0 ) is=-1
     truth  = ctruth * rad
     cotru  = ccotru * rad
     cenlon = corient * rad
     cenlon = mod(cenlon,twopi)
     if(cenlon.lt.0.e0) cenlon = twopi + cenlon
     if( truth.eq.cotru ) then
       cone = cos (hfpi-is*truth)
     else
       cone =(log(cos(truth))-log(cos(cotru)))/                                &
             (log(tan(qtpi-is*truth/2))-log(tan(qtpi-is*cotru/2)))
     endif
     r00 = rerth/cone*cos(truth)/(tan(qtpi-is*truth/2))**cone
     do i = 1,io
       r = sqrt(cx(i)*cx(i) + cy(i)*cy(i) )
       blon = cenlon + asin(cx(i)/r) / cone
       blon = mod(blon+twopi,twopi)
       blat = hfpi - 2 * is * atan ( (r/r00)**(1./cone) )
       clat(i) = blat / rad
       clon(i) = blon / rad
     enddo
   endif
!
   return
   end subroutine sfc_interp_xy2llr
!-------------------------------------------------------------------------------
