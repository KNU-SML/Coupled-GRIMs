#include <define.h>
   subroutine sfc_interp_setup
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [sfc_interp_setup]
!           |
!           |-- [sfc_interp_weight] *
!           |-- [sfc_interp_init] *
!
!-------------------------------------------------------------------------------
   end subroutine sfc_interp_setup
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_interp_weight(io,j,proj ,orient ,truth ,cotru ,              &
                  delx ,dely ,rlat1 ,rlon1 ,ii, jj,                            &
                  projo,oriento,trutho,cotruo,                                 &
                  delxo,delyo,rlat1o,rlon1o,glao,                              &
                  x00,y00,x00o,y00o,                                           &
                  n00,n10,n11,n01,d00,d10,d11,d01,                             &
                  glon,glat,gaul)
!-------------------------------------------------------------------------------
   real              ::  d00(io),d10(io),d11(io),d01(io)
   real              ::  glat(jj),glon(ii),gaul(jj)
   real, allocatable :: rlato(:),rlono(:)
   real, allocatable :: xlat (:), xlon(:)
   real, allocatable :: xo   (:),   yo(:)
   real, allocatable :: x    (:),    y(:)
!
   integer           ::  n00(io),n10(io),n11(io),n01(io)
!
   allocate (rlato(io),rlono(io))
   allocate (xlat (io),xlon (io))
   allocate (x    (io),y    (io))
!
!  get output latitude
!
   if(abs(projo).eq.3..or.abs(projo).eq.4.) then
     if(projo.eq.3.) then
       do i = 1,io
         rlato(i)=rlat1o+delyo*float(j-1)/1000.
       enddo
     elseif(abs(projo).eq.4.) then
       do i = 1,io
         rlato(i)=glao
       enddo
     elseif(projo.eq.-3.) then
       do i = 1,io
         rlato(i)=rlat1o-delyo*float(j-1)/1000.
       enddo
     endif
!
!  get output longitude
!
     do i = 1,io
       rlono(i)=mod(rlon1o+delxo*float(i-1)/1000.+360.,360.)
     enddo
     if(abs(proj).eq.3..or.abs(proj).eq.4.) then
       call sfc_interp_ll2gmp(io,proj,glat,glon,ii,jj,delx,dely,rlat1,rlon1,   &
                    rlato,rlono,x,y)
     else
       call sfc_interp_ll2xyr(io,proj,orient,truth,cotru,                      &
                                rlato,rlono,x,y)
     endif
!
   else
     allocate (xo(io),yo(io))
     do i = 1,io
       yo(i)=y00o+(j-1)*delyo
     enddo
     do i = 1,io
       xo(i)=x00o+(i-1)*delxo
     enddo
     call sfc_interp_xy2llr(io,projo,oriento,trutho,cotruo,xo,yo,              &
                             rlato,rlono)
     if(abs(proj).eq.3..or.abs(proj).eq.4.) then
       call sfc_interp_ll2gmp(io,proj,glat,glon,ii,jj,delx,dely,rlat1,rlon1,   &
                                rlato,rlono,x,y)
     else
       call sfc_interp_ll2xyr(io,proj,orient,truth,cotru,rlato,rlono,          &
                                x,y)
     endif
     deallocate (xo,yo)
   endif
!
!  input latitude and longitude
!
   if(abs(proj).eq.3..or.abs(proj).eq.4.) then
     do i = 1,io
       xlon(i)=x(i)
       xlat(i)=y(i)
     enddo
   else
     do i = 1,io
       xlon(i)=(x(i)-x00)/delx+1
       xlat(i)=(y(i)-y00)/dely+1
     enddo
   endif
!
!  computation of coefficients
!
   do i = 1,io
     lon=max(xlon(i),1.)
     lat=max(xlat(i),1.)
     lon=min(lon,ii)
     lat=min(lat,jj)
     d00(i)=(1.-(xlon(i)-lon)) * (1.-(xlat(i)-lat))
     d10(i)=(xlon(i)-lon) * (1.-(xlat(i)-lat))
     d11(i)=(xlon(i)-lon) * (xlat(i)-lat)
     d01(i)=(1.-(xlon(i)-lon)) * (xlat(i)-lat)
     n00(i)=lon   + (lat   -1)*ii
     n10(i)=lon+1 + (lat   -1)*ii
     n11(i)=lon+1 + (lat+1 -1)*ii
     n01(i)=lon   + (lat+1 -1)*ii
   enddo
!
   deallocate (rlato,rlono,xlat,xlon,x,y)
!
   return
   end subroutine sfc_interp_weight
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sfc_interp_init(proj ,orient ,truth ,cotru ,                     &
                     delx ,dely ,rlat1 ,rlon1 ,ii, jj,                         &
                     projo,oriento,trutho,cotruo,                              &
                     x00,y00,x00o,y00o,                                        &
                     rlat1o,rlon1o,jo,glao,glon,glat,gaul)
!-------------------------------------------------------------------------------
   logical*1 lflag
!
   real                 ::  glat(jj),glon(ii),gaul(jj)
   real                 ::  glao(jo)
   real                 ::  gauo(jo)
!
   if(abs(proj).eq.4.) then  ! gaussian grid lat/lon data
     call gaulat(gaul,jj)   ! n-->s ( 0 --> 180)
     do j = 1,jj
       if(proj.eq.-4.) then      ! n--> s
         glat(j) = 90.-gaul(j)
       elseif(proj.eq.4) then    ! s--> n
         glat(jj-j+1) = 90.-gaul(j)
       endif
     enddo
   endif
!
   if(abs(projo).eq.4.) then  ! gaussian grid lat/lon data
     call gaulat(gauo,jo)   ! n-->s ( 0 --> 180)
     do j = 1,jo
       if(projo.eq.-4.) then      ! n--> s
         glao(j) = 90.-gauo(j)
       elseif(projo.eq.4) then    ! s--> n
         glao(jo-j+1) = 90.-gauo(j)
       endif
     enddo
   endif
!
   if(abs(proj).eq.3.) then  ! lat/lon data
     do j = 1,jj
       if(proj.eq.3.) then
         glat(j) = rlat1 + dely/1000. * float(j-1)
       elseif(proj.eq.-3.) then
         glat(j) = rlat1 - dely/1000. * float(j-1)
       endif
     enddo
   endif
!
   if(abs(proj).eq.3..or.abs(proj).eq.4.) then
     do i = 1,ii
       glon(i) = rlon1 + delx/1000. * float(i-1)
       if(glon(i).lt.0.) glon(i) = glon(i) + 360.
     enddo
   else
     call sfc_interp_ll2xyr(1,proj,orient,truth,cotru,rlat1,rlon1,x00,y00)
#ifdef DBG
     print *, ' input grid rlat1 rlon1 x00 y00 ',                              &
                        rlat1,rlon1,x00,y00
#endif
   endif
#ifdef DBG
!
   print*,'glat ',(glat(j),j=1,jj)
   print*,'glon ',(glon(i),i=1,ii)
#endif
!
   if(abs(projo).ne.3..and.abs(projo).ne.4.) then
     call sfc_interp_ll2xyr(1,projo,oriento,trutho,cotruo,rlat1o,rlon1o,       &
                             x00o,y00o)
#ifdef DBG
     print *, ' output grid rlat1  rlon1  x00  y00 ',                          &
                         rlat1o,rlon1o,x00o,y00o
#endif
   endif
!
   return
   end subroutine sfc_interp_init
!-------------------------------------------------------------------------------
