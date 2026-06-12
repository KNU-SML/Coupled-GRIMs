#include <define.h> 
   module nislq
#ifdef NISLQ
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! Common variables
!
   integer                               ::  nx         ,& ! global longitude
                                             my         ,& ! global latitude
                                             my_max     ,& ! local  latitude
                                             lev        ,& ! global level
                                             ncld       ,& ! # of hydrometeor.
                                             nlevs      ,&
                                             nlevsp     ,&
!                                            jlistnum   ,& ! myrank lat. length
!                                            lat_s      ,& ! myrank lat. start
!                                            lat_e      ,& ! myrank lat. end
                                             lonfull    ,& ! nx
                                             latfull    ,& ! my*2
                                             lonhalf    ,& ! nx/2
                                             lathalf    ,& ! my
                                             lonpart    ,& ! lonhalf/nsize+1
                                             latpart    ,& ! lathalf/nsize+1
                                             mylonlen   ,& ! myrank lon. length
                                             mylatlen      ! myrank lat. length
   integer, allocatable, dimension(:)    ::  lonstr     ,& ! allpe lon. start
                                             lonlen     ,& ! allpe lon. length
                                             latstr     ,& ! allpe lat. start
                                             latlen     ,& ! allpe lat. length
                                             truej      ,&
                                             shflj      ,&
                                             jlist1
   real   , allocatable, dimension(:)    ::  glat       ,& ! gaussian lat.(my)
                                             gglat      ,& ! gaussian lat.(my*2)
                                             gglati     ,&
                                             gglon      ,& ! gaussian lon. (nx)
                                             ggloni     ,&
                                             racos      ,& ! 1./(Re*cos)
                                             racos2        ! 1./(Re.cos^2)
   real   , allocatable, dimension(:,:,:) :: slq_q1     ,& ! moisture at n-1
                                             slq_q2     ,& ! moisture at n
                                             slq_q3        ! moisture at n+1
   real   , allocatable, dimension(:,:)   :: slq_psfc2
   real   , allocatable, dimension(:,:,:) :: slq_u2,slq_v2,slq_w2
                                             
   contains
!-------------------------------------------------------------------------------
   subroutine nislq_init(nsize,myrank,colrad,rbs2)
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_,latg_,latg2_,levs_,ntotal_
   use paramodel, only : LEVSS
#ifdef DFS
   use dfsvar   , only : iope,ib,jbw,levh,latdef
#else
   use comio    , only : iope
   use paramodel, only : LONF2S,LATG2S,levh_
#ifdef MP
   use commpi   , only : latdef
#else
   use comfgrid , only : latdef
#endif
#endif /* DFS end */
   use constant , only : pi_,rrerth_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
   integer,intent(in)                    :: nsize, myrank
   real   ,intent(in), dimension(latg2_) :: colrad,rbs2
!
! local variables
!
   integer                               :: myh,my2,i,j,ii,jj,js,k
   real                                  :: hfpi,twopi,dlat,dlon
   integer                               :: n,nm,nr,lonp
   real, dimension(latg2_)               :: colatrad
   integer                               :: j1,j2
!
! assign parameter
!
   nx      = lonf_
   my      = latg_
   lev     = levs_
!sldp   ncld    = 1+ntotal_
   ncld    = ntotal_
   nlevs   = ncld*lev
   nlevsp  = ncld*LEVSS
!
   lonfull = nx
   lonhalf = nx / 2             ! nx has to be even
   lonpart = lonhalf/nsize+1
   latfull = my * 2
   lathalf = my 
!  latpart = lathalf/nsize+1
   latpart = 2*((lathalf/2)/nsize+1)
#ifdef MP
   my_max  = latpart
#else
   my_max  = latg_
#endif
!
! constant
!
   hfpi  = pi_ * 0.5
   twopi = pi_ * 2.0
   myh   = my/2
   my2   = my*2
!
! specific humidity for nislq
!
#ifdef DFS
   allocate(  slq_q1   (ib,jbw,levh)                                          ,&
              slq_q2   (ib,jbw,levh)                                          ,&
              slq_q3   (ib,jbw,levh)                                          ,&
              slq_psfc2(ib,jbw)                                               ,&
              slq_u2   (ib,jbw,levs_)                                         ,&
              slq_v2   (ib,jbw,levs_)                                         ,&
              slq_w2   (ib,jbw,levs_+1)                                        )
#else /* SPH */
   allocate(  slq_q1   (LONF2S,levh_  ,LATG2S)                                ,&
              slq_q2   (LONF2S,levh_  ,LATG2S)                                ,&
              slq_q3   (LONF2S,levh_  ,LATG2S)                                ,&
              slq_psfc2(LONF2S        ,LATG2S)                                ,&
              slq_u2   (LONF2S,levs_  ,LATG2S)                                ,&
              slq_v2   (LONF2S,levs_  ,LATG2S)                                ,&
              slq_w2   (LONF2S,levs_+1,LATG2S)                                 )
#endif /* DFS end */
!
! initialize slq_q
!
   slq_q1=0.
   slq_q2=0.
   slq_q3=0.
   slq_psfc2=0.
   slq_u2=0.
   slq_v2=0.
   slq_w2=0.
!
! informations for latitude
!
   allocate(jlist1(latg2_))
   allocate(racos(latg2_))
   allocate(racos2(latg2_))
!
! define latitude
!
   do j = 1,latg2_
#ifdef MP
     jlist1(j)=j
#else
     jlist1(j)=j
#endif
   enddo
!
! define colatitude in radian
!
   do j = 1,latg2_
     colatrad(j)=colrad(j)
   enddo
!
! ------------------- gaussian latitude 
!
   allocate ( glat(my) )
   allocate ( gglat(my2), gglati(my2+1) )
!
   do j = 1,myh
     glat(j)=hfpi-colatrad(j)
     glat(my+1-j)=-glat(j)
   enddo
!
   do j = 1,myh       ! co-latitude
     gglat(j) = hfpi - glat(j)
     gglat(my+1-j) = hfpi + glat(j)
   enddo
   do j = my+1,my2
     gglat(j) = twopi - gglat(my2+1-j)
   enddo
!
   gglati(myh+1) = hfpi
   dlat = (gglati(myh+1)-gglat(myh))*2.0
   do j = myh,2,-1
     gglati(j) = gglati(j+1) - dlat
     dlat      = (gglati(j) - gglat(j-1))*2.0
   enddo
   gglati(1) = 0.0
   do j = myh+2,my
     gglati(j) = pi_ - gglati(my+2-j)
   enddo
   gglati(my+1) = pi_
   do j = my+2,my2
     gglati(j) = twopi - gglati(my2+2-j)
   enddo
   gglati(my2+1) = twopi
#ifdef SLDBG
!
   if( iope ) then
     print *,' j interface =1    gglati=',gglati(1)
     do j = 1,my2
       print *,'               j mean =',j,'  gglat= ',gglat(j)
       print *,' j interface =',j+1,' gglati= ',gglati(j+1)
     enddo
   endif
#endif
!
!  determin the longitude with full grid
!
   allocate( gglon(nx), ggloni(nx+1) )
   dlon = twopi / nx
   do i = 1,nx
     gglon(i)=(i-1)*dlon
   enddo
   do i = 2,nx
     ggloni(i)=0.5*(gglon(i-1)+gglon(i))
   enddo
   ggloni(     1)=ggloni(   2)-dlon
   ggloni(nx+1)=ggloni(nx)+dlon
!
#ifdef SLDBG
   if( iope ) then
     print *,' ------ total edge nx=',nx,' -------'
     print *,' i edge number=1    ggloni=',ggloni(1)
     do i = 1,nx
       print *,'               i cell number=',i,'  gglon= ',gglon(i)
       print *,' i edge number=',i+1,' ggloni= ',ggloni(i+1)
     enddo
   endif
#endif
!
! --------------------- for parallel --------------------
!
! in nisl, we do great circle, so nx and my2 are full
! transpose will between (lonfull,lev,latpart)  by (nx  ,lev,my/nsize+1)
!                    and (latfull,lev,lonpart)  by (my*2,lev,nx/2/nsize+1)
!
   allocate( lonstr(nsize), lonlen(nsize) )
   allocate( latstr(nsize), latlen(nsize) )
!
! equally distribute len, no location is considered
!
   lonlen(1:nsize) = 0
   i=1
   do ii = 1,lonpart
     do n = 1,nsize
       if(i.le.lonhalf) then
         lonlen(n) = lonlen(n)+1
         i=i+1
       endif
     enddo
   enddo
!
! sequential location for longitude
!
   nm=1
   do n = 1,nsize
     lonstr(n) = nm
     nm = nm + lonlen(n)
   enddo
!
! check make_list  to have consistent latitude number for each pe here
!
! refer equdis
   latlen(1:nsize) = 0
   i=1
   n=1
   do jj = 1,myh
     latlen(n)=latlen(n)+1
     n=n+i
     if(n.eq.nsize+1) then
       i=-1
       n=n+i
     endif
     if(n.eq.0) then
       i=1
       n=n+i
     endif
   enddo
   latlen=latlen*2
#ifdef SLDBG
   if(iope) then
     do n = 1,nsize
       print *,'n,latlen=',n,latlen(n)
     enddo
   endif
#endif
!
! sequential location for latitude
!
   nm=1
   do n = 1,nsize
     latstr(n) = nm
     nm = nm + latlen(n)
   enddo
!
!   jlistnum=latlen(myrank+1)
!   lat_s=latstr(myrank+1)
!   lat_e=latstr(myrank+1)+latlen(myrank+1)-1
#ifdef SLDBG
!
! check by print
!
   if( iope ) then
     do n = 1,nsize
       print *,' pe lonstr lonlen ',n,lonstr(n),lonlen(n)
       print *,' pe latstr latlen ',n,latstr(n),latlen(n)
     enddo
   endif
#endif
!
! make true latitude index follow make_list
!
   allocate( truej(my), shflj(my) )
!!   j=1
!!   do jj=1,my/nsize+1
!!     js = jj
!!     do n=1,nsize
!!       if( j.le.my ) then
!!         truej (js) = j
!!         shflj (j ) = js
!!         j = j + 1
!!       endif
!!       js = js + latlen(n)
!!     enddo
!!   enddo
   do jj = 1,latg2_
     j1=jj*2-1
     j2=jj*2
     truej(j1)=latdef(jj)
     truej(j2)=latg_-latdef(jj)+1
   enddo
#ifdef SLDBG
!
! check by print
!
   if(iope) then
!     do j=1,my
!       print *,' true j=',j,' to shaffled j=',shflj(j),   &
!                 ' back to true j=',truej(shflj(j))
!     enddo
     do j = 1,my
       print *,' shfl j=',j,' to true j=',truej(j)                            ,&
                 ' back to shfl j=',shflj(truej(j))
     enddo
     j=0
     do n = 1,nsize
       print *,' --- start pe =',n-1
       do jj = 1,latlen(n)
         j=j+1
         print *,' shaffled j=',j,' to true j=',truej(j)
       enddo
     enddo
   endif
#endif
!
   mylonlen = lonlen(myrank+1)
   mylatlen = latlen(myrank+1)
! debug
!  if( n.eq.n ) then
!    call mpe_finalize
!    stop
!  endif
!
   return
   end subroutine nislq_init
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine cyclic_cell_massadvx(levs,ncld,delt,uc,qq,mass)
!-------------------------------------------------------------------------------
!
! compute local positive advection with mass conservation
! qq is advected by uc which is in radiance/sec from past to next position
!
! author: hann-ming henry juang 2008
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real   , parameter                    ::  fa1 = 9./16.                     ,&
                                             fa2 = 1./16.
   integer                               ::  levs,ncld,mass
   integer                               ::  i,k,im,n
   real                                  ::  delt, rm2, dist, sc
   real   , dimension(lonfull,levs)      ::  uc
   real   , dimension(lonfull,levs,ncld) ::  qq
   real   , dimension(lonfull)           ::  past,next,da,dxfact
   real   , dimension(lonfull+1)         ::  xpast,xnext,uint
!
! preparations ---------------------------
!
! x is equal grid spacing, so location can be specified by grid point number
!
   im = lonfull
!
   do k = 1,levs
!
! 4th order interpolation from mid point to cell interfaces
!
     do i = 3,im-1
       uint(i)=fa1*(uc(i,k)+uc(i-1,k))-fa2*(uc(i+1,k)+uc(i-2,k))
     enddo
     uint(2)=fa1*(uc(2,k)+uc(1 ,k))-fa2*(uc(3,k)+uc(im  ,k))
     uint(1)=fa1*(uc(1,k)+uc(im,k))-fa2*(uc(2,k)+uc(im-1,k))
     uint(im+1)=uint(1)
     uint(im  )=fa1*(uc(im,k)+uc(im-1,k)) -fa2*(uc(1,k)+uc(im-2,k))
!
! compute past and next positions of cell interfaces
!
     do i = 1,im+1
       dist     = uint(i) * delt
       xpast(i) = ggloni(i) - dist
       xnext(i) = ggloni(i) + dist
     enddo
!      
     if( mass.eq.1 ) then
       do i = 1,im
         dxfact(i) = (xpast(i+1)-xpast(i)) / (xnext(i+1)-xnext(i))
       enddo
     endif
!
!  mass positive advection
!
     sc=ggloni(im+1)-ggloni(1)
     do n = 1,ncld
       past(1:im) = qq(1:im,k,n)
       call cyclic_cell_ppm_intp(ggloni,past,xpast,da,im,im,im,sc)
       if(mass.eq.1) da(1:im) = da(1:im) * dxfact(1:im)
       call cyclic_cell_ppm_intp(xnext,da,ggloni,next,im,im,im,sc)
       qq(1:im,k,n) = next(1:im)
     enddo
!
   enddo
!
   return
   end subroutine cyclic_cell_massadvx
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine cyclic_cell_massadvy(levs,ncld,delt,vc,qq,mass)
!-------------------------------------------------------------------------------
!
! compute local positive advection with mass conserving
! qq will be advect by vc from past to next location with 2*delt
!
! author: hann-ming henry juang 2007
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real   , parameter                    ::  fa1 = 9./16.                     ,&
                                             fa2 = 1./16.
   integer                               ::  levs,ncld,mass                   ,&
                                             n,k,j,jm,jm2
   real                                  ::  delt,sc
   real   , dimension(latfull,levs)      ::  vc
   real   , dimension(latfull,levs,ncld) ::  qq
   real   , dimension(latfull)           ::  var,past,da,next,dyfact
   real   , dimension(latfull+1)         ::  ypast,ynext,dist
!
! preparations ---------------------------
!
   jm   = lathalf
   jm2  = latfull
!
   do k = 1,levs
!
     do j = 1,jm
       var(j)      =-vc(j   ,k) * delt
       var(j+jm)   = vc(j+jm,k) * delt
     enddo
!
     do j = 3,jm2-1
       dist(j)=fa1*(var(j)+var(j-1))-fa2*(var(j+1)+var(j-2))
     enddo
     dist(2)=fa1*(var(2)+var(1  ))-fa2*(var(3)+var(jm2  ))
     dist(1)=fa1*(var(1)+var(jm2))-fa2*(var(2)+var(jm2-1))
     dist(jm2+1)=dist(1)
     dist(jm2  )=fa1*(var(jm2)+var(jm2-1))-fa2*(var(1)+var(jm2-2))
!
     do j = 1,jm2+1
       ypast(j) = gglati(j) - dist(j)
       ynext(j) = gglati(j) + dist(j)
     enddo
!
     if( mass.eq.1 ) then
       do j = 1,jm2
         dyfact(j) = (ypast(j+1)-ypast(j)) / (ynext(j+1)-ynext(j))
       enddo
     endif
!
! advection all in y
!
     sc=gglati(jm2+1)-gglati(1)
     do n = 1,ncld
       past(1:jm2) = qq(1:jm2,k,n)
       call cyclic_cell_ppm_intp(gglati,past,ypast,da,jm2,jm2,jm2,sc)
       if( mass.eq.1 ) da(1:jm2) = da(1:jm2) * dyfact(1:jm2)
       call cyclic_cell_ppm_intp(ynext,da,gglati,next,jm2,jm2,jm2,sc)
       qq(1:jm2,k,n) = next(1:jm2)
     enddo
! 
   enddo
!
   return
   end subroutine cyclic_cell_massadvy
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine cyclic_cell_intpx(levs,imp,imf,qq)
!-------------------------------------------------------------------------------
!
! do  mass conserving interpolation from different grid at given latitude
!
! author: hann-ming henry juang 2008
!
!-------------------------------------------------------------------------------
   use constant, only : pi_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                          ::  levs, imp, imf, i,k,im
   real                             ::  two_pi,dxp,dxf,hfdxp,hfdxf,sc
   real   , dimension(lonfull,levs) ::  qq
   real   , dimension(lonfull)      ::  past,next
   real   , dimension(lonfull+1)    ::  xpast,xnext
!
   im = lonfull
! ..................................  
   if( imp.ne.imf ) then
! ..................................
     two_pi = 2.*pi_
     dxp = two_pi / imp
     dxf = two_pi / imf
     hfdxp = 0.5 * dxp
     hfdxf = 0.5 * dxf
!
     do i = 1,imp+1
       xpast(i) = (i-1) * dxp - hfdxp
     enddo
!
     do i = 1,imf+1
       xnext(i) = (i-1) * dxf - hfdxf
     enddo
!
     sc=two_pi
     do k = 1,levs
       do i = 1,imp
         past(i)=qq(i,k)
       enddo
       call cyclic_cell_ppm_intp(xpast,past,xnext,next,im,imp,imf,sc)
       do i = 1,imf
         qq(i,k)=next(i)
       enddo
     enddo
! .................       
   endif
! .................
   return
   end subroutine cyclic_cell_intpx
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine cyclic_cell_ppm_intp(pp,qq,pn,qn,lons,lonp,lonn,sc)
!-------------------------------------------------------------------------------
!
! mass conservation in cyclic bc interpolation: interpolate a group
! of grid point  coordiante call pp at interface with quantity qq at
! cell averaged to a group of new grid point coordinate call pn at
! interface with quantity qn at cell average with ppm spline.
! in horizontal with mass conservation is under the condition that
! pp(1)= pp(lons+1)=pn(lons+1)
!
! pp    location at interfac level as input
! qq    quantity at averaged-cell as input
! pn    location at interface of new grid structure as input
! qn    quantity at averaged-cell as output
! lons  numer of cells for dimension
! lonp  numer of cells for input
! lonn  numer of cells for output
! lev   number of vertical layers
! mono  monotonicity o:no, 1:yes
!
! author : henry.juang@noaa.gov
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                     ::  lons,lonp,lonn             ,&
                                   ik,le,kstr,kend            ,&
                                   i,k, kl, kh, kk, kkl, kkh
   real                        ::  length,sc                  ,&
                                   dqi,dqimax,dqimin          ,&
                                   tl,tl2,tl3,qql,dql         ,&
                                   th,th2,th3,qqh,dqh         ,&
                                   dpp,dpq,c1,c2              ,&
                                   dpph,dppl,dpqh,dpql
   real   , dimension(lons)    ::  qq,qn
   real   , dimension(lons+1)  ::  pp,pn
   real   , dimension(3*lonp)  ::  locs,mass,hh,dqmono,qmi,qpi
!
!  length = pp(lonp+1) - pp(1)
   length=sc
!
! arrange input array cover output location with cyclic boundary condition
!
   locs(lonp+1:2*lonp) = pp(1:lonp)
   do i = 1,lonp
     locs(i) = locs(i+lonp) - length
     locs(i+2*lonp) = locs(i+lonp) + length
   enddo
!
   find_kstr : do i = 1,3*lonp
     if( pn(1).ge.locs(i) .and. pn(1).lt.locs(i+1) ) then
       kstr = i
       exit find_kstr
     else
       cycle find_kstr
     endif
   enddo find_kstr
   kstr=max(1,kstr)
!
   mass(lonp+1:2*lonp) = qq(1:lonp)
   do i = 1,lonp
     mass(i) = mass(i+lonp)
     mass(i+2*lonp) = mass(i+lonp)
   enddo
!
! prepare grid spacing
!
   do i = lonp+1,2*lonp
     hh(i) = locs(i+1)-locs(i)
   enddo
   do i = 1,lonp
     hh(i) = hh(i+lonp)
     hh(i+2*lonp) = hh(i+lonp)
   enddo
!
! prepare location with monotonic concerns
!
   do i = lonp+1,2*lonp
     dqi = 0.25*(mass(i+1)-mass(i-1))
     dqimax = max(mass(i-1),mass(i),mass(i+1)) - mass(i)
     dqimin = mass(i) - min(mass(i-1),mass(i),mass(i+1))
     dqmono(i) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
   enddo
   do i = 1,lonp
     dqmono(i) = dqmono(i+lonp)
     dqmono(i+2*lonp) = dqmono(i+lonp)
   enddo
!
! compute value at interface with monotone
!
   do i = lonp+1,2*lonp
     qmi(i)=(mass(i-1)*hh(i)+mass(i)*hh(i-1))/(hh(i)+hh(i-1))                  &
             +(dqmono(i-1)-dqmono(i))/3.0
!    qmi(i)=(mass(i-1)*hh(i)+mass(i)*hh(i-1))/(hh(i)+hh(i-1))
   enddo
   qmi(2*lonp+1)=qmi(lonp+1)
   do i = lonp+1,2*lonp
     qpi(i)=qmi(i+1)
   enddo
!
! do less diffusive
!
   do i = lonp+1,2*lonp
     qmi(i)=mass(i)                                                            &
            -sign(min(abs(2.*dqmono(i)),abs(qmi(i)-mass(i))),2.*dqmono(i))
   enddo
   do i = lonp+1,2*lonp
     qpi(i)=mass(i)                                                            &
            +sign(min(abs(2.*dqmono(i)),abs(qpi(i)-mass(i))),2.*dqmono(i))
   enddo
!
! do monotonicity
!
!     if( mono.eq.1 ) then
   do i = lonp+1,2*lonp
     c1=qpi(i)-mass(i)
     c2=mass(i)-qmi(i)
     if( c1*c2.le.0.0 ) then
       qmi(i)=mass(i)
       qpi(i)=mass(i)
     endif
   enddo
   do i = lonp+1,2*lonp
     c1=(qpi(i)-qmi(i))*(mass(i)-0.5*(qpi(i)+qmi(i)))
     c2=(qpi(i)-qmi(i))*(qpi(i)-qmi(i))/6.
     if( c1.gt.c2 ) then
       qmi(i)=3.*mass(i)-2.*qpi(i)
     else if( c1.lt.-c2 ) then
       qpi(i)=3.*mass(i)-2.*qmi(i)
     endif
   enddo
!     endif
!
! extend array with cyclic condition
!
   do i = 1,lonp
     qmi(i)        = qmi(i+lonp)
     qmi(i+2*lonp) = qmi(i+lonp)
     qpi(i)        = qpi(i+lonp)
     qpi(i+2*lonp) = qpi(i+lonp)
   enddo
!
! start interpolation by integral of ppm spline
!
   kkl = kstr
   do i = 1,lonn
     kl = i
     kh = i + 1
! find kkh
     do kk = kkl+1,3*lonp
       if( pn(kh).lt.locs(kk) ) then
         kkh = kk-1
         go to 100
       endif
     enddo
100  continue
! mass interpolate
     if( kkh.eq.kkl ) then
!      print *,' condition 000000000000 ',kl,kh,kkl,kkh
!      print *,' pl ph ll lh ',pn(kl),pn(kh),locs(kkl),locs(kkl+1)
       tl=(pn(kl)-locs(kkl))/hh(kkl)
       tl2=tl*tl
       tl3=tl2*tl
       th=(pn(kh)-locs(kkl))/hh(kkl)
       th2=th*th
       th3=th2*th
       qqh=(th3-th2)*qpi(kkl)                                                  &
           +(th3-2.*th2+th)*qmi(kkl)+(-2.*th3+3.*th2)*mass(kkl)
       qql=(tl3-tl2)*qpi(kkl)                                                  &
           +(tl3-2.*tl2+tl)*qmi(kkl)+(-2.*tl3+3.*tl2)*mass(kkl)
       qn(i) = (qqh-qql)/(th-tl)
     else if( kkh.gt.kkl ) then
       tl=(pn(kl)-locs(kkl))/hh(kkl)
       tl2=tl*tl
       tl3=tl2*tl
       qql=(tl3-tl2)*qpi(kkl)                                                  &
           +(tl3-2.*tl2+tl)*qmi(kkl)+(-2.*tl3+3.*tl2)*mass(kkl)
       dql = mass(kkl)-qql
       th=(pn(kh)-locs(kkh))/hh(kkh)
       th2=th*th
       th3=th2*th
       dqh=(th3-th2)*qpi(kkh)                                                  &
           +(th3-2.*th2+th)*qmi(kkh)+(-2.*th3+3.*th2)*mass(kkh)
       dpp  = (1.-tl)*hh(kkl) + th*hh(kkh)
       dpq  = dql*hh(kkl) + dqh*hh(kkh)
       if( kkh-kkl.gt.1 ) then
!        print *,' condition 2222222222 ',kl,kh,kkl,kkh
!        print *,' pl ph ll lh ',pn(kl),pn(kh),locs(kkl),locs(kkh)
         do kk = kkl+1,kkh-1
           dpp = dpp + hh(kk)
           dpq = dpq + mass(kk)*hh(kk)
         enddo
       endif
       qn(i) = dpq / dpp
     else
       print *,' Error in cyclic_cell_ppm_intp location not found '
       call abort
     endif
! next one
     kkl = kkh
   enddo
!
   return
   end subroutine cyclic_cell_ppm_intp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vertical_cell_advect(lons,londim,levs,ncld,deltim,               &
                                   ppi,wwi,qql,mass)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                               ::  lons,londim,levs,ncld,i,k,n,mass
   real                                  ::  deltim
   real   , dimension(londim,levs+1)     ::  ppi,wwi
   real   , dimension(londim,levs,ncld)  ::  qql
   real   , dimension(lons,levs)         ::  dsfact,rqmm,rqnn,rqda
   real   , dimension(lons,levs+1)       ::  ppii,ppid,ppia
!
   do k = 1,levs+1
     do i = 1,lons
       ppii(i,k)=ppi(i,k)
       ppid(i,k)=ppi(i,k)-wwi(i,k)*deltim
       ppia(i,k)=ppi(i,k)+wwi(i,k)*deltim
     enddo
   enddo
!
   if( mass.eq.1) then
     do k = 1,levs
       do i = 1,lons
         dsfact(i,k)=(ppid(i,k)-ppid(i,k+1))/(ppia(i,k)-ppia(i,k+1))
       enddo
     enddo
   endif
!
   do n = 1,ncld                              !hmhj nisl
     do k = 1,levs
       do i = 1,lons
         rqmm(i,k) = qql(i,k,n)
       enddo
     enddo
     call vertical_cell_ppm_intp(ppii,rqmm,ppid,rqda,lons,levs)
     if( mass.eq.1 ) then
       do k = 1,levs
         do i = 1,lons
           rqda(i,k) = rqda(i,k) * dsfact(i,k)
         enddo
       enddo
     endif
     call vertical_cell_ppm_intp(ppia,rqda,ppii,rqnn,lons,levs)
     do k = 1,levs
       do i = 1,lons
         qql(i,k,n)=rqnn(i,k)
       enddo
     enddo
   enddo
!
   return
   end subroutine vertical_cell_advect
!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------
   subroutine vertical_cell_ppm_intp(pp,qq,pn,qn,lons,levs)
!-------------------------------------------------------------------------------
!
! mass conservation in vertical interpolation: interpolate a group
! of grid point  coordiante call pp at interface with quantity qq at
! cell averaged to a group of new grid point coordinate call pn at
! interface with quantity qn at cell average with ppm spline.
! in vertical with mass conservation is under the condition that
! pp(1)=pn(1), pp(lev+1)=pn(lev+1)
!
! pp    pressure at interfac level as input
! qq    quantity at layer as input
! pn    pressure at interface of new grid structure as input
! qn    quantity at layer as output
! lev  numer of verical layers
!
! author : henry.juang@noaa.gov
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                           ::  lons,levs                            ,&
                                         i,k, kl, kh, kk, kkl, kkh
   real   , dimension(lons,levs)   ::  qq,qn
   real   , dimension(lons,levs+1) ::  pp,pn
!
   real                              ::  dqi,dqimax,dqimin    ,&
                                         tl,tl2,tl3,qql,dql   ,&
                                         th,th2,th3,qqh,dqh   ,&
                                         dpp,dpq,c1,c2
   real   , dimension(levs)          ::  hh,dqmono,qmi,qpi
   real   , dimension(0:levs+1)      ::  mass
!
   do i = 1,lons
!
     if( pp(i,1).ne.pn(i,1) .or. pp(i,levs+1).ne.pn(i,levs+1) ) then
       print *,' Error in vertical_cell_ppm_intp for domain values '
       print *,' i,pp1 pn1 ppt pnt ',i,pp(i,1),pn(i,1),pp(i,levs+1),pn(i,levs+1)
       call abort
     endif
!
! prepare thickness for uniform grid
!
     do k = 1,levs
       hh(k) = pp(i,k+1)-pp(i,k)         ! (top to bottm) ??
     enddo
!
! prepare location with monotonic concerns
!
     mass(1:levs)=qq(i,1:levs)
!
     mass(0)=(3.*hh(1)+hh(2))*mass(1)-2.*hh(1)*mass(2)
     mass(0)=mass(0)/(hh(1)+hh(2))
     mass(levs+1)=(3.*hh(levs)+hh(levs-1))*mass(levs)-2.*hh(levs)*mass(levs-1)
     mass(levs+1)=mass(levs+1)/(hh(levs)+hh(levs-1))
     do k = 1,levs
       dqi = 0.25*(mass(k+1)-mass(k-1))
       dqimax = max(mass(k-1),mass(k),mass(k+1)) - mass(k)
       dqimin = mass(k) - min(mass(k-1),mass(k),mass(k+1))
       dqmono(k) = sign( min( abs(dqimin), dqimin, dqimax ), dqi)
     enddo
!
! compute value at interface with momotone
!
     do k = 2,levs
       qmi(k)=(mass(k-1)*hh(k)+mass(k)*hh(k-1))/(hh(k)+hh(k-1))                &
             +(dqmono(k-1)-dqmono(k))/3.0
     enddo
     do k = 1,levs-1
       qpi(k)=qmi(k+1)
     enddo
     qmi(1)=mass(1)
     qpi(1)=mass(1)
     qmi(levs)=mass(levs)
     qpi(levs)=mass(levs)
!
! do monotonicity
!
!     if( mono.eq.1 ) then
     do k = 1,levs
       c1=qpi(k)-mass(k)
       c2=mass(k)-qmi(k)
       if( c1*c2.le.0.0 ) then
         qmi(k)=mass(k)
         qpi(k)=mass(k)
       endif
     enddo
     do k = 1,levs
       c1=(qpi(k)-qmi(k))*(mass(k)-0.5*(qpi(k)+qmi(k)))
       c2=(qpi(k)-qmi(k))*(qpi(k)-qmi(k))/6.
       if( c1.gt.c2 ) then
         qmi(k)=3.*mass(k)-2.*qpi(k)
       else if( c1.lt.-c2 ) then
         qpi(k)=3.*mass(k)-2.*qmi(k)
       endif
     enddo
!     endif
!
! start interpolation by integral of ppm spline
!
     kkl = 1
     do k = 1,levs
       kl = k
       kh = k + 1
! find kkh
       do kk = kkl+1,levs+1
!         if( pn(i,kh).ge.pp(i,kk) ) then
         if( pn(i,kh).le.pp(i,kk) ) then        ! top to bottom (hwangso)
           kkh = kk-1
           go to 100
         endif
       enddo
! mass interpolate
 100   if( kkh.eq.kkl ) then
         tl=(pn(i,kl)-pp(i,kkl))/hh(kkl)
         tl2=tl*tl
         tl3=tl2*tl
         th=(pn(i,kh)-pp(i,kkl))/hh(kkl)
         th2=th*th
         th3=th2*th
         qqh=(th3-th2)*qpi(kkl)+(th3-2.*th2+th)*qmi(kkl)   &
             +(-2.*th3+3.*th2)*mass(kkl)
         qql=(tl3-tl2)*qpi(kkl)+(tl3-2.*tl2+tl)*qmi(kkl)   &
             +(-2.*tl3+3.*tl2)*mass(kkl)
         qn(i,k) = (qqh-qql)/(th-tl)
       else if( kkh.gt.kkl ) then
         tl=(pn(i,kl)-pp(i,kkl))/hh(kkl)
         tl2=tl*tl
         tl3=tl2*tl
         qql=(tl3-tl2)*qpi(kkl)+(tl3-2.*tl2+tl)*qmi(kkl)   &
             +(-2.*tl3+3.*tl2)*mass(kkl)
         dql = qq(i,kkl)-qql
         th=(pn(i,kh)-pp(i,kkh))/hh(kkh)
         th2=th*th
         th3=th2*th
         dqh=(th3-th2)*qpi(kkh)+(th3-2.*th2+th)*qmi(kkh)   &
             +(-2.*th3+3.*th2)*mass(kkh)
         dpp= (1.-tl)*hh(kkl) + th*hh(kkh)
         dpq= dql*hh(kkl) + dqh*hh(kkh)
         if( kkh-kkl.gt.1 ) then
           do kk=kkl+1,kkh-1
             dpp = dpp + hh(kk)
             dpq = dpq + qq(i,kk)*hh(kk)
           enddo
         endif
         qn(i,k) = dpq / dpp
       else
         print *,' Error in vertical_cell_ppm_intp for no lev found '
         print *,' i kh kl ',i,kh,kl
         print *,' pn ',(pn(i,kk),kk=1,levs+1)
         print *,' pp ',(pp(i,kk),kk=1,levs+1)
         call abort
       endif
! next one
       kkl = kkh
     enddo
!
   enddo
!
   return
   end subroutine vertical_cell_ppm_intp
!-------------------------------------------------------------------------------
#endif /* NISLQ end */
   end module nislq
