#include <define.h>
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!    [dyn_sph_driver] or [dyn_driver_onestep]
!             |
!             |--- [sph_nonlinear_tend] *
!             |--- [sph_nonlinear_tend_hybrid] *
!                      |--- [vcnhyb] *
!                              |--- [tridim_hyb] *
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2001-01-01  songyou hong           cvs version
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
#ifndef HYBRID
   subroutine sph_nonlinear_tend(lons2,lat,                                    &
#ifdef SMP
    dg,tg,rqg,qg,del,sl,rdel2,                                                 &
#ifdef DBG
    dtdt,drdt,vadt,adbc)
#else
    dtdt,drdt)
#endif
#else /* ~SMP */
    dg,tg,zg,ug,vg,rqg,dphi,dlam,                                              &
    rcl,del,rdel2,ci,p1,p2,h1,h2,tov,spdmax,                                   &
#ifdef CLM_CWF
    jrow,cgs,                                                                  &
#endif
    dtdf,dtdl,drdf,drdl,dudl,dvdl,dudf,dvdf,                                   &
#ifndef NISLQ
    dqdt,dtdt,drdt,dudt,dvdt)
#else
    dqdt,dtdt,drdt,dudt,dvdt,dot)
#endif
#endif /* SMP end */
!-------------------------------------------------------------------------------
   use constant, only  : cp_,rd_,rerth_,omega_
   use paramodel, only : LONF22S,LONF2S,LATG2S,levs_,levh_,levp1_,levm1_,      &
                         ntotal_,nwmass_,icloud_
#ifdef DCMIP
   use dcmip_grims, only : dcmip_rho=>rho,pwi
#endif
!-------------------------------------------------------------------------------
!
! input variables
!
   integer              ::  lons2,lat
   real                 ::  dg(LONF22S,levs_),tg(LONF22S,levs_)
   real                 ::  ug(LONF22S,levs_),vg(LONF22S,levs_)
   real                 ::  rqg(LONF22S,levs_,ntotal_),zg(LONF22S,levs_)
#ifdef SMP
#ifdef DBG
   real                 ::  vadt(LONF22S,levs_),adbc(LONF22S,levs_)
#endif
   real                 ::  qg(LONF22S), sl(levs_)
#else
   real                 ::  dphi(LONF22S),dlam(LONF22S)
#endif
!
#ifndef SMP
   real                 ::  dtdf(LONF22S,levs_),dtdl(LONF22S,levs_)
   real                 ::  drdf(LONF22S,levs_,ntotal_)
   real                 ::  drdl(LONF22S,levs_,ntotal_)
   real                 ::  dudl(LONF22S,levs_),dvdl(LONF22S,levs_)
   real                 ::  dudf(LONF22S,levs_),dvdf(LONF22S,levs_)
#endif
!
! output variables
!
   real                 ::  spdmax(levs_)
#ifndef SMP
   real                 ::  dudt(LONF22S,levs_),dvdt(LONF22S,levs_)
#endif
   real                 ::  dtdt(LONF22S,levs_),drdt(LONF22S,levs_,ntotal_)
#ifndef SMP
   real                 ::  dqdt(LONF22S)
#endif
#ifdef CLM_CWF
   integer              ::  jrow
   real                 ::  cgs(LONF2S*LATG2S,levs_)
   real                 ::  qcon(LONF2S,levs_), ucon(LONF2S,levs_)
#endif
!
! constant arrays
!
   real                 ::  del(levs_),rdel2(levs_)
   real                 ::  ci(levp1_),tov(levs_)
#ifndef SMP
   real                 ::  p1(levs_),p2(levs_),h1(levs_),h2(levs_)
#endif
!
! local variables
!
   real                 ::  cg (LONF22S,levs_), db(LONF22S,levs_)
   real                 ::  dot(LONF22S,levp1_),dup(LONF22S,levs_)
   real                 ::  dum(LONF22S,levs_ ),dvm(LONF22S,levs_)
   real                 ::  rmu(levs_ ),rnu(levs_),rho(levs_),si(levp1_)
   real                 ::  cb(LONF22S,levs_),dvp(LONF22S,levs_)
   real                 ::  ek(LONF22S,levs_)
!-------------------------------------------------------------------------------
   rk= rd_ /cp_
#ifndef SMP
   sinra=sqrt(1.-1./rcl)
   fnor=2.*omega_*sinra
   fsou=-fnor
   sinra=sinra/rerth_
#endif
!
   si(1)=1.0
   do k = 1,levs_
     si(k+1)=si(k)-del(k)
   enddo
!
#ifdef DCMIP
   do k = 1,levs_
     rho(k)=dcmip_rho(1,k,lat)
   enddo
#else
   do k = 1,levm1_
     rho(k)=alog(si(k)/si(k+1))
   enddo
   rho(levs_)=0.
#endif
!
   do k = 1,levs_
     rmu(k)=1.-si(k+1)*rho(k)/del(k)
   enddo
!
   do k = 1,levm1_
     rnu(k+1)=-1.+si(k)*rho(k)/del(k)
   enddo
   rnu(1)=0.
!
#ifndef SMP
   do k = 1,levs_
     spdmax(k)=0.
   enddo
   rcl2=.5e0*rcl
!
   do k = 1,levs_
     do j = 1,lons2
       ek(j,k)=(ug(j,k)*ug(j,k)+vg(j,k)*vg(j,k))*rcl
     enddo
   enddo
!
   do k = 1,levs_
     do j = 1,lons2
       if (ek(j,k) .gt. spdmax(k))  spdmax(k)=ek(j,k)
     enddo
   enddo
#ifdef CLM_CWF
   do n=icloud_,nwmass_
     do k = 1,levs_
       do i = 1,lons2
         ucon(i,k)=ug(i,k)*drdl(i,k,n)+vg(i,k)*drdf(i,k,n)
         qcon(i,k)=rqg(i,k,n)*dudl(i,k)+rqg(i,k,n)*dvdf(i,k)
       enddo
     enddo
     do k = 1,levs_
       do i = 1,lons2
         ij=(jrow-1)*LONF2S+i
         cgs(ij,k)=-(qcon(i,k)+ucon(i,k))
       enddo
     enddo
   enddo
#endif
!
!     compute c=v(true)*del(ln(ps)).divide by cos for del, cos for v
!
   do j = 1,lons2
     dphi(j)=dphi(j)*rcl
     dlam(j)=dlam(j)*rcl
   enddo
   do k = 1,levs_
     do j = 1,lons2
       cg(j,k)=ug(j,k)*dlam(j)+vg(j,k)*dphi(j)
     enddo
   enddo
!
   do j = 1,lons2
     db(j,1)=del(1)*dg(j,1)
     cb(j,1)=del(1)*cg(j,1)
   enddo
   do k = 1,levm1_
     do j = 1,lons2
       db(j,k+1)=db(j,k)+del(k+1)*dg(j,k+1)
       cb(j,k+1)=cb(j,k)+del(k+1)*cg(j,k+1)
     enddo
   enddo
!
!   store integral of cg in dlax
!
   do j = 1,lons2
     dqdt(j)= -cb(j,levs_)
   enddo
!
!   sigma dot computed only at interior interfaces.
!
   do j = 1,lons2
     dot(j,1)=0.e0
     dvm(j,1)=0.e0
     dum(j,1)=0.e0
     dot(j,levp1_)=0.e0
     dvp(j,levs_ )=0.e0
     dup(j,levs_ )=0.e0
   enddo

#ifdef DCMIP
   do k = 2,levs_
     do i = 1,lons2
       dot(i,k)=pwi(i,k,lat)*1.e-3/(exp(qg(j)-ptop)) * (-1.)
     enddo
   enddo
#else
   do k = 1,levm1_
     do j = 1,lons2
       dot(j,k+1)=dot(j,k)+                                                    &
                    del(k)*(db(j,levs_)+cb(j,levs_)-dg(j,k)-cg(j,k))
     enddo
   enddo
#endif
!
   do k = 1,levm1_
     do j = 1,lons2
       dvp(j,k  )=vg(j,k+1)-vg(j,k)
       dup(j,k  )=ug(j,k+1)-ug(j,k)
       dvm(j,k+1)=vg(j,k+1)-vg(j,k)
       dum(j,k+1)=ug(j,k+1)-ug(j,k)
     enddo
   enddo
!
   do j = 1,lons2
     dphi(j)=dphi(j)/rcl
     dlam(j)=dlam(j)/rcl
   enddo
!
   do k = 1,levs_
     do j = 1,lons2
       dudt(j,k)=dudt(j,k)-ug(j,k)*dudl(j,k)-vg(j,k)*dudf(j,k)                 &
                  -rdel2(k)*(dot(j,k+1)*dup(j,k)+dot(j,k)*dum(j,k))            &
                  -rd_*tg(j,k)*dlam(j)
!
       dvdt(j,k)=dvdt(j,k)-ug(j,k)*dvdl(j,k)-vg(j,k)*dvdf(j,k)                 &
                  -rdel2(k)*(dot(j,k+1)*dvp(j,k)+dot(j,k)*dvm(j,k))            &
                  -rd_*tg(j,k)*dphi(j)
     enddo
   enddo
!
   lons=lons2/2
   do k = 1,levs_
     do j = 1,lons
       dudt(j,k)=dudt(j,k)+vg(j,k)*fnor
       dudt(j+lons,k)=dudt(j+lons,k)+vg(j+lons,k)*fsou
!
       dvdt(j,k)=dvdt(j,k)-ug(j,k)*fnor-sinra*ek(j,k)
       dvdt(j+lons,k)=dvdt(j+lons,k)-ug(j+lons,k)*fsou+sinra*ek(j+lons,k)
     enddo
   enddo
!
   do k = 1,levs_
     do j = 1,lons2
       dudt(j,k)=dudt(j,k)*rcl
       dvdt(j,k)=dvdt(j,k)*rcl
     enddo
   enddo
!
#endif
!
#ifdef SMP
   do j = 1,lons2
     dot(j,1)=0.e0
     dot(j,levp1_)=0.e0
     do k = 1,levs_-1
       psfc = exp(qg(j))  ! cb
       dot(j,k+1) = (dg(j,k)+dg(j,k+1))/(2.*psfc) ! cb/sec / cb
     enddo
!
! ... dT/dt ...
!
     k = 1
     difm = 0.0
     sdotm = 0.0
     difp=tg(j,k+1) - tg(j,k)
     sdotp = dot(j,k+1)
     delk = si(k+1) - 1.0
     vadvT = -1.*(sdotp*difp + sdotm*difm)/delk
     pres_k = sl(k) * exp(qg(j))
     omega_k = dg(j,k)
     adiabatic = rk*tg(j,k)*omega_k/pres_k
     dtdt(j,k) = vadvT + adiabatic
#ifdef DBG
     vadt(j,k) = vadvT
     adbc(j,k) = adiabatic
     if (j.eq.1) print 1967, j,k,delk,sdotp,difp,sdotm,difm,                   &
                vadvT,pres_k,omega_k,adiabatic,dtdt(j,k)
#endif
!
     do k = 2,levs_-1
       difm = tg(j,k) - tg(j,k-1)
       sdotm = dot(j,k-1)
       difp = tg(j,k+1) - tg(j,k)
       sdotp = dot(j,k+1)
       delk = si(k+1) - si(k-1)
       vadvT = -1.*(sdotp*difp + sdotm*difm)/delk
       pres_k = sl(k) * exp(qg(j))
       omega_k = dg(j,k)
       adiabatic = rk*tg(j,k)*omega_k/pres_k
       dtdt(j,k) = vadvT + adiabatic
#ifdef DBG
       vadt(j,k) = vadvT
       adbc(j,k) = adiabatic
       if (j.eq.1) print 1967, j,k,delk,sdotp,difp,sdotm,difm,               &
                  vadvT,pres_k,omega_k,adiabatic,dtdt(j,k)
1967   format(2i5,10e13.5)
#endif
     enddo
!
     k = levs_
     difm = tg(j,k) - tg(j,k-1)
     sdotm = dot(j,k-1)
     difp =0.0
     sdotp = 0.0
     delk = 0.0 - si(k-1)
     vadvT = -1.*(sdotp*difp + sdotm*difm)/delk
     pres_k = sl(k) * exp(qg(j))
     omega_k = dg(j,k)
     adiabatic = rk*tg(j,k)*omega_k/pres_k
#ifdef DBG
     vadt(j,k) = vadvT
     adbc(j,k) = adiabatic
     dtdt(j,k) = vadvT + adiabatic
     if (j.eq.1) print 1967, j,k,delk,sdotp,difp,sdotm,difm,                   &
                vadvT,pres_k,omega_k,adiabatic,dtdt(j,k)
#endif
!
! ... dr/dt ...
!
     do n = 1, ntotal_
       k = 1
       difm = 0.0
       sdotm = 0.0
       difp = rqg(j,k+1,n) - rqg(j,k,n)
       sdotp = dot(j,k+1)
       delk = si(k+1) - 1.0
       drdt(j,k,n) = -1.*(sdotp*difp + sdotm*difm)/delk
!
       do k = 2, levs_-1
         difm = rqg(j,k,n) - rqg(j,k-1,n)
         sdotm = dot(j,k-1)
         difp = rqg(j,k+1,n) - rqg(j,k,n)
         sdotp = dot(j,k+1)
         delk = si(k+1) - si(k-1)
         drdt(j,k,n) = -1.*(sdotp*difp + sdotm*difm)/delk
#ifdef DBG
         if (j.eq.1) print 1966, j,k,delk,sdotp,difp,sdotm,difm,drdt(j,k,n)
  1966   format(2i5,6e13.5)
#endif
       enddo
!
       k = levs_
       difm = rqg(j,k,n) - rqg(j,k-1,n)
       sdotm = dot(j,k-1)
       difp = 0.0
       sdotp = 0.0
       delk = 0.0 - si(k-1)
       drdt(j,k,n) = -1.*(sdotp*difp + sdotm*difm)/delk
     enddo
   enddo
#else			/* SMP */
   do k = 1,levm1_
     do j = 1,lons2
       dup(j,k  )=tg(j,k+1)+tov(k+1)-tg(j,k)-tov(k)+2.*rk*rnu(k+1)*            &
                  (tg(j,k)+tov(k))
       dum(j,k+1)=tg(j,k+1)+tov(k+1)-tg(j,k)-tov(k)+2.*rk*rmu(k+1)*            &
                  (tg(j,k+1)+tov(k+1))
     enddo
   enddo
!
   do k = 1,levs_
     do j = 1,lons2
       dtdt(j,k)=dtdt(j,k)-ug(j,k)*dtdl(j,k)-vg(j,k)*dtdf(j,k)                 &
                  -rdel2(k)*(dot(j,k+1)*dup(j,k)+dot(j,k)*dum(j,k))
     enddo
   enddo
!
   do k = 1,levs_
     do j = 1,lons2
       dtdt(j,k)=dtdt(j,k)+rk*(tov(k)+tg(j,k))                                 &
                  *(cg(j,k)-cb(j,levs_)-db(j,levs_))
     enddo
   enddo
!
#ifndef NISLQ
   do n=1,ntotal_
     do k = 1,levm1_
       do j = 1,lons2
         dup(j,k  )=rqg(j,k+1,n)-rqg(j,k,n)
         dum(j,k+1)=rqg(j,k+1,n)-rqg(j,k,n)
       enddo
     enddo
     do j = 1,lons2
       dup(j,levs_)=0.e0
     enddo
     do k = 1,levs_
       do j = 1,lons2
         drdt(j,k,n)=drdt(j,k,n)-ug(j,k)*drdl(j,k,n)-vg(j,k)*drdf(j,k,n)       &
                        -rdel2(k)*(dot(j,k+1)*dup(j,k)+dot(j,k)*dum(j,k))
       enddo
     enddo
   enddo
#endif
#endif
!
   return
   end subroutine sph_nonlinear_tend
#endif /* ~HYBRID end */
!
#ifdef HYBRID
!------------------------------------------------------------------------------
   subroutine sph_nonlinear_tend_hybrid(lons_lat,lat,                          &
                        ak5,bk5,                                               &
                        dg,tg,zg,ug,vg,rqg,dphi,dlam,qg,                       &
                        rcl,spdmax,deltim,nvcn,xvcn,                           &
                        dtdf,dtdl,drdf,drdl,dudl,dvdl,dudf,dvdf,               &
#ifndef NISLQ
                        dqdt,dtdt,drdt,dudt,dvdt)
#else
                        dqdt,dtdt,drdt,dudt,dvdt,pdot)
#endif
!------------------------------------------------------------------------------
   use paramodel, only : LONF22S,levm1_,levp1_,levs_,lonf_,ntotal_
   use constant,  only : cp_,omega_,rd_,rerth_,akapa_,cvap_
   use comio   ,  only : iope
#ifdef MP
   use commpi
#endif
#ifdef DCMIP
   use dcmip_grims, only : icase,pwi
#endif
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
!
! add calculation of spdlat as in gfidi 
!
   integer  ::  lons_lat
   integer  ::  j,k,n,nvcn,ifirst,lat
   real     ::  coriol,rcl,sinra,deltim,xvcn,sinlat
   real     ::  dg(LONF22S,levs_), tg(LONF22S,levs_),  zg(LONF22S,levs_),     &
                ug(LONF22S,levs_), vg(LONF22S,levs_),                         &
                rqg(LONF22S,levs_,ntotal_),                                   &
                dphi(LONF22S), dlam(LONF22S), qg(LONF22S)
   real     ::  dtdf(LONF22S,levs_),       dtdl(LONF22S,levs_),               &
                dudf(LONF22S,levs_),       dudl(LONF22S,levs_),               &
                dvdf(LONF22S,levs_),       dvdl(LONF22S,levs_),               &
                drdf(LONF22S,levs_,ntotal_), drdl(LONF22S,levs_,ntotal_)
   real     ::  dudt(LONF22S,levs_),       dvdt(LONF22S,levs_),               &  
                dqdt(LONF22S),                                                &
                dtdt(LONF22S,levs_),                                          &
                drdt(LONF22S,levs_,ntotal_), spdmax(levs_)
!
   real     ::  pk5(lons_lat,levp1_), dpk(lons_lat,levs_)
   real     ::  zadv(lons_lat,levs_,3+ntotal_)
   real     ::  dot(lons_lat,levp1_), dotinv(lons_lat,levp1_),                &
                ek(lons_lat,levs_),      cg(lons_lat,levs_),                  &
                cb(lons_lat,levs_),      db(lons_lat,levs_),                  &
              zlam(lons_lat,levs_),    zphi(lons_lat,levs_),                  &
             worka(lons_lat,levs_),   workb(lons_lat,levs_),                  &
             workc(lons_lat,levs_),                                           &
              phiu(lons_lat,levs_),    phiv(lons_lat,levs_),                  &
              uprs(lons_lat,levs_),    vprs(lons_lat,levs_),                  &
              cofa(lons_lat,levs_),    cofb(lons_lat,levs_),                  &
              alfa(lons_lat,levs_),    rlnp(lons_lat,levs_),                  &
              px1u(lons_lat,levs_),    px1v(lons_lat,levs_),                  &
              px2u(lons_lat,levs_),    px2v(lons_lat,levs_),                  &
               px2(lons_lat,levs_),                                           &
              px3u(lons_lat,levs_),    px3v(lons_lat,levs_),                  &
              px4u(lons_lat,levs_),    px4v(lons_lat,levs_),                  &
              px5u(lons_lat,levs_),    px5v(lons_lat,levs_),                  &
              uphi(lons_lat,levs_),    vphi(lons_lat,levs_),                  &
              expq(lons_lat),                                                 &
              rdel(lons_lat,levs_),   rdel2(lons_lat,levs_),                  &
            sumdel(lons_lat),          del(lons_lat,levs_),                   &
                si(lons_lat,levp1_),     sl(lons_lat,levs_),                  &
                rk1,rkr,                 dif(lons_lat)
#if defined(NISLQ) || defined(DCMIP)
   real     ::  pdot(LONF22S,levp1_)
   real     ::  cgi(lons_lat,levp1_)
#endif
   real     ::  cons0,cons0p5,cons1,cons2,clog2   !constant
   real     ::  rmin,rmax,delta,delta1
   save clog2,ifirst,delta,delta1
   data ifirst /1/
! 
   real     ::  ak5(levp1_),bk5(levp1_),dbk(levs_),bkl(levs_),ck(levs_),diff
   integer  ::   lons,i
!------------------------------------------------------------------------------
   cons0   = 0.d0      !constant
   cons0p5 = 0.5d0     !constant
   cons1   = 1.d0      !constant
   cons2   = 2.d0      !constant
!
   sinra=sqrt(cons1-cons1/rcl)        !constant
   coriol=cons2*omega_*sinra          !constant
   sinra=sinra/rerth_
!------------------------------------------------------------------------------
!
! rcl = cons1/(cons1-sinlat*sinlat)  !constant
! 
   if(ifirst.eq.1)then
     clog2=log(cons2)     ! constant
     delta=cvap_/cp_      ! check these cpv cpd (at const p for vapor and dry
     delta1=delta-cons1
     rk1 = akapa_ + 1.e0
     rkr = 1.0/akapa_
     ifirst=0
   endif
!     
   do k = 1,levs_
     dbk(k) = bk5(k+1)-bk5(k)
     bkl(k) = (bk5(k+1)+bk5(k))*0.5
     ck(k)  = ak5(k+1)*bk5(k)-ak5(k)*bk5(k+1)
   enddo
!
   do j = 1,lons_lat
     expq(j)=exp(qg(j))
   enddo
!------------------------------------------------------------------------------
!
! get vertical coordinate for vcnhyb going bot. to top.
!
   do k = 1,levp1_
     do j = 1,lons_lat
       si(j,levs_+2-k)= ak5(k)+bk5(k)*expq(j) !ak(k) bk(k) go top to bottom
     enddo
   enddo
! 
   do  k=1,levs_
     do j = 1,lons_lat
       sl(j,k) = cons0p5*(si(j,k)+si(j,k+1))
     enddo
   enddo
!------------------------------------------------------------------------------
   do k = 1,levp1_
     do j = 1,lons_lat
       pk5(j,k)=ak5(k) + bk5(k)*expq(j)
     enddo
   enddo
! 
   do k = 1,levs_
     do j = 1,lons_lat
       dpk(j,k)=    pk5(j,k+1) - pk5(j,k)
       rdel(j,k)=    cons1/dpk(j,k)            ! constant
       rdel2(j,k)=  cons0p5/dpk(j,k)            ! constant
     enddo
   enddo
! 
   k=1
   do j = 1,lons_lat
     alfa(j,1)=clog2                          ! constant
   enddo
! 
   do j = 1,lons_lat
     rlnp(j,1)= 99999.99
   enddo
!
   do  k=2,levs_
     do j = 1,lons_lat
       rlnp(j,k)= log( pk5(j,k+1)/pk5(j,k) )
       alfa(j,k)= cons1-( pk5(j,k)/dpk(j,k) )*rlnp(j,k)
     enddo
   enddo
!
   spdmax=0.
   do  k=1,levs_
     do j = 1,lons_lat
       ek(j,k)=(ug(j,levp1_-k)*ug(j,levp1_-k)+                                 &
                vg(j,levp1_-k)*vg(j,levp1_-k))*rcl
       if (ek(j,k) .gt. spdmax(levp1_-k))  spdmax(levp1_-k)=ek(j,k)
     enddo
   enddo
!
   do k = 1,levs_
     do j = 1,lons_lat
       cg(j,k)=(ug(j,levp1_-k)*dlam(j)+vg(j,levp1_-k)*dphi(j))*rcl
     enddo
   enddo
#if defined(NISLQ) || defined(DCMIP)
!
! [t2b] cgi: cg*bk5 at half-level
!
   cgi(1:lons_lat,1)      =0.
   cgi(1:lons_lat,levs_+1)=0.
   do k = 2,levs_
     do j = 1,lons_lat
       cgi(j,k)=0.5*(cg(j,k)+cg(j,k-1))
     enddo
   enddo
#endif
!
   k=1
   do j = 1,lons_lat
     db(j,1)=dg(j,levs_)*dpk(j,1)
     cb(j,1)=cg(j,1)*dbk(1)
   enddo
! 
   do k = 1,levm1_
     do j = 1,lons_lat
       db(j,k+1)=db(j,k)+dg(j,levs_-k)*dpk(j,k+1)
       cb(j,k+1)=cb(j,k)+cg(j,k+1)*dbk(k+1)
     enddo
   enddo
! 
   do j = 1,lons_lat
     dqdt(j)= -db(j,levs_)/expq(j)-cb(j,levs_)
!    dqdt(j)= -cb(j,levs_)
     dot(j,    1)=cons0                    !constant
     dot(j,levp1_)=cons0                   !constant
   enddo
!
! [t2b] dot
!
#ifdef DCMIP
     do k = 2,levs_
       do j = 1,lons_lat
         if(icase.eq.11 .or. icase.eq.12) then
           dot(j,k)=pwi(j,levs_+2-k,lat)*1.e-3
         else if(icase.eq.13) then
           dot(j,k)=pwi(j,levs_+2-k,lat)*1.e-4
!          dot(j,k)=dot(i,j,k) - expq(j)*(dqdt(j)+cgi(j,k))*1.e-1
         else
           dot(j,k)=-expq(j)*(bk5(k)*dqdt(j)+cb(j,k-1)) -db(j,k-1)
         endif
       enddo
     enddo
     if(icase.eq.13) then
       do n = 1,ntotal_
         do j = 1,lons_lat
           rqg(j,1,n)=max(rqg(j,1,n),0.)
           rqg(j,2,n)=max(rqg(j,2,n),0.)
           rqg(j,3,n)=max(rqg(j,3,n),0.)
           rqg(j,4,n)=max(rqg(j,4,n),0.)
           rqg(j,5,n)=max(rqg(j,5,n),0.)
         enddo
       enddo
     endif
#else
   do k = 1,levm1_
     do j = 1,lons_lat
       dot(j,k+1)=-expq(j)*(bk5(k+1)*dqdt(j)+cb(j,k)) -db(j,k)
!      dot(j,k+1)=-expq(j)*(bk5(k+1)*(dqdt(j)-db(j,levs_)/expq(j))
!             +cb(j,k)) -db(j,k)
     enddo
   enddo
#endif /* DCMIP end */
!
! [b2t] dotinv
!
   do k = 1,levp1_
     do j = 1,lons_lat
       dotinv(j,k)=dot(j,levp1_+1-k)
     enddo
   enddo
#ifdef NISLQ
!
! [t2b] pdot: pressure velocity for nislq
!
   do k = 2,levs_
     do j = 1,lons_lat
#ifdef DCMIP
       if(icase.eq.11.or.icase.eq.12) then
         pdot(j,k)=pwi(j,levs_+2-k,lat)*1.e-3
       else if (icase.eq.13) then
         pdot(j,k)=pwi(j,levs_+2-k,lat)*1.e-4
       else
         !pdot(j,k)=dot(j,k)+expq(j)*bk5(k)*(dqdt(j)+cgi(j,k))
         pdot(j,k)=dot(j,k)
       endif
#else
       !pdot(j,k)=dot(j,k)+expq(j)*bk5(k)*(dqdt(j)+cgi(j,k))
       pdot(j,k)=dot(j,k)
#endif
     enddo
   enddo
   pdot(1:LONF22S,1)      =0.
   pdot(1:LONF22S,levs_+1)=0.
#endif /* NISLQ end */
!
! variables are in bottom to top order  !!!!!!!!!!!!!!!!!
! do horizontal advection.
! 
   k=1
   do j = 1,lons_lat
     dudt(j,levp1_-k)=-ug(j,levp1_-k)*dudl(j,levp1_-k)                         & 
                      -vg(j,levp1_-k)*dudf(j,levp1_-k)
     dvdt(j,levp1_-k)=-ug(j,levp1_-k)*dvdl(j,levp1_-k)                         &
                      -vg(j,levp1_-k)*dvdf(j,levp1_-k)
     dtdt(j,levp1_-k)=-ug(j,levp1_-k)*dtdl(j,levp1_-k)                         &
                      -vg(j,levp1_-k)*dtdf(j,levp1_-k)
   enddo
! 
   k=levs_
   do j = 1,lons_lat
     dudt(j,levp1_-k)=-ug(j,levp1_-k)*dudl(j,levp1_-k)                         &
                      -vg(j,levp1_-k)*dudf(j,levp1_-k)
     dvdt(j,levp1_-k)=-ug(j,levp1_-k)*dvdl(j,levp1_-k)                         &
                      -vg(j,levp1_-k)*dvdf(j,levp1_-k)
     dtdt(j,levp1_-k)=-ug(j,levp1_-k)*dtdl(j,levp1_-k)                         &
                      -vg(j,levp1_-k)*dtdf(j,levp1_-k)
   enddo
! 
   do k = 2,levm1_
     do j = 1,lons_lat
       dudt(j,levp1_-k)=-ug(j,levp1_-k)*dudl(j,levp1_-k)                       &
                        -vg(j,levp1_-k)*dudf(j,levp1_-k)
       dvdt(j,levp1_-k)=-ug(j,levp1_-k)*dvdl(j,levp1_-k)                       &
                        -vg(j,levp1_-k)*dvdf(j,levp1_-k)
       dtdt(j,levp1_-k)=-ug(j,levp1_-k)*dtdl(j,levp1_-k)                       &
                        -vg(j,levp1_-k)*dtdf(j,levp1_-k)
     enddo
   enddo
!
!  if (mype.eq.master) then
!    do k = 1,1
!      do i = 1,5
!        print'(a4,5e12.4,2i5)','dt1',                                         &
!            dudt(i,k),dvdt(i,k),dtdt(i,k),drdt(i,k,1),dqdt(i),i,lat
!      enddo
!    enddo
!  endif
!
! add coriolis,deformation note  coriolis sign for s.hemi
!
   lons=lons_lat/2
   do k = 1,levs_
     do j = 1,lons
       dudt(j,levp1_-k)=dudt(j,levp1_-k)                                       &
             +vg(j,levp1_-k)*coriol
 
       dudt(j+lons,levp1_-k)=dudt(j+lons,levp1_-k)                             & 
             -vg(j+lons,levp1_-k)*coriol
 
       dvdt(j,levp1_-k)=dvdt(j,levp1_-k)                                       &
             -ug(j,levp1_-k)*coriol                                            &
             -sinra*ek(j,k)

       dvdt(j+lons,levp1_-k)=dvdt(j+lons,levp1_-k)                             &
             +ug(j+lons,levp1_-k)*coriol                                       &
             +sinra*ek(j+lons,k)
     enddo
   enddo
!
!  if (mype.eq.master) then
!    do k = 1,1
!      do i = 1,5
!        print'(a4,5e12.4,2i5)','dt2',                                        &
!            dudt(i,k),dvdt(i,k),dtdt(i,k),drdt(i,k,1),dqdt(i),i,lat
!      enddo
!    enddo
!  endif
!------------------------------------------------------------------------------
!
! calculate pressure force:
!
   k=1
   do j = 1,lons_lat
     cofb(j,k)=-rdel(j,k)*(alfa(j,k)*dbk(k))
   enddo
! 
   do k = 2,levs_
     do j = 1,lons_lat
       cofb(j,k)=-rdel(j,k)*(bk5(k)*rlnp(j,k)+alfa(j,k)*dbk(k))
     enddo
   enddo
! 
   do k = 1,levs_
     do j = 1,lons_lat
       uprs(j,k)=cofb(j,k)*rd_*tg(j,levp1_-k)*expq(j)*dlam(j)
       vprs(j,k)=cofb(j,k)*rd_*tg(j,levp1_-k)*expq(j)*dphi(j)
     enddo
   enddo
!
   do k = 1,levs_
     do j = 1,lons_lat
       cofa(j,k)=-rdel(j,k)*(                                                  &
         bk5(k+1)*pk5(j,k)/pk5(j,k+1) - bk5(k)                                 &
        +rlnp(j,k)*( bk5(k)-pk5(j,k)*dbk(k)*rdel(j,k) )  )
     enddo
   enddo
!-------------------------------------------------------------------------------
   do k = 1,levs_
      do j = 1,lons_lat
        px1u(j,k)=cons0              ! grid topography =0 for testing
        px1v(j,k)=cons0              ! grid tpopgraphy =0 for testing
      enddo
   enddo
!
!  see programming notes for looping in  calculating  px2 and px3
!
   do j = 1,lons_lat
     px2(j,levs_)=cons0                             ! constant
     px2(j,levm1_)=                                                            &
        -rd_*( bk5(levp1_)/pk5(j,levp1_)-bk5(levs_)/pk5(j,levs_) )             &
             *tg(j,1)
   enddo
! 
   do k = 2,levm1_
     do j = 1,lons_lat
     px2(j,levs_-k)=px2(j,levp1_-k)                                            &
     -rd_*(bk5(levs_+2-k)/                                                     &
              pk5(j,levs_+2-k)-bk5(levp1_-k)/pk5(j,levp1_-k))*                 &
                                                        tg(j,k)
     enddo
   enddo
! 
   do k = 1,levs_
     do j = 1,lons_lat
       px2u(j,k)=px2(j,k)*expq(j)*dlam(j)
       px2v(j,k)=px2(j,k)*expq(j)*dphi(j)
     enddo
   enddo
!
   do j = 1,lons_lat
     px3u(j,levs_)=cons0 ! constant
     px3v(j,levs_)=cons0 ! constant
   enddo
! 
   do j = 1,lons_lat
     px3u(j,levm1_)=-rd_*rlnp(j,levs_)*dtdl(j,1)
     px3v(j,levm1_)=-rd_*rlnp(j,levs_)*dtdf(j,1)
   enddo
! 
   do k = 2,levm1_
     do j = 1,lons_lat
       px3u(j,levs_-k)=px3u(j,levp1_-k)-rd_*rlnp(j,levp1_-k)*dtdl(j,k)
       px3v(j,levs_-k)=px3v(j,levp1_-k)-rd_*rlnp(j,levp1_-k)*dtdf(j,k)
     enddo
   enddo
!
   do k = 1,levs_
#ifdef LINUX_PGI
!pgi$l novector
#endif
     do j = 1,lons_lat
       px3u(j,k)=px3u(j,k)/rcl
       px3v(j,k)=px3v(j,k)/rcl
     enddo
   enddo
!
   do k = 1,levs_
#ifdef LINUX_PGI
!pgi$l novector
#endif
     do j = 1,lons_lat
       px4u(j,k)=-rd_*alfa(j,k)*dtdl(j,levp1_-k)/rcl
       px4v(j,k)=-rd_*alfa(j,k)*dtdf(j,levp1_-k)/rcl
     enddo
   enddo
! 
   do k = 1,levs_
     do j = 1,lons_lat
       px5u(j,k)=-cofa(j,k)*rd_*tg(j,levp1_-k)*expq(j)*dlam(j)
       px5v(j,k)=-cofa(j,k)*rd_*tg(j,levp1_-k)*expq(j)*dphi(j)
     enddo
   enddo
! 
   do k = 1,levs_
     do j = 1,lons_lat
       uphi(j,k)=px1u(j,k)+px2u(j,k)+px3u(j,k)+px4u(j,k)+px5u(j,k)
       vphi(j,k)=px1v(j,k)+px2v(j,k)+px3v(j,k)+px4v(j,k)+px5v(j,k)
     enddo
   enddo
!
   do k = 1,levs_
     do j = 1,lons_lat
       dudt(j,levp1_-k)=dudt(j,levp1_-k)+uphi(j,k)+uprs(j,k)
       dvdt(j,levp1_-k)=dvdt(j,levp1_-k)+vphi(j,k)+vprs(j,k)
     enddo
   enddo
!------------------------------------------------------------------------------
!  if (mype.eq.master) then
!    do k = 1,1
!      do i = 1,5
!        print'(a4,5e12.4,2i5)','dt3',                                         &
!              dudt(i,k),dvdt(i,k),dtdt(i,k),drdt(i,k,1),dqdt(i),i,k           &
!              uphi(i,k),uprs(i,k),vphi(i,k),vprs(i,k),dqdt(i),i,k             &
!              px2(i,k),dlam(i),rlnp(i,k),dtdl(i,k),tg(i,k),i,lat
!      enddo
!    enddo
!  endif
!
   do k = 1,levs_
     do j = 1,lons_lat
       worka(j,k)=akapa_*tg(j,levp1_-k)/(cons1+(delta1)*                       &
         rqg(j,levp1_-k,1)) * rdel(j,k)
     enddo
   enddo
! 
   k=1
   do j = 1,lons_lat
     workb(j,1)=                                                               &
     alfa(j,1)*( dg(j,levs_)*dpk(j,1)+expq(j)*cb(j,1)*dbk(1) )
   enddo
! 
   do k = 2,levs_
     do j = 1,lons_lat
       workb(j,k)=rlnp(j,k)*( db(j,k-1)+expq(j)*cb(j,k-1) )                    &
       +alfa(j,k)*( dg(j,levp1_-k)*dpk(j,k)+expq(j)*cg(j,k)*dbk(k) )
     enddo
   enddo
! 
   k=1
   do j = 1,lons_lat
     workc(j,1)=expq(j)*cg(j,1)*dbk(1)
   enddo
! 
   do k = 2,levs_
     do j = 1,lons_lat
        workc(j,k)=expq(j)*cg(j,k)*( dbk(k)+ck(k)*rlnp(j,k)*rdel(j,k) )
     enddo
   enddo
! 
   do k = 1,levs_
     do j = 1,lons_lat
       dtdt(j,levp1_-k)=                                                       &
       dtdt(j,levp1_-k)+worka(j,k)*( -workb(j,k) + workc(j,k))
     enddo
   enddo
!
#ifndef NISLQ
   do 330 n=1,ntotal_
     k=1
     do j = 1,lons_lat
       drdt(j,levp1_-k,n)=-ug(j,levp1_-k)*drdl(j,levp1_-k,n)                   &
                           -vg(j,levp1_-k)*drdf(j,levp1_-k,n)
     enddo
     k=levs_
     do j = 1,lons_lat
       drdt(j,levp1_-k,n)=-ug(j,levp1_-k)*drdl(j,levp1_-k,n)                   & 
                           -vg(j,levp1_-k)*drdf(j,levp1_-k,n)
     enddo
     do k = 2,levm1_
       do j = 1,lons_lat
          drdt(j,levp1_-k,n)=-ug(j,levp1_-k)*drdl(j,levp1_-k,n)                &
                           -vg(j,levp1_-k)*drdf(j,levp1_-k,n)
       enddo
     enddo
330  continue
#endif
! 
! do vertical advection
!
   k=1
   do j = 1,lons_lat
     zadv(j,levp1_-k,1)=                                                       &
       -rdel2(j,k)*dot(j,k+1)*( ug(j,levs_-k)-ug(j,levp1_-k))
! 
     zadv(j,levp1_-k,2)=                                                       &
       -rdel2(j,k)*dot(j,k+1)*( vg(j,levs_-k)-vg(j,levp1_-k))
! 
     zadv(j,levp1_-k,3)=                                                       &
       -rdel2(j,k)*dot(j,k+1)*( tg(j,levs_-k)-tg(j,levp1_-k))
   enddo
! 
   k=levs_
   do j = 1,lons_lat
     zadv(j,levp1_-k,1)=                                                       &
       -rdel2(j,k)*dot(j,k)*( ug(j,levp1_-k)-ug(j,levs_+2-k) )
! 
     zadv(j,levp1_-k,2)=                                                       &
       -rdel2(j,k)*dot(j,k)*( vg(j,levp1_-k)-vg(j,levs_+2-k) )
! 
     zadv(j,levp1_-k,3)=                                                       &
       -rdel2(j,k)*dot(j,k)*( tg(j,levp1_-k)-tg(j,levs_+2-k) )
   enddo
! 
   do k = 2,levm1_
      do j = 1,lons_lat
        zadv(j,levp1_-k,1)=                                                    &
        -rdel2(j,k)*( dot(j,k+1)*( ug(j,levs_  -k)-ug(j,levp1_-k) ) +          &
                      dot(j,k  )*( ug(j,levp1_-k)-ug(j,levs_+2-k) ) )
 !
        zadv(j,levp1_-k,2)=                                                    &
        -rdel2(j,k)*( dot(j,k+1)*( vg(j,levs_  -k)-vg(j,levp1_-k) ) +          &
                      dot(j,k  )*( vg(j,levp1_-k)-vg(j,levs_+2-k) ) )
! 
        zadv(j,levp1_-k,3)=                                                    &
        -rdel2(j,k)*( dot(j,k+1)*( tg(j,levs_  -k)-tg(j,levp1_-k) ) +          &
                      dot(j,k  )*( tg(j,levp1_-k)-tg(j,levs_+2-k) ) )
      enddo
   enddo
!
!-------------------------------------------------------------------------------
!
   do 340 n=1,ntotal_
     k=1
     do j = 1,lons_lat
       zadv(j,levp1_-k,3+n)=                                                   &
       -rdel2(j,k)*dot(j,k+1)*( rqg(j,levs_-k,n)-rqg(j,levp1_-k,n) )
     enddo
! 
     k=levs_
     do j = 1,lons_lat
       zadv(j,levp1_-k,3+n)=                                                   &
       -rdel2(j,k)*dot(j,k)*( rqg(j,levp1_-k,n)-rqg(j,levs_+2-k,n) )
     enddo
! 
     do k = 2,levm1_
       do j = 1,lons_lat
         zadv(j,levp1_-k,3+n)=-rdel2(j,k)*                                     &
         ( dot(j,k+1)*( rqg(j,levs_  -k,n)-rqg(j,levp1_-k,n) ) +               &
           dot(j,k  )*( rqg(j,levp1_-k,n)-rqg(j,levs_+2-k,n) ) )
       enddo
     enddo
340  continue
!
!    if (mype.eq.master) then
!      do k = 1,1
!        do i = 1,5
!          print'(a5,4e12.4,2i5)','zad1',                                      &
!                zadv(i,k,1),zadv(i,k,2),zadv(i,k,3),zadv(i,k,4),i,lat
!        enddo
!      enddo
!    endif
!-------------------------------------------------------------------------------
!
#ifndef DCMIP
#ifdef NISLQ
   call vcnhyb(lons_lat,levs_,3        ,deltim,si,sl,dotinv,zadv,nvcn,xvcn)
#else
   call vcnhyb(lons_lat,levs_,3+ntotal_,deltim,si,sl,dotinv,zadv,nvcn,xvcn)
#endif
#endif /* ~DCMIP end */
!
!sela if(xvcn.ne.0.) print*,'xvcn=',xvcn,' nvcn=',nvcn
!-------------------------------------------------------------------------------
! add vertical filterd advection
!
!  if (mype.eq.master) then
!    do k = 1,1
!      do i = 1,5
!        print'(a5,4e12.4,2i5)','zad2',                                        &
!             zadv(i,k,1),zadv(i,k,2),zadv(i,k,3),zadv(i,k,4),i,lat
!      enddo
!    enddo
!  endif
!
   do k = 1,levs_
     do j = 1,lons_lat
       dudt(j,k)=dudt(j,k)+zadv(j,k,1)
       dvdt(j,k)=dvdt(j,k)+zadv(j,k,2)
       dtdt(j,k)=dtdt(j,k)+zadv(j,k,3)
     enddo
   enddo
#ifndef NISLQ
!
   do n = 1,ntotal_
     do k = 1,levs_
       do j = 1,lons_lat
         drdt(j,k,n)=drdt(j,k,n)+zadv(j,k,3+n)
       enddo
     enddo
   enddo
#endif
! 
!------------------------------------------------------------------------------
! this multiplication must be on  completed tendencies.
   do k = 1,levs_
     do j = 1,lons_lat
       dudt(j,levp1_-k)=dudt(j,levp1_-k)*rcl
       dvdt(j,levp1_-k)=dvdt(j,levp1_-k)*rcl
     enddo
   enddo
!
   return
   end subroutine sph_nonlinear_tend_hybrid
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
   subroutine vcnhyb(im,km,nm,dt,zint,zmid,zdot,zadv,nvcn,xvcn)
!------------------------------------------------------------------------------
!
! subprogram:    vcnhyb      vertical advection instability filter
!
! abstract: filters vertical advection tendencies
!   in the dynamics tendency equation in order to ensure stability
!   when the vertical velocity exceeds the cfl criterion.
!   the vertical velocity in this case is sigmadot.
!   for simple second-order centered eulerian advection,
!   filtering is needed when vcn=zdot*dt/dz>1.
!   the maximum eigenvalue of the linear advection equation
!   with second-order implicit filtering on the tendencies
!   is less than one for all resolvable wavenumbers (i.e. stable)
!   if the nondimensional filter parameter is nu=(vcn**2-1)/4.
!
! program history log:
!   1997-07-30  iredell
!   2000-01-01  hann-ming henry juang  mpi
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call vcnhyb(im,km,nm,dt,zint,zmid,zdot,zadv,nvcn,xvcn)
!
!   input argument list:
!     im       - integer number of gridpoints to filter
!     km       - integer number of vertical levels
!     nm       - integer number of fields
!     dt       - real timestep in seconds
!     zint     - real (im,km+1) interface vertical coordinate values
!     zmid     - real (im,km) midlayer vertical coordinate values
!     zdot     - real (im,km+1) vertical coordinate velocity
!     zadv     - real (im,km,nm) vertical advection tendencies
!
!   output argument list:
!     zadv     - real (im,km,nm) vertical advection tendencies
!     nvcn     - integer number of points requiring filtering
!     xvcn     - real maximum vertical courant number
!
!   subprograms called:
!     tridim_hyb   - tridiagonal matrix solver
!
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer,intent(in):: im,km,nm
   real,intent(in):: dt,zint(im,km+1),zmid(im,km),zdot(im,km+1)
   real,intent(inout):: zadv(im,km,nm)
   integer,intent(out):: nvcn
   real,intent(out):: xvcn
   integer i,j,k,n,ivcn(im)
   logical lvcn(im)
   real zdm,zda,zdb,vcn(im,km-1)
   real rnu,cm(im,km),cu(im,km-1),cl(im,km-1)
   real rr(im,km,nm)
!------------------------------------------------------------------------------
!
!  compute vertical courant number
!  increase by 10% for safety
!
   nvcn=0
   xvcn=0.
   lvcn=.false.
   do k = 1,km-1
     do i = 1,im
       zdm=zmid(i,k)-zmid(i,k+1)
       vcn(i,k)=abs(zdot(i,k+1)*dt/zdm)*1.1
       lvcn(i)=lvcn(i).or.vcn(i,k).gt.1
       xvcn=max(xvcn,vcn(i,k))
     enddo
   enddo
!------------------------------------------------------------------------------
!  determine points requiring filtering
   if(xvcn.gt.1) then
     do i = 1,im
       if(lvcn(i)) then
         ivcn(nvcn+1)=i
         nvcn=nvcn+1
       endif
     enddo
!------------------------------------------------------------------------------
!  compute tridiagonal matrim
     do j = 1,nvcn
       cm(j,1)=1
     enddo
     do k = 1,km-1
       do j = 1,nvcn
         i=ivcn(j)
         if(vcn(i,k).gt.1) then
          zdm=zmid(i,k)-zmid(i,k+1)
          zda=zint(i,k+1)-zint(i,k+2)
          zdb=zint(i,k)-zint(i,k+1)
           rnu=(vcn(i,k)**2-1)/4
           cu(j,k)=-rnu*zdm/zdb
           cl(j,k)=-rnu*zdm/zda
           cm(j,k)=cm(j,k)-cu(j,k)
           cm(j,k+1)=1-cl(j,k)
         else
           cu(j,k)=0
           cl(j,k)=0
           cm(j,k+1)=1
         endif
       enddo
     enddo
!------------------------------------------------------------------------------
!  fill fields to be filtered
     do n=1,nm
       do k = 1,km
         do j = 1,nvcn
           i=ivcn(j)
           rr(j,k,n)=zadv(i,k,n)
         enddo
       enddo
     enddo
!------------------------------------------------------------------------------
!  solve tridiagonal system
        call tridim_hyb(nvcn,im,km,km,nm,cl,cm,cu,rr,cu,rr)
!------------------------------------------------------------------------------
!  replace filtered fields
     do n=1,nm
       do k = 1,km
         do j = 1,nvcn
           i=ivcn(j)
           zadv(i,k,n)=rr(j,k,n)
         enddo
       enddo
     enddo
   endif
!
   end subroutine vcnhyb
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   subroutine tridim_hyb(l,lx,n,nx,m,cl,cm,cu,r,au,a)
!------------------------------------------------------------------------------
!
! subprogram:    tridim_hyb      solves tridiagonal matrix problems.
!
! abstract: this routine solves multiple tridiagonal matrix problems
!   with multiple right-hand-side and solution vectors for every matrix.
!   the solutions are found by eliminating off-diagonal coefficients,
!   marching first foreward then backward along the matrix diagonal.
!   the computations are vectorized around the number of matrices.
!   no checks are made for zeroes on the diagonal or singularity.
!
! program history log:
!   07-30-1997    mark iredell    development
!   04-30-2007    jung-eun kim    grims implementation
!
! usage:    call tridim_hyb(l,lx,n,nx,m,cl,cm,cu,r,au,a)
!
!   input argument list:
!     l        - integer number of tridiagonal matrices
!     lx       - integer first dimension (lx>=l)
!     n        - integer order of the matrices
!     nx       - integer second dimension (nx>=n)
!     m        - integer number of vectors for every matrix
!     cl       - real (lx,2:n) lower diagonal matrix elements
!     cm       - real (lx,n) main diagonal matrix elements
!     cu       - real (lx,n-1) upper diagonal matrix elements
!                (may be equivalent to au if no longer needed)
!     r        - real (lx,nx,m) right-hand-side vector elements
!                (may be equivalent to a if no longer needed)
!
!   output argument list:
!     au       - real (lx,n-1) work array
!     a        - real (lx,nx,m) solution vector elements
!
!------------------------------------------------------------------------------
   real  ::  cl(lx,2:n),cm(lx,n),cu(lx,n-1),r(lx,nx,m),                        &  
                                 au(lx,n-1),a(lx,nx,m)
!------------------------------------------------------------------------------
!  march up
   do i = 1,l
     fk=1./cm(i,1)
     au(i,1)=fk*cu(i,1)
   enddo
   do j = 1,m
     do i = 1,l
       fk=1./cm(i,1)
       a(i,1,j)=fk*r(i,1,j)
     enddo
   enddo
   do k = 2,n-1
     do i = 1,l
       fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
       au(i,k)=fk*cu(i,k)
     enddo
     do j = 1,m
       do i = 1,l
         fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
         a(i,k,j)=fk*(r(i,k,j)-cl(i,k)*a(i,k-1,j))
       enddo
     enddo
   enddo
!------------------------------------------------------------------------------
!  march down
   do j = 1,m
     do i = 1,l
       fk=1./(cm(i,n)-cl(i,n)*au(i,n-1))
       a(i,n,j)=fk*(r(i,n,j)-cl(i,n)*a(i,n-1,j))
     enddo
   enddo
   do k = n-1,1,-1
     do j = 1,m
       do i = 1,l
         a(i,k,j)=a(i,k,j)-au(i,k)*a(i,k+1,j)
       enddo
     enddo
   enddo
!------------------------------------------------------------------------------
   end subroutine tridim_hyb
#endif /* HYBRID end */
