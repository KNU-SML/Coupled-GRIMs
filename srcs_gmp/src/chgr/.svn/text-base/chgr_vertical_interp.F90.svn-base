#include "define.h"
   subroutine chgr_vertical_interp
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [chgr_vertical_interp]
!      |
!      |--- [chgr_vinterp_solver] *
!      |--- [chgr_sigma2sigma] *
!      |--- [chgr_moisture_extrap] *
!      |--- [chgr_moisture_extrap_old] *
!
!-------------------------------------------------------------------------------
   end subroutine chgr_vertical_interp
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_vinterp_solver(lnew,gzso                                    &
#ifdef DFS
                              ,tave                                            &
#endif
                          )
!-------------------------------------------------------------------------------
!
! subprogram:    chgr_vinterp_solver      vertically interpolate spectral data.
!
! abstract: this routine transforms spectral data to gridspace,
!           fixes moisture in the stratosphere if necessary,
!           interpolates surface pressure to the new orography,
!           vertically interpolates wind, temperature and moisture,
!           and transforms the data back to spectral space.
!
! program history log:
!   1991-03-15  mark iredell  docblock written (prehistorical program)
!   2000-01-01  song-you hong          usgs topo, kagwd options
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! usage:    call chgr_vinterp_solver(lnew)
!   input argument list:
!     lnew     - logical .true. to use newer moisture climatology
!
!   subprograms called:
!     chgr_divor2wind            - compute winds from vorticity and divergence
!     chgr_sph_funct             - compute legendre polynomials
!     chgr_sum_coeff             - transform scalar spectral to fourier
!     chgr_sum_coeff_wind        - transform vector spectral to fourier
!     sph_fft2grid               - fourier transform
!     chgr_moisture_extrap_old   - older moisture extrapolation
!     chgr_moisture_extrap       - newer moisture extrapolation
!     chgr_get_sfcp              - interpolate surface pressure
!     chgr_sigma2sigma           - interpolate upper air quantities
!     chgr_fourier_diff          - combine fourier coefficients
!     sph_coeff_lat              - transform scalar fourier to spectral
!     chgr_sph_derivative   - compute legendre derivatives
!     chgr_sph_operator   - transform vector fourier to spectral vorticity
!     chgr_div_coeff_lat   - transform vector fourier to spectral divergence
!
!-------------------------------------------------------------------------------
   use constant, only     :  pi_,g_,rd_,rv_,psat_,ttp_,hvap_,akapa_
#ifdef DFS
   use dfsvar,   only     :  mt,jl,jlg,coslat
#endif
   use paramter
   use parmchgr
   use comchgr, only      :  mdim, mdimv, jdimhf, idimt, kdimp, kdimpi,        &
                             eps, colrad, wgt, rcs2,                           &
#ifdef DFS
                             mdimd,                                            &
#else
                             pss, tts, qqs, uus, vvs,                          &
                             psa, tta, qqa, uua, vva,                          &
#endif
                             q , di,  ze,  te , rq ,  gz,                      &
                             qo, dio, zeo, teo, rqo, gzo,                      &
                             ps , tf , rqf , uf , vf,                          &
                             pso, tfo, rqfo, ufo, vfo,                         &
                             siin, slin, ak , bk,                              &
                             si  , sl  , ak5, bk5
   use module_trans, only :  dyn_trans2model_grid, dyn_trans2output_grid
!-------------------------------------------------------------------------------
!#define DBG
!#define BIN_DBG
   save
#ifdef DFS
   real                 ::  gzso(mdimd)
#else
   real                 ::  gzso(mdim)
#endif
#ifdef DFS
!
   real  ::  psa(idimt,jdimhf),tta(idimt,jdimhf,kdim),qqa(idimt,jdimhf,kdim),  &
             uua(idimt,jdimhf,kdim),vva(idimt,jdimhf,kdim),zsa(idimt,jdimhf)
   real  ::  grid(idim*jdim,kdim)
   real  ::  tave(kdim)
#endif
   real  ::  qlnt(mdim),qlnv(mdimv),                                           &
             qdert(mdim),qlnwct(mdim),                                         &
             dplnxn(mdim),dplnyn(mdim),                                        &
             plnw(mdim),plnd2w(mdim),dplnxw(mdim),dplnyw(mdim)
!
   real  ::  uln(mdimv,kdimi),vln(mdimv,kdimi)
   real  ::  zs(idimt),zso(idimt)
!
   real  ::  expps(idimt),rqfx(idimt,kdimi),rqfy(idimt,kdim)
   real  ::  exppso(idimt)
   real  ::  spin(idimt,kdimi+1),spln(idimt,kdimi)
   real  ::  spino(idimt,kdim+1),splno(idimt,kdim)
!
   logical  ::  lnew
#include <funvap.h>
!-------------------------------------------------------------------------------
!
! initialize local variables
!
#ifdef DFS
   psa=0. ; tta=0. ; qqa=0. ; uua=0. ; vva=0. ;  zsa=0. ; grid=0.
#endif
   qlnt=0. ; qlnv=0.   ; qdert=0.  ; qlnwct=0. ; dplnxn=0. ; dplnyn=0.
   plnw=0. ; plnd2w=0. ; dplnxw=0. ; dplnyw=0.
   uln=0.  ;  vln=0.   ; zs=0.     ; zso=0.    ; expps=0.
   rqfx=0. ;  rqfy=0.  ; exppso=0. ; spin=0.   ; spln=0.   ; spino=0. ; splno=0.
!
   in=0
   tensn=10.0
!
   rk1 = akapa_ + 1.
   rkinv=1./akapa_
!
! transfer from slin to sl
!
   print *,' transfer from slin to sl '
   print *,' slin ',slin
   print *,' sl   ',sl
!
#ifdef DFS
#define MDIM mdimd
#else
#define MDIM mdim
#endif
   do k = 1,kdimq
     do i = 1, MDIM
       rqo(i,k)=0.0
     enddo
   enddo
   do i = 1, MDIM
     qo(i)=0.0
   enddo
   do k = 1,kdim
     do i = 1, MDIM
       teo(i,k)=0.0
       dio(i,k)=0.0
       zeo(i,k)=0.0
     enddo
   enddo
#undef MDIM
   do k = 1,kdimi
     call chgr_divor2wind(di(1,k),ze(1,k),uln(1,k),vln(1,k),eps)
   enddo
#ifdef DFS
   call dfs_fft_driver(-1,zsa,idim,jdim,1,gzso,mt,jl,1,jlg,                    &
                       1,1,0,coslat,1)
   call dyn_trans2model_grid(zsa,1)
#endif
!
!  lat loop
!
   do lat = 1,jdimhf
     call chgr_sph_funct(qlnt,qlnv,colrad,lat)
!
     xlat=0.5*pi_-colrad(lat)
!
     call chgr_sum_coeff(gz,zs,qlnt,1)
     call sph_fft2grid(zs,zs,2,1)
!
#ifdef DFS
     zso(1:idimt)=zsa(1:idimt,lat)
#else
     call chgr_sum_coeff(gzso,zso,qlnt,1)
     call sph_fft2grid(zso,zso,2,1)
#endif
!
     call chgr_sum_coeff(q ,ps ,qlnt,     1)
     call chgr_sum_coeff(te,tf ,qlnt, kdimi)
     call chgr_sum_coeff(rq,rqf,qlnt,kdimqi)
!
     call chgr_sum_coeff_wind(uln,uf,qlnv,kdimi)
     call chgr_sum_coeff_wind(vln,vf,qlnv,kdimi)
!
     call sph_fft2grid(ps,ps,2,1)
     call sph_fft2grid(tf,tf,2*kdimi,1)
     call sph_fft2grid(rqf,rqf,2*kdimqi,1)
     call sph_fft2grid(uf,uf,2*kdimi,1)
     call sph_fft2grid(vf,vf,2*kdimi,1)
!
#ifdef DBG
     if (lat .eq. 1 ) then
       print *,' check grid point values before chgr_sigma2sigma.  lat=', lat
       call chgr_comp_maxmin(uln,mdimv,1,kdimi,' uln    ')
       call chgr_comp_maxmin(uf,idimt,1,kdimi,' uf     ')
       call chgr_comp_maxmin(vf,idimt,1,kdimi,' vf     ')
       call chgr_comp_maxmin(tf,idimt,1,kdimi,' tf     ')
       call chgr_comp_maxmin(rqf,idimt,1,kdimi,' rqf    ')
       call chgr_comp_maxmin(ps,idimt,1,1,' ps     ')
     endif
#endif
#ifdef AQUA_PLANET
     call chgr_symmetric(uf,idimt,1,kdimi)
     call chgr_symmetric(vf,idimt,1,kdimi)
     call chgr_symmetric(tf,idimt,1,kdimi)
     call chgr_symmetric(rqf,idimt,1,kdimi)
     call chgr_symmetric(ps,idimt,1,1)
     zs(:) = 0.0
     zso(:) = 0.0
#endif
!
!  at this point grid, values are available
!
!       tf ... temperature
!      rqf ... specific humidity
!       ps ... pai
!       uf ... wind at x direction
!       vf ... wind at y direction
!
     sinlat=cos(colrad(lat))
     rcl=rcs2(lat)
!
     do i = 1, idimt
       expps(i)=exp(ps(i))
     enddo
!
     rmaxsiin=siin(1)
     do k = 1,kdimpi
       rmaxsiin=max(rmaxsiin,siin(k))
     enddo
!
     if(rmaxsiin.gt.1.or.rmaxsiin.eq.0. )then
       do k = 1,kdimi+1  
         do i = 1,idimt
            spin(i,k) = ak(k)/expps(i) + bk(k)
         enddo
       enddo
       do k = 1,kdimi
         do i = 1,idimt
!          spln(i,k) = (spin(i,k)+spin(i,k+1))*0.5
           dif = spin(i,k)**rk1 - spin(i,k+1)**rk1
           dif = dif / (rk1*(spin(i,k)-spin(i,k+1)))
           spln(i,k) = dif**rkinv
         enddo
       enddo
     else
       do k = 1,kdimi+1
         do i = 1,idimt
           spin(i,k) = siin(k)
         enddo
       enddo
       do k = 1,kdimi
         do i = 1,idimt
           spln(i,k) = slin(k)
         enddo
       enddo
     endif
!
!  extrapolate humidity into stratosphere
! kdimqi --> kdimi for rqf --> rqfx
!
     if(lnew) then
       call chgr_moisture_extrap(expps,rqf,tf,rqfx,lat,xlat)
     else
       call chgr_moisture_extrap_old (expps,rqf,tf,rqfx)
     endif
!
! ---   end of grid point and chgr_vinterp_solver q to stratosphere at old lev
!
     call chgr_get_sfcp(tf,rqfx,ps,zs,zso,pso,spin)
!    do 300 i = 1, idimt
!    pso(i)=ps(i)
!300 continue
!
!
#ifndef HYBRID
     do k = 1,kdimp
       do i = 1,idimt
         spino(i,k)=si(k)
       enddo
     enddo
     do k = 1,kdim
       do i = 1,idimt
         splno(i,k)=sl(k)
       enddo
     enddo
#else
     do i = 1, idimt
       exppso(i)=exp(pso(i))
     enddo
!
     do k = 1,kdimp
       do i = 1,idimt
         spino(i,k) = ak5(k)/exppso(i) + bk5(k)
       enddo
     enddo
     do k = 1,kdim
       do i = 1,idimt
         dif = spino(i,k)**rk1 - spino(i,k+1)**rk1
         dif = dif / (rk1*(spino(i,k)-spino(i,k+1)))
         splno(i,k) = dif**rkinv
!        splno(i,k) = (spino(i,k)+spino(i,k+1))*0.5
       enddo
     enddo
#endif
!
! start chgr_sigma2sigma 
!
     call chgr_sigma2sigma(ps,spln,uf,pso,splno,ufo,idimt,1,                   &
                           kdimi,kdim,in,tensn,1)
     call chgr_sigma2sigma(ps,spln,vf,pso,splno,vfo,idimt,1,                   &
                           kdimi,kdim,in,tensn,1)
     call chgr_sigma2sigma(ps,spln,tf,pso,splno,tfo,idimt,1,                   &
                           kdimi,kdim,in,tensn,1)
     call chgr_sigma2sigma(ps,spln,rqfx,pso,splno,rqfy,idimt,1,                &
                           kdimi,kdim,in,tensn,1)
!
! chgr_sigma2sigma sets values to constant outside of input domain.
! fix lapse rate and relative humidity below input surface.
!
     gamma=6.5e-3
     do i = 1,idimt
       ps1=log(spln(i,1))+ps(i)
       rh1=rqfx(i,1)/fqs(tf(i,1),exp(ps1))
       do k = 1,kdim
         psk=log(splno(i,k))+pso(i)
!
         if(psk.lt.ps1) goto 299
!
         tfo(i,k)=tf(i,1)*exp(gamma*rd_/g_*(psk-ps1))
         rqfy(i,k)=rh1*fqs(tfo(i,k),exp(psk))
       enddo
       299 continue
     enddo
!
! kdim --> kdimq
!
     do kq=1,kdimq
       do i = 1, idimt
         rqfo(i,kq)=rqfy(i,kq)
       enddo
     enddo
!
! ======= check values ========
!
#ifdef DBG
     if(lat .eq. 1 ) then
       print *,' check grid point values after chgr_sigma2sigma.   lat=', lat
       call chgr_comp_maxmin(ufo,idimt,1,kdim,'    ufo ')
       call chgr_comp_maxmin(vfo,idimt,1,kdim,'    vfo ')
       call chgr_comp_maxmin(tfo,idimt,1,kdim,'    tfo ')
       call chgr_comp_maxmin(rqfo,idimt,1,kdimq,'   rqfo ')
       call chgr_comp_maxmin(pso,idimt,1,1,'    pso ')
     endif
#endif
!
#ifdef DFS
!
! save grid value
!
     zsa(1:idimt,lat)=zso(1:idimt)
     psa(1:idimt,lat)=pso(1:idimt)
     do k = 1,kdim
       uua(1:idimt,lat,k)=ufo(1:idimt,k)
       vva(1:idimt,lat,k)=vfo(1:idimt,k)
       tta(1:idimt,lat,k)=tfo(1:idimt,k)-tave(k)
     enddo
     do k = 1,kdimq
       qqa(1:idimt,lat,k)=rqfo(1:idimt,k)
     enddo
   enddo
   call dyn_trans2output_grid(zsa,1)
   call dyn_trans2output_grid(psa,1)
   call dyn_trans2output_grid(uua,kdim)
   call dyn_trans2output_grid(vva,kdim)
   call dyn_trans2output_grid(tta,kdim)
   call dyn_trans2output_grid(qqa,kdimq)
#ifdef BIN_DBG
   call file_write_bin(112,zsa,idim,jdim,1,0)
   call file_write_bin(112,psa,idim,jdim,1,0)
   call file_write_bin(112,uua,idim,jdim,kdim,0)
   call file_write_bin(112,vva,idim,jdim,kdim,0)
   call file_write_bin(112,tta,idim,jdim,kdim,0)
   call file_write_bin(112,qqa,idim,jdim,kdimq,0)
#endif
!
!  dfs transform
!
! not use gzo variable
!      call dfs_fft_driver(1,zsa,idim,jdim,1,gzo,mt,jl,1,jl,
!     &                  1,1,0,coslat,1)
   call dfs_fft_driver(1,psa,idim,jdim,1,qo ,mt,jl,1,jl,                       &
                    1,1,0,coslat,1)
   call dfs_fft_driver(1,tta,idim,jdim,kdim,teo,mt,jl,kdim,jl,                 &
                    kdim,kdim,1,coslat,1)
   call dfs_fft_driver(1,qqa,idim,jdim,kdimq,rqo,mt,jl,kdimq,jl,               &
                    kdim,kdim,kdimq/kdim,coslat,1)
   call dfs_wind2divor(uua,vva,idim,jdim,kdim,zeo,dio,mt,jl,kdim,coslat,1)
#ifdef BIN_DBG
!
! wave
!
   call file_write_bin(119,dio,mt,jl,kdim,0)
   call file_write_bin(119,zeo,mt,jl,kdim,0)
   call file_write_bin(119,teo,mt,jl,kdim,0)
   call file_write_bin(119,qo,mt,jl,kdimq,0)
   call dfs_fft_driver(-1,uua,idim,jdim,kdim,dio,mt,jl,kdim,jl,                &
                    kdim,kdim,1,coslat,1)
   call dfs_fft_driver(-1,vva,idim,jdim,kdim,zeo,mt,jl,kdim,jl,                &
                    kdim,kdim,1,coslat,1)
!
! grid
!
   call file_write_bin(112,uua,idim,jdim,kdim,0)
   call file_write_bin(112,vva,idim,jdim,kdim,0)
#endif /* end BIN_DBG */
!
#else /* not DFS */
!  spherical transform
!
     call sph_fft2grid(pso,pss,2,-1)
     call sph_fft2grid(tfo,tts,2*kdim,-1)
     call sph_fft2grid(rqfo,qqs,2*kdimq,-1)
     call sph_fft2grid(ufo,uus,2*kdim,-1)
     call sph_fft2grid(vfo,vvs,2*kdim,-1)
!
     do i = 1, mdim
       qdert(i) = qlnt(i) * wgt(lat)
     enddo
!
     call chgr_fourier_diff(pss,psa,1)
     call sph_coeff_lat(pss,psa,qo,qdert,1)
!
     call chgr_fourier_diff(tts,tta,kdim)
     call sph_coeff_lat(tts,tta,teo,qdert,kdim)
!
     call chgr_fourier_diff(qqs,qqa,kdimq)
     call sph_coeff_lat(qqs,qqa,rqo,qdert,kdimq)
!
     call chgr_sph_derivative(qlnt,qlnv,qdert,eps,lat,qlnwct,rcs2,wgt(lat))
!
     call chgr_fourier_diff(uus,uua,kdim)
     call chgr_fourier_diff(vvs,vva,kdim)
     call chgr_sph_operator(uua,uus,vva,vvs,zeo,qlnwct,qdert,kdim)
     call chgr_div_coeff_lat(uua,uus,vva,vvs,dio,qlnwct,qdert,kdim)
!
   enddo
!
   do k = 1,kdim
     do m = 1,mdim
       dio(m,k)=-dio(m,k)
     enddo
   enddo
#endif
#ifdef DBG
!
! ======= check values ========
!
#ifdef DFS
#define MDIM mdimd
#else
#define MDIM mdim
#endif
   print *,'check coef(before noah_soilm_transpiration) after chgr_sigma2sigma'
   call chgr_comp_maxmin(zeo,MDIM,1,kdim, '   zeo  ')
   call chgr_comp_maxmin(dio,MDIM,1,kdim, '   dio  ')
   call chgr_comp_maxmin(teo,MDIM,1,kdim, '   teo  ')
   call chgr_comp_maxmin(rqo,MDIM,1,kdimq,'   rqo  ')
   call chgr_comp_maxmin(qo,MDIM,1,1,'   qo   ')
#undef MDIM
#endif
!
   return
   end subroutine chgr_vinterp_solver
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_sigma2sigma(ps1,s1,var1,ps2,s2,var2,im,jm,km1,km2,          &
                               in,tensn,iuv)               
!-------------------------------------------------------------------------------
!
! abstract: this subroutine transfer from one defined sigma coordinate (s1) to
!   another sigma coordinate (s2). 
!
! program history log:
!
! usage:    call chgr_sigma2sigma(ps1,s1,var1,ps2,s2,var2,im,jm,km1,km2,       &
!                           in,tensn,iuv)  
!  input argument list:                                              
!    ps1      - primary ground pressure   alog(psfc)                         
!    ps2      - secondary gound pressure alog(psfc(cb))                      
!    s1       - primary sigma coordinate                                     
!    var1     - primary 3 dimensional variable                               
!    s2       - secondary sigma coordinate                                   
!    var2     - secondary 3 dimensional variable                             
!    im       - dimension in x                                               
!    jm       - dimension in y                                               
!    km1      - dimension in z for primary coordinate and variable           
!    km2      - dimension in z for secondary coordinate and variable         
!    in       - control index:   0 for initial transformation                
!                                1 for transformation as previous            
!    tensn    - factor of tension 0 for cubic spline                         
!                                 50 for linear interpolation                 
!    iuv      - index for wind avoid extrapolation                           
!                                                                               
!  local dimensions: st1 st2 sv1 sv2 q vinc  
!                                                                               
!  common block /spl/ ovh sh iflag theta    
!                                                                               
!  routines related: chgr_cubic_coeff chgr_cubic_compute         
!                                                                               
!-------------------------------------------------------------------------------
   save
!
   common /spl/ ovh(100), sh(100), iflag, jflag, theta
   real   ::  s1(im,jm,km1),s2(im,jm,km2)                 
   real   ::  var1(im,jm,km1),var2(im,jm,km2)                 
   real   ::  st1(100),st2(100),st0(100),ps1(im,jm),ps2(im,jm)                
   real   ::  sv1(100),sv2(100),q(100),vinc(100)                              
!-------------------------------------------------------------------------------
!
   theta=tensn      ! theta is tension factor      
                                                                             
   if( in .eq. 0 ) then                                                      
     in=1                                                                      
     print *, ' initial values for transformation.'                            
     iflag = 0                                                                 
     print *,'  iflag=',iflag,' tension factor=',theta                         
!
     do k = 1,km1
       do j = 1,jm
         do i = 1,im
           kk=km1-k+1
           if( s1(i,j,kk) .gt. 0.0  .and.  s1(i,j,kk) .lt. 1.0 ) then
             st0(k)=log(s1(i,j,kk))
           else
             print *, ' error in s1(k) values at k=',k,' sigma=',s1(i,j,kk)
             call abort
           endif
         enddo
       enddo
     enddo
!
     do k = 1,km2
       do j = 1,jm
         do i = 1,im
           kk=km2-k+1
           if( s2(i,j,kk) .gt. 0.0  .and.  s2(i,j,kk) .lt. 1.0 ) then
             st2(k)=log(s2(i,j,kk))
           else
             print *, ' error in s2(i,j,k) values at k=',k,'                   &
                          sigma=',s2(i,j,kk)
             call abort
           endif
         enddo
       enddo
     enddo
!
   endif                                                                     
!                                                                             
   do i = 1,im                                                              
     do j = 1,jm                                                              
       do k = 1,km1            
         kk=km1-k+1
         st1(k)=log(s1(i,j,kk))+ps1(i,j)-ps2(i,j)     
         sv1(k)=var1(i,j,km1-k+1)           
       enddo
!
       do k = 1,km2
         kk=km2-k+1
         st2(k)=log(s2(i,j,kk))
       enddo
       vst=st2(1)
       do k = 1,km2-1
         vinc(k)=st2(k+1)-st2(k)
       enddo
!
       call chgr_cubic_coeff(km1,st1,q,sv1)        
       call chgr_cubic_compute(q,vst,vinc,st1,sv1,km2,km1,sv2) 
!
       kmm=km2/2 + 1 
!
       do k = kmm,1,-1                                                          
         if( sv2(k) .eq. 99999.9 ) then                
!                                                                            
!    set to the lowest primary sigma level value         
!                                                    
           do kk = 1,km1                             
             if(st1(kk).gt.st2(k)) then          
               diff1=st1(kk)-st2(k)            
               if(kk.gt.1) then               
                 diff2=st1(kk-1)-st2(k)     
               else                         
                 diff2=999.              
               endif                       
               if(diff1.lt.diff2) then   
                 sv2(k)=sv1(kk)        
               else                    
                 sv2(k)=sv1(kk-1)    
               endif                 
               exit
             endif              
           enddo                
         endif                  
       enddo
!
       do k = kmm,km2                                                           
         if( sv2(k) .eq. 99999.9 ) then        
!                                                                            
!    set to the lowest primary sigma level value 
!                                            
         do kk = 1,km1                     
           if(st1(kk).gt.st2(k)) then  
             diff1=st1(kk)-st2(k)    
             if(kk.gt.1) then       
               diff2=st1(kk-1)-st2(k)             
             else                                 
               diff2=999.                       
             endif                              
             if(diff1.lt.diff2) then           
               sv2(k)=sv1(kk)                
             else                            
               sv2(k)=sv1(kk-1)            
             endif                         
             exit
           else                   
             sv2(k)=sv1(km1)    
           endif                
         enddo                  
       endif                                      
     enddo
!                                                                             
     do k = 1,km2                                                             
       var2(i,j,k)=sv2(km2-k+1)    
     enddo
!
     enddo
   enddo
!
   return                                                                    
   end subroutine chgr_sigma2sigma
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_moisture_extrap(pstar,qin,t,chgrp,lat,xlat)    
!-------------------------------------------------------------------------------
!
! abstract: this code extrapolates moisture up into model dry layers
!   this code extrapolates moisture up into model dry layers                  
!   using exponential decrease up to a specified q at 80 mb,                    
!   then use lims profile above(see middle atmos program handbook,              
!   vol 22,sep 86,ed. j.m.russell ---> j alpert has copy)....                   
!
! usage:    call chgr_moisture_extrap(pstar,qin,t,chgrp,lat,xlat)              
!  input argument list:
!    pstar    - sfc pres (cb)                                              
!    qin      - moisture - levh lyrs only                                 
!    t        - temperature - model virtual temperature                      
!    xlat     - lat in radians 
!                                                                               
!  output argument list:                                                      
!    chgrp    - levs moisture lyrs (use lims data abuv 8 cb)                
!             - thermodynamic temperature back into t {revised to            
!             - restore temp. to virtual temp.                            
!
!-------------------------------------------------------------------------------
   use constant, only : pi_,rd_,rv_,t0c_,psat_
   use paramter, only : idim
   use parmchgr, only : kdimi, kdimqi
   use comchgr, only : idimt, slin, ak, bk
!-------------------------------------------------------------------------------
   save
!
   real                 ::  xlamb,expon
   real                 ::  qin(idimt,kdimqi),t(idimt,kdimi),qout(kdimi)
   real                 ::  qsat(kdimi),prs(kdimi),pstar(idimt)
   real                 ::  chgrp(idimt,kdimi)
!
!  qq stores lims (satellite-1986-see j.alpert) observed h2o above              
!  80 mb...(k,j)-k=1,8 for 80,70,..,10 mb                                       
!  j=1,19 for every 10 deg lat (starting at n.pole)                             
!                                                                               
   real                 ::  qq(8,19)                                           
   data qq /32*5.  ,  4.8,4.55,4.5,4.6,4.7,4.75,4.75,4.7,                      &
      4.,3.85,3.9,4.,4.25,4.5,4.5,4.6, 3.,2.75,3.,3.3,3.7,4.2,4.3,4.5,         &
      2.1,2.,2.,2.5,3.,3.9,4.,4. ,                                             &
      2.,2.,2.,2.3,2.8,3.5,3.7,3.75 , 2.,2.,2.,2.3,2.7,3.25,3.45,3.75,         &
      2.,2.,2.,2.3,2.8,3.5,3.7,3.75 , 2.1,2.,2.,2.5,3.,3.9,4.,4. ,             &
      3.,2.75,3.,3.3,3.7,4.2,4.3,4.5,                                          &
      4.,3.85,3.9,4.,4.25,4.5,4.5,4.6 ,                                        &
      4.8,4.55,4.5,4.6,4.7,4.75,4.75,4.7 , 32*5.  /                    
   data itst / 0 /                                                           
!-------------------------------------------------------------------------------
!                                                                               
!  first convert mdl virtual temp to thermodynmic temp                          
!                                                                               
   if (itst.gt.0) go to 1                                                    
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
   itst = 1                                                                  
!
   1 continue                                                                  
!                                                                               
!  loop for hemisphere                                                          
!                                                                               
   do ihem = 1,2                                                          
!                                                                               
!  get latitude in degrees..                                                    
!                                                                               
     dlat = xlat * 180.0/pi_                                                   
     if (ihem.eq.2) dlat =-dlat                                                
!                                                                               
!  lon loop                                                                     
!                                                                               
     do i = 1,idim                                                            
       ix=(ihem-1)*idim+i 
!                                                                               
!  model temp is virtual (since hydrostatic from hgts initially),               
!                                                                               
       do 2 k=1,kdimqi 
!
         if (qin(ix,k).le.0.0) go to 2 
!
         t(ix,k)=t(ix,k)*(1.0+qin(ix,k))/(1.0+rv_/rd_*qin(ix,k)) 
       2  continue
!                                                                               
!  get lyr pressure and then saturated moisture                                 
!                                                                               
       do k = 1,kdimi 
#ifndef HYBRID
         prs(k) = pstar(ix) * slin(k)
#else
         prs(k) = pstar(ix) * bk(k) + ak(k)
#endif
       enddo
!                                                                               
!  compute saturation specific humidity(dimensionless) from                     
!  temperature t (deg k) and pressure (cb)                                      
!  conversion to specific humidity follows from definition                      
!                                                                               
       do k = 1,kdimi
         expon = 7.50*(t(ix,k)-t0c_)/((t(ix,k)-t0c_)+237.30) 
         es = psat_*1.e-2 * 10.0**expon
         qsat(k)=rd_/rv_*es/(prs(k)*10.0-(1.0-rd_/rv_)*es)
       enddo
!                                                                               
!  limit moisture in lowest levh layers--rh le 1. but ge .15                    
!          the latter avoids the negative q problem...                          
!                                                                               
       do k = 1,kdimqi
         qout(k) = qin(ix,k)
         rh = qin(ix,k) / qsat(k) 
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
         xlnqm =  log(q8) 
!                                                                               
!  extrapolate moisture to min valu(q8) at pressure pmin(80 mb)                 
!   use exponential decrease from layer levh----                                
!   i.e.  q=q(levh)*(p/p(levh)) ** xlamb                                        
!   where xlamb is computed to fit q8 at pmin and q,p at levh.                  
!                                                                               
         qqout = qin(ix,kdimqi)
         rh = qqout / qsat(kdimqi)
         if (rh.le..15) qqout = .15 * qsat(kdimqi)
         if (rh.gt. 1.0) qqout = qsat(kdimqi) 
         xlnpq =  log(prs(kdimqi)) 
         xlnqq =  log(qqout)
         xlamb = (xlnqm-xlnqq) / (xlnpm-xlnpq) 
         do k = ldry,kdimi
           if (prs(k).ge.pmin) go to 21 
           if (prs(k).ge.pmm) go to 22 
!                                                                               
!  above pmm(10 mb) use constant value from table                               
!                                                                               
           qout(k) = q1 
!
           go to 13 
           21 continue 
!
           qout(k) = qout(kdimqi)*(prs(k)/prs(kdimqi)) ** xlamb 
!
           go to 13
           22 continue
!                                                                               
!  above 8 cb so complete linear interp from table..                            
!                                                                               
           kdy = 9. - prs(k)
           dy = 9. - prs(k) - kdy 
           qout(k) = qq(kdy,jdx  )*(1.-dy)*(1.-dx)+qq(kdy+1,jdx+1)*dy*dx       &
                      +qq(kdy,jdx+1)*(1.-dy)*dx+qq(kdy+1,jdx)*dy*(1.-dx) 
!
           13 rh = qout(k) / qsat(k)
!
           if(rh.gt.1.) qout(k)=qsat(k) 
         enddo
! 
       endif  
!                                                                               
!  store extrapolated moisture                                                  
!                                                                               
       do k = 1,kdimi 
         chgrp(ix,k) = qout(k) 
       enddo
!
! ---- temp. is restored to virtual temp.                                       
!
       do 25 k = 1,kdimqi  
         if (qin(ix,k).le.0.) go to 25 
         t(ix,k) = t(ix,k) * (1. + rv_/rd_ * qin(ix,k))                        &
                  /(1. + qin(ix,k))
       25 continue 
!                                                                               
     enddo
!                                                                               
   enddo
!                                                                               
   return                                                                    
   end subroutine chgr_moisture_extrap 
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine chgr_moisture_extrap_old(pstar,qin,t,chgrp)
!-------------------------------------------------------------------------------
!
! abstract: extrapolates moisture up into model dry layers
!    using exponential decrease up to a minimum valu of q
!    first convert mdl virtual temp to thermodynmic temp
!
! usage:    call chgr_moisture_extrap_old(pstar,qin,t,chgrp)
!   input argument list:
!     pstar    - sfc pres (cb)                                          
!     qin      - moisture - levh lyrs only                           
!     t        - temperature - model virtual temperature                
!                                                                               
!   output argument list:                                                       
!     chgrp    - levs moisture lyrs                                     
!              - thermodynamic temperature back into t                         
!              - virtual temperature should be restored                 
!
!-------------------------------------------------------------------------------
   use constant, only : rd_,rv_,t0c_,psat_
   use paramter, only : idim
   use parmchgr, only : kdimi, kdimqi
   use comchgr, only  : slin, ak, bk
!-------------------------------------------------------------------------------
   save                                                                      
!
   real                 ::  qin(idim,kdimqi),t(idim,kdimi),qout(kdimi)
   real                 ::  qsat(kdimi),prs(kdimi),pstar(idim)
   real                 ::  chgrp(idim,kdimi)
!
   data itst / 0 /                                                           
!-------------------------------------------------------------------------------
!                                                                               
   if (itst.gt.0) go to 1                                                    
!
      pmin = 5.                                                                 
      qmin = 3.e-6                                                              
      xlnpm = log(pmin)                                                         
      xlnqm = log(qmin)                                                         
      itst = 1                                                                  
   1 continue    
!
   do i = 1,idim                                                            
!
!----    model temp is virtual (since hydrostatic from hgts initially),         
!          so we convert to thermodynamic temp for radiation calculation        
!
     do 2 k = 1,kdimqi 
!
       if (qin(i,k).le.0.0) go to 2
!
       t(i,k) = t(i,k) * (1.0 + qin(i,k))                                      &
                /(1.0 + rv_/rd_ * qin(i,k))                            
     2  continue 
!                                                                               
!----    get lyr pressure and then saturated moisture                           
!
     do k = 1,kdimi 
#ifndef HYBRID
       prs(k) = pstar(i) * slin(k)
#else
       prs(k) = pstar(i) * bk(k) + ak(k)
#endif
     enddo
!
!----    compute saturation specific humidity(dimensionless) from               
!          temperature t (deg k) and pressure (cb)                              
!        method is tetens formula for saturation vapor pressure (mb) ,          
!          given by haurwitz p. 9 - dynamic meteorology , 1941.                 
!----    conversion to specific humidity follows from definition                
!          ibid. , p.10                                                         
!                                                                               
     do k = 1,kdimi                                                            
       expon = 7.50*(t(i,k)-t0c_)/((t(i,k)-t0c_)+237.30) 
       es = psat_*1.e-2 * 10.0**expon                   
       qsat(k)=rd_/rv_*es/(prs(k)*10.0-(1.0-rd_/rv_)*es)
     enddo
!
!----    limit moisture in lowest levh layers--rh le 1. but ge .15              
!          the latter avoids the negative q problem...                          
!
     do k = 1,kdimqi                                    
       qout(k) = qin(i,k)                             
       rh = qin(i,k) / qsat(k)                       
     enddo
     ldry = kdimqi + 1                               
!                                                                        
     if( ldry .lt. kdimi ) then                                                
!                                                                               
!----    extrapolate moisture to min valu(qmin) at pressure pmin(cb)            
!         use exponential decrease from layer levh----                          
!         i.e.  q=q(levh)*(p/p(levh)) ** xlamb                                  
!        where xlamb is computed to fit qmin at pmin and q,p at levh.           
!
       xlnpq = log(prs(kdimqi))                    
       xlnqq = log(qout(kdimqi))                  
       xlamb = (xlnqm-xlnqq) / (xlnpm-xlnpq)     
       do k = ldry,kdimi                        
         qout(k) = qout(kdimqi)*(prs(k)/prs(kdimqi)) ** xlamb
         if (prs(k).le.pmin) qout(k) = qmin                 
         if (qout(k).le.qmin) qout(k) = qmin               
         rh = qout(k) / qsat(k)                           
         if (rh.gt.1.0) qout(k) = qsat(k)                
       enddo
!                                                                        
     endif                                                                     
!
!----    store extrapolated moisture                                            
!
     do k = 1,kdimi                                       
       chgrp(i,k) = qout(k)                             
     enddo
!
! ---- temp. is restored to virtual temp.                                       
!
     do 25 k=1,kdimqi                                                          
       if (qin(i,k).le.0.0) go to 25   
       t(i,k) = t(i,k) * (1.0 + rv_/rd_ * qin(i,k))                            &
                /(1.0 + qin(i,k))                                     
     25 continue 
   enddo
!
   return                                                                    
   end subroutine chgr_moisture_extrap_old
!-------------------------------------------------------------------------------
!
#ifdef AQUA_PLANET
!-------------------------------------------------------------------------------
   subroutine chgr_symmetric(a,idimt,one,kdim)
!-------------------------------------------------------------------------------
!
!  set the same value for norther hemisphere to that in the southern
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                 :: idim,idimt,kdim,i,k,one
   real                    :: a(idimt,kdim)
!
   idim = idimt/2
   do k = 1,kdim
     do i = 1,idim
       a(i,k) = a(i+idim,k)
     enddo
   enddo
!
   return                                                                    
   end subroutine chgr_symmetric
!-------------------------------------------------------------------------------
#endif   
