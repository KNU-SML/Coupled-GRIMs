#include <define.h>
   subroutine phys_gwd_alpert(dudt,dvdt,u1,v1,t1,q1,                           &
                    delprsi,prsi,prsl,prslk,zl,rcl,                            &
                    deltim,lat,kdt,hprime,dusfc,dvsfc,                         &
                    g,cp,rd,rv,fv,                                             &
                    ids,ide, jds,jde, kds,kde,                                 &
                    ims,ime, jms,jme, kms,kme,                                 &
                    its,ite, jts,jte, kts,kte)
!-------------------------------------------------------------------------------
!
!  subprogram:    phys_gwd_alpert       includes gravity wave drag.
!
!  abstract: using the gwd parameterizations of ps-glas and ph-
!    gfdl technique, the time tendencies of u v
!    are altered to include the effect of mountain induced
!    gravity wave drag from sub-grid scale orography including
!    convective breaking, shear breaking and the presence of
!    critical levels.
!
! program history log:
!   1987-06-03  jordan c. alpert       development 
!   2000-01-01  song-you hong          physcis options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!  input argument list:
!     dudt     - negative non-lin tendency for u wind component.
!     dvdt     - negative non-lin tendency for v wind component.
!     u1       - zonal      wind component *cos(lat)  m/sec at t0-dt.
!     v1       - meridional wind component *cos(lat)  m/sec at t0-dt.
!     t1       - temperature deg k at t0-dt.
!     q1       - specific humidity at t0-dt.
!     prsi(n)  - p/psfc at base of layer n.
!     del(n)   - positive increment of p/psfc across layer n.
!     prsl(n)  - p/psfc at middle of layer n.
!     rcl      - reciprocal of square of cos(lat).
!     deltim   - time step  secs.
!     lat      - latitude  number.
!     kdt      - time step number.
!     hprime   - topographic standard deviation  (m).
!
!  output argument list:
!     dudt     - as augmented by tendency due to migwd.
!     dvdt     - as augmented by tendency due to migwd.
!
!  output files:
!     ft06f001 - printout file.
!
!  ::: structure :::
!
!    [phys_gwd_alpert] --- [phys_gwd_alpert_sub] *
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#undef DBG
!
   integer              ::  lat, kdt,                                          & 
                            ids,ide, jds,jde, kds,kde,                         &
                            ims,ime, jms,jme, kms,kme,                         &
                            its,ite, jts,jte, kts,kte
   real                 ::  g,rd,rv,cp,fv,rcl,deltim,rdz
   real                 ::  dudt(ims:ime,kms:kme),dvdt(ims:ime,kms:kme),       &
                            u1(ims:ime,kms:kme),v1(ims:ime,kms:kme),           &
                            t1(ims:ime,kms:kme),q1(ims:ime,kms:kme),           &
                            prsi(ims:ime,kms:kme+1),prsl(ims:ime,kms:kme),     &
                            delprsi(ims:ime,kms:kme),slk(kms:kme),             &
                            prslk(ims:ime,kms:kme),zl(ims:ime,kms:kme),        &
                            hprime(ims:ime),                                   &
                            dusfc(ims:ime),dvsfc(ims:ime)
!
!  constants for migwd *j!
!  only do constants first time through monnin
!
   real,parameter       ::  dw2min = 1.
   real,parameter       ::  rimin  = -100.
   real,parameter       ::  rlowlv = 0.7
   real,parameter       ::  xl     = 4.0e4
!
!  local variables
!
   logical              ::  ldrag(ims:ime)
   integer              ::  i,k,kbj,kbjbeg,kbjp1,kpblmax,lcap,ksm,ksmm1
   real                 ::  rcs,cs,value,dw2,bvf2,rdz2,xlinv,aj,ti,shr2,       &
                            gmax,gr2,savem,factop,wtkbj,rcsks,rdelks,          &
                            sira,veleps
   real                 ::  velco(ims:ime,kms:kme-1)
   real                 ::  rdzt(ims:ime,kms:kme-1),                           &
                            delks(ims:ime),delks1(ims:ime)
!
!  mountain induced gravity wave drag
!  unit14  - subgrid scale mountain variance height input
!  common to be added to smf,gloo for migwd        *j!
!  ite = 256 and kte = 18 for example...
!
   real                 ::  taub(ims:ime),                                     &
                            vtj(ims:ime,kms:kme),                              &
                            usqj(ims:ime,kms:kme),                             &
                            bnv2(ims:ime,kms:kme),                             &
                            ro(ims:ime,kms:kme),                               &
                            taud(ims:ime,kms:kme),                             &
                            dtfac(ims:ime,kms:kme),                            &
                            taup(ims:ime,kms:kme+1)
   real                 ::  xn(its:ite),yn(its:ite),                           &
                            fr(its:ite),gf(its:ite),                           &
                            ubar(its:ite),vbar(its:ite),                       &
                            ulow(its:ite),bnv(its:ite),                        &
                            dtaux(its:ite,kms:kme),                            &
                            dtauy(its:ite,kms:kme),                            &
                            vtk(its:ite,kms:kme),                              &
                            roll(its:ite)
#ifdef REIGH_FRICTION
   real                 ::  sigma, damping
   real, parameter      ::  sigma_cr = 0.05, r_day_scale = 1./(5.*86400.)
#endif
!-------------------------------------------------------------------------------
!
!  kbj is the bottom of the low 1/3 level usually = 1
!
   kbj = 1
!
   do k = kbj,kte
     do i = its,ite
       if (prsi(i,k)/prsi(i,1) .lt. rlowlv) then
         ksm = k
         !
         go to 16
         !
       endif
     enddo
   enddo
!
16 continue
!
!  ksm -1 intervals in the lower third of atm (sigma < .667)
!
   ksmm1  = ksm - 1
!
   do i = its,ite
     delks(i)  = (prsi(i,kbj)-prsi(i,ksm))/prsi(i,1)
     delks1(i) = (prsl(i,kbj)-prsl(i,ksm))/prsi(i,1)
   enddo
!
!  above, the low layer delta sigma
!  below the starting sigma level for ps stress calc defaults to 2
!
   kpblmax = 2
   lcap    = kte
   factop  = 0.5
!
   gr2    = 2.0 * g * g / rd
   gmax   = 1.
   aj     = 1.
   xlinv  = 1.0 / xl
   veleps = 1.0
   rcs    = sqrt(rcl)
   cs     = 1. / rcs
!
!  saving richardson number in usqj for migwd        *j!
!
   do k = kts,kte-1
     do i = its,ite
       rdzt(i,k) = g/rd*prsi(i,k+1)/(prsl(i,k)-prsl(i,k+1))
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       vtj(i,k) = t1(i,k)*(1.+fv*q1(i,k))
       vtk(i,k) = vtj(i,k) / prslk(i,k)
     enddo
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       ti        = 0.5*(t1(i,k)+t1(i,k+1))
       rdz       = rdzt(i,k)/ti
       dw2       = rcl*((u1(i,k)-u1(i,k+1))**2+(v1(i,k)-v1(i,k+1))**2)
       shr2      = max(dw2,dw2min)*rdz**2
       bvf2      = g*(g/cp+rdz*(vtj(i,k+1)-vtj(i,k)))/ti
       usqj(i,k) = max(bvf2/shr2,rimin)
     enddo
   enddo
!
!  the linear mountain induced gravity mode p&s prameterization explitly done
!  this routine computes the deceleration of the zonal wind and
!  meridional wind due to mountain gravity drag.
!
!  code variables          description
!
!  xn,yn    projections of "low-level" wind in zonal & meridional directions
!
!  ulow             "low-level" wind magnitude -        (= u)
!                   averaged up to 2km above surface
!
!  bnv2             bnv2 = n**2
!
!  hprime           sub-grid scale mountain height      (= h)
!                   from navy tape, averaged,'envelope'std. va
!                   read in in smf,common-ed to gloo
!
!  taub             base momentum flux
!                   = -(ro * u**3/(n*xl)*gf(fr) for n**2> 0
!                   = 0.                        for n**2 < 0
!
!  fr               froude    =   n*hprime / u
!  g                gmax*fr**2/(fr**2+aj**2)
!  gmax             = 1.0
!  aj               = 1.0
!
!  ksm is defined as the number of levels up 1/3 from the lowest used
!  to calculate the "low-level" averages.
!
!  initialize arrays     (on cyber)
!
   do i = its,ite
     xn(i)   = 0.0
     yn(i)   = 0.0
     ubar(i) = 0.0
     vbar(i) = 0.0
     roll(i) = 0.0
     taub(i) = 0.0
     ulow(i) = 0.0
     taup(i,kte+1) = 0.0
   enddo
!
   do k = kts,kte
     do i = its,ite
       taup(i,k) = 0.0
       ro(i,k)   = prsl(i,k) / ( rd * vtj(i,k) )
     enddo
   enddo
!
!  density   tons/meter**3
!
!  compute low level averages
!  (u,v)*cos(lat)  use uv=(u1,v1) which is wind at t0-1
!  use rcs=1/cos(lat) to get wind field
!  ksm   the top of the lowest 1/3 layer "the low level" is 6
!
   do k = kbj,ksmm1
     do i = its,ite
       rcsks   = rcs * delprsi(i,k) / prsi(i,1) / delks(i)
       ubar(i) = ubar(i) + rcsks * u1(i,k)
       vbar(i) = vbar(i) + rcsks * v1(i,k)
     enddo
   enddo
!
!  compute the "low level" or 1/3 wind magnitude (m/s)
!
   do i = its,ite
     ulow(i) = sqrt( ubar(i) * ubar(i) + vbar(i) * vbar(i) )
   enddo
!
   do i = its,ite
     value   = 1.0
     ulow(i) = max( ulow(i), value )
   enddo
!
!  calculate squared low level brunt vaisala frequency over the
!  first ksm levels then average
!  rdelks (del(k)/delks) vert ave factor so we can * instead of /
!
   do i = its,ite
     bnv2(i,1) = 0.
   enddo
!
#ifdef DBG
   print*,'--------------------------------------------------------------------'
   print*,'lat kdt',lat,kdt
   print*,'--------------------------------------------------------------------'
#endif
!
#ifdef DBG
print *,'bnv2 in phys_gwd_alpert'
#endif
!
   do k = kbj,ksmm1
     do i = its,ite
       rdz2      = 1./(zl(i,k+1) - zl(i,k))
       bnv2(i,k) = 2*g*rdz2*(vtk(i,k+1)-vtk(i,k))/(vtk(i,k+1)+vtk(i,k))
!
#ifdef DBG
if(i.eq.ite) print*,'i,k',i,k,bnv2(i,k)
#endif
!
     enddo
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       velco(i,k) = (0.5*rcs)*( (u1(i,k) + u1(i,k+1)) * ubar(i) +              &
                                (v1(i,k) + v1(i,k+1)) * vbar(i))
       velco(i,k) = velco(i,k)/ulow(i)
       if ((velco(i,k).lt.veleps).and.(velco(i,k).ge.0.)) then
         velco(i,k) = veleps
       endif
     enddo
   enddo
!
!  no drag when critical level in the base layer
!
   do i = its,ite
     ldrag(i) = velco(i,1).le.0.
   enddo
!
   do k = 2,ksmm1
     do i = its,ite
       ldrag(i) = ldrag(i).or. velco(i,k).le.0.
     enddo
   enddo
!
!  no drag when bnv2.lt.0
!
   do k = 1,ksmm1
     do i = its,ite
       ldrag(i) = ldrag(i).or. bnv2(i,k).lt.0.
     enddo
   enddo
!
!  the low level weighted average ri is stored in usqj(1,1; ite)
!  the low level weighted average n**2 is stored in bnv2(1,1; ite)
!  this is called bnvl2 in phys_gwd_alpert_sub not bnv2
!
   kbjp1     = kbj + 1
!
   do i = its,ite
     wtkbj     = (prsl(i,kbj)-prsl(i,kbjp1))/prsi(i,1)/delks1(i)
     usqj(i,1) = wtkbj * usqj(i,kbj)
     bnv2(i,1) = wtkbj * bnv2(i,kbj)
   enddo    
!
   do k = kbjp1,ksmm1
     do i = its,ite
       rdelks    = (prsl(i,k)-prsl(i,k+1))/prsi(i,1)/delks1(i)
       bnv2(i,1) = bnv2(i,1) + bnv2(i,k) * rdelks
       usqj(i,1) = usqj(i,1) + usqj(i,k) * rdelks
     enddo   
   enddo
!
   do i = its,ite
     ldrag(i) = ldrag(i).or. bnv2(i,1).le.0.0
     ldrag(i) = ldrag(i).or. ulow(i).eq.1.0
   enddo  
!
!  set all ri low level values to the low level value
!
   kbjbeg = kbj
   if(kbj .eq. 1) kbjbeg = 2
   do k = kbjbeg,ksmm1
     do i = its,ite
       usqj(i,k) = usqj(i,1)
     enddo   
   enddo
!
!  low level density
!
   do k = kbj,ksmm1
     do i = its,ite
       rdelks = delprsi(i,k) / prsi(i,1) / delks(i)
       roll(i) = roll(i) + ro(i,k) * rdelks
     enddo   
   enddo
!
   do i = its,ite
     if (.not.ldrag(i) ) then
!
!  vector square root function - vsqrt -  used to compute bnv
!
       bnv(i) = sqrt( bnv2(i,1) )
!
!  calculate fr  froude    ---- n*hprime / u
!
       fr(i) = bnv(i) * hprime(i) / ulow(i)
!
!  continue w/ where block
!
!  calculate g   the universal flux function
!
       gf(i) = gmax * fr(i) * fr(i) / ( fr(i) * fr(i) + aj * aj )
!
!  calculate taub - (the base flux)
!  remember - the low level n is in bnv2(1,1;ite) = bnv = bnvl2
!
       taub(i) =  -xlinv * roll(i) * ulow(i) * ulow(i) *                       &
                                     ulow(i) * gf(i) / bnv(i)
!
!  calculate xn, yn
!
       xn(i) = ubar(i) / ulow(i)
       yn(i) = vbar(i) / ulow(i)
     else
       taub(i) = 0.0
       xn(i)   = 0.0
       yn(i)   = 0.0
     endif
   enddo
!
!  the call to phys_gwd_alpert_sub:  
!   taup are returned other parameters from monn
!
   call phys_gwd_alpert_sub(u1,v1,t1,vtj,usqj,kbj,kpblmax,g,                   &
               velco,bnv2,ro,taub,prslk,zl,rcs,                                &
               lat,kdt,hprime,xlinv,taup,                                      &
               ids,ide, jds,jde, kds,kde,                                      &
               ims,ime, jms,jme, kms,kme,                                      &
               its,ite, jts,jte, kts,kte)
!
#ifdef DBG
   call print_maxmin_seven(taup,ite,ime,kte,kts,kte,                           &
                          'taup after phys_gwd_alpert_sub')
#endif
!
   if(lcap.lt.kte) then
     do k = lcap+1, kte
       do i = its,ite
         sira          = prsi(i,k) / prsi(i,lcap)
         taup(i,k) = sira * taup(i,lcap)
       enddo
     enddo
   endif
!
!  fix up the level 1 (or more) stress to be linear with level
!  kpblmax and the stress at the bottom taub (if kpblmax is .gt. 1)
!
   if (kpblmax .gt. kbj) then
     do k = kbj,kpblmax-1
       do i = its,ite
         savem       = ((prsi(i,k+1) - prsi(i,kpblmax+1))/prsi(i,1)            &
                         / (prsi(i,kpblmax+1) - prsi(i,1)))
         taup(i,k+1) = taup(i,kpblmax+1)-savem*(taup(i,kbj)-taup(i,kpblmax+1))
       enddo      
     enddo
   endif
!
#ifdef DBG
   call print_maxmin_seven(taup,ite,ime,kte,kts,kte,                           &
                          'taup in phys_gwd_alpert')
#endif
!
!  keep in mind that taup is zero-ed out before each call
!  vertically difference stress for d tau / d sigma  from si to sl
!
   do k = kts,kte
     do i = its,ite
!
!  the stress in gwsdrag has been calc using -taub which now must be
!  returned to -(amplitude) below the old way for tau as -(ro*u**3/nl
!  instead of ro*uamp*k*n*hprime**2
!
       taud(i,k) = (taup(i,k+1) - taup(i,k) ) / delprsi(i,k) * prsi(i,1)
!
!  where del= prsi(i,k)-prsi(i,k+1)  (sign 'switched' in sela code)
!
     enddo
   enddo
!
#ifdef DBG
   call print_maxmin_seven(taud,ite,ime,kte,kts,kte,                           &
                          'taud in phys_gwd_alpert -- step 1')
#endif
!
!  calculate deceleration terms - dtaux,dtauy
!
!org
   do k = kts,kte
     do i = its,ite
       taud(i,k) = taud(i,k) / prsi(i,1)
     enddo
   enddo
!org
!
!  limit de-acceleration (momentum deposition ) at top to 1/2 value
!  the idea is some stuff must go out the 'top'
!
!  limit de-acceleration (momentum deposition ) at top to 1/2 value
!  the idea is some stuff must go out the 'top'
!
   do k = lcap, kte
     do i = its,ite
       taud(i,k) = taud(i,k) * factop
     enddo
   enddo 
!
#ifdef DBG
   call print_maxmin_seven(taud,ite,ime,kte,kts,kte,                           &
                          'taud in phys_gwd_alpert -- step 2')
#endif
!
!
!  *g and * by cos(lat) for mrf tendencies
!
   do k = kts,kte
     do i = its,ite
       taud(i,k) = taud(i,k) * cs * g
     enddo 
   enddo
!
#ifdef DBG
   call print_maxmin_seven(taud,ite,ime,kte,kts,kte,                           &
                          'taud in phys_gwd_alpert -- step 3')
#endif
!
!
!  if the gravity wave drag would force a critical line
!  in the lower ksmm1 layers during the next 2*deltim timestep,
!  then only apply drag until that critical line is reached.
!
   do k = kts,kte
     do i = its,ite
       dtfac(i,k) = 1.
     enddo
   enddo  
!
   do k = kts,ksmm1
     do i = its,ite
       if(taud(i,k).ne.0.)                                                     &
         dtfac(i,k)=min(dtfac(i,k),abs(velco(i,k)/(2.*deltim*rcs*taud(i,k))))
     enddo  
   enddo
!
   do k = kts,kte
     do i = its,ite
       taud(i,k) = taud(i,k)*dtfac(i,k)
     enddo  
   enddo
!
#ifdef DBG
   call print_maxmin_seven(taud,ite,ime,kte,kts,kte,                           &
                          'taud in phys_gwd_alpert -- step 4')
#endif
!
   do k = kts,kte
     do i = its,ite
       dtaux(i,k) = xn(i) * taud(i,k)
       dtauy(i,k) = yn(i) * taud(i,k)
     enddo      
   enddo
!
!  done with calculation - add it to old a and old b
!  a corresponds to dtauy term and b to dtaux
!
   do i = its,ite
     dusfc(i) = 0.
     dvsfc(i) = 0.
   enddo 
!
   do k = kts,kte
     do i = its,ite
       dudt(i,k) = dtaux(i,k) + dudt(i,k)
       dvdt(i,k) = dtauy(i,k) + dvdt(i,k)
       dusfc(i)  = dusfc(i)+dtaux(i,k)*delprsi(i,k)/prsi(i,1)
       dvsfc(i)  = dvsfc(i)+dtauy(i,k)*delprsi(i,k)/prsi(i,1)
     enddo 
   enddo
#ifdef REIGH_FRICTION
!
   do k = kte-10,kte
     do i = its,ite
       sigma = prsl(i,k)/prsi(i,1)
       damping = r_day_scale*log(sigma_cr/sigma)
       if(sigma.le.sigma_cr) then
         dudt(i,k)  = dudt(i,k) - damping * u1(i,k)
         dvdt(i,k)  = dvdt(i,k) - damping * v1(i,k)
       endif
     enddo
   enddo
#endif
!
   do i = its,ite
     dusfc(i) = -1.e3/g*rcs*prsi(i,1)*dusfc(i)
     dvsfc(i) = -1.e3/g*rcs*prsi(i,1)*dvsfc(i)
   enddo
!
!  diagnostic flag .ne.0 on, otherwise off, num(715)=output unit (=6)
!
   return
   end subroutine phys_gwd_alpert
!-------------------------------------------------------------------------------
   subroutine phys_gwd_alpert_sub(u1,v1,t1,vtj,usqj,kbj,kpblmax,g,             &
                     velco,bnvl2,ro,taub,prslk,zl,rcs,                         &
                     lat,kdt,hprimx,akwnmb,tensio,                             &
                     ids,ide, jds,jde, kds,kde,                                &
                     ims,ime, jms,jme, kms,kme,                                &
                     its,ite, jts,jte, kts,kte)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
!  subprogram:    phys_gwd_alpert_sub   performs gity wave drag computations.
!
!  abstract: performs gity wave drag computations.
!
!  program history log:
!   1987-06-03  jordan c. alpert       development
!   1989-02-01  hann-ming henry juang  change to fortran 77.
!   1991-03-03  sela-rozwodoski        cray guard code and constants -
!   2000-01-01  song-you hong          cvs version
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!  usage:    call phys_gwd_alpert_sub(u1,v1,t1,pstar,vtj,usqj,ksm,kbj,kpblmax,
!                       velco,bnvl2,roll,ro,taub,si,del,sl,sigk,rcl,
!                       lat,kdt,hprimx,akwnmb,tensio)
!  input argument list:
!     u1       - zonal      wind component *cos(lat)  m/sec at t0-dt.
!     v1       - meridional wind component *cos(lat)  m/sec at t0-dt.
!     t1       - temperature deg k at t0-dt.
!     vtj      - virtual temperature.
!     usqj     - richardson number.
!     kbj      - bottom of low level layer for phys_gwd_alpert_sub set to 1.
!     kpblmax  - bottom starting sigma level for p&s stress calculation
!     rcs      - square reciprocal of square of cos(lat).
!     lat      - latitude  number.  used only as diagnostic.
!     kdt      - time step number.  used only as diagnostic.
!     hprimx   - topographic standard deviation  (m).
!
!  output argument list:
!     velco    - comp of wind along the direction of low level layer.
!     bnvl2    - brunt-viasila freq as (n) and also (n**2).
!     roll     - low level ro.
!     ro       - density  mts.
!     taub     - surface stress.
!     prslk    - single dimension array of length kte which holds
!              - the constants of inverse sigma values raised to
!              - r/cp power from subroutine ggwdps.
!     akwnmb   - length scale.
!     tensio   - stress.
!
!  remarks: list caveats, other helpful hints or information
!
!  g l a s   m i g w d  (phys_gwd_alpert_sub)
!
   integer              ::  kbj,kpblmax,lat,kdt,                               &
                            ids,ide, jds,jde, kds,kde,                         &
                            ims,ime, jms,jme, kms,kme,                         &
                            its,ite, jts,jte, kts,kte
   real                 ::  g,rcs,akwnmb
   real                 ::  u1(ims:ime,kms:kme),v1(ims:ime,kms:kme),           &
                            t1(ims:ime,kms:kme),                               &
                            vtj(ims:ime,kms:kme),usqj(ims:ime,kms:kme),        &
                            velco(ims:ime,kms:kme-1),bnvl2(ims:ime),           &
                            ro(ims:ime,kms:kme),taub(ims:ime),                 &
                            hprimx(ims:ime),                                   &
                            prslk(ims:ime,kms:kme),zl(ims:ime,kms:kme),        &
                            tensio(ims:ime,kms:kme)
!
!  local variables
!
   logical              ::  icrilv(its:ite)
   integer              ::  i,k
   real                 ::  frocut,value,rdz
   real                 ::  hprime(its:ite),hco(its:ite),                      &
                            hsi(its:ite),crif2(its:ite),                       &
                            fro2(its:ite),bnv2(its:ite), ulow(its:ite),        &
                            vtk(its:ite,kts:kte),                              &
                            cl(kts:kte), zmean(kts:kte,2) 
!-------------------------------------------------------------------------------
!
!  only do constants first time through monnin
!  gsfc to nmc bridge constants
!
   frocut = 0.85 * 0.85
!
!  the variance of topography (comes in as std dev) on a lat pair
!
   do i = its,ite
     hprime(i) = hprimx(i) * hprimx(i)
   enddo  
!
!  constrain variance to be not greater than 160000 m**2
!     only if p&s low level stress is used - comment for gfcl low lev
!
   do i = its,ite
     value     = 1.6e5
     hprime(i) = min(hprime(i), value)
   enddo  
!
!  initialize critical level control vector bits all to zero
!
   do i = its,ite
     icrilv(i) = .false.
   enddo
!
!  bnvl2, the low level brunt-viasla frequency is not a fct of k
!  sqrt (n**2) nb: i use bnvl2 for n**2 and now n itself
!
   do i = its,ite
     value    = 0.0
     bnvl2(i) = max(bnvl2(i),value)
     bnvl2(i) = sqrt(bnvl2(i))
   enddo
!
!  set initial values for stress
!
   do k = kts,kte
     do i = its,ite
       tensio(i,k) = 0.0
       vtk(i,k)    = vtj(i,k) / prslk(i,k)
     enddo
   enddo
!
!  level loop
!  set up bottom values of stress if we are not starting from
!  from level 1 - nb (-) here by convention
!
   do k = kts,kpblmax
     do i = its,ite
       tensio(i,k) = -1. * taub(i)
     enddo
   enddo
!
   do k = kpblmax,kte-1
     do i = its,ite
       fro2(i) = 0.0
!
!  calculate squared brunt vaisala frequency at level k
!  n**2 as function of k - branch on low level <ksm - tv is used (vtj)
!  or just allow this type of operation on usqj  (ri)
!
       rdz     = 1./(zl(i,k+1) - zl(i,k))
       bnv2(i) = 2*g*rdz*(vtk(i,k+1)-vtk(i,k))/(vtk(i,k+1)+vtk(i,k))
!
#ifdef DBG
       if(i.eq.ite) print*,'i,k',i,k,bnv2(i) 
#endif
!
!  unstable layer if ri < 0 - - using ri in place of n**2
!
       icrilv(i) = icrilv(i) .or. ( usqj(i,k) .le. 0.0 )
       icrilv(i) = icrilv(i) .or. ( usqj(i,k) .lt. 0.25 )
!
!  compute critical froude    ---> means stable (icrilv=0)
!
       if ( .not. icrilv(i) ) then
         crif2(i) = 1 - .25 / usqj(i,k)
         crif2(i) = crif2(i) * crif2(i)
       else
         crif2(i) = 0.0
       endif
!
!  unstable layer if upper air vel comp along surf vel <=0 (crit lay)
!  at (u-c)=0. crit layer exists and bit vector should be set (.le.)
!
       icrilv(i) = icrilv(i) .or. velco(i,k) .le. 0.0
!
!  sqrt (n**2)  nb: i use bnv2 for n**2 and now n  itself
!
       value   = 0.
       bnv2(i) = max(bnv2(i),value)
       bnv2(i) = sqrt(bnv2(i))
!
     enddo
!
!  computing stress at surface  and 1 level up & limit max value
!
!  using taub at kbj level
!
     if(k .eq. kbj)  then
       do i = its,ite
         crif2(i) = min(crif2(i),frocut)
!
         if ( .not. icrilv(i) )  then
!
!  gfdl low level surface stress is negitive with respect to glas
!  becasue glas subtracts tendency while pbl routine adds we change taub
!  to -taub (nb but change back for pbl tendency
!
           tensio(i,k) = -1. * taub(i)
!
!  fr**2 at surface only
!
           fro2(i) = bnvl2(i) * bnvl2(i) * hprime(i)                           &
                                         / (velco(i,k) * velco(i,k))
!
         endif
       enddo
     else
!
!  in glas version there is a calculation at the boundary layer
!  which is not a "part" of the model std levels and a calc at the
!  first model layer.  the value for tensio at the bndy layer is set
!  to the value at the first model layer.  in nmc the bndy layer can
!  be ignored all together because the stress (tensio) is on sis,
!  interfaces while the deacceration is on layers sls.
!
!  compute the local froude   for the stable case - make sure
!  that by chance the projection of the local wind (u1,v1)
!  should not be smaller then 1.m/s, since it is cubed in the
!  denominator
!
       do i = its,ite
         if ( .not.  icrilv(i) ) then
           fro2(i) = bnv2(i) / ( (akwnmb * 0.5) * ( ro(i,k) + ro(i,k+1) ) *    &
                     velco(i,k) * velco(i,k) * velco(i,k) ) * tensio(i,k)
         endif
       enddo
!
     endif
!
!  compute stress at level in question for stable case
!
     do i = its,ite
       if( .not.icrilv(i) .and. fro2(i) .gt. crif2(i) ) then
!
!  fro2 changed to> from .ge.
!
         tensio(i,k+1) = tensio(i,k) * crif2(i) / fro2(i)
       endif
!
!  constant stress if crit froude   not met (.le. 6/1)
!
       if( .not. icrilv(i) .and. fro2(i) .le. crif2(i) )  then
         tensio(i,k+1) = tensio(i,k)
       endif
     enddo    
!
!  all done - pass back stress profile and vertically diff
!
   enddo
!
   return
   end subroutine phys_gwd_alpert_sub
!
!-------------------------------------------------------------------------------
