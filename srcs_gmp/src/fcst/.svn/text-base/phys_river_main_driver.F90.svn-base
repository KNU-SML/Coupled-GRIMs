#include <define.h>
   subroutine phys_river_main_driver
!------------------------------------------------------------------------------
!
! subroutine: phys_river_main_driver
!
!  program history log:
!   2007-02-09  kei yoshimura          development
!   2008-07-01  jung-eun kim           initial implementation
!   2012-02-24  suryun ham             rmp implementation
!
!  compute river discharge by TRIP (Oki and Sud, 1996)
!     S(t+1) = S(t)*exp(-dc*dt)+(1-exp(-dc*dt))*I/dc
!     O = (S(t+1)-S(t))/dt+I
!     S:river storage 
!     O:river discharge
!     I:input to the flow (inflow from upstream and total runoff 
!       accumulated in a grid)
!     dc = v/l (v: flow velocity, l: flowpath length)
!
!     input runoff are taken from LSM and output river discharge etc.
!     are written to riv.ftXX in grib format.
!
!  variables:
!     gdriv: river storage [kg/m/m]
!     rflow: river discharge [kg/s]
!     imap : flow direction (0:sea, 9:rivermouth)
!           8  1  2
!           7  *  3
!           6  5  4
!     roff : converted runoff [kg/m/m/s]
!
!  initial and restart file:
!     rivi (if not exsited, all river channel storage is assumed
!           to have 10mm water)
!
!  Note1: 
!     Calculation is in a single processor
!  Note2: 
!     Map(s) of riverflow direction is very important.
!     Currently, 1deg x 1deg version (Oki and Sud, 1996) and
!     0.5deg x 0.5deg version are available.
!     The map should be placed/linked as fort.41
!  Note3: 
!     io_ and jo_ in paramodel.h are used as io2_ and jo2_.
!     t126 (io_ = 360,jo_ = 181) and t248 (io_ = 720,jo_ = 361) are
!     usable in default, but other resolution io2_ and jo2_ 
!     needs to be selected manually. => in comfriv.h and post_gauss2latlon_river.F
!
!     Grid coordination (1deg) is centered, such as 
!     (0.5E, 89.5N)  - (0.5W, 89.5N) ; i = 1,io_
!                    :               ; j = 1,jo_-1
!     (0.5E, 89.5S)  - (0.5W, 89.5S) 
!     (        dummy column        ) ; j = jo_
!
!  included h file:
!     comfriv.h (in include directory)
!  subroutine called: 
!     post_gauss2latlon: gausian to latlon-centered cordination
!
!  influenced subroutine:
!     gmp_integrate        : call this subroutine (phys_river_main_driver) and
!                            write restart file (rivi)
!     phys_main_solver     : calculate total runoff (surface+subsurface)
!     file_write_flux      : writing to riv.ftXX in grib format.
!     idsdef               : (in share directory) define output variable IDs.
!                            gdriv(199), rflow(200), imap(202), roff(203)
!
!------------------------------------------------------------------------------
#ifdef RIVER
   use paramodel, only        :  LONF2S, LATG2S, LONF2F, LATG2F
   use paramodel, only        :  LONFD, LATGD
   use paramodel, only        :  io2_,jo2_
#ifdef RMP
   use paramodel, only        :  delxo2_, delyo2_
#endif
#ifdef MP
   use commpi
#endif
   use constant, only         :  pi_, rerth_, g_
#ifndef RMP
   use module_trans, only     :  dyn_trans2output_grid
#else
   use module_trans, only     :  rmp_trans2output_grid
#endif
   use comfver
   use comfrivh
   use comcon
!
   integer                ::  len, ioerr
   real, allocatable      ::  work(:)
   real                   ::  runof(io2_,jo2_)
   real                   ::  otflw(io2_,jo2_)
   real                   ::  inflw(io2_,jo2_)
   real                   ::  gdrivo(io2_,jo2_) 
   real,save,allocatable  ::  rdest(:,:)
   real,save,allocatable  ::  area(:,:)
   real                   ::  alon(io2_),dlon(io2_)
   real                   ::  alat(jo2_),dlat(jo2_)
!
   integer,save,allocatable  ::  idest(:,:)
   integer,save,allocatable  ::  jdest(:,:)
   integer                   ::  iofs(0:9)
   integer                   ::  jofs(0:9)
   integer                   ::  idx,idy
   integer                   ::  itmp(io2_)
!                     0  1  2  3  4  5  6  7  8  9
   DATA       IOFS  / 0, 0, 1, 1, 1, 0,-1,-1,-1, 0 /
   DATA       JOFS  / 0,-1,-1, 0, 1, 1, 1, 0,-1, 0 /
!
   character c3*3
   integer  ::  ifirst
   data ifirst/0/
   integer  ::  nrvmap,nrvini
   data nrvmap,nrvini/41,42/
   real     ::  vriver
   data vriver/0.385e0/       !! v = 0.5 m/s, meandering ratio = 1.3
   save ifirst,vriver
   real     ::  gdrivall,gdrivallo,rflowall,inflwall,rivall,runofall
   real     ::  areaall,rbud
!
   len = LONFD*LATGD
   allocate(work(len))
   if(.not.allocated(rdest)) allocate(rdest(io2_,jo2_))
   if(.not.allocated(area))  allocate(area(io2_,jo2_))
   if(.not.allocated(idest)) allocate(idest(io2_,jo2_))
   if(.not.allocated(jdest)) allocate(jdest(io2_,jo2_))
!
!  Initial setting
!    
#ifdef RMP
   deltim = con(11)
#endif
!
   if (ifirst.eq.0) then
#ifdef MP
     if (mype.eq.master) then
#endif
       ifirst = 1
!
!   reading river direction map
!
       rewind (nrvmap)
       write(c3,'(i3.3)') io2_
       do j = 1,jo2_-1
         read (nrvmap,'('//c3//'i1)') (itmp(i),i = 1,io2_)
         do i = 1,io2_
            imap(i,j) = real(itmp(i))
         enddo
       enddo
!
!   identify downstream grid
!
       do j = 1,jo2_-1
         do i = 1,io2_
           idest(i,j) = i+iofs(int(imap(i,j)))
           if (idest(i,j).eq.  0  ) idest(i,j) = io2_
           if (idest(i,j).eq.io2_+1) idest(i,j) = 1
           jdest(i,j) = j+jofs(int(imap(i,j)))
         enddo
       enddo
!
!   reading initial
!
       rewind (nrvini)
       read (nrvini,iostat=ioerr) gdriv
       if(ioerr .eq. 1 ) then
         goto 101
       else
         goto 102
       endif
102    continue 
       print*,'no rivi'
       do j = 1,jo2_-1
         do i = 1,io2_
           gdriv(i,j) = 10.
         enddo
       enddo
101    continue
! 
!   identify river mouth
!
       do j = 1,jo2_-1
         do i = 1,io2_
           idx = idest(i,j)
           idy = jdest(i,j)
           if (imap(idx,idy).eq.0 .and. (imap(i,j).ne.0)) then
             imap(idx,idy) = 9
             print*,'new rivermouth',idx,idy
           endif
         enddo
       enddo
!
!   polarious location
!
       do i = 1,io2_
         alon(i) = (real(i)-0.5)/real(io2_)*2.*pi_
         dlon(i) = 1./real(io2_)
       enddo
       do j = 1,jo2_-1
         alat(j) = (0.5-(real(j)-0.5)/real(jo2_-1))*pi_
         dlat(j) = 1./real(jo2_-1)
       enddo
!        
!   calculate rdest
!
       do j = 1,jo2_-1
         do i = 1,io2_
           rdest(i,j) = 0.
           idx = idest(i,j)
           idy = jdest(i,j)
           if (imap(i,j).ge.1 .and. imap(i,j).le.8) then
             cosx = sin(alat(j))*sin(alat(idy))                                &
                   +cos(alat(j))*cos(alat(idy))                                &
                   *cos(alon(i)-alon(idx))
              rdest(i,j) = vriver/(acos(cosx)*rerth_)
           elseif (imap(i,j).eq.9) then !! assume flows toward east
             idxx = mod(i,io2_)+1
             cosx = sin(alat(j))*sin(alat(j))                                  &
                   +cos(alat(j))*cos(alat(j))                                  &
                   *cos(alon(i)-alon(idxx))
             rdest(i,j) = vriver/(acos(cosx)*rerth_)
           endif
         enddo
       enddo
!
!  calculate area
!
       do j = 1,jo2_-1
         do i = 1,io2_
           area(i,j) = 2.*rerth_*rerth_*pi_*pi_*dlon(i)*dlat(j)                &
                      *cos(alat(j))
         enddo
       enddo
!
#ifdef MP
     endif
#endif
   endif !! end of initial setting
 100  format('phys_river_main_driver first step: ',4i4,f4.0,2e10.2)
!
#ifdef MP
!
!  MPI -> single
!
#ifdef RMP
#define MPGP2F rmpgp2f
#else
#define MPGP2F mpgp2f
#endif
   call MPGP2F(trunof,LONF2S,LATG2S,work,LONF2F,LATG2F,1)
#else
   n = 0
   do j = 1,LATG2F
     do i = 1,LONF2F
       n = n+1
       work(n) = trunof(i,j)
     enddo
   enddo
#endif
!
!  single
!
#ifdef MP
   if (mype.eq.master) then
#endif
#ifdef RMP
     call rmp_trans2output_grid(work,1)
!
!  interpolate rmp to center-lonlat
!
     call post_rmp2latlon_river(work,LONFD,LATGD,runof,io2_,jo2_,delxo2_,delyo2_)
#else
     call dyn_trans2output_grid(work,1)
!
!  interpolate gaussian to center-lonlat
!
     call post_gauss2latlon(work,LONFD,LATGD,                                  &
                               180./float(io2_),90.-90./float(jo2_-1),         &
                               360./float(io2_),180./float(jo2_-1),            &
                               runof,io2_,jo2_)
#endif
!
!  reset
!        
     do j = 1,jo2_-1
       do i = 1,io2_
         gdrivo(i,j) = gdriv(i,j)
         inflw(i,j) = 0.
       enddo
     enddo
!            
!  trip
!
     do j = 1,jo2_-1
       do i = 1,io2_
         inflw(i,j) = inflw(i,j)+runof(i,j)
         idx = idest(i,j)
         idy = jdest(i,j)
         if (imap(i,j).ge.1) then
           if (rdest(i,j).eq.0.) then
             print*,'RDEST',i,j,rdest(i,j)
           endif
           gdriv(i,j) = gdriv(i,j)*exp(-rdest(i,j)*deltim)                     &
                       + (1.e0-exp(-(rdest(i,j)*deltim)))                      &
                       * runof(i,j)/rdest(i,j)
           otflw(i,j) = ((gdrivo(i,j)-gdriv(i,j))/deltim                       &
                       + inflw(i,j))*area(i,j)
           if (imap(i,j).ne.9) then
             gdriv(idx,idy) = gdriv(idx,idy)                                   &
                             + otflw(i,j)/area(idx,idy)*deltim
             inflw(idx,idy) = inflw(idx,idy)                                   &
                             + otflw(i,j)/area(idx,idy)
           endif
         else !! at seas or big lakes
           otflw(i,j) = 0.
           gdriv(i,j) = 0.
         endif
         rflow(i,j) = rflow(i,j)+otflw(i,j)*deltim !! average for output
         roff(i,j) = roff(i,j)+runof(i,j)*deltim   !! averege for output
       enddo
     enddo
!
!  budget check    
!
     gdrivall = 0.
     gdrivallo = 0.
     rflowall = 0.
     inflwall = 0.
     rivall = 0.
     runofall = 0.
     areaall = 0.
     do j = 1,jo2_-1
       do i = 1,io2_
         if (imap(i,j).ge.1) then
           areaall = areaall+area(i,j)
           gdrivall = gdrivall+gdriv(i,j)*area(i,j)
           gdrivallo = gdrivallo+gdrivo(i,j)*area(i,j)
           rflowall = rflowall+otflw(i,j)*deltim
           inflwall = inflwall+inflw(i,j)*area(i,j)*deltim
           runofall = runofall+runof(i,j)*area(i,j)*deltim
           if (imap(i,j).le.8) then
             rivall = rivall+otflw(i,j)*deltim
           endif
           rbud = gdriv(i,j)-gdrivo(i,j)-                                      &
                  (inflw(i,j)-otflw(i,j)/area(i,j))*deltim
#ifdef DBG
           if (abs(rbud).gt.0.1) then
             print'(a15,2i4,f4.0,5e15.7)',                                     &
                      'Riv Inbalance',i,j,imap(i,j),rbud,                      &
                      gdriv(i,j),gdrivo(i,j),inflw(i,j)*deltim,                &
                      otflw(i,j)/area(i,j)*deltim
           endif
#endif
         endif
       enddo
     enddo
#ifdef DBG
     print '(a20,2f20.7)','*** phys_river_main_driver_trip bud ',              &
             (gdrivall-gdrivallo-(inflwall-rflowall))/areaall,                 &
             (inflwall-runofall-rivall)/areaall
     print '(a20,e20.7)','*** trip check A   ',areaall
     print '(a20,f20.7)','*** trip check S1  ',gdrivallo/areaall
     print '(a20,f20.7)','*** trip check S2  ',gdrivall/areaall
     print '(a20,f20.7)','*** trip check O   ',rflowall/areaall
     print '(a20,f20.7)','*** trip check I   ',inflwall/areaall
     print '(a20,f20.7)','*** trip check R   ',runofall/areaall
     print '(a20,f20.7)','*** trip check Riv ',rivall/areaall
#endif
!
#ifdef MP
    endif
#endif
    deallocate(work)
#endif
!
    return
    end subroutine phys_river_main_driver
!
