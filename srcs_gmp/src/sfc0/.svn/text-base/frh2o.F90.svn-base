!
   function frh2o (tkelv,smc,sh2o,smcmax,bexp,psis)
!-------------------------------------------------------------------------------
!
! function: frh2o
!
! abstract: calculate amount of supercooled liquid soil water content if
! temperature is below 273.15k (t0).  requires newton-type iteration to
! solve the nonlinear implicit equation given in eqn 17 of koren et al
! (1999, jgr, vol 104(d16), 19569-19585).
!
! new version (june 2001): much faster and more accurate newton
! iteration achieved by first taking log of eqn cited above -- less than
! 4 (typically 1 or 2) iterations achieves convergence.  also, explicit
! 1-step solution option for special case of parameter ck=0, which
! reduces the original implicit equation to a simpler explicit form,
! known as the "flerchinger eqn". improved handling of solution in the
! limit of freezing point temperature t0.
!
! usage:    call frh2o(tkelv,smc,sh2o,smcmax,bexp,psis)
!   input argument list:
!     tkelv    - temperature (kelvin)
!     smc      - total soil moisture content (volumetric)
!     sh2o     - liquid soil moisture content (volumetric)
!     smcmax   - saturation soil moisture content (from noah_read_parameter)
!     b        - soil type "b" parameter (from noah_read_parameter)
!     psis     - saturated soil matric potential (from noah_read_parameter)
!
!   output argument list:
!   frh2o      - supercooled liquid water content
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real, parameter   ::  ck = 8.0
!      real, parameter  ::  ck = 0.0
   real, parameter   ::  blim = 5.5
   real, parameter   ::  error = 0.005
   real, parameter   ::  hlice = 3.335e5
   real, parameter   ::  gs = 9.81
   real, parameter   ::  dice = 920.0
   real, parameter   ::  dh2o = 1000.0
   real, parameter   ::  t0 = 273.15
!
   real              ::  bexp
   real              ::  bx
   real              ::  denom
   real              ::  df
   real              ::  dswl
   real              ::  fk
   real              ::  frh2o
   real              ::  psis
   real              ::  sh2o
   real              ::  smc
   real              ::  smcmax
   real              ::  swl
   real              ::  swlk
   real              ::  tkelv
   integer           ::  nlog
   integer           ::  kcount
!-------------------------------------------------------------------------------
!
! limits on parameter b: b < 5.5  (use parameter blim)
! simulations showed if b > 5.5 unfrozen water content is
! non-realistically high at very low temperatures.
!
   bx = bexp
   if (bexp .gt. blim) bx = blim
! 
! initializing iterations counter and iterative solution flag.
!
   nlog=0
   kcount=0
! 
!  if temperature not significantly below freezing (t0), sh2o = smc
!
   if (tkelv .gt. (t0 - 1.e-3)) then                                           
     frh2o = smc                                                              
   else
     if (ck .ne. 0.0) then
!-------------------------------------------------------------------------------
! option 1: iterated solution for nonzero ck
! in koren et al, jgr, 1999, eqn 17
!-------------------------------------------------------------------------------
! initial guess for swl (frozen content)
!
       swl = smc-sh2o
!
! keep within bounds.
!
       if (swl .gt. (smc-0.02)) swl = smc-0.02
       if (swl .lt. 0.) swl = 0.
! 
!  start of iterations
!
       do while ( (nlog .lt. 10) .and. (kcount .eq. 0) )
         nlog = nlog+1
         df = alog(( psis*gs/hlice ) * ( ( 1.+ck*swl )**2. ) *                 &
         ( smcmax/(smc-swl) )**bx) - alog(-(tkelv-t0)/tkelv)
         denom = 2. * ck / ( 1.+ck*swl ) + bx / ( smc - swl )
         swlk = swl - df/denom
! 
! bounds useful for mathematical solution.
!
         if (swlk .gt. (smc-0.02)) swlk = smc - 0.02
         if (swlk .lt. 0.) swlk = 0.
! 
! mathematical solution bounds applied.
! 
         dswl = abs(swlk-swl)
         swl = swlk
! 
! if more than 10 iterations, use explicit method (ck=0 approx.)
! when dswl less or eq. error, no more iterations required.
!
         if ( dswl .le. error )  then
           kcount = kcount+1
         endif
       end do
! 
!  end of iterations
! 
! bounds applied within do-block are valid for physical solution.
! 
       frh2o = smc - swl
! 
! end option 1
! 
     endif
!-------------------------------------------------------------------------------
! option 2: explicit solution for flerchinger eq. i.e. ck=0
! in koren et al., jgr, 1999, eqn 17
! apply physical bounds to flerchinger solution
!-------------------------------------------------------------------------------
     if (kcount .eq. 0) then
!lu.....comment out the following line to shorten the standard output
!clu      print*,'flerchinger used in new version. iterations=',nlog
       fk = (((hlice/(gs*(-psis)))*                                            &
              ((tkelv-t0)/tkelv))**(-1/bx))*smcmax
       if (fk .lt. 0.02) fk = 0.02
       frh2o = min (fk, smc)
! 
! end option 2
! 
     endif
   endif
! 
! end function frh2o
! 
   return
   end
!-------------------------------------------------------------------------------
