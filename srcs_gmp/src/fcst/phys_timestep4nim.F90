#include "define.h"
!-------------------------------------------------------------------------------
#ifdef NIM_BIN
   subroutine phys_timestep4nim(idate,deltat,initial_step,u_p,v_p,t_p)
#else
   subroutine phys_timestep4nim(idate,deltat,initial_step)
#endif
!-------------------------------------------------------------------------------
#ifdef NIM  /* NIM */
   use paramodel, only    :  itsbeg,LATG2S,LONF2S,levs_
   use comfgsm,   only    :  solsec,dthr,hdthr
   use comio,     only    :  iope
   use comfphys
   use comfver
#ifdef GFS_SFC
   use comsfc
#endif
#ifdef NIM_BIN
   use comnim
#endif
!-------------------------------------------------------------------------------
   real, intent(in)      ::  deltat
   integer, intent(in)   ::  idate(4)
   logical, intent(in)   ::  initial_step
#ifdef NIM_BIN
   real*8, intent(in)    ::  u_p(LONF2S,levs_),v_p(LONF2S,levs_)
   real*8, intent(in)    ::  t_p(LONF2S,levs_)
#endif
!-------------------------------------------------------------------------------
!
! set hours
!
   limlow=1
   inistp=1
   deltim=deltat
   stepone=.false.        ! for leap-frog in GRIMs   
!
   dtpost=0.0             ! not used, but common var. in comfver
   maxstp=48              ! maxstp=num(7)
   dthr = deltat/3600.e0  ! dthr=con(1)/3600.e0
   hdthr = 0.5 * dthr
!
!  set initial solhr
! 
   if (initial_step) then    ! first step
     if(itsbeg.eq.1) then
       shour=0.
       fhour=0.
     else                       ! itsbeg.ne.1 --> restart
       shour=deltim*(itsbeg-1.) ! shour (second)
       fhour=AINT(shour/3600.)  ! fhour (hr)
     endif
     thour=fhour
     solhr=fhour+idate(1)
     iday=solhr/24.e0
     solhr=solhr-iday*24.e0
     solsec=solhr*3600.
   endif
!
! reset times
!
   if ( .not. initial_step) then
     shour=shour+deltim
     if(inistp.eq.1) inistp=0
!
     if( .not. stepone ) then
       thour=ifix(shour/3600.)
       solsec=solsec+deltim
       fhour=thour
       solhr=solsec/3600.e0
       iday=solhr/24.e0
       solhr=solhr-iday*24.e0
     endif
     if( stepone ) then
       stepone=.false.
       solsec=solsec+deltim
       solhr=solsec/3600.0
     endif
   endif
#ifdef NIM_DBG2
   if( iope ) then
      write(6,*) 'fhour, shour, solhr, solsec = ',fhour, shour, solhr, solsec
      write(6,*) 'idate(1),iday = ',idate(1),iday
   endif
#endif
#ifdef LANDSEA
!
! set up for LADNSEA contrast exp
!
   if (initial_step) then    ! first step
     do j=1,LATG2S
       do i=1,LONF2S
         smc(i,j,1)=0.315
         stype(i,j)=7.
         vfrac(i,j)=0.7
         vtype(i,j)=7.
       enddo
     enddo
   endif
#endif
#ifdef NIM_BIN
   if (initial_step) then    ! first step
     do i = 1,LONF2S
       wor1d(i,1)=tsea(i,1)
       wor1d(i,2)=smc(i,1,1)
       wor1d(i,3)=snoweq(i,1)
       wor1d(i,4)=stc(i,1,1)
       wor1d(i,5)=tg3(i,1)
       wor1d(i,6)=z0cm(i,1)
       wor1d(i,7)=albedo(i,1,1)
       wor1d(i,8)=albedo(i,1,2)
       wor1d(i,9)=albedo(i,1,3)
       wor1d(i,10)=albedo(i,1,4)
       wor1d(i,11)=slmsk(i,1)
       wor1d(i,12)=vfrac(i,1)
       wor1d(i,13)=canopy(i,1)
       wor1d(i,14)=f10m(i,1)
       wor1d(i,15)=vtype(i,1)
       wor1d(i,16)=stype(i,1)
       wor1d(i,17)=snwdph(i,1)
       wor1d(i,18)=snoalb(i,1)
     do k=1,levs_
       wor2d(i,k,1)=u_p(i,k)
       wor2d(i,k,2)=v_p(i,k)
       wor2d(i,k,3)=t_p(i,k)
     enddo
     enddo
   endif
#endif
!
   return
#endif  /* NIM */
   end subroutine phys_timestep4nim
