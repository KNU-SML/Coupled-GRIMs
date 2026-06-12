   subroutine dfs_standard_temp(iunit,sl,to,kmo,psref)
!-------------------------------------------------------------------------------
!
! subprogram: dfs_standard_temp          
!
! abstract: interpolate standard atmospheric temperature on
!   model sigma level.
!
! program history log:
!   2006-0101     hoon park
!   2000-01-01  song-you hong          usgs topo, kagwd options
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
! 
! Assumed fine format
!   ASCII file contains
!      rec1 : number of level(ascending order)
!      rec2 : description1
!      rec3 : description2
!      rec4 : description3
!      rec5 : description4
!      rec6 : first level data (Altitude(km) Temperature(K) Pressure(pa)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                           :: iunit,kmo
   real                              :: to(kmo),sl(kmo),psref
!
   integer                           :: kmi,i,idx
   real, allocatable, dimension(:)   :: pi,ti,coef
   real                              :: po(kmo),hgt,den,dp(1)
!
   data idx /1/
!-------------------------------------------------------------------------------
!
!  open(iunit,file='standard_atm.txt')
   read(iunit,*)kmi
   allocate(pi(kmi),ti(kmi),coef(kmi))
!
   do i = 1,4
     read(iunit,*)
   enddo
!
! for ascending order
!
   do i = 1,kmi
     read(iunit,*) hgt,ti(i),pi(i),den
   enddo
!  close(iunit)

!
! pressure(hPa)
!
   do i = 1,kmi
     pi(i)=log(pi(i)*0.01)
   enddo
!
   do i = 1,kmo
     po(i)=log(sl(i)*psref)
   enddo

!
! interpolate
!
   if(IDX.NE.0) then
     call chgr_cubic_spline(to,po,kmo,ti,pi,kmi,IDX)
   else
     call chgr_linear_interp(to,po,kmo,ti,pi,kmi)
   endif
!
   deallocate(pi,ti,coef)
!
   return
   end
!-------------------------------------------------------------------------------
