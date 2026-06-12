#include "define.h"
   subroutine sfc_read_file_driver(nsfc,fn,sfcftyp,                            &
                    labs,iy,im,id,ih,fh,                                       &
                    sfcfcs,idim,jdim,ioflag)
!-------------------------------------------------------------------------------
!
! subroutine: sfc_read_file_driver
!
! abstract: read sfc file
!
!  nsfc   : integer unit number for io (input)
!  fn     : input sfc file name
!  sfcftyp: ['osu1','osu2','noa1' 'vic1'] (input)
!  iy     : integer year of initial state (input/output)
!  im     : integer month of intiial state (input/output)
!  id     : integer day of intiial state (input/output)
!  ih     : integer hour of intiial state (input/output)
!  fh     : real forecast hour (input/output)
!  sfcfcs : real output/input array (input/output)
!  idim   : integer x-dimension of sfcfcs. partial dim for mp (input)
!  jdim   : integer y-dimension of sfcfcs. partial dim for mp (input)
!  ioflag : integer  =0 read, =1 write (input)
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use varsfc
#if defined(MP)
   use commpi          ! mype,master
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#include "abort.h"
   character(len=128)               ::  fn
   character(len=4)                 ::  sfcftyp
   character(len=32)                ::  labs
   character(len=8), allocatable    ::  gvar(:)
   real                             ::  fh
   real                             ::  fhour
   integer                          ::  idim,jdim
   integer                          ::  nsfc
   integer                          ::  ioflag
   integer                          ::  idate(4)
   real                             ::  sfcfcs(idim*jdim,*)
!
   integer                          ::  iy,im,id,ih
   integer                          ::  nrecs,maxlev
   integer                          ::  m,n,nchar
   integer, allocatable             ::  lev(:)
!
   if(fn(1:6).ne.'linked') then
#ifdef MP
     if( mype.eq.master ) then
#endif
       call sfc_grib_numchar(fn,nchar)
       do n = nchar+1,128
         fn(n:n)=' '
       enddo
       close(nsfc)
       open(unit=nsfc,file=fn(1:nchar),form='unformatted',err=811)
       !
       go to 810
       !
811    continue
       write(6,*) ' error in opening file ',fn(1:nchar)
       call MPABORT
810    continue
#ifdef MP
     endif
#endif
   endif
!
   if(ioflag.eq.0) then
#ifdef MP
     if( mype.eq.master ) then
#endif
       rewind nsfc
       read(nsfc) labs
       read(nsfc) fhour,idate
       if ( fhour .ne. fh ) then
         print*,'fcst hour does not match:sfc=',fhour,'sig=',fh
         call MPABORT
       endif
#ifdef MP
     endif
#endif
#ifdef MP
#ifndef RMP
     call mpbcastc(labs,32)
     call mpbcastr(fh,1)
     call mpbcasti(idate,4)
#else
     call rmpbcastc(labs,32)
     call rmpbcastr(fh,1)
     call rmpbcasti(idate,4)
#endif
#endif
     iy=idate(4)
     im=idate(2)
     id=idate(3)
     ih=idate(1)
#ifndef NOPRINT
#ifdef MP
     if(mype.eq.master) then
#endif
99     format(1h ,'fh, idate=',f10.1,2x,4(1x,i4))
       print *,'sfc field read in from unit=',nsfc
       print 99,fh,idate
#ifdef MP
     endif
#endif
#endif
   else
#ifdef MP
     if( mype.eq.master ) then
#endif
       write(nsfc) labs
#ifdef MP
     endif
#endif
     idate(4)=iy
     idate(2)=im
     idate(3)=id
     idate(1)=ih
#ifdef MP
     if( mype.eq.master ) then
#endif
       write(nsfc) fh, idate
#ifdef MP
     endif
#endif
#ifdef DBG
#ifdef MP
     if(mype.eq.master) then
#endif
       print *,'sfc field write to unit=',nsfc
       print 99,fh,idate
#ifdef MP
     endif
#endif
#endif
   endif
!
!  get sfc file properties
!
   call sfcfld(sfcftyp,0,nrecs,lev,gvar,maxlev)
   allocate (lev(nrecs),gvar(nrecs))
   call sfcfld(sfcftyp,1,nrecs,lev,gvar,maxlev)
!
   m=1
   do n = 1,nrecs
#ifdef SMP_RA2SFC
     call sfc_read_file(nsfc,gvar(n),lev(n),sfcfcs(1,m),idim,jdim,ioflag)
#else
     call sfc_read_file(nsfc,gvar(n),lev(n),sfcfcs(1,m),ioflag)
#endif
     m=m+lev(n)
   enddo
!
   deallocate (lev,gvar)
!
#ifdef MP
   if( mype.eq.master ) then
#endif
#ifndef NOPRINT
     print *,' sfc file io completed.'
#endif
#ifdef MP
   endif
#endif
!
   return
   end subroutine sfc_read_file_driver
!-------------------------------------------------------------------------------
