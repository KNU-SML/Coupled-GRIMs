!
   program fgb
!-------------------------------------------------------------------------------
!
! program: fgb
!
! program history log:
!   2005-01-01  song-you hong          development
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatble with namelist input
!
!-------------------------------------------------------------------------------
!
!  set defaults and read namelist
!
   use paramodel, only : para_init
!-------------------------------------------------------------------------------
   namelist/namfgb/ fhs,fhe,fhinc
   read(5,namfgb)
!
   call para_init
!
   nflx = 21
   nfgb = 51
!
!  read first sigma header record
!
   fh = fhs
   do while (fh.le.fhe)
     call fgb1(nflx,nfgb,fh)
     close(nflx)
     close(nfgb)
     print *,' finish fgb for fh = ',fh
     fh = fh+nint(fhinc)
   enddo
!
   stop
   end
!-------------------------------------------------------------------------------
