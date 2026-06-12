#include <define.h>
   subroutine sfc_land_type(vegtype,soiltype,slimsk,idim,jdim)
!-------------------------------------------------------------------------------
   real                 ::  vegtype(idim*jdim),soiltype(idim*jdim),            &
                            slimsk(idim*jdim)
   real, parameter      ::  vegdflt=1.,soidflt=1.
!
   vegocn = 0.
#ifndef USGS
   soiocn = 0.
   soiice = 9.
#else
   soiocn = 14.
   soiice = 16.
#endif
!
!  correction of soiltype and soiltype over ice
!
   do j = 1,jdim
     do i = 1,idim
       ij = (j-1)*idim+i
#ifdef USGS
       if(soiltype(ij).eq.0.) then
         if(slimsk(ij).eq.0..or.slimsk(ij).eq.1.) then
           soiltype(ij) = soiocn
         endif
       endif
#endif
       if(slimsk(ij).eq.2.) then
         soiltype(ij) = soiice
       endif
     enddo
   enddo
!
!  make sure that the soil type and veg type are non-water over land
!
   do j = 1,jdim
     do i = 1,idim
       ij = (j-1)*idim+i
!
!  check for vegitation type
!
       if(slimsk(ij).eq.1.and.vegtype(ij).eq.vegocn) then
         do ii = i+1,idim
           iij = (j-1)*idim+ii
           if(vegtype(iij).ne.vegocn) then
             vegr = vegtype(iij)
             iincr = ii-i
#ifdef DBG
             print *,'vegr j=',j,' i=',i,' ii=',ii,' vegr=',vegr,         &
                          ' iincr=',iincr
#endif
             !
             go to 10
             !
           endif
         enddo
         iincr = 0
10       continue
         do ii = i-1,1,-1
           iij = (j-1)*idim+ii
           if(vegtype(iij).ne.vegocn) then
             vegl = vegtype(iij)
             iincl = i-ii
#ifdef DBG
             print *,'vegl j=',j,' i=',i,' ii=',ii,' vegl=',vegl,         &
                          ' iincl=',iincr
#endif
             !
             go to 20
             !
           endif
         enddo
         iincl = 0
20       continue
!
         if(iincr.eq.0.and.iincl.eq.0) then
           do jj = j+1,jdim
             iij = (jj-1)*idim+i
             if(vegtype(iij).ne.vegocn) then
               vegu = vegtype(iij)
               jincu = jj-j
#ifdef DBG
               print *,'vegu j=',j,' i=',i,' jj=',jj,' vegu=',vegu,      &
                             ' jincu=',jincu
#endif
               !
               go to 15
               !
             endif
           enddo
           jincu = 0
15         continue
           do jj = j-1,1,-1
             iij = (jj-1)*idim+i
             if(vegtype(iij).ne.vegocn) then
               vegl = vegtype(iij)
               jincl = j-jj
#ifdef DBG
               print *,'vegl j=',j,' i=',i,' jj=',jj,' vegl=',vegl,      &
                            ' jincl=',jincl
#endif
               !
               go to 25
               !
             endif
           enddo
           jincl = 0
25         continue
           if(jincu.ne.0.and.jincl.ne.0.and.jincu.le.jincl) then
             vegtype(ij) = vegu
#ifdef DBG
             print *,'vegtype on land-water filled with vegu=',vegu
             print *,'at i=',i,' j=',j
#endif
           elseif(jincu.ne.0.and.jincl.ne.0.and.jincu.gt.jincl) then
             vegtype(ij) = vegl
#ifdef DBG
             print *,'vegtype on land-water filled with vegl=',vegl
             print *,'at i=',i,' j=',j
#endif
           elseif(jincu.ne.0.and.jincl.eq.0) then
             vegtype(ij) = vegu
#ifdef DBG
             print *,'vegtype on land-water filled with vegu=',vegu
             print *,'at i=',i,' j=',j
#endif
           elseif(jincu.eq.0.and.jincl.ne.0) then
             vegtype(ij) = vegl
#ifdef DBG
             print *,'vegtype on land-water filled with vegl=',vegl
             print *,'at i=',i,' j=',j
#endif
           elseif(jincu.eq.0.and.jincl.eq.0) then
             print *,'WARNING: no nearby veg land point found'
             vegtype(ij) = vegdflt
             print *,'leave it alone and continue'
           endif
         elseif(iincr.ne.0.and.iincl.ne.0.and.iincr.le.iincl) then
           vegtype(ij) = vegr
#ifdef DBG
           print *,'vegtype on land-water filled with vegr=',vegr
           print *,'at i=',i,' j=',j
#endif
         elseif(iincr.ne.0.and.iincl.ne.0.and.iincr.gt.iincl) then
           vegtype(ij) = vegl
#ifdef DBG
           print *,'vegtype on land-water filled with vegl=',vegl
           print *,'at i=',i,' j=',j
#endif
         elseif(iincr.ne.0.and.iincl.eq.0) then
           vegtype(ij) = vegr
#ifdef DBG
           print *,'vegtype on land-water filled with vegr=',vegr
           print *,'at i=',i,' j=',j
#endif
         elseif(iincr.eq.0.and.iincl.ne.0) then
           vegtype(ij) = vegl
#ifdef DBG
           print *,'vegtype on land-water filled with vegl=',vegl
           print *,'at i=',i,' j=',j
#endif
         endif
       endif
!
!  check for soil type
!
       if(slimsk(ij).eq.1.and.soiltype(ij).eq.soiocn) then
         do ii = i+1,idim
           iij = (j-1)*idim+ii
           if(soiltype(iij).ne.soiocn) then
             soir = soiltype(iij)
             iincr = ii-i
             !
             go to 30
             !
           endif
         enddo
         iincr = 0
30       continue
         do ii = i-1,1,-1
           iij = (j-1)*idim+ii
           if(soiltype(iij).ne.soiocn) then
             soil = soiltype(iij)
             iincl = i-ii
             !
             go to 40
             !
           endif
         enddo
         iincl = 0
40       continue
         if(iincr.eq.0.and.iincl.eq.0) then
           do jj = j+1,jdim
             iij = (jj-1)*idim+i
             if(soiltype(iij).ne.soiocn) then
               soiu = soiltype(iij)
               jincu = jj-j
#ifdef DBG
               print *,'soiu j=',j,' i=',i,' jj=',jj,' soiu=',soiu,      &
                             ' jincu=',jincu
#endif
               !
               go to 35
               !
             endif
           enddo
           jincu = 0
35         continue
           do jj = j-1,1,-1
             iij = (jj-1)*idim+i
             if(soiltype(iij).ne.soiocn) then
               soil = soiltype(iij)
               jincl = j-jj
#ifdef DBG
               print *,'soil j=',j,' i=',i,' jj=',jj,' soil=',soil,      &
                             ' jincl=',jincl
#endif
               !
               go to 45
               !
             endif
           enddo
           jincl = 0
45         continue
           if(jincu.ne.0.and.jincl.ne.0.and.jincu.le.jincl) then
             soiltype(ij) = soiu
#ifdef DBG
             print *,'soiltype on land-water filled with soiu=',soiu
             print *,'at i=',i,' j=',j
#endif
           elseif(jincu.ne.0.and.jincl.ne.0.and.jincu.gt.jincl) then
             soiltype(ij) = soil
#ifdef DBG
             print *,'soiltype on land-water filled with soil=',soil
             print *,'at i=',i,' j=',j
#endif
           elseif(jincu.ne.0.and.jincl.eq.0) then
             soiltype(ij) = soiu
#ifdef DBG
             print *,'soiltype on land-water filled with soiu=',soiu
             print *,'at i=',i,' j=',j
#endif
           elseif(jincu.eq.0.and.jincl.ne.0) then
             soiltype(ij) = soil
#ifdef DBG
             print *,'soiltype on land-water filled with soil=',soil
             print *,'at i=',i,' j=',j
#endif
           elseif(jincu.eq.0.and.jincl.eq.0) then
             print *,'WARNING: no nearby soil land point found'
             soiltype(ij) = soidflt
             print *,'leave it alone and continue'
           endif
         elseif(iincr.ne.0.and.iincl.ne.0.and.iincr.le.iincl) then
           soiltype(ij) = soir
#ifdef DBG
           print *,'soiltype on land-water filled with ',soir
           print *,'at i=',i,' j=',j
#endif
         elseif(iincr.ne.0.and.iincl.ne.0.and.iincr.gt.iincl) then
           soiltype(ij) = soil
#ifdef DBG
           print *,'soiltype on land-water filled with soil=',soil
           print *,'at i=',i,' j=',j
#endif
         elseif(iincr.ne.0.and.iincl.eq.0) then
           soiltype(ij) = soir
#ifdef DBG
           print *,'soiltype on land-water filled with soir=',soir
           print *,'at i=',i,' j=',j
#endif
         elseif(iincr.eq.0.and.iincl.ne.0) then
           soiltype(ij) = soil
#ifdef DBG
           print *,'soiltype on land-water filled with soil=',soil
           print *,'at i=',i,' j=',j
#endif
         endif
       endif
     enddo
   enddo
!
   return
   end subroutine sfc_land_type
!-------------------------------------------------------------------------------
