#include <define.h>
   module funct_rad_cloud
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!      [cvintfx] - [cvintpx]
!      [ggavet]  - [gintp]
!
!-------------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------------
   subroutine cvintfx(cvin,cvtin,cvbin,iin,jtwidl,jin,                         &
                     cvout,cvtout,cvbout,iout,jpout,jout,                      &
                     xx,wgt,tt,bb,sum,nn,                                      &
                     ltwidl,latrd1,latinb)                                   
!-------------------------------------------------------------------------------
!
!       code bilinearly interpolates cld amt between gaussian grids--         
!       clone of ggintp for interpolation of convective cld amt (cv).         
!         special interp procedure for tops(cvt) and bots(cvb)...             
!       j = 1 is just belo n.pole, i = 1 is greenwich (then go east).         
!      iin,jin are i,j dimensions of input grid--iout,jout for output         
!      jin2,jout2=jin/2,jout/2                                                
!                                          --k.campana - june 1988            
!
!-------------------------------------------------------------------------------
   integer              ::  iin,jtwidl,jin,iout,jpout,jout,ltwidl,latrd1,latinb
   real                 ::  cvin(iin,jtwidl),cvtin(iin,jtwidl)
   real                 ::  cvbin(iin,jtwidl)            
   real                 ::  cvout(iout,jpout)                                               
   real                 ::  cvtout(iout,jpout),cvbout(iout,jpout)                           
   real                 ::  xx(iout,4),wgt(iout,4),tt(iout,4)
   real                 ::  bb(iout,4),sum(iout,4)        
   integer              ::  nn(iout)                                                        
!
   iii = iin                                                                 
   jbb = jtwidl                                                              
   jjj = jin                                                                 
   iiiout = iout                                                             
   lbb = ltwidl                                                              
   lr1 = latrd1                                                              
   do latout = 1,jpout                                                      
     lat=latout+latinb-1                                                      
     if(lat.eq.1) then
       inslat=-1  
       wgtlat=0.
     else
       inslat=lat-1
       wgtlat=1.
     endif
!
!        if output lat is poleward of input lat=1 ,then simpl average           
!          (small region and cld amt wouldn t extrapolate well)                 
!
     call cintpx(iii,jbb,jjj,iiiout,inslat,wgtlat,                            &
              cvin,cvtin,cvbin,cvout(1,latout),                               &
              cvtout(1,latout),cvbout(1,latout),                              &
               xx,wgt,tt,bb,sum,nn,lbb,lr1)                                  
   enddo
!
   return                                                                    
   end subroutine cvintfx                                                                      
!
!-------------------------------------------------------------------------------
   subroutine cintpx(iin,jtwidl,jin,iout,                                      &
                    inslat,wgtlat,                                             &
                    cv,cvt,cvb,camt,ctop,cbot,                                 &
                    xx,wgt,tt,bb,sum,nn,ltwidl,latrd1)
!-------------------------------------------------------------------------------
!
!        simpl linear interpolation of cldamt, unless only 1,2 of the
!         surrounding pts has cv. then,if output gridpt not close enuf
!         do not interpolate to it(prevents spreading of cv clds)..
!           for 1 pt convection-intrp wgt ge (.7)**2 ...
!           for 2 pt convection-sum of intrp wgt ge .45...
!              .45 used rather than .5 to give better result for
!              diagonally opposed pts...
!===>    for tops(cvt) and bots(cvb) just take average of surrounding
!         non-zero cv points.....
!         nn will be number of surrounding pts with cld (gt zero)
!---     nhsh = 1,-1 for northern,southern hemisphere
!         here instead of an extrapolation,just do a simple mean....
!
!-------------------------------------------------------------------------------
   integer              ::  iin,jtwidl,jin,iout,inslat,ltwidl,latrd1
   real                 ::  wgtlat
   real                 ::  cv(iin,jtwidl),cvt(iin,jtwidl),cvb(iin,jtwidl)
   real                 ::  camt(iout),ctop(iout),cbot(iout)
   real                 ::  xx(iout,4),wgt(iout,4),tt(iout,4)
   real                 ::  bb(iout,4),sum(iout,4)
   integer              ::  nn(iout)
!
! local 
!
   integer              ::  i,j,lonf,inth,inth1,ileft,irght,kpt,ltop,lbot
   integer              ::  iout2
   real                 ::  wgtlon
!-------------------------------------------------------------------------------
   lonf=iin/2
!
   if (inslat.lt.0) go to 600
!
   inth = mod(ltwidl + inslat + jtwidl - latrd1 - 1,jtwidl) + 1
   inth1 = mod(inth,jtwidl) + 1
!
   if (inslat.eq.jin) go to 105
!
   do i = 1,iout
     ileft=i
     irght=i+1
     if(mod(ileft,lonf).eq.0) irght=irght-lonf
     wgtlon=0.0
     wgtlat=1.0
!
!----   normalized distance from upper lat to gaussian lat
!
     xx(i,1) = cv(ileft,inth)
     xx(i,2) = cv(ileft,inth1)
     xx(i,3) = cv(irght,inth)
     xx(i,4) = cv(irght,inth1)
     wgt(i,1) = (1.e0-wgtlon)*(1.e0-wgtlat)
     wgt(i,2) = (1.e0-wgtlon)*wgtlat
     wgt(i,3) = wgtlon*(1.e0-wgtlat)
     wgt(i,4) = wgtlon*wgtlat
     tt(i,1) = cvt(ileft,inth)
     tt(i,2) = cvt(ileft,inth1)
     tt(i,3) = cvt(irght,inth)
     tt(i,4) = cvt(irght,inth1)
     bb(i,1) = cvb(ileft,inth)
     bb(i,2) = cvb(ileft,inth1)
     bb(i,3) = cvb(irght,inth)
     bb(i,4) = cvb(irght,inth1)
   enddo
!
   go to 130
!
105 continue
!
   do i = 1,iout
     ileft=i
     irght=i+1
     if(mod(ileft,lonf).eq.0) irght=irght-lonf
     wgtlon=0.0
     wgtlat=1.0
!
!----   normalized distance from upper lat to gaussian lat
!
     xx(i,1) = cv(ileft,inth)
     xx(i,3) = cv(irght,inth)
     wgt(i,1) = (1.e0-wgtlon)*(1.e0-wgtlat)
     wgt(i,2) = (1.e0-wgtlon)*wgtlat
     wgt(i,3) = wgtlon*(1.e0-wgtlat)
     wgt(i,4) = wgtlon*wgtlat
     tt(i,1) = cvt(ileft,inth)
     tt(i,3) = cvt(irght,inth)
     bb(i,1) = cvb(ileft,inth)
     bb(i,3) = cvb(irght,inth)
   enddo
   iout2 = iout / 2
   do i = 1,iout2
     ileft=i
     irght=ileft+1
     if(mod(ileft,lonf).eq.0) irght=irght-lonf
     ileft2=i+iout2
     irght2=ileft2+1
     if(mod(ileft2,lonf).eq.0) irght2=irght2-lonf
     xx(i      ,2) = cv(ileft2,inth)
     xx(i+iout2,2) = cv(ileft ,inth)
     xx(i      ,4) = cv(irght2,inth)
     xx(i+iout2,4) = cv(irght ,inth)
     bb(i      ,2) = cvb(ileft2,inth)
     bb(i+iout2,2) = cvb(ileft ,inth)
     bb(i      ,4) = cvb(irght2,inth)
     bb(i+iout2,4) = cvb(irght ,inth)
     tt(i      ,2) = cvt(ileft2,inth)
     tt(i+iout2,2) = cvt(ileft ,inth)
     tt(i      ,4) = cvt(irght2,inth)
     tt(i+iout2,4) = cvt(irght ,inth)
   enddo
!
!---      nn will be number of surrounding pts with cld (gt zero)
!
130 continue
!
   do i = 1,iout
     nn(i) = 0
   enddo
   do j = 1,4
     do i = 1,iout
       sum(i,j) = 0.e0
     enddo
   enddo
   do kpt = 1,4
     do i=1,iout
       if (xx(i,kpt).gt.0.e0) then
         nn(i) = nn(i) + 1
         sum(i,1) = sum(i,1) + wgt(i,kpt)
         sum(i,2) = sum(i,2) + tt(i,kpt)
         sum(i,3) = sum(i,3) + bb(i,kpt)
       endif
     enddo
     do i = 1,iout
       sum(i,4) = sum(i,4) + wgt(i,kpt) * xx(i,kpt)
     enddo
   enddo
   do i = 1,iout
     if (nn(i).eq.1.and.sum(i,1).gt.0.49e0) go to 17
     if (nn(i).eq.2.and.sum(i,1).ge.0.45e0) go to 17
     if (nn(i).ge.3) go to 17
!
     ctop(i) = 0.e0
     cbot(i) = 100.e0
     camt(i) = 0.e0
!
     cycle
!
17   continue
!
     ltop = sum(i,2)/nn(i) + 0.5e0
     lbot = sum(i,3)/nn(i) + 0.5e0
     lamt = sum(i,4) + 0.5e0
     ctop(i) = ltop
     cbot(i) = lbot
     camt(i) = lamt
   enddo
!
   return
!
!--- polar region-no extrapolation
!
600 continue
!
   ja = iabs(inslat)
   do i = 1,iout
     ileft=i
     irght=ileft+1
     if(mod(ileft,lonf).eq.0) irght=irght-lonf
     wgtlon=0.0
!
!----    get left point on nearest latitude
!
     xx(i,1) = cv(ileft,ja)
     xx(i,2) = cv(irght,ja)
     wgt(i,1) = 1.e0-wgtlon
     wgt(i,2) = wgtlon
     tt(i,1) = cvt(ileft,ja)
     tt(i,2) = cvt(irght,ja)
     bb(i,1) = cvb(ileft,ja)
     bb(i,2) = cvb(irght,ja)
   enddo
!
!---      nn will be number of surrounding pts with cld (gt zero)
!
   do i = 1,iout
     nn(i) = 0
   enddo
   do j = 1,4
     do i = 1,iout
       sum(i,j) = 0.e0
     enddo
   enddo
   do kpt = 1,2
     do i = 1,iout
       if (xx(i,kpt).gt.0.e0) then
         nn(i) = nn(i) + 1
         sum(i,1) = sum(i,1) + wgt(i,kpt)
         sum(i,2) = sum(i,2) + tt(i,kpt)
         sum(i,3) = sum(i,3) + bb(i,kpt)
       endif
     enddo
     do i = 1,iout
       sum(i,4) = sum(i,4) + wgt(i,kpt) * xx(i,kpt)
     enddo
   enddo
   do i = 1,iout
     if (nn(i).eq.1.and.sum(i,1).gt.0.7e0) go to 27
     if (nn(i).eq.2) go to 27
!
     ctop(i) = 0.e0
     cbot(i) = 100.e0
     camt(i) = 0.e0
     cycle
!
27   continue
!
     ltop = sum(i,2)/nn(i) + 0.5e0
     lbot = sum(i,3)/nn(i) + 0.5e0
     lamt = sum(i,4) + 0.5e0
     ctop(i) = ltop
     cbot(i) = lbot
     camt(i) = lamt
   enddo
!
   return
   end subroutine cintpx
!
!-------------------------------------------------------------------------------
   subroutine ggavet(cttin,iin,jtwidl,jin,cttout,iout,jpout,jout,              &
                     tt,sum,nn,ltwidl,latrd1,latinb)                         
!-------------------------------------------------------------------------------
!    put cloud top temperature onto fcst model grid......          *         
!       only average those points which have cld (ie temp nonzero) *         
!     j = 1 is just belo n.pole, i = 1 is greenwich (then go east).*         
!    iin,jin are i,j dimensions of input grid--iout,jout for output*         
!    jin2,jout2=jin/2,jout/2                                       *         
!-------------------------------------------------------------------------------
   real     ::  cttin(iin,jtwidl)                                               
   real     ::  cttout(iout,jpout)                                              
   real     ::  tt(iout,4),sum(iout)                                            
   integer  ::  nn(iout)                                                        
!-------------------------------------------------------------------------------
   iii = iin                                                                 
   jbb = jtwidl                                                              
   jjj = jin                                                                 
   iiiout = iout                                                             
   lbb = ltwidl                                                              
   lr1 = latrd1                                                              
   do latout = 1,jpout                                                      
     lat=latout+latinb-1                                                      
     if(lat.eq.1) then
       inslat=-1  
       wgtlat=0.
     else
       inslat=lat-1
       wgtlat=1.
     endif
!
!===>    if output lat is poleward of input lat=1 ,then simpl average           
!          (small region and cld amt wouldn t extrapolate well)                 
     call gintp(iii,jbb,jjj,iiiout,                                            &
                inslat,wgtlat,                                                 &
                cttin,cttout(1,latout),tt,sum,nn,lbb,lr1)                     
   enddo
   return                                                                    
   end subroutine ggavet     
!
!-------------------------------------------------------------------------------
   subroutine gintp(iin,jtwidl,jin,iout,                                       &
                    inslat,wgtlat,                                             &
                    ctt,cldt,tt,sum,nn,ltwidl,latrd1)                        
!-------------------------------------------------------------------------------
!    for top temp just take average of surrounding                          
!       non-zero points (these are the cloud-filled ones)....                 
!       nn will be number of surrounding pts with cld (gt zero)               
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  iin,jtwidl,jin,iout,inslat,ltwidl,latrd1
   real                 ::  wgtlat
   real                 ::  ctt(iin,jtwidl)                                                 
   real                 ::  cldt(iout)                                                      
   real                 ::  tt(iout,4),sum(iout)                                            
   integer              ::  nn(iout)                                                        
!
! local
!
   integer              ::  lonf,inth,inth1,i,ileft,irght,iout2,ileft2,irght2
   integer              ::  kpt,ja
!-------------------------------------------------------------------------------
   lonf=iin/2
!
   if (inslat.lt.0) go to 600                                                
!
   inth = mod(ltwidl + inslat + jtwidl - latrd1 - 1,jtwidl) + 1              
   inth1 = mod(inth,jtwidl) + 1                                              
!
   if (inslat.eq.jin) go to 105                                              
!
   do i = 1,iout                                                           
     ileft=i
     irght=i+1
     if(mod(ileft,lonf).eq.0) irght=irght-lonf
     tt(i,1) = ctt(ileft,inth)                                            
     tt(i,2) = ctt(ileft,inth1)                                           
     tt(i,3) = ctt(irght,inth)                                            
     tt(i,4) = ctt(irght,inth1)                                           
   enddo
!
   go to 130                                                                 
!
105 continue
!
   do i = 1,iout                                                           
     ileft=i
     irght=i+1
     if(mod(ileft,lonf).eq.0) irght=irght-lonf
     tt(i,1) = ctt(ileft,inth)                                            
     tt(i,3) = ctt(irght,inth)                                            
   enddo
   iout2 = iout / 2                                                          
   do i = 1,iout2                                                          
     ileft=i
     irght=ileft+1
     if(mod(ileft,lonf).eq.0) irght=irght-lonf
     ileft2=i+iout2
     irght2=ileft2+1
     if(mod(ileft2,lonf).eq.0) irght2=irght2-lonf
     tt(i      ,2) = ctt(ileft2,inth)                                      
     tt(i+iout2,2) = ctt(ileft ,inth)                                      
     tt(i      ,4) = ctt(irght2,inth)                                      
     tt(i+iout2,4) = ctt(irght ,inth)                                      
   enddo
!
!---      nn will be number of surrounding pts with cld (gt zero)               
!
130 continue
!
   do i = 1,iout                                                            
      nn(i) = 0                                                               
   enddo
   do i = 1,iout                                                            
      sum(i) = 0.e0                                                          
   enddo
   do kpt = 1,4                                                            
     do i = 1,iout                                                           
       if (tt(i,kpt).gt.0.e0) then                                            
         nn(i) = nn(i) + 1                                                     
         sum(i) = sum(i) + tt(i,kpt)                                           
       endif                                                                   
     enddo
   enddo
   do i = 1,iout                                                            
     if (nn(i).lt.1) then                                                    
       cldt(i) = 0.e0                                                       
     else                                                                    
       cldt(i) = sum(i) / nn(i)                                              
     end if                                                                  
   enddo
!
   return                                                                    
!
!--- polar region                                                               
!
600 continue                                                                  
!
   ja = iabs(inslat)                                                         
   do i = 1,iout                                                           
     ileft=i
     irght=ileft+1
     if(mod(ileft,lonf).eq.0) irght=irght-lonf
     tt(i,1) = ctt(ileft,ja)                                              
     tt(i,2) = ctt(irght,ja)                                              
   enddo
!
!---      nn will be number of surrounding pts with cld (gt zero)               
!
   do i = 1,iout                                                            
     nn(i) = 0                                                               
   enddo
   do i = 1,iout                                                            
     sum(i) = 0.e0                                                           
   enddo
   do kpt = 1,2                                                            
     do i = 1,iout                                                          
       if (tt(i,kpt).gt.0.e0) then                                          
         nn(i) = nn(i) + 1                                                   
         sum(i) = sum(i) + tt(i,kpt)                                         
       endif                                                                 
     enddo
   enddo
   do i = 1,iout                                                            
     if (nn(i).lt.1) then                                                    
       cldt(i) = 0.e0                                                       
     else                                                                    
       cldt(i) = sum(i) / nn(i)                                              
     end if                                                                  
   enddo
!
   return                                                                    
   end subroutine gintp
!-------------------------------------------------------------------------------
   end module funct_rad_cloud
