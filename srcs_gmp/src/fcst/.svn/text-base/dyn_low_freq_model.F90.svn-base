#include <define.h>
   subroutine dyn_low_freq_model(ifstep,thour)                                          
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains subprograms related to LFM.
!
!      [gmp_start]   --- [dyn_lfm_initialize] *
!    [gmp_integrate] --- [dyn_low_freq_model] *
!
! program history log:
!   2000-01-01  masao kanamitsu        low frequecy model option
!
!-------------------------------------------------------------------------------
   use paramodel, only : latg2_,levh_,levs_,lnt2_,lonf2_
   use varsfc, only    : lsoil_,lalbd_
   use comsfc
   use comfspec_vr
   use comfphys
   use comlfm
!-------------------------------------------------------------------------------
#ifndef NOPRINT
   write(6,*) 'doing dyn_low_freq_model for ifstep=',ifstep                             
!
#endif
   do l = 1,lnt2_                                                              
     fq(l)=fq(l)+q(l)*weight(ifstep)                                         
   enddo                                                                     
!
   do k = 1,levs_                                                              
     do l = 1,lnt2_                                                            
       fte(l,k)=fte(l,k)+te(l,k)*weight(ifstep)                              
       fdi(l,k)=fdi(l,k)+di(l,k)*weight(ifstep)                              
       fze(l,k)=fze(l,k)+ze(l,k)*weight(ifstep)                              
     enddo                                                                   
   enddo                                                                     
!
   do k = 1,levh_                                                              
     do l = 1,lnt2_                                                            
       frq(l,k)=frq(l,k)+rq(l,k)*weight(ifstep)                              
     enddo                                                                   
   enddo                                                                     
!                                                                               
   do j = 1,latg2_                                                             
     do i = 1,lonf2_                                                           
       ftsea(i,j)=ftsea(i,j)+tsea(i,j)*weight(ifstep)                        
       fsnoweq(i,j)=fsnoweq(i,j)+snoweq(i,j)*weight(ifstep)                  
       ftg3(i,j)=ftg3(i,j)+tg3(i,j)*weight(ifstep)                           
       fz0cm(i,j)=fz0cm(i,j)+z0cm(i,j)*weight(ifstep)                        
       fplantr(i,j)=fplantr(i,j)+plantr(i,j)*weight(ifstep)                  
       fcv(i,j)=fcv(i,j)+cv(i,j)*weight(ifstep)                              
       do il = 1, 4
       falbedo(i,j,il)=falbedo(i,j,il)+albedo(i,j,il)*weight(ifstep)                  
       enddo
       ff10m(i,j)=ff10m(i,j)+f10m(i,j)*weight(ifstep)                        
       fcanopy(i,j)=fcanopy(i,j)+canopy(i,j)*weight(ifstep)                  
       isl=nint(slmsk(i,j))+1                                                
       islmsk(i,j,isl)=islmsk(i,j,isl)+1                                     
       if(cvb(i,j).ne.cvb0) then                                             
         fcvb(i,j)=fcvb(i,j)+cvb(i,j)*weight(ifstep)                         
         wcvb(i,j)=wcvb(i,j)+weight(ifstep)                                  
       endif                                                                 
       if(cvt(i,j).ne.cvt0) then                                             
         fcvt(i,j)=fcvt(i,j)+cvt(i,j)*weight(ifstep)                         
         wcvt(i,j)=wcvt(i,j)+weight(ifstep)                                  
       endif                                                                 
     enddo                                                                   
   enddo                                                                     
!
   do k = 1,lsoil_                                                             
     do j = 1,latg2_                                                           
       do i = 1,lonf2_                                                         
         fsmc(i,j,k)=fsmc(i,j,k)+smc(i,j,k)*weight(ifstep)                   
         fstc(i,j,k)=fstc(i,j,k)+stc(i,j,k)*weight(ifstep)                   
       enddo                                                                 
     enddo                                                                   
   enddo                                                                     
!
   return                                                                    
   end subroutine dyn_low_freq_model
!
!-------------------------------------------------------------------------------
   subroutine dyn_lfm_initialize(ipstep,fhour)
!-------------------------------------------------------------------------------
   use paramodel, only : latg2_,levh_,levs_,lnt2_,lonf2_
   use varsfc, only : lsoil_
   use comsfc
   use comfspec_vr
   use comfphys
   use comlfm
   use comio
!-------------------------------------------------------------------------------
#include "abort.h"
!
   real,parameter       ::  cvb0=100.,cvt0=0.
!-------------------------------------------------------------------------------
!
!     dimension idate(4)
!
#ifndef NOPRINT
   write(6,*) 'dyn_lfm_initialize called.  fhour=',fhour
#endif
!
!  non-existent nlfmsgi is used for start of filtering integration
!
   read(nlfmsgi,end=999,err=999) lab
   go to 888
!
!  start of new filtering integration
!
999   continue
   ipstep=1
!
!  zero out output surface array
!
   do j = 1,latg2_
     do i = 1,lonf2_
       islmsk(i,j,1)=0
       islmsk(i,j,2)=0
       islmsk(i,j,3)=0
     enddo
   enddo
!
   do l = 1,lnt2_
     fq(l)=fq(l)+qm(l)*weight(1)
   enddo
!
   do k = 1,levs_
     do l = 1,lnt2_
       fte(l,k)=tem(l,k)*weight(1)
       fdi(l,k)=dim(l,k)*weight(1)
       fze(l,k)=zem(l,k)*weight(1)
     enddo
   enddo
!
   do k = 1,levh_
     do l = 1,lnt2_
       frq(l,k)=rm(l,k)*weight(1)
     enddo
   enddo
!
   do j = 1,latg2_
     do i = 1,lonf2_
       ftsea(i,j)=tsea(i,j)*weight(1)
       fsnoweq(i,j)=snoweq(i,j)*weight(1)
       ftg3(i,j)=tg3(i,j)*weight(1)
       fz0cm(i,j)=z0cm(i,j)*weight(1)
       fplantr(i,j)=plantr(i,j)*weight(1)
       fcv(i,j)=cv(i,j)*weight(1)
       do il = 1, 4
         falbedo(i,j,il)=albedo(i,j,il)*weight(1)
       enddo
       ff10m(i,j)=f10m(i,j)*weight(1)
       fcanopy(i,j)=fcanopy(i,j)*weight(1)
       isl=nint(slmsk(i,j))+1
       islmsk(i,j,isl)=islmsk(i,j,isl)+1
       if(cvb(i,j).ne.cvb0) then
         fcvb(i,j)=cvb(i,j)*weight(1)
         wcvb(i,j)=weight(1)
       else
         fcvb(i,j)=0.
         wcvb(i,j)=0.
       endif
       if(cvt(i,j).ne.cvt0) then
         fcvt(i,j)=cvt(i,j)*weight(1)
         wcvt(i,j)=weight(1)
       else
         fcvt(i,j)=0.
         wcvt(i,j)=0.
       endif
     enddo
   enddo
!
   do k = 1,lsoil_
     do j = 1,latg2_
       do i = 1,lonf2_
         fsmc(i,j,k)=smc(i,j,k)*weight(1)
         fstc(i,j,k)=stc(i,j,k)*weight(1)
       enddo
     enddo
   enddo
!
   return
!
!  continuation of filtering integration
!
  888 continue
   read(nlfmsgi) ipstep,idate
#ifndef NOPRINT
   write(6,*) 'ipstep,idate of filtered sig=',ipstep,idate
#endif
   read(nlfmsgi)
   read(nlfmsgi)(fq(i),i=1,lnt2_)
!
   do k = 1,levs_
     read(nlfmsgi) (fte(i,k),i=1,lnt2_)
   enddo
!
   do k = 1,levs_
     read(nlfmsgi) (fdi(i,k),i=1,lnt2_)
     read(nlfmsgi) (fze(i,k),i=1,lnt2_)
   enddo
!
   do k = 1,levh_
     read(nlfmsgi) (frq(i,k),i=1,lnt2_)
   enddo
!
   read(nlfmsfi) lab
   read(nlfmsfi) igstep,idate
#ifndef NOPRINT
   write(6,*) 'igstep,idate of filtered sfc=',igstep,idate
#endif
   if(igstep.ne.ipstep) then
     write(6,*) 'no. of steps on sig and sfc does not match'
     call MPABORT
   endif
   read(nlfmsfi) ftsea
   read(nlfmsfi) fsmc
   read(nlfmsfi) fsnoweq
   read(nlfmsfi) fstc
   read(nlfmsfi) ftg3
   read(nlfmsfi) fz0cm
   read(nlfmsfi) fcv
   read(nlfmsfi) fcvb,wcvb
   read(nlfmsfi) fcvt,wcvt
   read(nlfmsfi) falbedo
   read(nlfmsfi) islmsk
   read(nlfmsfi) fplantr
   read(nlfmsfi) fcanopy
   read(nlfmsfi) ff10m
!
   return
   end subroutine dyn_lfm_initialize
!
