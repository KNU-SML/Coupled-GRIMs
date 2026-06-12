#include <define.h>
   subroutine chgr_read_sigma(n,fhour,idate,gzi,qi,tei,dii,zei,rqi             &
      ,waves,xlayers,trun,order,realform,gencode                               &
      ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,gases                               &
      ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid                    &
      ,pdryini,dummy2,water,iret)
!-------------------------------------------------------------------------------
!
! subprogram:    chgr_read_sigma      read input sigma file for chgr.
!                                                                               
! abstract: reads a sigma file.                                                 
!                                                                               
! program history log:                                                          
!   1991-03-15  mark iredell  docblock written (prehistorical program)
!   2000-01-01  song-you hong          usgs topo, kagwd options
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call chgr_read_sigma(n,fhour,idate,gzi,qi,tei,dii,zei,rqi,iret)   
!   input argument list:                                                        
!     n        - logical unit number to read                                    
!                                                                               
!   output argument list:                                                       
!     fhour    - forecast hour                                                  
!     idate    - initial hour,month,day,year                                    
!     gzi      - orography                                                      
!     qi       - ln(psfc)                                                       
!     tei      - temperature                                                    
!     dii      - divergence                                                     
!     zei      - vorticity                                                      
!     rqi      - moisture                                                       
!     iret     - return code (0: normal, 1: end of file, 2: i/o error)          
!                                                                               
!   subprograms called:                                                         
!     dyn_sigma_setup   - set sigma values
!                                                                               
!-------------------------------------------------------------------------------
   use paramodel, only : igen_,ilonf_,ilatg_,ilevs_,ijcap_
   use comchgr, only   : mdimi,kdimi, kdimpi,lab,                              &
                         siin,slin,delin,ciin,clin,rpiin,ak,bk
   use parmchgr, only  : kdimqi
!-------------------------------------------------------------------------------
   save
!-------------------------------------------------------------------------------
   integer, parameter   ::   kdum2=21,kens=2
   integer              ::   kdum
   real, allocatable    ::   dummy(:)
   real                 ::   dummy2(kdum2),ensemble(kens)
!                                                                               
   integer              ::   idate(4)
   real                 ::   gzi(mdimi),qi(mdimi),tei(mdimi,kdimi),            &
                             dii(mdimi,kdimi),zei(mdimi,kdimi),                &
                             rqi(mdimi,kdimqi)
!-------------------------------------------------------------------------------
!                                                                               
!     spectral input data file format                                           
!          lab                                                                  
!          hour,idate(4),siin(kdimpi),slin(kdimi)                               
!          zlni qi tei dii zei rqi                                              
!                                                                               
!-------------------------------------------------------------------------------
   kdum=201-ilevs_-1-ilevs_
   allocate(dummy(kdum))
!
   rewind n                                                                  
   read(n,end=100,err=200) lab                                               
   print *,' end read lab.'                                                  
!  read(n) fhour,idate,siin,slin                                             
!
   read(n,err=201)fhour,idate,siin,slin                                        &
          ,dummy,waves,xlayers,trun,order,realform,gencode                     &
          ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,water                           &
          ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid                &
          ,pdryini,dummy2,gases
   if(water.eq.0.) then
     water=1.
     print *,'water reset to 1.'
   endif
   igases=nint(gases)
   iwater=nint(water)
   print *,'tread unit,fhour,idate=',n,fhour,idate
   print *,' number of water input = ',iwater
   print *,' number of gases input = ',igases
!
   goto 202
!
   201   continue
!
   rewind n
   read(n) lab
   read(n,err=201)fhour,idate,siin,slin
!
   do i = 1,kdum
     dummy(i)=0.
   enddo
   waves=ijcap_
   xlayers=ilevs_
   trun=1.
   order=2.
   realform=1.
   gencode=igen_
   rlond=ilonf_
   rlatd=ilatg_
   rlonp=ilonf_
   rlatp=ilatg_
   rlonr=ilonf_
   rlatr=ilatg_
   gases=0.
   water=1.
   pdryini=0.
   subcen=0.
   do i = 1,kens
     ensemble(i)=0.
   enddo
   ppid=0.
   slid=0.
   vcid=0.
   vmid=0.
   vtid=0.
   do k = 1,kdum2
     dummy2(k)=0.
   enddo
!
   igases=nint(gases)
   iwater=nint(water)
!
   print *,'tread old format unit,fhour,idate=',n,fhour,idate
   print *,' number of water input  = ',iwater
   print *,' number of gases input = ',igases
!
   202   continue
!
   print *,' end read siin'                                                  
!
   rmaxsiin=siin(1)
   do k = 1,kdimpi
     rmaxsiin=max(rmaxsiin,siin(k))
   enddo
!
   if(rmaxsiin.gt.1. .or. rmaxsiin.eq.0.) then
     call chgr_hybrid_coeff(kdimi,siin,slin,ak,bk)
     print *,' end call chgr_hybrid_coeff' 
   else
     call dyn_sigma_setup(ciin,siin,delin,slin,clin,rpiin)
     print *,' end call dyn_sigma_setup'  
   endif
!
   read(n) gzi                                                               
   print *,' end read gzi.'                                                  
   read(n) qi                                                                
   print *,' end read qi .'                                                  
!                                                                               
   do k = 1,kdimi                                                           
     read(n)(tei(i,k),i=1,mdimi)                                               
   enddo
   do k = 1,kdimi                                                           
     read(n)(dii(i,k),i=1,mdimi)                                               
     read(n)(zei(i,k),i=1,mdimi)                                               
   enddo
!                                                                               
   do k = 1,kdimqi                                                          
     read(n)(rqi(i,k),i=1,mdimi)                                               
   enddo
!                                                                               
   iret=0
   return 
!                                                                               
   100   continue 
!
   iret=1                                                                    
   return                                                                    
!
200   continue
!
   iret=2
!
   deallocate(dummy)
   return
!                                                                               
   end
!-------------------------------------------------------------------------------
