#include <define.h>
   subroutine chgr_make_header(lab,idate,fhour,ifin)        
!-------------------------------------------------------------------------------
!
! subprogram:    chgr_make_header   sets beginning of office note 85 label.   
!                                                                               
! abstract: sets first 8 bytes of office note 85 label.                         
!                                                                               
! program history log:                                                          
!   1988-01-01  sela                   initial mrf
!   2000-01-01  song-you hong          cvs version
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call chgr_make_header (lab, idate, fhour, ifin)          
!   input argument list:                                                        
!     lab      - office note 85 label.                                          
!     idate    - idate(1)=initial hour (gmt) of forecast.                       
!                idate(2)=month (1-12).                                         
!                idate(3)=day of the month.                                     
!                idate(4)=year.                                                 
!     fhour    - forecast hour.                                                 
!     ifin     - integer switch.                                                
!                when ifin.gt.0,                                                
!                first 8 bytes of lab are set to                                
!                'siginib2', 'sigin6b2', 'sigge6b2', 'siggesb2',                
!                depending on values of fhour and idate(1).                     
!                when ifin.le.0,                                                
!                first 8 bytes of lab are set to                                
!                'smsxxxb2', where xxx is the forecast hour, fhour.             
!                                                                               
!   output argument list:                                                       
!     lab      - office note 85 label with first 8 bytes set.                   
!                                                                               
!   output files:                                                               
!     output   - print file.                                                    
!
!-------------------------------------------------------------------------------
   save
!
   integer              ::  idigts(3) 
   integer              ::  idate (4)
!                                                                               
   character(len=4)     ::  lab (8) 
   character(len=4)     ::  lchars   
!-------------------------------------------------------------------------------
   ihour=fhour+0.5                                                           
   write(lchars,103) ihour                                                   
      103  format (i4)
   read(lchars,105) idigts                                                   
      105  format (1x, 3i1)
   write(lab(1),107) idigts(1)
      107  format (3x,i1)
   write(lab(2),108) idigts(2),idigts(3)
      108  format (2i1,'b2')
!                                                                               
   return                                                                    
   end 
!-------------------------------------------------------------------------------
