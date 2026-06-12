#define NCPUS ncpus
   subroutine dyn_cpu_number(ncpus)                                                  
!-------------------------------------------------------------------------------
   use paramodel, only : ncpus_
!-------------------------------------------------------------------------------
!
! subprogram: dyn_cpu_number     gets environment number of cpus                    
!                                                                               
! abstract: gets and returns the environment variable ncpus,                    
!   designating the number of processors over which to parallelize.             
!                                                                               
! program history log:
!   1994-08-19  iredell                                                           
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call dyn_cpu_number(ncpus)  
!   output arguments:                                                           
!     ncpus        integer number of cpus                                       
!                                                                               
! subprograms called:                                                           
!   getenv       get environment variable                                       
!                                                                               
!-------------------------------------------------------------------------------
   integer              ::  getenv
   character(len=8)     ::  cncpus                                                        
!
   ncpus=ncpus_                                                         
#ifdef GETENV
   iret=getenv('ncpus',cncpus)                                        
   if(iret.eq.1) then                                                   
     read(cncpus,'(bn,i8)',iostat=ios) ncpus                            
     ncpus=max(ncpus,1)                                                 
#ifndef NOPRINT
     print *,'ncpus=',ncpus                                             
#endif
   endif                                                                
#endif
#ifdef PXFGETENV
   iret=0                                                             
   call pxfgetenv('ncpus',5,cncpus,linval,jret)                       
   if(jret.eq.0) iret=1                                               
#endif
!
   return                                                                    
   end                                                                       
