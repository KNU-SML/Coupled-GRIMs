#include <define.h>
   subroutine print_maxmin_seven(a,len,lenx,k,k1,k2,ch) 
!-------------------------------------------------------------------------------
! subprogram documentation block
!                .      .    .
! abstract:  do print maximum and minimum of a given array.
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
! usage:    call  print_maxmin_seven(a,len,k,k1,k2,ch)
!   input argument list:
!     a    - array for computing max and min (len,k)
!     len   - the first dimension of a
!     k    - the second dimension of a
!     k1   - lower limit of second dimension to print
!     k2   - upper limit to print
!     ch    - charcter string to print
!                 fpl and fml
!
!   output argument list:
!
!   input files: none
!
!   output files:
!     standard output
!
!   subprograms called:
!     intrinsic functions: amax1 amin1
!
!   remark: none
!
!-------------------------------------------------------------------------------
   real       ::  a(lenx,k)                                                        
   character  ::  ch*(*)                                                          
!                                                                               
   print   *,ch
!
   do j = k1,k2                                                            
#ifdef SMP
     aa1 = a(1,j)                                                            
     aa2 = a(2,j)                                                            
     print 100,aa1,aa2,aa2-aa1,j 
100  format('          (1)= ',e20.10,'     (2)= ',e20.10,                     &
             '     (D)= ',e20.10,'    k =',i3)
#else
     aamax = a(1,j)                                                            
     aamin = a(1,j)                                                            
     do m = 1,len                                                             
       aamax = max( aamax, a(m,j) )                                              
       aamin = min( aamin, a(m,j) )                                              
     enddo
!    write(0,*)ch,' has max=',aamax,' min=',aamin,' at k=',j                   
     print   100,aamax,aamin,j                   
!    print*,'     max=',aamax,' min=',aamin,' at k=',j                   
100  format('          max= ',e20.10,'     min= ',e20.10,'    k =',i3)
#endif
   enddo
!
   return                                                                    
   end subroutine print_maxmin_seven
!-------------------------------------------------------------------------------
