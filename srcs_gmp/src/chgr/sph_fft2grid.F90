#include <define.h>
   subroutine sph_fft2grid (a, b, lot, isign)   
!-------------------------------------------------------------------------------
   use paramter, only : idim
   use comchgr, only : idimt, mwvp2
!-------------------------------------------------------------------------------
   save                                                                      
   integer, parameter   ::  lotmin=32,lotmax=64,ncpu=8
   real                 ::  a( idim ,lot)
   real                 ::  b( idim ,lot)
!                                                                               
   integer, save           ::  ifax(100)
   real, allocatable,save  ::  trigs(:,:)                                       
   integer, save           ::  itest
!                                                                               
   real                    ::  work(idim,lotmax,4), al(idimt,lotmax)
!                                                                               
   data itest/0/                                                            
!
   if(.not.allocated(trigs)) allocate(trigs(idim,2))
!
   if (itest.eq.0)  then                                                   
     itest=1                                                             
     jump = idim *2                                                      
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
     call fftfax (idim,ifax,trigs)                                       
#endif
#ifdef ASLES
#undef DEFAULT
!    print *,'new fft initializaiton starts'
     call ldfrmfb(idim,0,0,0,0,0,ifax,trigs,0,ierr)
#endif
#ifdef DEFAULT
     call    fax (ifax,idim,3)                                           
     call fftrig (trigs,idim,3)                                          
#endif
#ifdef ASLES
     if (ierr .ge. 3000)  print 120
     if (ierr .ge. 3000)  stop
#else
     if (ifax(1) .eq. -99)  print 120                                    
     if (ifax(1) .eq. -99)  stop     
#endif                                    
120  format (' error in fft idim .   idim  not factorable. ')            
     print 140, idim                                                     
140  format (' fftfax called in fft idim .  lonf = ',i4)                 
   endif                                                                     
!                                                                               
   if (isign .eq.  1)  then                                                  
!                                                                               
!  multiple fast fourier transform - synthesis.  isign=1                        
!  good for zonal wave number  80 .                                             
!                                                                               
!     dimension a( idim ,lot)                                                   
!                                                                               
!   input - lot sets of complex coefficients in                                 
!           a(1,j), a(2,j), ..., a( mwavep *2,j), j=1,...,lot.                  
!           a( mwavep *2+1,j), ..., a( idim ,j), j=1,...,lot are not set        
!           before call fft idim .                                              
!                                                                               
!  output - lot sets of grid values in                                          
!           a(1,j), a(2,j), ..., a( idim ,j), j=1,...,lot.                      
!                                                                               
     nlot=max0(lot/ncpu,lotmin)                                                
     nlot=min0(nlot    ,lotmax)                                                
!
     do i = 1,lot,nlot                                                       
       lots = min0(nlot, lot-i+1) 
       do j = i,i+lots-1                                                       
         do l=1, mwvp2                                                         
           al(l,j-i+1) = a(l,j)                                               
         enddo
         do l= mwvp2+1 , idimt                                                 
           al(l,j-i+1) = 0.0                                                  
         enddo
       enddo
!                                                                               
!     call fft for systhesis.                                              
!                                                                               
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
         call rfftmlt(al,work,trigs,ifax,1,jump, idim ,lots,1)                   
#endif
#ifdef ASLES
#undef DEFAULT
!      print *,'new fft, fft starts'
       call ldfrmbf(idim,lots,al,1,jump,-1,ifax,trigs,work,ierr)
       do jasl=1,lots
         do iasl=1,idim+2
           al(1+(iasl-1)*1+(jasl-1)*jump)=                                 &
           al(1+(iasl-1)*1+(jasl-1)*jump)/dble(idim)
         enddo
       enddo
#endif
#ifdef DEFAULT
         call fft99m (al,work,trigs,ifax,1,jump, idim ,lots,1)
#endif

#ifdef ASLES
       if (ierr .ge. 3000)  print 150
       if (ierr .ge. 3000)  stop
150    format (' error in asles FFT ')
#endif
       do j = i,i+lots-1                                                       
         do l=1, idim                                                          
               a(l,j) = al(l,j-i+1) 
         enddo
       enddo
     enddo
     do j = 1,lot                                                            
       do i = 1,idim                                                           
            b(i,j)=a(i,j)                                                             
       enddo
     enddo
!                                                                               
   endif                                                                     
!                                                                               
   if (isign .eq. -1)  then                                                  
!                                                                               
!  multiple fast fourier transform - analysis.  isign=-1                        
!  good for zonal wave number  80 .                                             
!                                                                               
!     dimension a( idim ,lot), b( idim ,lot)                                    
!                                                                               
!   input - lot sets of grid values in                                          
!           a(1,j), a(2,j), ..., a( idim ,j), j=1,...,lot.                      
!           a array is not changed by subroutine fft idim .                     
!                                                                               
!  output - lot sets of complex coefficients in                                 
!           b(1,j), b(2,j), ..., b( mwavep *2,j), j=1,...,lot.                  
!           b( mwavep *2+1,j), ..., b( idim ,j), j=1,...,lot are not set        
!                                                                               
      nlot=max0(lot/ncpu,lotmin)                                                
      nlot=min0(nlot    ,lotmax)                                                
!
      do i = 1,lot,nlot                                                       
         lots = min0(nlot, lot-i+1) 
         do j = i,i+lots-1         
            do l=1, idim          
               al(l,j-i+1) = a(l,j)
            enddo
            do l= idim+1 , idimt   
               al(l,j-i+1) = 0.0  
            enddo
         enddo
!                                                                               
!       call fft for analysis.                                               
!                                                                               
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
         call rfftmlt(al,work,trigs,ifax,1,jump, idim ,lots,-1)  
#endif
#ifdef ASLES
#undef DEFAULT
!        print *,'new fft another fft starts'
         call ldfrmbf(idim,lots,al,1,jump,1,ifax,trigs,work,ierr)
         do jasl=1,lots
           do iasl=1,idim+2
!          do iasl=1,idim
             al(1+(iasl-1)*1+(jasl-1)*jump)=                                   &
             al(1+(iasl-1)*1+(jasl-1)*jump)/dble(idim)
           enddo
         enddo
#endif
#ifdef DEFAULT
         call fft99m (al,work,trigs,ifax,1,jump, idim ,lots,-1)    
#endif
#ifdef ASLES
         if (ierr .ge. 3000)  print 150
         if (ierr .ge. 3000)  stop
#endif
         do j = i,i+lots-1        
            do l=1, mwvp2        
               b(l,j) = al(l,j-i+1) 
            enddo
         enddo
      enddo
   endif                                                                     
!                                                                               
   return                                                                    
   end                                                           
!-------------------------------------------------------------------------------
