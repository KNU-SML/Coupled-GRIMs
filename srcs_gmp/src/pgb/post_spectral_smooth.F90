#include <define.h>
   subroutine post_spectral_smooth(im2,jm2,km,m,ncpus,                         &
                                   clat,slat,wlat,trig,ifax,                   &
                                   u,v,o,z,t,a)                   
!-------------------------------------------------------------------------------
   use paramodel, only : io_,jo_,ko_
!-------------------------------------------------------------------------------
!                                                                               
! subprogram: post_spectral_smooth          perform spectral functions on latlon fields        
!                                                                               
! abstract: this subprogram performs spectral functions on latlon fields        
!   winds, vertical velocities, heights and temperatures are smoothed           
!   and the absolute vorticity is computed from the winds.                      
!   all fields are on a global latlon grid that includes the poles              
!   with northern and southern hemispheres paired.                              
!                                                                               
! program history log:                                                          
!    1992-10-31  iredell                development
!    2000-03-09  songyou hong           cvs verion setup, prognostic clouds
!    2009-10-01  jung-eun kim           f90 format with standard physics modules
!    2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:    call post_spectral_smooth(im2,jm2,km,m,ncpus,clat,slat,wlat,trig,ifax,             
!    &                 u,v,o,z,t,a)                                             
!   input arguments:                                                            
!     im2          integer twice the number of longitude points                 
!     jm2          integer half the number of latitude points                   
!     km           integer number of levels                                     
!     m            integer spectral truncation (jm2-2 recommended)              
!     ncpus        integer number of cpus over which to distribute work         
!     clat         real (jm2) cosines of latitude                               
!     slat         real (jm2) positive sines of latitude                        
!     wlat         real (jm2) gaussian weights                                  
!     trig         real (im2) trigonometric quantities for the fft              
!     ifax         integer (100) factors for the fft                             
!     u            real (im2,jm2,km) zonal wind in m/s                          
!     v            real (im2,jm2,km) meridional wind in m/s                     
!     o            real (im2,jm2,km) vertical velocity in pa/s                  
!     z            real (im2,jm2,km) heights in m                               
!     t            real (im2,jm2,km) temperature in k                           
!   output arguments:                                                           
!     u            real (im2,jm2,km) smoothed zonal wind in m/s                 
!     v            real (im2,jm2,km) smoothed meridional wind in m/s            
!     o            real (im2,jm2,km) smoothed vertical velocity in pa/s         
!     z            real (im2,jm2,km) smoothed heights in m                      
!     t            real (im2,jm2,km) smoothed temperature in k                  
!     a            real (im2,jm2,km) absolute vorticity in 1/s                  
!                                                                               
! subprograms called:                                                           
!   post_spectral_field    compute spectral constants         
!   rfftmlt                compute fft                                
!   post_sph_poly          compute legendre functions               
!   post_sph_fft_coeff     compute spectral from fourier     
!   post_wind2divor        compute divergence and vorticity in spectral space  
!   post_divor2wind        compute wind components in spectral space
!   post_synth_coeff1      compute fourier from spectral                
!                                                                               
!-------------------------------------------------------------------------------
   real, parameter  ::  omega= 7.2921e-5
   real             ::  slat(jm2),clat(jm2),wlat(jm2)
   real             ::  trig(im2)
   integer          ::  ifax(100)
   real             ::  u(im2,jm2,km),v(im2,jm2,km),o(im2,jm2,km)
   real             ::  z(im2,jm2,km),t(im2,jm2,km),a(im2,jm2,km)
!                                                                               
#define DEFAULT
#ifdef DEFAULT
#undef DYNAMIC_ALLOC
#endif
#ifdef DYNAMIC_ALLOC
#undef DEFAULT
   real             ::  eps((m+1)*(m+2)/2),epstop(m+1)
   real             ::  enn1((m+1)*(m+2)/2),elonn1((m+1)*(m+2)/2)
   real             ::  eon((m+1)*(m+2)/2),eontop(m+1)
   real             ::  wfft(im2,2*6*km)
   real             ::  pln((m+1)*(m+2)/2),plntop(m+1)
   real             ::  f1(im2/2+3,2,5*km,ncpus),f2(im2/2+3,2,6*km)
   real             ::  s1((m+1)*(m+2)+1,6*km),s1top(2*(m+1),6*km)
   real             ::  sd((m+1)*(m+2))
   integer          ::  mp(6*km)
#endif
   integer             ::  mfil
#ifdef DEFAULT
   integer, parameter  ::  mcpus=1
   integer             ::  mp(6*ko_)
   real                ::  wfft(2*io_,2*6*ko_)
   real                ::  f1(2*io_/2+3,2,5*ko_,mcpus),f2(2*io_/2+3,2,6*ko_)
   real, allocatable   ::  eps(:),epstop(:)
   real, allocatable   ::  enn1(:),elonn1(:)
   real, allocatable   ::  eon(:),eontop(:)
   real, allocatable   ::  pln(:),plntop(:)
   real, allocatable   ::  s1(:,:),s1top(:,:)
   real, allocatable   ::  sd(:)
#endif
!-------------------------------------------------------------------------------
   mfil=(jo_+1)/2-2
#ifdef DEFAULT
   allocate(eps((mfil+1)*(mfil+2)/2),epstop(mfil+1) )
   allocate(enn1((mfil+1)*(mfil+2)/2),elonn1((mfil+1)*(mfil+2)/2) )
   allocate(eon((mfil+1)*(mfil+2)/2),eontop(mfil+1) )
   allocate(pln((mfil+1)*(mfil+2)/2),plntop(mfil+1) )
   allocate(s1((mfil+1)*(mfil+2)+1,6*ko_),s1top(2*(mfil+1),6*ko_) )
   allocate(sd((mfil+1)*(mfil+2)) )
#endif
!
!  set transform constants                                                      
!                                                                               
   im=im2/2                                                                  
   ix=im+3                                                                   
!                                                                               
   do k = 1,km*2                                                               
     mp(k)=1                                                                 
   enddo                                                                     
!
   do k = 2*km+1,6*km                                                          
     mp(k)=0                                                                 
   enddo                                                                     
!
   if(m.ne.mfil) then                                                        
     write(6,*) 'm.ne.nfil in post_spectral_smooth'                      
     call abort                                                              
   endif                                                                     
!                                                                               
   nc=(m+1)*(m+2)+1                                                          
   nctop=2*(m+1)                                                             
!                                                                               
   call post_spectral_field(m,eps,epstop,enn1,elonn1,eon,eontop)        
!                                                                               
!  transform to spectral space                                                  
!                                                                               
   do k = 1,km                                                                 
     do i = 1,nc                                                               
       s1(i,k)=0.                                                            
       s1(i,k+km)=0.                                                         
       s1(i,k+2*km)=0.                                                       
       s1(i,k+3*km)=0.                                                       
       s1(i,k+4*km)=0.                                                       
     enddo                                                                   
     do i=1,nctop                                                            
       s1top(i,k)=0.                                                         
       s1top(i,k+km)=0.                                                      
       s1top(i,k+2*km)=0.                                                    
       s1top(i,k+3*km)=0.                                                    
       s1top(i,k+4*km)=0.                                                    
       s1top(i,k+5*km)=0.                                                    
     enddo                                                                   
   enddo                                                                     
!                                                                               
!     warning!! wrong in original code 'do j1=2,jm2,ncpus'                      
!                                                                               
   do j1 = 1,jm2,ncpus                                                         
     j2=min(j1+ncpus-1,jm2)                                                  
     do j = j1,j2                                                              
       jc=j-j1+1                                                             
       do k = 1,km                                                             
         do i = 1,im                                                           
           if(clat(j).ne.0.) then                                            
             f1(i,1,k,jc)=u(i,j,k)/clat(j)**2                                
             f1(i,2,k,jc)=u(i+im,j,k)/clat(j)**2                             
             f1(i,1,k+km,jc)=v(i,j,k)/clat(j)**2                             
             f1(i,2,k+km,jc)=v(i+im,j,k)/clat(j)**2                          
           else                                                              
             f1(i,1,k,jc)=0.                                                 
             f1(i,2,k,jc)=0.                                                 
             f1(i,1,k+km,jc)=0.                                              
             f1(i,2,k+km,jc)=0.                                              
           endif                                                             
           f1(i,1,k+2*km,jc)=o(i,j,k)                                        
           f1(i,2,k+2*km,jc)=o(i+im,j,k)                                     
           f1(i,1,k+3*km,jc)=z(i,j,k)                                        
           f1(i,2,k+3*km,jc)=z(i+im,j,k)                                     
           f1(i,1,k+4*km,jc)=t(i,j,k)                                        
           f1(i,2,k+4*km,jc)=t(i+im,j,k)                                     
         enddo                                                               
       enddo                                                                 
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
       call rfftmlt(f1(1,1,1,jc),wfft,trig,ifax,1,ix,im,2*5*km,-1)         
#endif
#ifdef DEFAULT
       call fft99m (f1(1,1,1,jc),wfft,trig,ifax,1,ix,im,2*5*km,-1)           
#endif
     enddo ! j                                                                  
     do j = j1,j2                                                              
       jc=j-j1+1                                                             
       call post_sph_poly(m,slat(j),clat(j),eps,epstop,pln,plntop)  
       wa=wlat(j)                                                            
       call post_sph_fft_coeff(m,im,ix,nc,nctop,5*km,wa,clat(j),               &
                               pln,plntop,mp,                                  &
                               f1(1,1,1,jc),s1,s1top)             
     enddo                                                                   
   enddo ! j1                                                                    
!                                                                               
!  spectrally compute vorticity                                                 
!                                                                               
   do k = 1,km                                                                 
     call post_wind2divor(m,enn1,elonn1,eon,eontop,                            &
                s1(1,k),s1(1,k+km),s1top(1,k),s1top(1,k+km),                   &
                sd,s1(1,k+5*km))                                             
     call post_divor2wind(m,enn1,elonn1,eon,eontop,                            &
                sd,s1(1,k+5*km),                                               &
                s1(1,k),s1(1,k+km),s1top(1,k),s1top(1,k+km))                 
   enddo                                                                     
!                                                                               
!  transform to grid                                                            
!                                                                               
   do j = 1,jm2                                                                
     call post_sph_poly(m,slat(j),clat(j),eps,epstop,pln,plntop)        
     call post_synth_coeff2(m,im,ix,nc,nctop,6*km,clat(j),pln,plntop,mp,       &
                  s1,s1top,f2)                                                
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
     call rfftmlt(f2,wfft,trig,ifax,1,ix,im,2*6*km,1)                      
#endif
#ifdef DEFAULT
     call fft99m (f2,wfft,trig,ifax,1,ix,im,2*6*km,1)                        
#endif
     do k = 1,km                                                               
       do i = 1,im                                                             
         u(i,j,k)=f2(i,1,k)                                                  
         u(i+im,j,k)=f2(i,2,k)                                               
         v(i,j,k)=f2(i,1,k+km)                                               
         v(i+im,j,k)=f2(i,2,k+km)                                            
         o(i,j,k)=f2(i,1,k+2*km)                                             
         o(i+im,j,k)=f2(i,2,k+2*km)                                          
         z(i,j,k)=f2(i,1,k+3*km)                                             
         z(i+im,j,k)=f2(i,2,k+3*km)                                          
         t(i,j,k)=f2(i,1,k+4*km)                                             
         t(i+im,j,k)=f2(i,2,k+4*km)                                          
         a(i,j,k)=f2(i,1,k+5*km)+2*omega*slat(j)                             
         a(i+im,j,k)=f2(i,2,k+5*km)-2*omega*slat(j)                          
       enddo                                                                 
     enddo                                                                   
   enddo                                                                     
!
#ifdef DEFAULT
   deallocate(eps,epstop)
   deallocate(enn1,elonn1 )
   deallocate(eon,eontop)
   deallocate(pln,plntop)
   deallocate(s1,s1top)
   deallocate(sd)
#endif
   return                                                                    
   end                                                                       
