#include <define.h>
   subroutine mtn_trans_wave_grid(idir,grid,wave,mlt,fac,                      &
                                  imax,jmax,maxwv,iromb)           
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:  mtn_trans_wave_grid     spherical transform                
!                                                                               
! abstract: transforms a field between grid and spectral domains.               
!   this versatile routine will:                                                
!     ...transform grid to spectral or transform spectral to grid;              
!     ...pass a gaussian grid or an equally-space grid;                         
!     ...pass a triangular truncation or a rhomboidal truncation;               
!     ...optionally transform with derivatives of legendre functions;           
!     ...optionally transform from grid dividing by coslat**2;                  
!     ...optionally multiply spectral field by complex factors.                 
!                                                                               
! program history log:                                                          
!   1991-01-01  iredell                initial mrf
!   2000-01-01  hann-ming henry juang  mpi
!   2000-01-01  song-you hong          usgs topo, kagwd options
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:  call mtn_trans_wave_grid(idir,grid,wave,mlt,fac,imax,jmax,maxwv,iromb) 
!                                                                               
!   input argument list:                                                        
!     idir     - must be one of 1,2,3,4,101,102,103,104,-1,-2,-101,-102,        
!                idir> 0 to transform grid to wave,                            
!                idir < 0 to transform wave to grid,                            
!                abs(idir) < 100 for gaussian grid,                             
!                abs(idir)> 100 for equally-spaced grid,                       
!                idir is odd for normal transform,                              
!                idir is even for legendre derivative transform,                
!                idir last digit> 2 to transform with cosine latitude.         
!     grid     - if idir> 0, real (imax,jmax) field to transform.              
!                grid starts at north pole and greenwich meridian.              
!     wave     - if idir < 0, complex (kmax) field to transform,                
!                where kmax=(maxwv+1)*(iromb+1)*maxwv+2)/2.                     
!                wave starts at the global mean component and then              
!                contains the zonally symmetric components.                     
!     mlt      - multiplication option on wave values.                          
!                mlt = 0 for no multiplication,                                 
!                mlt = 1 to multiply wave by fac                                
!                mlt = -1 to multiply wave by sqrt(-1)*fac                      
!     fac      - if mlt.ne.0, real (kmax) field to multiply wave.               
!     imax     - longitudinal dimension of the grid                             
!     jmax     - latitudinal dimension of the grid                              
!     maxwv    - spectral truncation of the wave                                
!     iromb    - iromb = 0 for triangular truncation                            
!                iromb = 1 for rhomboidal truncation                            
!                                                                               
!   output argument list:                                                       
!     grid     - if idir < 0, real (imax,jmax) field output.                    
!                grid starts at north pole and greenwich meridian.              
!     wave     - if idir> 0, complex (kmax) field output,                      
!                where kmax=(maxwv+1)*(iromb+1)*maxwv+2)/2.                     
!                wave starts at the global mean component and then              
!                contains the zonally symmetric components.                     
!                                                                               
!   subprograms called:                                                         
!     unique:                                                                   
!     mtn_gaussian_lat   - compute gaussian latitudes         
!     mtn_sph_lat        - compute equally-spaced latitudes            
!     (mtn_legendre)     - compute legendre polynomials   
!     fftfax     - fft (library call can be substituted)                        
!     rfftmlt    - fft (library call can be substituted)                        
!                                                                               
!-------------------------------------------------------------------------------
#ifdef SMP_RA2SFC
   use paramodel, only : jcap_=>jcapscm_,latg_=>latgscm_,lonf_=>lonfscm_
#else
   use paramodel, only : jcap_,latg_,lonf_
#endif
!-------------------------------------------------------------------------------
   save
   complex              ::  wave((maxwv+1)*((iromb+1)*maxwv+2)/2)
   real                 ::  grid(imax,jmax)
   real                 ::  fac((maxwv+1)*((iromb+1)*maxwv+2)/2)
!                                                                               
   integer              ::  ngg
   integer              ::  npnm
   integer              ::  ntrigs
   complex, allocatable ::  ww(:)
   complex              ::  w2(-1:1), ws
   real, allocatable    ::  gg(:,:)
!
!---we  real trigs(ntrigs),ifax(100)                                               
!
   real, allocatable    ::  trigs(:)
   integer              ::  ifax(100),indir
!                                                     
   integer              ::  nwrkft
   integer              ::  nepx
   real, allocatable    ::  cosclt(:),wgtclt(:) 
   real, allocatable    ::  wrkfft(:)               
   real, allocatable    ::  pnm  (:)                   
   real, allocatable    ::  ex(:),px(:)         
   integer, allocatable ::  is(:)                        
!   integer              ::  ngg=2*(1+(lonf_+1)/2)
!   integer              ::  npnm=(jcap_+1)*((0+1)*jcap_+2)/2
!   integer              ::  ntrigs=2*lonf_
!   integer              ::  nwrkft=2*lonf_*latg_
!   integer              ::  nepx=jcap_+1
!   complex              ::  ww(npnm)
!   real                 ::  gg(ngg,latg_),trigs(ntrigs)
!   real                 ::  cosclt(latg_),wgtclt(latg_)
!   real                 ::  wrkfft(nwrkft),pnm(npmn)
!   real                 ::  ex(0:nepx),px(-1:nepx)
!   integer              ::  is(npnm)
!
   ngg=2*(1+(lonf_+1)/2)
   npnm=(jcap_+1)*((0+1)*jcap_+2)/2
   ntrigs=2*lonf_
   nwrkft=2*lonf_*latg_     
   nepx=jcap_+1               
   if(.not.allocated(ww)) allocate(ww(npnm))
   if(.not.allocated(gg)) allocate(gg(ngg,latg_))
   if(.not.allocated(trigs)) allocate(trigs(ntrigs))
   if(.not.allocated(cosclt)) allocate(cosclt(latg_))
   if(.not.allocated(wgtclt)) allocate(wgtclt(latg_))
   if(.not.allocated(wrkfft)) allocate(wrkfft(nwrkft) )
   if(.not.allocated(pnm)) allocate(pnm  (npnm) )
   if(.not.allocated(ex)) allocate(ex(0:nepx))
   if(.not.allocated(px)) allocate(px(-1:nepx) )
   if(.not.allocated(is)) allocate(is(npnm)  )
!
   if (idir.gt.0) then
     indir=1
   else
     indir=-1
   endif
!
   kmax=(maxwv+1)*((iromb+1)*maxwv+2)/2                                      
   ipd=1-mod(abs(indir),2)                                                    
   isd=1-2*ipd                                                               
   icd=(mod(abs(indir),10)-1)/2                                               
   jump=2*(1+(imax+1)/2)                                                     
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
   call fftfax (lonf_,ifax,trigs)                                      
#endif
#ifdef DEFAULT
   call    fax (ifax, lonf_,3)                                         
   call fftrig (trigs,lonf_,3)                                         
#endif
   if(abs(idir).lt.100) then                                                 
     call mtn_gaussian_lat(jmax,cosclt,wgtclt)                      
   else                                                                      
#ifdef DFS
     call mtn_dfs_lat(jmax,cosclt,wgtclt)                                       
#else
     call mtn_sph_lat(jmax,cosclt,wgtclt)                                       
#endif
   endif                                                                     
!
   k=0                                                                       
   do m = 0,maxwv                                                              
     do n = m,iromb*m+maxwv                                                    
       k=k+1                                                                 
       is(k)=isd*(1-2*mod(n-m,2))                                            
     enddo                                                                   
   enddo                                                                     
!
   if(indir.gt.0) then                                                        
     do j = 1,jmax                                                             
       do i = 1,imax                                                           
         gg(i,j)=grid(i,j)                                                   
       enddo                                                                 
     enddo                                                                   
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
     call rfftmlt(gg,wrkfft,trigs,ifax,1,jump,imax,jmax,-1)                  
#endif
#ifdef DEFAULT
     call fft99m (gg,wrkfft,trigs,ifax,1,jump,imax,jmax,-1)                  
#endif
     do k = 1,kmax                                                             
       wave(k)=0.                                                            
     enddo                                                                   
     do j = 1,(jmax+1)/2                                                       
       jr=jmax+1-j                                                           
       call mtn_legendre(ipd,cosclt(j),maxwv,iromb,ex,px,pnm)                    
       wj=wgtclt(j)                                                          
       if(icd.ne.0.and.cosclt(j).lt.1.) wj=wj/(1.-cosclt(j)**2)              
       k=0                                                                   
       do m = 0,maxwv                                                          
         w2(1)=wj*cmplx(gg(2*m+1,j)+gg(2*m+1,jr),                              &
                       gg(2*m+2,j)+gg(2*m+2,jr))                            
         w2(-1)=wj*cmplx(gg(2*m+1,j)-gg(2*m+1,jr),                             &
                       gg(2*m+2,j)-gg(2*m+2,jr))                           
         do n = m,iromb*m+maxwv                                                
           k=k+1                                                             
           wave(k)=wave(k)+w2(is(k))*pnm(k)                                  
         enddo                                                               
       enddo                                                                 
     enddo                                                                   
     if(mlt.ne.0) then                                                       
        ws=cmplx(1.,0.)                                                       
        if(mlt.lt.0) ws=cmplx(0.,1.)                                          
        do k = 1,kmax                                                           
          wave(k)=wave(k)*ws*fac(k)                                           
        enddo                                                                 
     endif                                                                   
   else                                                                      
     if(mlt.ne.0) then                                                       
       ws=cmplx(1.,0.)                                                       
       if(mlt.lt.0) ws=cmplx(0.,1.)                                          
       do k = 1,kmax                                                           
         ww(k)=wave(k)*ws*fac(k)                                             
       enddo                                                                 
     else                                                                    
       do k = 1,kmax                                                           
         ww(k)=wave(k)                                                       
       enddo                                                                 
     endif                                                                   
     do j = 1,(jmax+1)/2                                                       
       call mtn_legendre(ipd,cosclt(j),maxwv,iromb,ex,px,pnm)                    
       k=0                                                                   
       do m = 0,maxwv                                                          
         w2(1)=0.                                                            
         w2(-1)=0.                                                           
         do n = m,iromb*m+maxwv                                                
           k=k+1                                                             
           w2(1)=w2(1)+ww(k)*pnm(k)                                          
           w2(-1)=w2(-1)+ww(k)*pnm(k)*is(k)                                  
         enddo                                                               
         gg(2*m+1,jmax+1-j)=real(w2(-1))                                     
         gg(2*m+2,jmax+1-j)=aimag(w2(-1))                                    
         gg(2*m+1,j)=real(w2(1))                                             
         gg(2*m+2,j)=aimag(w2(1))                                            
       enddo                                                                 
     enddo                                                                   
     do i = 2*maxwv+3,jump                                                     
       do j = 1,jmax                                                           
         gg(i,j)=0.                                                          
       enddo                                                                 
     enddo                                                                   
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
     call rfftmlt(gg,wrkfft,trigs,ifax,1,jump,imax,jmax,1)                   
#endif
#ifdef DEFAULT
     call fft99m (gg,wrkfft,trigs,ifax,1,jump,imax,jmax,1)                   
#endif
     do j = 1,jmax                                                             
       do i = 1,imax                                                           
         grid(i,j)=gg(i,j)                                                   
       enddo                                                                 
     enddo                                                                   
   endif                                                                     
!                                                                             
   return                                                                    
   end subroutine mtn_trans_wave_grid
!-------------------------------------------------------------------------------
