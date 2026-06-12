#include <define.h>
   subroutine sph_trans_wave2grid
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!    [sph_trans_wave2grid]
!        |
!        |--- [sph_wave2grid] *
!        |--- [sph_grid2wave] *
!        |--- [difference_wave] *
!        |--- [difference_grid] *
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   end subroutine sph_trans_wave2grid
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_wave2grid (a, b, lot, lcapf, lonff, lat, isign)
!-------------------------------------------------------------------------------
#ifndef DFS
   use paramodel, only : latg_,lonf_,lonf22_
#ifdef MP
#ifdef RMP
   use paramodel, only : igrd1p_
#else
   use paramodel, only : lonfp_
#endif
#else /* ~MP */
#ifdef RMP
   use paramodel, only : igrd_,igrd1_
#endif
#endif /* MP end */
#include "abort.h"
!
#define DEFAULT
#ifdef DCRFT
#undef DEFAULT
#ifndef RMP
   use comfcst, only : crscale, rcscale
#endif
#endif
#ifdef FFT_SGIMATH
#undef DEFAULT
   use comfcst, only  : scale,trig
#endif
#ifdef FFTW
#undef DEFAULT
   use comfcst, only : scale, iplan_c_to_r,iplan_r_to_c
#endif
#ifdef DEFAULT
   use comfcst, only : trigs, ifax
#endif
!-------------------------------------------------------------------------------
   integer, parameter   ::  ncpu=1
   integer              ::  lotmin,lotmax
!                                                                               
#ifdef DCRFT
   integer              ::  jumph,jump
   real                 ::  table(42000)
#endif
#ifdef FFT_SGIMATH
   integer              ::  jump
#endif
#ifdef FFTW
   integer              ::  jump
   real, allocatable    ::  ak(:,:)
#endif
#ifdef DEFAULT
   integer              ::  jump
   real, allocatable    ::  work(:,:,:)
#endif
   real, allocatable    ::  al(:,:)
   real                 ::  a(lonf22_,lot)                    
   real                 ::  b(lonf22_,lot)                   
!-------------------------------------------------------------------------------
#ifdef REDUCE_GRID
#define LATX lat
#else
#define LATX 1
#endif
!
#ifdef MP
#ifdef RMP
#define NVECT (igrd1p_)*2
#else
#define NVECT lonfp_*2
#endif
#else
#ifdef RMP
#define NVECT (igrd_+1)*2
#else
#define NVECT 64
#endif
#endif  /* MP */
   lotmin=NVECT
   lotmax=NVECT
#ifdef DCRFT
   jumph=(lonf_+2)/2
   jump=jumph*2
#endif
#ifdef FFT_SGIMATH
   jump=lonf_ + 3
#endif
#ifdef FFTW
   jump=lonf_ + 3
   allocate(ak(jump,lotmax))
#endif
#ifdef DEFAULT
   jump=lonf_ + 3
   allocate(work(lonf_,lotmax,2))
#endif
   allocate(al(jump,lotmax))
!                                                        
   jcapf=lcapf-1
   twoj1=2*(jcapf+1)
!
   if (isign .eq.  1)  then                                                  
!                                                                               
!  multiple fast fourier transform - synthesis.  isign=1                        
!  good for zonal wave number jcap_.                                            
!                                                                               
!     dimension a(lonf_,lot)                                                    
!                                                                               
!   input - lot sets of complex coefficients in                                 
!           a(1,j), a(2,j), ..., a(jcap1_*2,j), j=1,...,lot.                    
!           a(jcap1_*2+1,j), ..., a(lonf_,j), j=1,...,lot are not set           
!           before call fftlonf.                                                
!                                                                               
!  output - lot sets of grid values in                                          
!           a(1,j), a(2,j), ..., a(lonf_,j), j=1,...,lot.                       
!
     nlot=max0(lot/ncpu,lotmin)                                                
     nlot=min0(nlot    ,lotmax)                                                
!
     do i = 1,lot,nlot                                                       
       lots = min0(nlot, lot-i+1)                        
       do j = i,i+lots-1                                                       
         jj=j-i+1
         k=(j+1)/2                                     
         loff=mod(j+1,2)*lonf_                        
#ifdef FFTW
         al(1,jj)=a(1+loff,k)
         do l = 2,jcapf+1
           al(l        ,jj)=a(2*l-1+loff,k)
           al(lonff+2-l,jj)=a(2*l  +loff,k)
         enddo
         do l = jcapf+2,lonff-jcapf
           al(l,jj)=0.0
         enddo
         do l = lonff+1,jump
           al(l,jj)=0.0
         enddo
#else
         do l = 1,twoj1                                                         
           al(l,jj) = a(l+loff,k)                                          
         enddo
         do l = twoj1+1,jump
           al(l,jj) = 0.0                                                  
         enddo
#endif
       enddo ! j
!                                              
!     call  fft for systhesis.                
!                                            
#define DEFAULT
#ifdef DCRFT
#undef DEFAULT
#ifdef RMP
       scale=1.0
       call dcrft(1,al,jumph,al,jump,lonf_,lots,                               &
            -1,scale,table,22000,table(22001),20000)
       call dcrft(0,al,jumph,al,jump,lonf_,lots,                               &
            -1,scale,table,22000,table(22001),20000)
#else
       call dcrft(1,al,jumph,al,jump,lonff,lots,                               &
            -1,crscale(LATX),table,22000,table(22001),20000)
       call dcrft(0,al,jumph,al,jump,lonff,lots,                               &
            -1,crscale(LATX),table,22000,table(22001),20000)
#endif
#endif
#ifdef FFT_SGIMATH
#undef DEFAULT
       call zdfftm1du(-1,lonff,lots,al,1,jump,trig(1,LATX))  
#endif
#ifdef FFTW
#undef DEFAULT
       call rfftw_f77(iplan_c_to_r(LATX),lots,al,1,jump,ak,1,jump)
#endif
#ifdef RFFTMLT
#undef DEFAULT
       call rfftmlt(al,work,trigs(1,1,LATX),ifax(1,LATX),1,jump,lonff,lots,1)
#endif
#ifdef ASLES
#undef DEFAULT
!        print *,'new fft, fft starts'
       call ldfrmbf(lonff,lots,al,1,jump,-1,ifax(1,LATX),                      &
            trigs(1,1,LATX),work,ierr)
#endif
#ifdef DEFAULT
       call fft99m (al,work,trigs(1,1,LATX),ifax(1,LATX),                      &
                                   1,jump,lonff,lots,1)
#endif
#ifdef ASLES
       if (ierr .ge. 3000)  write(6,*) ' error in asles FFT '
       if (ierr .ge. 3000)  call MPABORT
#endif
!
#ifdef FFTW
#define ALK ak
#else
#define ALK al
#endif
       do j = i,i+lots-1                                                       
         k=(j+1)/2    
         loff=mod(j+1,2)*lonff     
         do l = 1,lonff             
           b(l+loff,k) = ALK(l,j-i+1)
         enddo
       enddo
!
     enddo ! i
!
   endif ! isign.eq.1                                                                    
!                                                                               
   if (isign .eq. -1)  then                                                  
!                                                                               
!  multiple fast fourier transform - analysis.  isign=-1                        
!  good for zonal wave number jcap_.                                            
!                                                                               
!     dimension a(lonf_,lot), b(lonf_,lot)                                      
!                                                                               
!   input - lot sets of grid values in                                          
!           a(1,j), a(2,j), ..., a(lonf_,j), j=1,...,lot.                       
!           a array is not changed by sr fftlonf.                               
!                                                                               
!  output - lot sets of complex coefficients in                                 
!           b(1,j), b(2,j), ..., b(jcap1_*2,j), j=1,...,lot.                    
!           b(jcap1_*2+1,j), ..., b(lonf_,j), j=1,...,lot are not set.          
!                                                                               
     nlot=max0(lot/ncpu,lotmin)                                                
     nlot=min0(nlot    ,lotmax)                                                
     do i = 1,lot,nlot                                                       
       lots = min0(nlot, lot-i+1)                                                
       do j = i,i+lots-1                                                       
         k=(j+1)/2                   
         loff=mod(j+1,2)*lonff      
         do l = 1,lonff              
           al(l,j-i+1) = a(l+loff,k) 
         enddo
       enddo
!                                                                               
!     call fft for analysis.                                               
!
#define DEFAULT
#ifdef DCRFT
#undef DEFAULT
#ifdef RMP
       scale=1./float(lonf_)
       call drcft(1,al,jump,al,jumph,lonf_,lots,                               &
              1,scale,table,22000,table(22001),20000)
       call drcft(0,al,jump,al,jumph,lonf_,lots,                               &
              1,scale,table,22000,table(22001),20000)
#else
       call drcft(1,al,jump,al,jumph,lonff,lots,                               &
              1,rcscale(LATX),table,22000,table(22001),20000)
       call drcft(0,al,jump,al,jumph,lonff,lots,                               &
              1,rcscale(LATX),table,22000,table(22001),20000)
#endif
#endif
#ifdef FFT_SGIMATH
#undef DEFAULT
       call dzfftm1du(1,lonff,lots,al,1,jump,trig(1,LATX))                    
       call dscalm1d(lonff,lots,scale(LATX),al,1,jump)
#endif
#ifdef FFTW
#undef DEFAULT
       call rfftw_f77(iplan_r_to_c(LATX),lots,al,1,jump,ak,1,jump)
#endif
#ifdef RFFTMLT
#undef DEFAULT
       call rfftmlt(al,work,trigs(1,1,LATX),ifax(1,LATX),                      &
                                   1,jump,lonff,lots,-1)                   
#endif
#ifdef ASLES
#undef DEFAULT
!        print *,'new fft, fft starts'
       call ldfrmbf(lonff,lots,al,1,jump,1,ifax(1,LATX),                       &
                               trigs(1,1,LATX),work,ierr)
!
       do jasl = 1,lots
         do iasl = 1,lonff+2
           al(1+(iasl-1)*1+(jasl-1)*jump)=                                     &
              al(1+(iasl-1)*1+(jasl-1)*jump)/dble(lonff)
         enddo
       enddo
#endif
#ifdef DEFAULT
       call fft99m (al,work,trigs(1,1,LATX),ifax(1,LATX),                      &
                                   1,jump,lonff,lots,-1)                     
#endif
#ifdef ASLES
       if (ierr .ge. 3000)  write(6,*)' error in asles FFT '
       if (ierr .ge. 3000)  call MPABORT
#endif
!                                                                               
       do j = i,i+lots-1                                                       
         k=(j+1)/2     
         loff=mod(j+1,2)*lonf_      
#ifdef FFTW
         jj=j-i+1
         b(1+loff,k)      =ak(1,jj)*scale
         b(2+loff,k)      =0.0
         do l = 2,jcapf+1
           b(2*l-1+loff,k)=ak(l        ,jj)*scale
           b(2*l  +loff,k)=ak(lonff+2-l,jj)*scale
         enddo
#else
         do l = 1,twoj1              
           b(l+loff,k) = al(l,j-i+1)    
         enddo
#endif
       enddo
!
     enddo ! i
!
   endif ! isign.eq.-1                                                                    
!                                                                               
#ifdef FFTW
   deallocate(ak)
#endif
#ifdef DEFAULT
   deallocate(work)
#endif
   deallocate(al)
#endif  /* not DFS */
!
   return                                                                    
   end subroutine sph_wave2grid
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine sph_grid2wave(flp,flm,anl,lwvlen,lota)
!-------------------------------------------------------------------------------
#ifndef RMP
   use paramodel, only : JCAP1S,LCAPS,LCAP22S
!-------------------------------------------------------------------------------
!
!  for loopa                                                                    
!                                                                               
   integer              ::  lwvlen,lota
   real                 ::  flp(2,JCAP1S,lota),flm(2,JCAP1S,lota)
   real                 ::  anl(LCAP22S,lota)
!-------------------------------------------------------------------------------
!
   do k = 1,lota                                                            
     do  ll = 1,lwvlen
       !
       ! do n.hemi                    
       !
       flp(1,ll,k)=anl(2*(ll-1)+1,k)+anl(2*(ll-1)+1+LCAPS,k)
       flp(2,ll,k)=anl(2*(ll-1)+2,k)+anl(2*(ll-1)+2+LCAPS,k)
       !                             
       ! do s.hemi                  
       !                           
       flm(1,ll,k)=anl(2*(ll-1)+1,k)-anl(2*(ll-1)+1+LCAPS,k) 
       flm(2,ll,k)=anl(2*(ll-1)+2,k)-anl(2*(ll-1)+2+LCAPS,k)
     enddo                                                                     
   enddo                                                                     
#endif
!
   return                                                                    
   end subroutine sph_grid2wave
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine difference_wave(iunit,gdata,mt,jl,kmx,str)
!-------------------------------------------------------------------------------
#ifdef DFS
   use dfsvar, only : mtg,jlg,iope
!-------------------------------------------------------------------------------
   integer              ::  iunit,mt,jl,kmx
   character(len=10)    ::  str
   real                 ::  gdata (mt,jl,kmx)
   real                 ::  fdata (mtg,jlg,kmx)
#ifdef MP
   real                 ::  full(mtg,jlg,kmx)
#endif
!-------------------------------------------------------------------------------
#ifdef MP
   call mpsp2f(gdata,mt,jl,full,mtg,jlg,kmx)
#define GDATA full
#endif
!
   if (str(1:4)=='SAVE') then
     if (iope) then
       write(iunit)GDATA(:,:,:)
       call flush(iunit)
     endif
     return
   endif
!
   if (iope) then
!
     read(iunit)fdata(:,:,:)
     dmax=-99999999999.
     dmin= 99999999999.
     idif=0
     icnt=0
     do k = 1,kmx
       do j = 1,jlg
         do i = 1,mtg
           if ( GDATA(i,j,k) /= fdata(i,j,k) ) then
             diff=abs(gdata(i,j,k)-fdata(i,j,k))
             icnt=icnt+1
             if (diff>dmax) dmax=diff
             if (diff<dmin) dmin=diff
             if (idif==0)                                                      &
               write(*,*)str,' data diff:gdata,fdata,i,j,k=',                  &
               GDATA(i,j,k),fdata(i,j,k),i,j,k
             idif=1
           endif
         enddo
       enddo
     enddo
!
     if (idif==0) then
       write(*,*)str,' is same'
     else
       write(*,*)str,' total ',icnt,' data diff, max,min=',dmax,dmin
     endif
!
     if (str(1:4)=='COPY') then
       do k = 1,kmx
         do j = 1,jlg
           GDATA(:,j,k)=fdata(:,j,k)
         enddo
       enddo
       write(*,*)'data was COPIED'
     else
       write(*,*)str
     endif
!
   endif ! iope
!
#ifdef MP
#undef GDATA
   if (str(1:4)=='COPY') then
     call mpsf2p(full,mtg,jlg,gdata,mt,jl,kmx)
   endif
#endif
#endif
   return
   end subroutine difference_wave
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine difference_grid(iunit,gdata,imx,jmx,kmx,str)
!-------------------------------------------------------------------------------
#ifdef DFS
   use dfsvar, only : iba,jbwa,iope
!-------------------------------------------------------------------------------
   integer              ::  iunit,imx,jmx,kmx
   character(len=10)    ::  str
   real                 ::  gdata (imx,jmx,kmx)
   real                 ::  fdata (iba,jbwa,kmx)
#ifdef MP
   real                 ::  full(iba,jbwa,kmx)
#endif
!-------------------------------------------------------------------------------
!
#ifdef MP
   call mpgp2f(gdata,imx,jmx,full,iba,jbwa,kmx)
#define GDATA full
#endif
!
   if (str(1:4)=='SAVE') then
     if (iope) then
       write(iunit)GDATA(:,:,:)
       call flush(iunit)
     endif
     return
   endif
!
   if (iope) then
     read(iunit)fdata(:,:,:)
     dmax=-99999999999.
     dmin= 99999999999.
     idif=0
     icnt=0
     do k = 1,kmx
       do j = 1,jbwa
         do i = 1,iba
           if ( GDATA(i,j,k) /= fdata(i,j,k) ) then
             diff=abs(gdata(i,j,k)-fdata(i,j,k))
             icnt=icnt+1
             if (diff>dmax) dmax=diff
             if (diff<dmin) dmin=diff
             if (idif==0)                                                      &
                write(*,*)str,' data diff:gdata,fdata,i,j,k=',                 &
                GDATA(i,j,k),fdata(i,j,k),i,j,k
             idif=1
           endif
         enddo
       enddo
     enddo
!
     if (idif==0) then
       write(*,*)str,' is same'
     else 
       write(*,*)str,' total ',icnt,' data diff, max,min=',dmax,dmin
     endif
!
     if (str(1:4)=='COPY') then
       do k = 1,kmx
         do j = 1,jbwa
           GDATA(:,j,k)=fdata(:,j,k)
         enddo
       enddo
       write(*,*)'data was COPIED'
     else
       write(*,*)str
     endif
!
   endif ! iope
!
#ifdef MP
#undef GDATA
   if (str(1:4)=='COPY') then
     call mpgf2p(full,iba,jbwa,gdata,imx,jmx,kmx)
   endif
#endif
#endif
!
   return
   end subroutine difference_grid
!-------------------------------------------------------------------------------
