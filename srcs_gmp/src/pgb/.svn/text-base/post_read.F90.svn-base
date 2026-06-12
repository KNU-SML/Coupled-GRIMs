#include <define.h>
   subroutine post_read
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [post_read]
!           |
!           |-- [post_read_sigma] *
!           |-- [post_read_coeff] *
!                    |-- [post_kinetic_energy] * 
!                             |-- [post_wind_speed] * 
!
!-------------------------------------------------------------------------------
   end subroutine post_read
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_read_sigma(nsig,fh,fhour,idate,si,sl,iret,igases,iwater)
!-------------------------------------------------------------------------------
   use paramodel, only : levp1_,levs_
!-------------------------------------------------------------------------------
!                                                                               
! subprogram:    post_read_sigma       read sigma file header record                      
!                                                                               
! abstract: reads the header record from the sigma file.                        
!                                                                               
! program history log:                                                          
!    1992-10-31  iredell                development
!    2000-03-09  songyou hong           cvs verion setup, prognostic water
!    2009-10-01  jung-eun kim           f90 format with standard physics modules
!    2010-07-01  myung-seo koo          dimension allocatable with namelist input
!                                                                               
! usage:      call post_read_sigma(nsig,fhour,idate,si,sl,iret)       
!                                                                               
!   input argument list:                                                        
!     nsig     - integer unit from which to read header                         
!                                                                               
!   output argument list:                                                       
!     fhour    - real forecast hour                                             
!     idate    - integer (4) date                                               
!     si       - real (levs+1) sigma interfaces                                 
!     sl       - real (levs) sigma levels                                       
!                                                                               
!   input files:                                                                
!     nsig     - sigma file                                                     
!                                                                               
! subprograms called:                                                           
!   maxfac       return maximum prime factor                                    
!                                                                               
!-------------------------------------------------------------------------------
   character(LEN=32)            ::  clabe                                                        
   integer                      ::  idate(4)                                                        
   real                         ::  si(levs_+1)
   real                         ::  sl(levs_)
   real                         ::  dummy(201-levp1_-levs_)
   real                         ::  ensemble(2),dummy2(21)
   character(LEN=3), parameter  ::  fni='sig'
   integer, parameter           ::  nchi=3
!
   character(LEN=80)            ::  fno
#ifdef ASSIGN
   character(LEN=80)            ::  asgnstr
#endif
!-------------------------------------------------------------------------------
   call file_name(fni,nchi,fh,fno,ncho)
#ifdef ASSIGN
   write(asgnstr,'(9hassign u:,I2,)') nsig
   call assign('assign -R')
   call assign(asgnstr)
#endif
   open(unit=nsig,file=fno(1:ncho),form='unformatted',err=900)
   go to 901
900 continue
   write(6,*) ' error in opening file ',fno(1:ncho)
   call abort
901 continue
   write(6,*) ' file ',fno(1:ncho),' opened. unit=',nsig
!
!  read and extract header record                                               
!  read sigma spectral file header and determine gaussian grid                  
!
!     print *,'reading lab'                                                     
!     call flush(6)                                                             
   read(nsig,end=91,err=92) clabe                                            
!     print *,'reading fhour,....'                                              
!     call flush(6)                                                             
   read(nsig,err=201) fhour,idate,si,sl                                        &
                      ,dummy,waves,xlayers,trun,order,realform,gencode         &
                      ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,water               &
                      ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid    &
                      ,pdryini,dummy2,gases
   iwater=nint(water)
   igases=nint(gases)
   if(fh.ne.fhour) then
     print *,'fh and fhour does not match'
     call abort
   endif
   goto 202
201   continue
   rewind nsig
   read(nsig,end=91,err=92) clabe                                            
   read(nsig,err=201)fhour,idate,si,sl
   iwater=1
   igases=0
202   continue
#ifndef NOPRINT
   print *,'tread unit,fhour,idate=',nsig,fhour,idate
   print *,' number of water input = ',iwater
   print *,' number of gases input = ',igases
#endif
   iret=0
   return                                                                    
!
!  end of file encountered                                                      
!
91 iret=1                                                                    
   return                                                                    
!
!  i/o error encountered                                                        
!
92 iret=2                                                                    
!
   return                                                                    
   end subroutine post_read_sigma
!-------------------------------------------------------------------------------
!
!#define BIN_DBG
#ifdef DFS
!-------------------------------------------------------------------------------
   subroutine post_read_coeff(nss,keid,nflds,igases,iwater,nc,jcap,            &
                              fhour,si,fxys,io,jo,io2,io22,johf)
!-------------------------------------------------------------------------------
   use dfsvar, only : mta,jla
#else
!-------------------------------------------------------------------------------
   subroutine post_read_coeff(nss,keid,nflds,igases,iwater,nc,nctop,jcap,      &
                              fhour,si,clat,slat,wlat,trig,ifax,eps,epstop,    &
                              ss,sstop,io,jo,io2,io22,johf,gauss_lat)
!-------------------------------------------------------------------------------
#endif
   use paramodel, only : nwater_,ntotal_,ngases_,levs_,jcap_,lonf_,latg_,      &
                         levs => levs_
   use constant, only  : rd  => rd_ , g    => g_   , rv   => rv_
   use module_file_write, only : file_write_bin
   use module_parm_ptr
!-------------------------------------------------------------------------------
!
! subprogram:    post_read_coeff        read data from a sigma spectral file
!
! abstract: reads the records of orography, surface pressure,
!           divergence and vorticity, temperature and humidity
!           from a sigma spectral file.  it is assumed that the first
!           two header records of the file have already been read.
!           the gradients of orography and log surface pressure
!           and the wind components are also computed in spectral space.
!           the geopotential of the pressure gradient is computed too.
!           also, some spectral transform utility fields are computed.
!           subprogram post_wave2grid should be used to transform to grid
!           as well as compute dry temperature and surface pressure
!           and winds and gradients without a cosine latitude factor.
!
! program history log:
!   91-10-31  mark iredell
!
! usage:    call post_read_coeff(nss,jcap,nc,nctop,johf,io2,levs,sl,
!   &                 clat,slat,trig,ifax,eps,epstop,ss,sstop)
!
!   input argument list:
!     nss      - integer unit from which to read file
!     jcap     - integer spectral truncation
!     nc       - integer number of spectral coefficients
!     nctop    - integer number of spectral coefficients over top
!     johf    - integer number of latitude pairs in gaussian grid
!     io2    - integer number of valid data points per latitude pair
!     levs     - integer number of levels
!     sl       - real (levs) sigma full level values
!
!   output argument list:
!     clat     - real (johf) cosines of latitude
!     slat     - real (johf) sines of latitude
!     trig     - real (io2) trigonometric quantities for the fft
!     ifax     - integer (100) factors for the fft
!     eps      - real ((jcap+1)*(jcap+2)/2) sqrt((n**2-l**2)/(4*n**2-1))
!     epstop   - real (jcap+1) sqrt((n**2-l**2)/(4*n**2-1)) over top
!     ss       - real (nc,6*levs+6) spectral coefs
!     sstop    - real (nctop,6*levs+6) spectral coefs over top
!                (:,1)                             orography
!                (:,2)                             d(orog)/dx
!                (:,3)                             d(orog)/dy
!                (:,4)                             surface pressure
!                (:,5)                             d(lnps)/dx
!                (:,6)                             d(lnps)/dy
!                (:,7:levs+6)                      divergence
!                (:,levs+7:2*levs+6)               vorticity
!                (:,2*levs+7:3*levs+6)             zonal wind
!                (:,3*levs+7:4*levs+6)             meridional wind
!                (:,4*levs+7:5*levs+6)             temperature
!                (:,5*levs+7:5*levs+iwater*levs+6)  water
!                (:,5*levs+7:5*levs+iwater*levs+6+igases*levs)  gases
!
!   input files:
!     nss      - sigma spectral file
!
! subprograms called:
!   post_equall_lat         compute latitudes
!   fftfax                  compute utility fields for fft
!   post_spectral_field     compute utility fields for spectral transform
!   post_grad_scalar        compute gradient in spectral space
!   post_divor2wind         compute vector components in spectral space
!
!-------------------------------------------------------------------------------
   integer                   ::  io,jo,io2,io22,johf
   integer                   ::  nwav
#ifdef DFS
   real,parameter            ::  fv=rv/rd -1.
   integer                   ::  jlha,jlg,mt
   integer                   ::  keid
   real                      ::  fhour
   real                      ::  fxys(io,jo,nflds)
   real, allocatable         ::  dumm(:,:,:)
   real, allocatable,save    ::  AMATm (:,:,:,:),DMATm (:,:,:,:)
   real, allocatable,save    ::  DMATm_(:,:,:,:),AMATm_(:,:,:,:)
   real, allocatable,save    ::  aMSQUAR(:)
   real, allocatable,save    ::  clatd(:,:)
   real                      ::  slatd(latg_)
   real, allocatable         ::  ss(:,:,:)
   real                      ::  tave(levs)
   integer, allocatable,save ::  jcol2js(:,:)
   data ist/1/
!keid
   real, dimension((mta+1)*(mta+2),levs)  ::  div,vor
#else  /* SPH */
   real                      ::  clat(johf),slat(johf)
   real                      ::  eps((jcap+1)*(jcap+2)/2),epstop(jcap+1)
   real                      ::  ss(nc,nflds)
   real                      ::  sstop(nctop,nflds)
!
   real                      ::  enn1((jcap+1)*(jcap+2))
   real                      ::  elonn1((jcap+1)*(jcap+2)/2)
   real                      ::  eon((jcap+1)*(jcap+2)/2),eontop(jcap+1)
   real                      ::  pln((jcap+1)*(jcap+2)/2),plntop(jcap+1)
   real                      ::  plndx((jcap+1)*(jcap+2)/2)
   real                      ::  plndy((jcap+1)*(jcap+2)/2)
   real                      ::  f(io2/2+3,2,3),wfft(io2,2*3)
   real                      ::  cosclt(jo),wgt(jo)
   logical                   ::  gauss_lat
#endif
   real                      ::  ken(0:jcap,levs+1)     ! kinetic energy
   real                      ::  dke(0:jcap,levs+1)     ! divergent
   real                      ::  vke(0:jcap,levs+1)     ! vorticity
   character*100 fmt
   real                      ::  si(levs+1),del(levs)
   real                      ::  wlat(johf)
   real                      ::  trig(io2)
   integer                   ::  ifax(100)
   character*200 ename
!-------------------------------------------------------------------------------
#ifdef DFS
   jlha=jla/2
   jlg=jla+1
   mt=2*mta+1
   nwav=mt*jlg
!
   if(.not.allocated(AMATm))   allocate( AMATm (-mta:mta,3,0:JLHA,2) )
   if(.not.allocated(DMATm))   allocate( DMATm (-mta:mta,3,0:JLHA,2) )
   if(.not.allocated(DMATm_))  allocate( DMATm_(-mta:mta,3,0:JLHA,2) )
   if(.not.allocated(AMATm_))  allocate( AMATm_(-mta:mta,3,0:JLHA,2) )
   if(.not.allocated(aMSQUAR)) allocate( aMSQUAR(-mta:mta) )
   if(.not.allocated(jcol2js)) allocate( jcol2js(0:JLHa,2) )
   if(.not.allocated(clatd))   allocate( clatd(latg_,3) )
   allocate( ss(nwav,levs,max(ntotal_,2)) )
   allocate( dumm(lonf_,latg_,levs*max(ntotal_,4)) )
#else
   nwav=nc
#endif
!
   if (keid .ne.0 ) then
     do k = 1,levs
       del(k)=si(k)-si(k+1)
     enddo
     if (levs.le.98) then
       write(fmt,'(A4,I2,A6)')'(I5,',levs+1,'F13.7)'
     else
       write(fmt,'(A4,I3,A6)')'(I5,',levs+1,'F13.7)'
     endif
   endif
!
   ktotal = igases*levs+iwater*levs
!
#ifdef DFS
   if (ist==1) then
     ist=0
     call post_equall_lat(latg_,slatd,clatd,1)
     call dfs_semi_matrix(AMATm,DMATm,aMSQUAR,DMATm_,AMATm_,jcol2js,mt,jlha)
   endif
!
   read(nss) ss(1:nwav,1,1)                  ! topo
   read(nss) ss(1:nwav,1,2),tave(1:levs)     ! ps,tavexy
!
   call print_maxmin_seven(ss(1,1,1),nwav,nwav,1,1,1,'topo wave')
   call print_maxmin_seven(ss(1,1,2),nwav,nwav,1,1,1,'lnps wave')
   call post_dfs_sfcp(ss(1,1,1),ss(1,1,2),mta,jla,                             &
                      dumm(1,1,1),dumm(1,1,2),dumm(1,1,3),                     &
                      dumm(1,1,4),dumm(1,1,5),dumm(1,1,6),lonf_,latg_,clatd)
   ! kszs,kszsx,kszsy,ksps,kspsx,kspsy (6)
   call post_dfs_latlon(dumm(1,1,1),lonf_,latg_,fxys(1,1,kszs),io,jo,6)
#ifdef BIN_DBG
   call file_write_bin(119,fxys(1,1,kszs),io,jo,1,0)
   call file_write_bin(119,fxys(1,1,ksps),io,jo,1,0)
#endif
!
! Virtual Temperature
!
   do k = 1,levs
     read(nss) ss(1:nwav,k,1)
   enddo
!
   call dfs_fft_driver(-1,dumm,lonf_,latg_,levs,                               &
             ss,mt,jlg,levs,jlg,levs,levs,1,clatd,+1 )
   call post_dfs_latlon(dumm,lonf_,latg_,fxys(1,1,kst),io,jo,levs)
!
   do k = 0,levs-1
     fxys(:,:,kst+k)=fxys(:,:,kst+k)+tave(k+1)
   enddo
!
#ifdef BIN_DBG
   call file_write_bin(119,fxys(1,1,kst),io,jo,levs,0)
!
#endif
   do k = 0,levs-1
     read(nss) ss(1:nwav,k+1,1)               ! Divergence
     read(nss) ss(1:nwav,k+1,2)               ! Vorticity
   enddo
!
   call dfs_divor2wind_kinetic( ss(1,1,1),ss(1,1,2),mt,jlg,levs,               &
                        AMATm,DMATm,AMATm_,DMATm_,aMSQUAR,jcol2js,             &
                        dumm(1,1,1),dumm(1,1,levs+1),                          &
                        dumm(1,1,2*levs+1),dumm(1,1,3*levs+1),                 &
                        lonf_,latg_,clatd)
   ! ksd,ksz,ksu,ksv
   call post_dfs_latlon(dumm(1,1,1),lonf_,latg_,fxys(1,1,ksd),io,jo,4*levs)
!
   if (keid.ne.0) then
     ! Hoon Park
!     call post_kinetic_energy(dumm(1,1,2*levs+1),dumm(1,1,3*levs+1),          &
!                               lonf_,latg_,levs,                               &
!                               ken,mta,del)
     ! grid->wave &  Koshyk (2001)
     do k = 1,levs
       call post_trans_wave_grid(101,dumm(1,1,k     ),div(1,k),0,0.,lonf_,latg_,mta,0)
       call post_trans_wave_grid(101,dumm(1,1,levs+k),vor(1,k),0,0.,lonf_,latg_,mta,0)
     enddo
     call post_kinetic_energy2(div,vor,dke,vke,ken,(mta+1)*(mta+2),mta,levs)
!
     call file_name('total_ke',8,fhour,ename,nchout)
     open(keid,file=ename(1:nchout),form='formatted',status='unknown')
     call file_name('diver_ke',8,fhour,ename,nchout)
     open(keid+1,file=ename(1:nchout),form='formatted',status='unknown')
     call file_name('vorti_ke',8,fhour,ename,nchout)
     open(keid+2,file=ename(1:nchout),form='formatted',status='unknown')
     do i = 0,mta
       write(keid  ,trim(fmt))i,(ken(i,k),k=1,levs+1)
       write(keid+1,trim(fmt))i,(dke(i,k),k=1,levs+1)
       write(keid+2,trim(fmt))i,(vke(i,k),k=1,levs+1)
     enddo
     close(keid)
     close(keid+1)
     close(keid+2)
   endif
#ifdef BIN_DBG
   call file_write_bin(119,fxys(1,1,ksd),io,jo,levs,0)
   call file_write_bin(119,fxys(1,1,ksz),io,jo,levs,0)
   call file_write_bin(119,fxys(1,1,ksu),io,jo,levs,0)
   call file_write_bin(119,fxys(1,1,ksv),io,jo,levs,0)
#endif
   lot=igases*levs+iwater*levs
   nvar=lot/levs
!
   do n = 1,nvar
     do k = 1,levs
       read(nss) ss(1:nwav,k,n)
     enddo
   enddo
!
   call dfs_fft_driver(-1,dumm,lonf_,latg_,lot,                                &
                      ss,mt,jlg,lot,jlg,levs,levs,nvar,clatd,1)
   call post_dfs_latlon(dumm,lonf_,latg_,fxys(1,1,ksq),io,jo,lot)
!
! virtual temp to temp
!
   do k = 0,levs-1
     do j = 1,jo
       do i = 1,io
         fxys(i,j,kst+k)= fxys(i,j,kst+k)/(1.+fv*fxys(i,j,ksq+k))
       enddo
     enddo
   enddo
!
#ifdef BIN_DBG
   call file_write_bin(119,fxys(1,1,kst),io,jo,levs,0)
   call file_write_bin(119,fxys(1,1,ksq),io,jo,levs,0)
#endif
   call print_maxmin_seven(fxys(1,1,ksq),io*jo,io*jo,levs,1,levs,'sph')
   call print_maxmin_seven(fxys(1,1,kst),io*jo,io*jo,levs,1,levs,'temp')
!
#else                   /* not DFS */
!
#ifdef SMP
   do nf = 1,nflds
     do i = 1,nc
       ss(i,nf) = 0.0
     enddo
   enddo
#else
!
!  compute utility fields
!
   if(gauss_lat) then
     print *,' gaussian latitude'
     call post_gaussian_lat(jo,cosclt,wgt)
     do j = 1,johf
       slat(j)=cosclt(j)
       clat(j)=sqrt(1.-cosclt(j)**2)
       wlat(j)=wgt(j)
     enddo
   else
     print *,' equall latitude'
     call post_equall_lat(johf,slat,clat,wlat)
   endif
#define DEFAULT
#ifdef RFFTMLT
#undef DEFAULT
   call fftfax (io2/2,ifax,trig)
#endif
#ifdef DEFAULT
   call    fax (ifax,io2/2,3)
   call fftrig (trig,io2/2,3)
#endif
   call post_spectral_field(jcap,eps,epstop,enn1,elonn1,eon,eontop)
#endif
!
!  read sigma spectral data
!
   nr=(jcap+1)*(jcap+2)
   read(nss) (ss(i,1),i=1,nr)
#ifdef DBG
   print*,' read spectral coefficients :: '
  call print_maxmin_seven(ss(1,1),nr,nr,1,1,1,'ss record 1')
#endif
   read(nss) (ss(i,4),i=1,nr)
#ifdef DBG
   call print_maxmin_seven(ss(1,4),nr,nr,1,1,1,'ss record 2')
#endif
!
   do k = 1,levs
     read(nss) (ss(i,4*levs+6+k),i=1,nr)
   enddo
!
#ifdef DBG
   call print_maxmin_seven(ss(1,4*levs+7),nr,nr,levs,1,levs,'ss record 3')
#endif
!
   do k = 1,levs
     read(nss) (ss(i,6+k),i=1,nr)
     read(nss) (ss(i,levs+6+k),i=1,nr)
   enddo
!
#ifdef DBG
   call print_maxmin_seven(ss(1,7),nr,nr,levs,1,levs,'ss record 4')
   call print_maxmin_seven(ss(1,7+levs),nr,nr,levs,1,levs,'ss record 5')
#endif
!
   do k = 1,igases*levs+iwater*levs
     read(nss) (ss(i,5*levs+6+k),i=1,nr)
   enddo
!
#ifdef DBG
   call print_maxmin_seven(ss(1,5*levs+7),nr,nr,ktotal,1,ktotal,'ss record 6')
#endif
#ifdef SMP
!
   do k = 1,levs
     read(nss) (ss(i,2*levs+6+k),i=1,nr)  ! work space for SMP
   enddo
!
#endif
#ifndef SMP
!
   do k = 1,nflds
     do l = 0,jcap
       sstop(2*l+1,k)=0.
       sstop(2*l+2,k)=0.
     enddo
   enddo
! 
!  compute gradients and winds
!
   call post_grad_scalar(jcap,enn1,elonn1,eon,eontop,ss(1,1),                  &
                         ss(1,2),ss(1,3),sstop(1,3))
   call post_grad_scalar(jcap,enn1,elonn1,eon,eontop,ss(1,4),                  &
                         ss(1,5),ss(1,6),sstop(1,6))
!
   do k = 1,levs
     call post_divor2wind(jcap,enn1,elonn1,eon,eontop,ss(1,6+k),ss(1,levs+6+k),&
                ss(1,2*levs+6+k),ss(1,3*levs+6+k),                             &
                sstop(1,2*levs+6+k),sstop(1,3*levs+6+k))
   enddo
!
   if (keid.ne.0) then
     ! Hoon Park
!     call post_kinetic_energy(ss(1,2*levs+7),sstop(1,2*levs+7),               &
!                               lonf_,latg_,levs,                               &
!                               ken,jcap,del,eps,epstop)
     ! Koshyk (2001)
     call post_kinetic_energy2(ss(1,7),ss(1,levs+7),dke,vke,ken,nc,jcap,levs)
!
     call file_name('total_ke',8,fhour,ename,nchout)
     open(keid,file=ename(1:nchout),form='formatted',status='unknown')
     call file_name('diver_ke',8,fhour,ename,nchout)
     open(keid+1,file=ename(1:nchout),form='formatted',status='unknown')
     call file_name('vorti_ke',8,fhour,ename,nchout)
     open(keid+2,file=ename(1:nchout),form='formatted',status='unknown')
     do i = 0,jcap
       write(keid  ,trim(fmt))i,(ken(i,k),k=1,levs+1)
       write(keid+1,trim(fmt))i,(dke(i,k),k=1,levs+1)
       write(keid+2,trim(fmt))i,(vke(i,k),k=1,levs+1)
     enddo
     close(keid)
     close(keid+1)
     close(keid+2)
   endif
! 
#endif
#endif			/* DFS */
!
#ifdef DFS
   deallocate( ss )
#endif
#undef DEFAULT
!
   return
   end subroutine post_read_coeff
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_kinetic_energy2(div,vor,dke,vke,ken,nc,jcap,levs)
   use constant, only : rerth_
!-------------------------------------------------------------------------------
! Abstract : Calculate kinetic energy spectrum
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                                     ::  nc,jcap,levs,jc
   integer                                     ::  n,l,k,m,j
   real                                        ::  const
   real   , dimension(nc,levs)                 ::  div, vor
   real   , dimension(0:jcap,levs+1)           ::  dke,vke,ken
!
! offset for triangular coefficient
!
   jc(n,l) = (jcap+1)*(jcap+2)-(jcap+1-l)*(jcap+2-l)+2*(n-l)
!
! constant
!
   const=0.25*rerth_*rerth_
!
   do k = 1,levs
     dke(:,k)=0.
     vke(:,k)=0.
     ken(:,k)=0.
     do n = 1,jcap   ! total
       do m = 0,n    ! zonal
         j=jc(n,m)+1
         dke(n,k)=dke(n,k)+( div(j,k)**2+div(j+1,k)**2 )  ! divergent
         vke(n,k)=vke(n,k)+( vor(j,k)**2+vor(j+1,k)**2 )  ! rotational
       enddo
       dke(n,k)=dke(n,k)*const/real((n*(n+1)))
       vke(n,k)=vke(n,k)*const/real((n*(n+1)))
       ken(n,k)=dke(n,k)+vke(n,k)
     enddo
   enddo
   dke(:,levs+1)=0.
   vke(:,levs+1)=0.
   ken(:,levs+1)=0.
!
   return
!-------------------------------------------------------------------------------
   end subroutine post_kinetic_energy2
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
#define NWAV (mt+1)*(mt+2)
#ifdef DFS
   subroutine post_kinetic_energy(u,v,ib,jbw,km,ke,mt,del)
#else
   subroutine post_kinetic_energy(u,utop,ib,jbw,km,ke,mt,del,eps,epstop)
#endif
!-------------------------------------------------------------------------------
!#define BIN_DBG
   use constant, only : pi_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                         ::  mt,km,ib,jbw
#ifdef DFS
   real, dimension(ib,jbw,km)      ::  u,v
#else
   real, dimension(NWAV+1,2*km)    ::  u
   real, dimension((mt+1)*2,2*km)  ::  utop
   real, dimension(NWAV/2)         ::  eps
   real, dimension((mt+1))         ::  epstop
#endif
   real                            ::  del(km)
   real                            ::  ke(0:mt,km+1)
!
! local variables
! 
   integer                         ::  i,j,k,n,idir
!
! sperical harmonic coefficients
! 
   real                            ::  ske(NWAV,km+1)
#ifndef DFS
   real, dimension(NWAV/2   )      ::  pln
   real, dimension((mt+1)    )     ::  plntop
   real, dimension(ib*2+4 ,2*km)   ::  ug,duma
#ifdef BIN_DBG
   real, dimension(ib,jbw,2*km)    ::  wd
#endif
   real                            ::  gke(ib,jbw,km),gak(ib,jbw)
   real, allocatable, save         ::  glat(:),triga(:)
   real, save                      ::  ifaxa(20)
   integer                         ::  nctop,nc
   real                            ::  sinlat,coslat,rcs2,rlat,torad,rcos
   logical                         ::  first
   data first/.true./
#else
   real                            ::  gke(ib,jbw),gak(ib,jbw)
#endif
!-------------------------------------------------------------------------------
#ifndef DFS
   idir=1
   if (first) then
     first=.false.
     allocate(glat(jbw))
     allocate(triga(ib*2))
     call    fax(ifaxa,ib,3)
     call fftrig(triga,ib,3)
     call  chgr_gaussian_lat(glat,jbw)
   endif
! 
! wave to grid
!
   nc=NWAV+1
   nctop=(mt+1)*2
   torad=pi_/180.
!
   do j = 1,jbw/2
     ug(:,:)=0.0
     rlat=glat(j)*torad
     sinlat=sin(rlat)
     coslat=cos(rlat)
     rcos=1./coslat
     rcs2=1./(coslat*coslat)
     call post_sph_poly(mt,sinlat,coslat,eps,epstop,pln,plntop)
     call post_synth_coeff1(mt,ib+2,nc,nctop,2*km,pln,plntop,u,utop,ug)
     call fft99m (ug,duma,triga,ifaxa,1,ib+2,ib,4*km,1)
     !
     do i = 1,ib
       gak(i,j)=0.0
       gak(i,jbw-j+1)=0.0
     enddo
     do k = 1,km
       do i = 1,ib
         gke(i,j,k)=(ug(i,k)**2+ug(i,k+km)**2)*0.5*rcs2
         gak(i,j)=del(k)*gke(i,j,k)+gak(i,j)
         gke(i,jbw-j+1,k)=                                                     &
                (ug(ib+2+i,k)**2+ug(ib+2+i,k+km)**2)*0.5*rcs2
         gak(i,jbw-j+1)=del(k)*gke(i,jbw-j+1,k)+gak(i,jbw-j+1)
#ifdef BIN_DBG
         ! U
         wd(i,j,k)=ug(i,k)*rcos
         wd(i,jbw-j+1,k)=ug(i+ib+2,k)*rcos
         ! V
         wd(i,j,km+k)=ug(i,km+k)*rcos
         wd(i,jbw-j+1,km+k)=ug(i+ib+2,km+k)*rcos
#endif
       enddo
     enddo
   enddo
!
#ifdef BIN_DBG
   call file_write_bin(222,wd,ib,jbw,2*km,0)
   call file_write_bin(222,gke,ib,jbw,km,0)
   call file_write_bin(222,gak,ib,jbw,1 ,0)
#endif
!
! to wave
!
   do k = 1,km
     call post_trans_wave_grid(idir,gke(1,1,k),ske(1,k),0,0.,ib,jbw,mt,0)
   enddo
   call post_trans_wave_grid(idir,gak,ske(1,km+1),0,0.,ib,jbw,mt,0)
#else        /* DFS */
#ifdef BIN_DBG
   call file_write_bin(222,u,ib,jbw,km,0)
   call file_write_bin(222,v,ib,jbw,km,0)
#endif
   idir=101
   gak(:,:)=0.
   do k = 1,km
     do j = 1,jbw
       do i = 1,ib
         gke(i,j)= (u(i,j,k)**2+v(i,j,k)**2)*0.5
         !
         ! vertical average
         !
         gak(i,j)=del(k)*gke(i,j)+gak(i,j)
       end do
     end do
#ifdef BIN_DBG
     call file_write_bin(222,gke,ib,jbw,1,0)
#endif
     call post_trans_wave_grid(idir,gke,ske(1,k),0,0.,ib,jbw,mt,0)
   enddo        ! k
#ifdef BIN_DBG
   call file_write_bin(222,gak,ib,jbw,1,0)
#endif
   call post_trans_wave_grid(idir,gak,ske(1,km+1),0,0.,ib,jbw,mt,0)
#endif
!
   call post_wind_speed(ske,NWAV,km+1,ke,mt)
!
#ifdef BIN_DBG
   do k = 1,km+1
     call post_trans_wave_grid(-idir,gak,ske(1,k),0,0.,ib,jbw,mt,0)
     call file_write_bin(222,gak,ib,jbw,1,0)
   enddo
#endif
#undef NWAV
!
   return
   end subroutine post_kinetic_energy
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine post_wind_speed(var,nwav,km,vabs,mt)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                   ::  nwav,km,mt
   integer                   ::  j,k,n,l,jc
   real, dimension(nwav,km)  ::  var
   real, dimension(0:mt,km)  ::  vabs
!-------------------------------------------------------------------------------
   jc(n,l) = (mt+1)*(mt+2)-(mt+1-l)*(mt+2-l)+2*(n-l)
!
   do k = 1,km
     vabs(:,k)=0.
     do l = 0,mt
       do n = l,mt
         j=jc(n,l)+1
         vabs(n,k)=vabs(n,k)+sqrt(var(j,k)**2+var(j+1,k)**2)
       enddo
     enddo
   enddo
!
   return
   end subroutine post_wind_speed
!-------------------------------------------------------------------------------
