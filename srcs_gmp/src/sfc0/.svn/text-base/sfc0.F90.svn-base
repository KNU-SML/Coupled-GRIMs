#include <define.h>
   program sfc0
!-------------------------------------------------------------------------------
!
! conversion of various types of sfc files, then
! reads grib sfc files and generates initial sfc file
!
!-------------------------------------------------------------------------------
   use paramodel, only               :  para_init
#ifdef SMP_RA2SFC
   use paramodel, only               :  ILOTS=>lonf2scm_,JLOTS=>latg2scm_
   use paramodel, only               :  GXLON=>gxlon_,GYLAT=>gylat_
#else
   use paramodel, only               :  ILOTS=>LONF2F,JLOTS=>LATG2F
#endif
   use vargrb
   use varsfc
   use comsfc                        ! init_comsfc,sfcfcs,sfcftyp
#ifdef AQUA_PLANET
   use constant, only                : t0c_, pi_,degrad_
#endif
#ifdef RMP
   use rscomloc
#endif
   use module_trans_reduce, only     :  dyn_trans2model_grid,                  &
                                        dyn_trans2output_grid 
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                          ::  idim,jdim,ijdim
!
! initialize variables
!
   character(len=128)               ::  fnin,fnout
   character(len=4)                 ::  sfcftypin
   character(len=8), allocatable    ::  svar(:)     
   character(len=8)                 ::  labs(4)
   real,             allocatable    ::  sfcfcsin(:,:)     
   real                             ::  fh
   integer                          ::  nsfc
   integer,          allocatable    ::  lev(:)
   integer                          ::  numsfcvin,maxlev,numchi
   integer                          ::  iy,im,id,ih,jy,jm,jd,jh
#ifdef RMP
!  real rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,                      &
!           rlftgrd,rbtmgrd
!  real cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,                      &
!           clftgrd,cbtmgrd
   namelist /namloc/                                                           &
   rproj,rtruth,rcotru,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd,    &
   cproj,ctruth,ccotru,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd
#endif
   integer                          ::  n,i,j,numsfcsin,ij
   real                             ::  dlat,dlon
#ifdef SMP_RA2SFC
   real                             ::  scmsfcfcs(2,numsfcs)
   real, allocatable                ::  ylin(:),work(:)
   real, allocatable                ::  rlon(:),rlat(:)
   real                             ::  gval(numsfcs)
   real*8, allocatable              ::  coscolat(:),wrk(:)
   real                             ::  xlon,ylat,gvl
   integer,parameter                ::  sidim=1, sjdim=2
   integer                          ::  i,k,ij,idx
!
   real             ::  smcfracf(lsoil_)
   real,allocatable ::  slmsks(:),stypes(:),stcs(:,:)
   real,allocatable ::  smcs(:,:),slcs(:,:)
#endif
#ifdef AQUA_PLANET
   real, allocatable                ::  rlon(:),rlat(:),glat(:),glon(:)
#endif
   data     nsfc,fh /19,0./
!
!  namelist variables
!    nsfc     : integer   .. unit number for i/o
!    fnin     : character .. input sfc file name
!    fnout    : character .. output sfc file name
!    sfcftypin: character .. input sfc file type [osu1,osu2,noa1,vic1]
!
   namelist/namsfc0/nsfc,fnin,fnout,sfcftypin,fh
!
   call para_init
!
   idim=ILOTS/2
   jdim=JLOTS*2
   ijdim=idim*jdim
#ifdef SMP_RA2SFC
   allocate(ylin(jdim),work(ijdim))
   allocate(coscolat(jdim),wrk(jdim))
   allocate(slmsks(ijdim),stypes(ijdim),stcs(ijdim,lsoil_))
   allocate(smcs(ijdim,lsoil_),slcs(ijdim,lsoil_))
#endif
#ifdef AQUA_PLANET
   allocate(rlat(jdim),rlon(idim),glat(ijdim),glon(ijdim))
#endif
!-------------------------------------------------------------------------------
   read(95,namsfc0)
#ifdef RMP
   read(95,namloc)
#endif
!
! call by org grid size 
!
   call init_comsfc(idim,jdim)
!
! find this file type's record structure and allocate space
! for reading variables.
!
!    numsfcv is number of surface file variables
!    numsfcs is number of surface file 2-D recs.
!
   call sfcfld(sfcftypin,0,numsfcvin,lev,svar,maxlev)
   allocate (lev (numsfcvin))
   allocate (svar(numsfcvin))
   call sfcfld(sfcftypin,1,numsfcvin,lev,svar,maxlev)
   numsfcsin=0
   do n = 1,numsfcvin
     numsfcsin=numsfcsin+lev(n)
   enddo
   allocate (sfcfcsin(ijdim,numsfcsin))
!
!  read input sfc file. the file type is sfcftypin and may
!  not be the same as the sfc file type of the model you are running.
!
   print*, '1st sfc_read_file_driver, nsfc,fnin,sfcftypin,labs,idim,jdim',     &
                        nsfc,fnin,sfcftypin,labs,idim,jdim
   call sfc_read_file_driver(nsfc,fnin,sfcftypin,                              &
              labs,iy,im,id,ih,fh,                                             &
              sfcfcsin,idim,jdim,0)
#ifdef GDAS
   call incdte(iy,im,id,ih,                                                    &
               jy,jm,jd,jh,nint(fh))
   iy=jy
   im=jm
   id=jd
   ih=jh
   fh=0.
#else
   if(fh.ne.0.) then
     print *,'fh.ne.0. reset to 0.'
     fh=0.
   endif
#endif
!
!  convert input sfc file type to current model surface type
!
   if(sfcftypin.ne.sfcftyp) then 
     call sfcfcnv(sfcftypin,sfcfcsin,idim,jdim,numsfcsin)
   else
     do n = 1,numsfcsin
       do j = 1,ijdim
         sfcfcs(j,n)=sfcfcsin(j,n)
       enddo
     enddo
   endif
!
!  create initial states by reading grib files and
!  merging it to the sfc file. output is sfcfcs in common
!
   call sfc_cycle_driver(sfcftypin,nsfc,iy,im,id,ih,fh)
#ifdef VICLSM1
!
!  initialize some specific parameters and variables used by VIC1
!  one-tile
!
   print *,' in sfc0 before vic1ini'
!
   if(sfcftypin.ne.sfcftyp) then
     call vic1ini(ijdim,im)
   endif
#endif
#ifdef VICLSM2
!
!  initialize some specific parameters and variables used by VIC2
!  multi-tile
!
   print *,' in sfc0 before vic2ini'
!
   if(sfcftypin.ne.sfcftyp) then
     call vic2ini(ijdim,im)
   endif
#endif
!
!  write out
!
#ifdef SMP_RA2SFC
!
!  interpolate for smp obs site
!
   xlon = GXLON
   ylat = GYLAT
   call sfc_gausst62_smp(ylin,jdim,coscolat,wrk)
   print*,'sfc0 xlon ylat ',xlon,ylat
!
   do n = 1, numsfcs
     call dyn_trans2output_grid(sfcfcs(ij,n),1)
     do ij = 1, ijdim
       work(ij)=sfcfcs(ij,n)
     enddo
     call sfc_hintp_smp(idim,jdim,ylin,work,xlon,ylat,gvl)
     if (n.eq.jsli) then
       gvl=real(ifix(gvl+.5))
       print *,'slimsk ',gvl-1.
     endif
#ifdef VICLSM1
     if (n.eq.jvet) then
       gvl = 11.0
     endif
#endif
     gval(n) = gvl
     write(*,'(a,i3,2x,f12.7)') 'sfc0 n gval ',n,gval(n)
     call dyn_trans2model_grid(sfcfcs(ij,n),1)
   enddo
!
! use observed sm
!
   idx = jsmc
   write(6,*) 'sfc0 idx jsmc',idx,jsmc
#ifdef OSULSM2
!
! 10 cm SM for osu2
!
   smcfracf(1)=0.345
   gval(jsmc)=smcfracf(1)*(0.387-0.119)+0.119
!
! 200 cm SM for osu2
!
   smcfracf(1)=(0.16092-0.1)/(0.47-0.1)
   gval(jsmc+1)=smcfracf(1)*(0.387-0.119)+0.119
!
!   gval(jsmc)=0.345
!
   write(6,*) 'use the observed SM n=',jsmc,gval(jsmc)
   write(6,*) 'use the observed SM n=',jsmc+1,gval(jsmc+1)
#endif
#ifdef NOALSM1
!
! 10,40,100 cm SM for noah
!
   smcfracf(1)=0.345
   smcfracf(2)=0.384
   smcfracf(3)=(0.16092-0.1)/(0.47-0.1)
   smcfracf(4)=(0.16092-0.1)/(0.47-0.1)
   gval(jsmc)  =smcfracf(1)*(0.387-0.119)+0.119
   gval(jsmc+1)=smcfracf(2)*(0.387-0.119)+0.119
   gval(jsmc+2)=smcfracf(3)*(0.387-0.119)+0.119
   gval(jsmc+3)=smcfracf(4)*(0.387-0.119)+0.119
   write(6,*) 'use the observed SM n=',jsmc,gval(jsmc)
   write(6,*) 'use the observed SM n=',jsmc+1,gval(jsmc+1)
   write(6,*) 'use the observed SM n=',jsmc+2,gval(jsmc+2)
   write(6,*) 'use the R2       SM n=',jsmc+3,gval(jsmc+3)
!
! update slc depending on the soil moisutre
!
  do i = 1,sidim*sjdim
    slmsks(i)=gval(jsli)
    stypes(i)=gval(jsot)
!
    do k = 1,lsoil_
      stcs(i,k)=gval(jstc+k-1)
      smcs(i,k)=gval(jsmc+k-1)
      print*,'ready stc smc',stcs(i,k),smcs(i,k)
    enddo
   enddo
!
   call getslc(slmsks,stypes,stcs,smcs,ijdim,lsoil_,slcs)
!
   do k = 1,lsoil_
     gval(jslc+k-1)=slcs(1,k)
     print*,'update slc n=',jslc+k-1,gval(jslc+k-1)
   enddo
#endif
#ifdef VICLSM1
!
! 30 cm SM for vic
! soil moisture spinup
!
! 1st layer
!
   smcfracf(1)=0.363
   gval(jsmc)=smcfracf(1)*(0.329-0.139)+0.139
!
! 2nd layer
!
   smcfracf(2)=(0.16092-0.1)/(0.47-0.1)
   gval(jsmc+1)=smcfracf(2)*(0.329-0.139)+0.139
!
! 3nd layer
!
   smcfracf(3)=(0.16092-0.1)/(0.47-0.1)
   gval(jsmc+2)=smcfracf(3)*(0.329-0.139)+0.139
!
   write(6,*) 'use the observed SM n=',jsmc,  gval(jsmc)
   write(6,*) 'use the observed SM n=',jsmc+1,gval(jsmc+1)
   write(6,*) 'use the observed SM n=',jsmc+2,gval(jsmc+2)
#endif
!
#ifdef WET_SOIL
   idx = jsmc
   do n = 1,lsoil_
     if (gval(jsli).eq.1) then
       gval(idx)=0.47
       print *,'wet idx=',idx,gval(idx)
       idx = idx + 1
     endif
   enddo
#endif
#ifdef DRY_SOIL
!
   idx = jsmc
   do n = 1,lsoil_
     if (gval(jsli).eq.1) then
       gval(idx)=0.1
       print *,'dry idx=',idx,gval(idx)
       idx = idx + 1
     endif
   enddo
#endif
!
!  write sfco for smp
!
   do n = 1,numsfcs
     do j = 1,sidim
       do i = 1,sjdim
         scmsfcfcs(i*j,n)=gval(n)
       enddo
     enddo
   enddo
#endif
#ifdef AQUA_PLANET
#ifndef DFS
   call gaulat(rlat,jdim)
!
   do j = 1,jdim
     rlat(j)=90.-rlat(j)
   enddo
#else
   dlat = 180./jdim
   do j = 1,jdim
     rlat(j)=90.-(j-0.5)*dlat
   enddo
#endif
!
   dlon = 360./float(idim)
   do i = 1,idim
     rlon(i)=dlon*0.5 + float(i-1)*dlon
   enddo
!
   do ij = 1,ijdim
     i = mod(ij-1,idim) + 1
     j = int((ij-1)/(idim*2)) + 1
     glon(ij) = rlon(i)/degrad_
     if(i.gt.idim) then 
       glat(ij) = -rlat(j)/degrad_
     else
       glat(ij) = rlat(j)/degrad_
     endif
   enddo
!
   do ij = 1,ijdim
     if(abs(glat(ij)).lt.pi_/3.) then
       sfcfcs(ij,1) = t0c_ + 27.*(1.-sin(1.5*glat(ij))**2.)
     else
       sfcfcs(ij,1) = t0c_
     endif
#ifdef AQUA_3KW1
     if(abs(glat(ij)).lt.pi_/6.) then
       sfcfcs(ij,1) =  sfcfcs(ij,1)                                            &
                     + 3.*cos(glon(ij)-pi_/2.)*cos(pi_/2.*glat(ij)/(pi_/6.))**2
     endif
#endif
     print*,' ij glat glon sst ',ij,glat(ij)*degrad_,glon(ij)*degrad_,         &
            sfcfcs(ij,1)
   enddo
   smc(:,:,:) = 0.47
   snoweq(:,:) = 0.0
   stc(:,:,:) = 0.0
   tg3(:,:) = 0.0
   cv(:,:) = 0.7
   cvb(:,:) = 1
   cvt(:,:) = 1
   slmsk(:,:) = 0.0
   call print_maxmin_seven(sfcfcs(1,1),ijdim,ijdim,1,1,1,'sst in sfc0')
#endif
   print*,'write iy,im,id,ih,fh=',iy,im,id,ih,fh
   print*, '1st sfc_read_file_driver, nsfc,fnin,sfcftypin,labs,idim,jdim',     &
                        nsfc,fnin,sfcftypin,labs,idim,jdim
   call sfc_read_file_driver(nsfc,fnout,sfcftyp,                               &
              labs,iy,im,id,ih,fh,                                             &
#ifdef SMP_RA2SFC
              scmsfcfcs,sidim,sjdim,1)
#else
              sfcfcs,idim,jdim,1)
#endif
!
#ifdef SMP_RA2SFC
   deallocate(ylin,work)
   deallocate(coscolat,wrk)
   deallocate(slmsks,stypes,stcs)
   deallocate(smcs,slcs)
#endif
   deallocate (lev )
   deallocate (svar)
   deallocate (sfcfcsin)
   stop
   end
!-------------------------------------------------------------------------------
