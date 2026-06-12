#include "define.h"
   program chgr_main
!-------------------------------------------------------------------------------
!
! main program documentation block
!                .      .    .                                       .
! main program:  chgr        change resolution of mrf sig and sfc files.
!
! abstract: changes the horizontal and vertical resolution
!           of the spectral model sigma and surface files.
!           for the sigma file, first the input data
!           are truncated or extended (with zero coefficients)
!           to the specified output horizontal resolution.
!           then the data are transformed to gridspace.
!           there the surface pressure is interpolated to the
!           specified output orography.  the temperature lapse rate
!           in this calculation is computed from the data
!           except below ground where it is set to -6.5k/km.
!           if necessary, the input moisture profile is extended
!           to the full input model layers using climatology.
!           then the wind, temperature and moisture profiles
!           are interpolated to the pressure locations determined
!           from the output surface pressure and the specified
!           output sigma layer values.  the profiles are never
!           extrapolated but held constant outside the input domain,
!           except below the input surface the temperature lapse rate
!           is again fixed at -6.5k/km and the relative humidity
!           is held fixed instead of the specific humidity.
!           the data are then transformed back to spectral space.
!           for the surface file, the input data are interpolated
!           horizontally to the output grid.  the output sea-land mask
!           is not interpolated but specified.
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   1991-03-15  mark iredell  docblock written (prehistorical program)
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatble with namelist input
!
! usage:
!   input files:
!     unit05   - namelist namchgr
!                ltrn     = .true. for ibm order files
!                lnew     = .true. for newer moisture climatology
!                lrep     = .true. to use orography from unit12
!     unit11   - input sigma file
!                (if empty, sigma file interpolation is skipped)
!     unit12   - output spectral orography
!     unit13   - output sea-land mask
!     unit14   - input surface file
!                (if empty, surface file interpolation is skipped)
!
!   output files:
!     unit06   - printed output and diagnostics
!     unit51   - output sigma file
!     unit52   - output surface file
!
!   subprograms called:
!     chgr_new_sigma           - get output sigma values
!     chgr_sph_gaussian        - get gaussian latitudes
!     sph_poly_epsilon1        - compute spectral constants
!     chgr_read_sigma          - read input sigma file
!     chgr_horizontal_wavecut  - change horizontal spectral resolution
!     sph_comp_index           - compute indices for transposing
!     dyn_transpose_coeff      - transpose spectral coefficients
!     rmsgt                    - print statistics
!     chgr_vinterp_solver      - transform to grid and interpolate vertically
!     dyn_transpose_coeff      - transpose spectral coefficients
!     chgr_write_sigma         - write output sigma file
!     chgr_interp_gaussian_sph - gaussian grid to gaussian grid interpolation
!
!-------------------------------------------------------------------------------
   use paramodel, only : para_init
#ifdef DFS
   use paramodel, only : jcap_,latg_,lonf_,levs_,ngases_,nwater_
   use dfsvar, only    : get_dfs_dim,mt,jl,jla,coslat
#endif
   use constant, only : pi_
   use paramter
   use parmchgr
   use comchgr,  only : comchgr_init, mdim                                    ,&
#ifdef DFS
                        mdimd                                                 ,&
#endif
                        jdimhf, kdimp, ijdim, ijdimi                          ,&
                        lab, idate, delin                                     ,&
                        si, sl, del, ci, cl, rpi, ak5, bk5                    ,&
                        eps, colrad, wgt, wgtcs, rcs2                         ,&
                        q , te , di , ze , rq , gz                            ,&
                        qi, tei, dii, zei, rqi, gzi                           ,&
                        qo, teo, dio, zeo, rqo, gzo
   use module_sph_legendre, only : sph_poly_epsilon1
!-------------------------------------------------------------------------------
   real, allocatable  ::   gzso(:)
#ifdef DFS
   real, allocatable  ::   tave(:)
   real, allocatable  ::   zsa(:,:)
#endif
   logical ltrn,lnew,lrep
!
   character(len=4)     ::  sfcftyp
   character(len=16)    ::  svar(1000)
   real, allocatable    ::  ggridi(:,:),ggrido(:,:)
   real, allocatable    ::  slmsko(:)
   real                 ::  ensemble(2),dummy2(21)
   integer              ::  lev(1000)
   integer              ::  issfc,issig
   data issfc,issig/1,1/
!
   namelist/namchgr/ nin,nou,ltrn,lnew,lrep,sfcftyp,iunit,issfc,issig
!
   data nin,norog,nslmsk,nsfcin/11,12,13,14/
   data nou,nsfcou/51,52/
   data ltrn/.true./
   data lnew/.true./
   data lrep/.true./
   data iunit/18/
!-------------------------------------------------------------------------------
   call para_init
   call comchgr_init
#ifdef DFS
   allocate(gzso(mdimd),tave(kdim),zsa(idim,jdim))
#else
   allocate(gzso(mdim))
#endif
   allocate(slmsko(ijdim))
!
! initialize local variables
!
   gzso=0.  ; slmsko=0.
#ifdef DFS
   tave=0.  ; zsa=0.
#endif
!
   read (5,namchgr)
   write(6,namchgr)
!
#ifdef DFS
   call get_dfs_dim(jcap_,levs_,ngases_,nwater_,lonf_,latg_)
#ifdef ALIASED
#ifdef ALIASED2
   if (jdim-2+1.ne.jl) then
#else
   if ((mw-mod(mw,2)+1).ne.jl) then
#endif
     print*,'jl.ne.nw',jl,mw-mod(mw,2)+1
     call exit(1)
   endif
#endif /* ALIASED end */
!
! calculate coslat
!
   call dfs_sincos_lat
   rlat=pi_/jdim
   do j = 1,jdim/2
     colrad(j)=(j-0.5)*rlat
     rcs2(j)=1./(sin(colrad(j))*sin(colrad(j)))
   enddo
#else /* not DFS */
   call chgr_sph_gaussian(jdimhf, colrad, wgt, wgtcs, rcs2)
#endif /* DFS end */
   call sph_poly_epsilon1(eps,mwave)
!
!  read in original initial data
!
   if (issig.eq.1) then
!
   call chgr_read_sigma(nin,fhour,idate,gzi,qi,tei,dii,zei,rqi                 &
        ,waves,xlayers,trun,order,realform,gencode                             &
        ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,water                             &
        ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid                  &
        ,pdryini,dummy2,gases                                                  & 
        ,iret)
   if(iret.ne.0) then
     print *,'hit end of file reading first record of sig file'
     call abort
   endif
!
   print 1101,fhour,idate
   1101 format(1h0,1x,'original initial data read in',/,                       &
              3x,'fhour=',f6.1,' ihour=',i2,' month=',i2,' day=',i2,           &
              ' year=',i4)
   print 1102,mwavei,idimi,jdimi,kdimi,kdimqi
   1102 format(1x,'input  file resolution =',5i10)
#ifndef HYBRID
   call chgr_new_sigma(ci,si,del,sl,cl,rpi)
#else
   call chgr_hybrid_setup(ak5,bk5,ci,si,del,sl,cl,rpi)
   do k = 1,kdimp
     print 1002, k,ak5(k),bk5(k)
   enddo
   1002 format (1x,'ak bk',i4,2f15.7)
#endif
   print 1103,mwave ,idim ,jdim ,kdim ,kdimq
   1103 format(1x,'output file resolution =',5i10)
!
!  change horizontal resolution by padding 0 or truncation
!
   call chgr_horizontal_wavecut(mwavei, qi,mwave, q)
   call chgr_horizontal_wavecut(mwavei,gzi,mwave,gz)
   do k = 1,kdimi
     call chgr_horizontal_wavecut(mwavei,tei(1,k),mwave,te(1,k))
     call chgr_horizontal_wavecut(mwavei,dii(1,k),mwave,di(1,k))
     call chgr_horizontal_wavecut(mwavei,zei(1,k),mwave,ze(1,k))
   enddo
   do k = 1,kdimqi
     call chgr_horizontal_wavecut(mwavei,rqi(1,k),mwave,rq(1,k))
   enddo
!
#ifdef DFS
#undef CHGR_SIMPLE
#define ISO_THERMAL
#ifndef ISO_THERMAL
   call dfs_standard_temp(iunit,sl,tave,kdim,1000.)
#else
   tave(:)=300.
#endif
#endif /* DFS end */
#ifdef CHGR_SIMPLE
   do m = 1,mdim
     gzso(m)=gz(m)
     qo (m)=q (m)
   enddo
   if(lrep) then
     read(norog) gzso
   endif
   do k = 1,kdimi
     do m = 1,mdim
       teo(m,k)=te(m,k)
       dio(m,k)=di(m,k)
       zeo(m,k)=ze(m,k)
       rqo(m,k)=rq(m,k)
     enddo
   enddo
   ltrn=.false.
#else /* not CHGR_SIMPLE */
   call sph_comp_index
   if(ltrn) then
     call dyn_transpose_coeff(gz,1)
     call dyn_transpose_coeff( q,1)
     call dyn_transpose_coeff(te,kdimi)
     call dyn_transpose_coeff(di,kdimi)
     call dyn_transpose_coeff(ze,kdimi)
     call dyn_transpose_coeff(rq,kdimqi)
   endif
   if(lrep) then
     read(norog) gzso
     print *,'gz of output resolution read in'
#ifndef DFS
     if(ltrn) then
       call dyn_transpose_coeff(gzso,1)
     endif
#endif
   else
#ifdef DFS
     print*,'The lrep must .TRUE. for DFS'
     call exit(1)
#endif
     do m = 1,mdim
       gzso(m)=gz(m)
     enddo
   endif
   call chgr_rms_input(q,di,te,ze,delin,rq)
#ifdef DFS
   call chgr_vinterp_solver(lnew,gzso,tave)
#else
   call chgr_vinterp_solver(lnew,gzso)
#endif
#endif			/* not CHGR_SIMPLE end */
!
   call chgr_rms_output(qo,dio,teo,zeo,del,rqo)
#ifdef ALIASED
   call dfs_cut_alias(qo,mt,0,jla,1   )
   call dfs_cut_alias(dio,mt,0,jla,kdim)
   call dfs_cut_alias(teo,mt,0,jla,kdim)
   call dfs_cut_alias(zeo,mt,0,jla,kdim)
   call dfs_cut_alias(rqo,mt,0,jla,kdim)
   call dfs_cut_alias(gzso,mt,0,jla,1   )
#endif
#ifndef HYBRID
   call chgr_write_sigma(nou,fhour,idate,qo,teo,dio,zeo,rqo,sl,si,gzso         &
#else
   call chgr_write_sigma(nou,fhour,idate,qo,teo,dio,zeo,rqo,ak5,bk5,gzso       &
#endif
#ifdef DFS
        ,tave                                                                  &
#endif
        ,ltrn                                                                  &
        ,waves,xlayers,trun,order,realform,gencode                             &
        ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,water                             &
        ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid                  &
        ,pdryini,dummy2,gases                                                  &
        ,iret)
!
! issig.eq.1
!
   endif
!
!  addition for the change resolution of surface file
!  note that this program uses only bilinear interpolation, disregarding
!  whether is it land-ocean mask related or it is index of types (e.g. 
!  vegetaion index).  It is assumed that sfc program is run to overwrite
!  these inproperly interpolated fields from higher resolution climatology.
!  In therms of scripts, SFC0=yes should be specified.
!
   if (issfc.eq.1) then
!
     call sfcfld(sfcftyp,1,nrecs,lev,svar,mxlv)
     allocate (ggridi(ijdimi,mxlv),ggrido(ijdim,mxlv))
!
     rewind nsfcin
     read (nsfcin,end=800,err=800) lab
     print *,' end read lab of sfc.',lab,' nrec=',nrecs
     write(nsfcou) lab
     read (nsfcin) fhour,idate
     print *,' read fhour,idate of sfc.'
     print *,'fhour,idate=',fhour,idate
     write(nsfcou) fhour,idate
!
     do nfld = 1,nrecs
       print*,'chgr_main:nsfcin,nfld,lev=',nsfcin,nfld,lev(nfld)
       read(nsfcin,end=909,err=909)                                            &
          ((ggridi(ij,l),ij=1,ijdimi),l=1,lev(nfld))
       print *,' read gaussian grid sfc. rec =',nfld
!
!  special handling of land sea mask
!   interpolate sea ice mask first.  Other areas, use input mask (from mtn).
!
#ifdef DFS
#define GL2GL chgr_interp_gaussian_dfs
#else
#define GL2GL chgr_interp_gaussian_sph
#endif
       if(nfld.eq.11) then
         do i = 1,ijdimi
           if(ggridi(i,1).eq.1.) then
             ggridi(i,1)=0.
           endif
         enddo
#ifdef DBG
         print *,'GL2GL (1)'
         call chgr_comp_maxmin(ggridi,idimi,jdimi,1,' ggrdi  ')
#endif
         call GL2GL(ggridi,idimi,jdimi,ggrido,idim,jdim)
         read(nslmsk) slmsko
         do i = 1,ijdim
           if(ggrido(i,1).gt.1..and.slmsko(i).eq.0.) then
             ggrido(i,1)=2.
           else
             ggrido(i,1)=slmsko(i)
           endif
         enddo
         call chgr_comp_maxmin(ggrido,idim,jdim,1,' ggrdo  ')
       else
         do l = 1,lev(nfld)
           call GL2GL(ggridi(1,l),idimi,jdimi,ggrido(1,l),idim,jdim)
           call chgr_comp_maxmin(ggrido(1,l),idim,jdim,1,' ggrido ')
         enddo
       endif
       write(nsfcou) ((ggrido(ij,l),ij=1,ijdim),l=1,lev(nfld))
     enddo
     stop
     909 continue
     print *,'hit eof or error reading sfc at nfld = ',nfld, ' .. but continued'
     write(nsfcou) ((ggrido(ij,l),ij=1,ijdim),l=1,lev(nfld))
!    call abort
     800 continue
     print *,'warning warning warning -- surface file empty'
     deallocate (ggridi,ggrido)
!
! issfc.eq.1
!
   endif
#ifdef DFS
   deallocate(tave,zsa)
#endif
   deallocate(gzso,slmsko)
!
   stop
   end program chgr_main
!-------------------------------------------------------------------------------
