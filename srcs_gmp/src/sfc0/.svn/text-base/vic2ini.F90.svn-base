!
   subroutine vic2ini(ijdim,im)
!-------------------------------------------------------------------------------
!  
! subprogram: vic2ini
! 
! abstract: initialize the rest variables and parameters after merging
!   grib data, analysis data, and other land surface initial dat
!
! program history log:
!   2003-03-09  ji chen                cvs verion setup
!   2008-03-09  kyeng hee seol         debugging
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
! 
! The following parameters and variables should be initialized after 
! calling sfc
!    
!   jtsf,jsmc,jsno,jstc,jtg3,jzor,jcv ,jcvb,jcvt,jalb,
!   jsli,jveg,jcpy,jf10,jvet,jrot,jalf,just,jffm,jffh,
!   jbif,jds ,jdsm,jws ,jcef,jexp,jkst,jdph,jbub,jqrt,
!   jbkd,jsld,jwcr,jwpw,jsmr,jslz,jsnz,jnve,jlai
!
! We need initialize the following parameters and variables here:
!
!     Name   level        Comments
!  1. jsmx:  lsoil        maximum soil moisture (mm)
!  2. jdpn:  nsoil        thickness of soil node (m)
!  3. jsxn:  nsoil        maximum sm at soil node (m3/m3)
!  4. jepn:  nsoil        Para the vari of Ksat at soil node (N/A)
!  5. jbbn:  nsoil        bubbling pressure at soil node (cm)
!  6. japn:  nsoil        para alpha at soil node
!  7. jbtn:  nsoil        para beta at soil node
!  8. jgmn:  nsoil        para gamma at soil node
!  9. jmsn:  msub         snow depth for tile (m)
! 10. jmsm:  kslmb        soil moisture for tile (mm)
! 11. jmsi:  kslmb        soil ice for tile (mm)
! 12. jmst:  nslmb        soil temperature for tile (k)
! 13. jcsn:  msub         canopy snow  (m h2o)
! 14. jrsn:  msub         snow density (kg/m^3)
! 15. jtsn:  msub         snow surface temperature (K)
! 16. jtpk:  msub         snow pack temperature (K)
! 17. jsfw:  msub         surface snow water equivalent (m h2o)
! 18. jpkw:  msub         snow pack snow water equivalent (m h2o)
! 19. jlst:  msub         time step since last snow fall
!-------------------------------------------------------------------------------
   use varsfc, only : lsoil=>lsoil_,nsoil=>nsoil_, msub=>msub_
!  use comsfc, only : sfcfcs
   use comsfc
#ifdef VICLSM2
   use vic_veglib
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer ijdim,im
!
#ifdef VICLSM2
!
   integer,parameter    ::  kslmb=lsoil*msub,nslmb=nsoil*msub

   integer ij, m, n

   real                 ::  vdep(lsoil), vdzn(nsoil), sumd, sumdp, totdep
   real                 ::  smx(lsoil), bub(lsoil), expt(lsoil)
   real                 ::  smmaxnd(nsoil),expnd(nsoil),bubnd(nsoil)
   real                 ::  alpha(nsoil),beta(nsoil),gamma(nsoil)
!
   real                 ::  tmp_T, tmp_smc, smmax, buble, exple, sfc_max_unfwat
!
   integer,parameter    ::  t0c = 273.15
!
!-------------------------------------------------------------------------------
!
   do ij = 1,ijdim
     if(sfcfcs(ij,jsli).eq.1.) then
       do m = 1,msub
         sfcfcs(ij,jmsn+m-1) = sfcfcs(ij,jtsf)
         sfcfcs(ij,jcsn+m-1) =   0.0
         sfcfcs(ij,jrsn+m-1) =  50.0
         if(sfcfcs(ij,jsno).gt.0.) then
           sfcfcs(ij,jlst+m-1) =  0.0
         else
           sfcfcs(ij,jlst+m-1) = -1.0
         endif
         sfcfcs(ij,jtsn+m-1) = min(sfcfcs(ij,jtsf),t0c)
         sfcfcs(ij,jtpk+m-1) = min(sfcfcs(ij,jtsf),t0c)
         sfcfcs(ij,jsfw+m-1) = 0.0
         sfcfcs(ij,jpkw+m-1) = 0.0
       enddo

       do m = 1,lsoil
         if(sfcfcs(ij,jsld+m-1).gt.0) then
           sfcfcs(ij,jsmx+m-1) = (sfcfcs(ij,jsld+m-1)-                         &
                    sfcfcs(ij,jbkd+m-1))/sfcfcs(ij,jsld+m-1) ! fraction
         else
           sfcfcs(ij,jsmx+m-1) = 0.0
         end if
         if(sfcfcs(ij,jsmc+m-1).gt.sfcfcs(ij,jsmx+m-1)) then
           sfcfcs(ij,jsmc+m-1) = sfcfcs(ij,jsmx+m-1)
         endif
         sfcfcs(ij,jsmx+m-1)=sfcfcs(ij,jsmx+m-1)*                              &
                    sfcfcs(ij,jdph+m-1)*1000.0 ! maximum sm (mm)
       end do
       do m = 1,lsoil
         tmp_T= sfcfcs(ij,jstc+m)
         smmax= sfcfcs(ij,jsmx+m-1)
         buble= sfcfcs(ij,jbub+m-1)
         exple = sfcfcs(ij,jexp+m-1)
         if(tmp_T.lt.t0c)then
           tmp_ice=sfcfcs(ij,jsmc+m-1)*                                        &
                              sfcfcs(ij,jdph+m-1)*1000.0-                      &
                              sfc_max_unfwat(tmp_T,smmax,buble,exple)
           if(tmp_ice.lt.0.0) tmp_ice=0.0
         else
           tmp_ice = 0.0
         endif
         tmp_smc=sfcfcs(ij,jsmc+m-1)*sfcfcs(ij,jdph+m-1)*1000.0
         if(tmp_ice.gt.tmp_smc) tmp_ice = tmp_smc
         do nv = 1, msub
           sfcfcs(ij,jmsm+(m-1)*msub+nv-1)=tmp_smc
           sfcfcs(ij,jmsi+(m-1)*msub+nv-1)=tmp_ice
         enddo
       enddo
       do n = 1,nsoil
         do nv = 1, msub
           sfcfcs(ij,jmst+(m-1)*msub+nv-1)=sfcfcs(ij,jstc+m-1)
         enddo
       enddo
     endif
   end do
! 
!  define the depth of soil nodes, the number of soil nodes >= 3
!  currently, we are using a very simple way to find the soil temperture
!  profile for each soil node (Ji 2003)
!
   print *,'in vic1ini for computing node parameters'
!
   do ij = 1,ijdim
     if(sfcfcs(ij,jsli).eq.1.) then
       do m = 1,lsoil
         vdep(m) = sfcfcs(ij,jdph+m-1)
       end do
       if(vdep(1) .gt.0.0) then
         vdzn(1) = 0.0
         vdzn(2) = vdep(1)
         sumd = 0
         do n = 3,nsoil
           vdzn(n) = (n-1)*0.1
           sumd = sumd + vdzn(n)
         end do
!
         sumdp = 0
         do m = 2,lsoil
           sumdp = sumdp + vdep(m)
         end do
!
         do n = 3,nsoil
           if(sumd.gt.0) then
             vdzn(n) = vdzn(n)/sumd*sumdp
           else
             vdzn(n) = 0
           end if
         end do
       else
         do n = 1,nsoil
           vdzn(n) = 0.0
         end do
       end if
       do n = 1,nsoil
         sfcfcs(ij,jdpn+n-1) = vdzn(n)
       end do
     endif
   end do
!
!-- compute thermal node soil parameters ---
!
   do ij = 1,ijdim
     if(sfcfcs(ij,jsli).eq.1.) then
       totdep = 0
       do m = 1,lsoil
         vdep(m) =  sfcfcs(ij,jdph+m-1)
         totdep  = totdep + vdep(m)
       end do
       do n = 1,nsoil
         vdzn(n) = sfcfcs(ij,jdpn+n-1)
       end do
       if(totdep.gt.0.0) then
         do m = 1,lsoil
           smx(m)  = sfcfcs(ij,jsmx+m-1)
           bub(m)  = sfcfcs(ij,jbub+m-1)
           expt(m) = sfcfcs(ij,jexp+m-1)
         end do
         call nodepara(lsoil, nsoil,  vdep,    vdzn,                           &
                            smx,   bub,  expt, smmaxnd,                        &
                          expnd, bubnd, alpha,    beta,                        &
                          gamma)
         do n = 1,nsoil
           sfcfcs(ij,jsxn+n-1) = smmaxnd(n)
           sfcfcs(ij,jepn+n-1) = expnd(n)
           sfcfcs(ij,jbbn+n-1) = bubnd(n)
           sfcfcs(ij,japn+n-1) = alpha(n)
           sfcfcs(ij,jbtn+n-1) = beta(n)
           sfcfcs(ij,jgmn+n-1) = gamma(n)
         end do
       else
         do n = 1,nsoil
           sfcfcs(ij,jsxn+n-1) = 0.0
           sfcfcs(ij,jepn+n-1) = 0.0
           sfcfcs(ij,jbbn+n-1) = 0.0
           sfcfcs(ij,japn+n-1) = 0.0
           sfcfcs(ij,jbtn+n-1) = 0.0
           sfcfcs(ij,jgmn+n-1) = 0.0
         end do
       end if
     endif
   end do
!
#endif
   return
   end subroutine vic2ini
!-------------------------------------------------------------------------------
