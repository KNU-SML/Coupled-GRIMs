   program sph_reduced_grid_main
!-------------------------------------------------------------------------------
!
! program: sph_reduced_grid_main
!
! program history log:
!   1988-04-25  joseph sela
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  hann-ming henry juang  reduced grid 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only               : jcap_,lnt2_,lnut2_,latg2_,lonf_
   use module_sph_legendre, only     : sph_gaussian_lat,                       &
                                       sph_ploy_funct_init, sph_ploy_funct
   use module_sph_reduced_grid, only : sph_reduced_grid
!-------------------------------------------------------------------------------
   real                 ::  tmpqtt(lnt2_)
   real                 ::  tmpqvv(lnut2_)
   real                 ::  colrad(latg2_)
   real                 ::  wgt(latg2_)
   real                 ::  wgtcs(latg2_)
   real                 ::  rcs2(latg2_)
!-------------------------------------------------------------------------------
   read(5,*) ndigit
   print *,' ndigit = ',ndigit
!
   call sph_gaussian_lat(latg2_, colrad, wgt, wgtcs, rcs2)
   print *,' colrad ',colrad
!
   call sph_ploy_funct_init
!
   qmaxall=0.0
   do lat = 1,latg2_
     call sph_poly_funct(tmpqtt,tmpqvv,colrad,lat)
     do n = 1,lnt2_
       qmaxall=max(qmaxall,abs(tmpqtt(n)))
     enddo
   enddo
!
   qttcut=qmaxall/(10.**ndigit)
   print *,' qttcut = ',qttcut
!
   sum_reduce=0.0
   sum_total=0.0
   do lat = 1,latg2_
     call sph_poly_funct(tmpqtt,tmpqvv,colrad,lat)
! ----- dyn_trans2model_grid, dyn_trans2output_grid depend on  sph_reduced_grid
     call sph_reduced_grid(tmpqtt,lnt2_,jcap_,qttcut,                          &
                   lcapd,lonfd)
     print *,' lat lcapd lonfd ',lat,lcapd,lonfd
     sum_reduce=sum_reduce+lonfd
     sum_total=sum_total+lonf_
   enddo
   reduce=(1.0-sum_reduce/sum_total)*100.
   print *,' percentage of reduce ',reduce
!
   stop
   end program sph_reduced_grid_main
