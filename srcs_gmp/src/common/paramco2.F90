!
   module paramco2
!-------------------------------------------------------------------------------
   use paramodel, only  :  idim  => lonf_                                     ,&
                           jdim  => latg_                                     ,&
                           mwave => jcap_                                     ,&
                           kdim  => levs_
   use varsfc, only     :  lsoil => lsoil_
!
   end module paramco2
