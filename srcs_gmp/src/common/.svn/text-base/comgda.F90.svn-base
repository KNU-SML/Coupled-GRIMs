!
   module comgda
!-------------------------------------------------------------------------------
   use paramodel, only   : LONF2S,LATG2S,levs_
!-------------------------------------------------------------------------------
   private              :: LONF2S,LATG2S,levs_
!
!  diagnostic indexes and flags
!
   integer, parameter   ::  kdgda=13
   integer, parameter   ::  kdtlarg=1,kdtconv=2,kdqconv=3,kdtshal=4,kdqshal=5, &
                            kdtvrdf=6,kduvrdf=7,kdvvrdf=8,kdqvrdf=9,           &
                            kdthsw=10,kdthlw=11,kdtcld=12,kdtccv=13
   integer, parameter   ::  ntgda=92
   character(len=8)     ::  cnmgda(kdgda)
!
   integer              ::  ipugda(kdgda),ibmgda(kdgda)
   integer              ::  nrgda, nwgda
   real,allocatable     ::  gdd(:)
!
   contains
!-------------------------------------------------------------------------------
   subroutine comgda_init
!-------------------------------------------------------------------------------
   nrgda=LATG2S
   nwgda=((LONF2S*levs_-1)/512+1)*512
   allocate(gdd(nwgda*kdgda*nrgda))
!
   end subroutine comgda_init
!-------------------------------------------------------------------------------
   end module comgda
