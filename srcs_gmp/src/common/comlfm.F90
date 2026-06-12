!
   module comlfm
!-------------------------------------------------------------------------------
   use paramodel, only   :  LNT22S,LONF2S,LATG2S,levs_,levh_
   use varsfc, only      :  lsoil_,lalbd_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   private              ::  LNT22S,LONF2S,LATG2S,levs_,levh_
   private              ::  lsoil_,lalbd_
!
   integer              ::  klenp,nlfmsgi,nlfmsfi,nlfmsgo,nlfmsfo
   real                 ::  weight(1000)
   real                 ::  filtwin
   integer, allocatable ::  islmsk(:,:,:)
   real, allocatable    ::  fq (:)                                            ,&
                            fte(:,:)                                          ,&
                            fdi(:,:)                                          ,&
                            fze(:,:)                                          ,&
                            frq(:,:)                                          ,&
                            ftsea(:,:)                                        ,&
                            fsmc(:,:,:)                                       ,&
                            fsnoweq(:,:)                                      ,&
                            fstc(:,:,:)                                       ,&
                            ftg3(:,:)                                         ,&
                            fz0cm(:,:)                                        ,&
                            fplantr(:,:)                                      ,&
                            fcv(:,:)                                          ,&
                            fcvb(:,:)                                         ,&
                            fcvt(:,:)                                         ,&
                            falbedo(:,:,:)                                    ,&
                            fslmsk(:,:)                                       ,&
                            ff10m(:,:)                                        ,&
                            fcanopy(:,:)                                      ,&
                            wcvb(:,:)                                         ,&
                            wcvt(:,:)
!
   contains
!-------------------------------------------------------------------------------
   subroutine comlfm_init
!-------------------------------------------------------------------------------
   allocate(                islmsk(LONF2S,LATG2S,3)                            )
   allocate(                fq (LNT22S)                                       ,&
                            fte(LNT22S,levs_)                                 ,&
                            fdi(LNT22S,levs_)                                 ,&
                            fze(LNT22S,levs_)                                 ,&
                            frq(LNT22S,levh_)                                 ,&
                            ftsea(LONF2S,LATG2S)                              ,&
                            fsmc(LONF2S,LATG2S,lsoil_)                        ,&
                            fsnoweq(LONF2S,LATG2S)                            ,&
                            fstc(LONF2S,LATG2S,lsoil_)                        ,&
                            ftg3(LONF2S,LATG2S)                               ,&
                            fz0cm(LONF2S,LATG2S)                              ,&
                            fplantr(LONF2S,LATG2S)                            ,&
                            fcv(LONF2S,LATG2S)                                ,&
                            fcvb(LONF2S,LATG2S)                               ,&
                            fcvt(LONF2S,LATG2S)                               ,&
                            falbedo(LONF2S,LATG2S,lalbd_)                     ,&
                            fslmsk(LONF2S,LATG2S)                             ,&
                            ff10m(LONF2S,LATG2S)                              ,&
                            fcanopy(LONF2S,LATG2S)                            ,&
                            wcvb(LONF2S,LATG2S)                               ,&
                            wcvt(LONF2S,LATG2S)                                )
!
   end subroutine comlfm_init
!-------------------------------------------------------------------------------
   end module comlfm
