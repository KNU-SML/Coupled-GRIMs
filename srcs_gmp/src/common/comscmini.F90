!
   module comscmini
!-------------------------------------------------------------------------------
   use paramodel, only  :  ilonf_,ilatg_, ilevs_
!-------------------------------------------------------------------------------
   integer             ::  imax, jmax, kmax, ijmax
   real, allocatable   ::  tsea  (:,:)                                        ,&
                           smc   (:,:,:)                                      ,&
                           snoweq(:,:)                                        ,&
                           stc   (:,:,:)                                      ,&
                           tg3   (:,:)                                        ,&
                           z0cm  (:,:)                                        ,&
                           cv    (:,:)                                        ,&
                           cvb   (:,:)                                        ,&
                           cvt   (:,:)                                        ,&
                           albedo(:,:,:)                                      ,&
                           slmsk (:,:)                                        ,&
                           plantr(:,:)                                        ,&
                           canopy(:,:)                                        ,&
                           f10m  (:,:)                                        ,&
                           z0cmt (:,:)
!
   contains
!-------------------------------------------------------------------------------
   subroutine comscmini_init
!-------------------------------------------------------------------------------
   imax=ilonf_
   jmax=ilatg_
   kmax=ilevs_
   ijmax=imax*jmax+1
   allocate(               tsea  (imax,jmax)                                  ,&
                           smc   (imax,jmax,2)                                ,&
                           snoweq(imax,jmax)                                  ,&
                           stc   (imax,jmax,2)                                ,&
                           tg3   (imax,jmax)                                  ,&
                           z0cm  (imax,jmax)                                  ,&
                           cv    (imax,jmax)                                  ,&
                           cvb   (imax,jmax)                                  ,&
                           cvt   (imax,jmax)                                  ,&
                           albedo(imax,jmax,1)                                ,&
                           slmsk (imax,jmax)                                  ,&
                           plantr(imax,jmax)                                  ,&
                           canopy(imax,jmax)                                  ,&
                           f10m  (imax,jmax)                                  ,&
                           z0cmt (imax,jmax)                                   )
!
   end subroutine comscmini_init
!-------------------------------------------------------------------------------
   end module comscmini
