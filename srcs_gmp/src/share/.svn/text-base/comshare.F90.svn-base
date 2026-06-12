!
module comshare
!-------------------------------------------------------------------------------
   use paramodel, only : lnt2_, jcap1_,jcap2_
!-------------------------------------------------------------------------------
!
! common gozcom
!
   integer                                 ::  jfir
   real   , allocatable, dimension(:)      ::  dxa, dxb, dxc, dxd
!
! common scrtch
!
   real   , allocatable, dimension(:,:,:)  ::  d, z, u, v
!
! common comind
!                                                                               
!    indxnn(mdim) :  1-d index of converting input form spher coeff array         
!                    to transposed form array                                     
!    indxmm(mdim) :  1-d index of converting transposed form spher coeff          
!                    array to input form spherical coeff array                    
!  
   real   , allocatable, dimension(:)      ::  indxnn, indxmm
!
   contains
!-------------------------------------------------------------------------------
   subroutine comshare_init
!-------------------------------------------------------------------------------
   allocate(  dxa(lnt2_)                                                      ,&
              dxb(lnt2_)                                                      ,&
              dxc(lnt2_)                                                      ,&
              dxd(lnt2_)              )
   allocate(  d(2,jcap1_,jcap1_)                                              ,&
              z(2,jcap1_,jcap1_)                                              ,&
              u(2,jcap1_,jcap2_)                                              ,&
              v(2,jcap1_,jcap2_)      )
   allocate(  indxnn(lnt2_)                                                   ,&
              indxmm(lnt2_)           )

   end subroutine comshare_init
!-------------------------------------------------------------------------------
!
end module comshare
