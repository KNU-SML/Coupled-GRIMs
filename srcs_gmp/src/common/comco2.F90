!
   module comco2
!-------------------------------------------------------------------------------
   use paramodel, only : levp1_, levp2_
!-------------------------------------------------------------------------------
!
!  common input
!
   real     ::  p1, p2, trnslo
   integer  ::  ia, ja, n
!
!  common press
!
   real     ::  pa(109)
!
!  common tran
!
   real     ::  transa(109,109)
!
!  common output
!
   real, allocatable, dimension(:,:)  ::  trns
!
!  common inputp
!
   real, allocatable, dimension(:)    ::  p, pd
!
!  common coefs
   real               ::  xa(109), ca(109), eta(109), sexpv(109)
   real               ::  core, uexp, sexp
!
   contains
!-------------------------------------------------------------------------------
   subroutine comco2_init
!-------------------------------------------------------------------------------
   allocate(trns(levp1_,levp1_))
   allocate(p(levp1_),pd(levp2_))
!
   end subroutine comco2_init
!-------------------------------------------------------------------------------
   end module comco2


