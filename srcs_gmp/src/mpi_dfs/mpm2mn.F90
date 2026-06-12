!
   subroutine mpm2mn(a,jld,b,jl,mt,nvar)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpm2mn
!            
!-------------------------------------------------------------------------------
   use commpi, only  :  ncol
   use dfsvar, only  :  ngs,nge,nls
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  mt,jld,jl,nvar
   real  ::  a(mt,jld,nvar),b(mt,jl,nvar)
!
! local
!
   integer  ::  n,j,jj
!
! option for 1-d decomposition
!
   if( ncol.eq.1 ) then
     do n = 1,nvar
       do j = 1,jl
         b(1:mt,j,n) = a(1:mt,j,n)
       enddo
     enddo
     return
   endif
!
! otherwise
! cut the part and through away the rest.
!
   do n = 1,nvar
     do j = 1,jl
       jj=nls-ngs+j
       b(1:mt,j,n)=a(1:mt,jj,n)
     enddo
   enddo
!
   return
   end subroutine mpm2mn
!-------------------------------------------------------------------------------
