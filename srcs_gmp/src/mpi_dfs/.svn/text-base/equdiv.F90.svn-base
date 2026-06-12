!
   subroutine equdiv(len,ncut,lenarr)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    equdiv
!            
! abstract: cut len into ncut pieces with load balancing
!
! usage:   equdiv(len,ncut,lenarr)
!
!    input argument lists:
!   len   - integer total length 
!   ncut   - integer number of subgroup
!
!    output argument list:
!   lenarr   - integer (ncut) length of each subgroup
! 
! subprograms called: none
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  len,ncut,n0,n1,n
   integer  ::  lenarr(ncut)
!
   n0=len/ncut
   n1=mod(len,ncut)
!
   do n = 0,n1-1
     lenarr(ncut-n)=n0+1
   enddo
   do n = 1,ncut-n1
     lenarr(n)=n0
   enddo
!
   return
   end subroutine equdiv
!-------------------------------------------------------------------------------

