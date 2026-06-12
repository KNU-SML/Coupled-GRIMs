#include <define.h>
   SUBROUTINE dfs_inverse_matrix( A,N,NP,INDX,Y )  ! N=NP; Inverse of A= Y
                                                   ! A is destroyed when returned.
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [dfs_inverse_matrix] *
!           |
!           |-- [dfs_lu_backsub] *
!           |-- [dfs_lu_decomp] *
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
   DIMENSION A(NP,NP),Y(NP,NP),INDX(NP)
!-------------------------------------------------------------------------------
!
   do i = 1,N
     do j = 1,N
       Y(i,j)= 0.
     end do
     Y(i,i)= 1.D0
   end do
!
   CALL dfs_lu_decomp( A,N,NP,INDX,D )
!
   do j = 1,N
     CALL dfs_lu_backsub( A,N,NP,INDX,Y(1,j) )    ! Y= Inverse Matrix of A
   end do                                         ! A is destroyed.
!
   RETURN
   END SUBROUTINE dfs_inverse_matrix
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_lu_backsub( A,N,NP,INDX,B )
!-------------------------------------------------------------------------------
   DIMENSION A(NP,NP),INDX(N),B(N)
!-------------------------------------------------------------------------------
   ii= 0
!
   do i = 1,N ! 12
     LL= indx(i)
     sum= b(LL)
     b(LL)= b(i)
!
     if(ii.ne.0) then
       do j = ii,i-1
         sum= sum-a(i,j)*b(j)
       end do
     elseif(sum.ne.0.) then
       ii= i
     endif
!
     b(i)= sum
   end do ! 12
!
   do i = N,1,-1
     sum= b(i)
     do j = i+1,N
       sum= sum-a(i,j)*b(j)
     end do
     b(i)= sum/a(i,i)
   end do
!
   RETURN
   END SUBROUTINE dfs_lu_backsub
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_lu_decomp( A,N,NP,INDX,D )
!-------------------------------------------------------------------------------
   IMPLICIT  REAL*8 (A-H,O-Z)
!
   DIMENSION A(NP,NP),INDX(N),VV(999)
   DATA TINY/1.D-20/
!-------------------------------------------------------------------------------
!
   if(N.GT.999) stop ' Too small work-array size ! '
   D= 1.D0
!
   do i = 1,N
     aamax= 0.
!
     do j = 1,N
#ifdef IBMSP
       if(abs(dble(a(i,j))).gt.aamax) aamax= abs(dble(a(i,j)))
#else
       if(dabs(dble(a(i,j))).gt.aamax) aamax= dabs(dble(a(i,j)))
#endif
     end do
     if(aamax.eq.0) then
       write(*,*) 'error in LUDCMP : Singular Matrix in ludcmp'
       call exit(11)
     endif
     vv(i)= 1.D0/aamax
   end do
!
   do j = 1,N ! 19
!
     do i = 1,j-1
       sum= a(i,j)
       do k = 1,i-1
         sum= sum-a(i,k)*a(k,j)
       end do
       a(i,j)= sum
     end do
!
     aamax= 0.
     do i = j,N ! 16
       sum= a(i,j)
       do k = 1,j-1
         sum= sum-a(i,k)*a(k,j)
       end do
       a(i,j)= sum
#ifdef IBMSP
       dum= vv(i)*abs(dble(sum))
#else
       dum= vv(i)*dabs(dble(sum))
#endif
       if(dum.ge.aamax) then
         imax= i
         aamax= dum
       endif
     end do ! 16
!
     if(j.ne.imax) then
       do k = 1,N
         dum= a(imax,k)
         a(imax,k)= a(j,k)
         a(j,k)= dum
       end do
       D= -D
       vv(imax)= vv(j)
     endif
!
     indx(j)= imax
     if(a(j,j).eq.0) a(j,j)= TINY
!
     if(j.ne.N) then
       dum= 1.D0/a(j,j)
       do i = j+1,N
         a(i,j)= a(i,j)*dum
       end do
     endif
   end do ! 19
!
   RETURN
   END SUBROUTINE dfs_lu_decomp
!-------------------------------------------------------------------------------
