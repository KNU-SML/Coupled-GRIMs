#include "define.h"
   subroutine equdis(ind,ls,le,ncut,lenarr,lendef)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    equdis
!            
! abstract: cut len into ncut pieces with load balancing by 
!           symmetric distribution
!
! usage:   equdiv(len,ncut,lenarr)
!
!    input argument lists:
!   ind   - integer spread direction:       1 for regular,
!                                          -1 for reverse
!   len   - integer total length 
!   ncut   - integer number of subgroup
!
!    output argument list:
!   lenarr   - integer (ncut) length of each subgroup
!   lendef   - integer (len) redefine the index 
! 
! subprograms called: none
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  ind,len,ncut,nn,n,i,lens,lene,ls,le
   integer  ::  lenarr(ncut),lendef(le-ls+1)
   integer, allocatable  ::  lentmp(:)
!
   len=le-ls+1
   allocate(lentmp(len))

   do i = 1,ncut
     lenarr(i)=0
   enddo
   if( ind.eq.1 ) then
     lens=ls
     lene=le
   else
     lens=le
     lene=ls
   endif
   i=1
   n=1
   do nn = lens,lene,ind
     lenarr(n)=lenarr(n)+1
     lentmp(nn)=n
     n=n+i
     if(n.eq.ncut+1) then
       i=-1
       n=n+i
     endif
     if(n.eq.0) then
       i=1
       n=n+i
     endif
   enddo
!
   n=0
   do i = 1,ncut
     do nn = lens,lene,ind
       if(lentmp(nn).eq.i) then
         n=n+1
         lendef(n)=nn
       endif
     enddo
   enddo
!
   deallocate(lentmp)
!
   return
   end subroutine equdis
!-------------------------------------------------------------------------------
