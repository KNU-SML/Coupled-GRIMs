!
   subroutine equlw12(jcap,nrow,lwvdis,lnpdis,jcaprm)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    equlw12
!            
! abstract: about equal division of all triangular truncation waves
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   equlw12(jcap,nrow,lwvdis,lnpdis,jcaprm)
!
!    input argument lists:
!   jcap   - integer spectral wavenumber
!   nrow   - integer number of nrow
!
!    output argument list:
!   lwvdis   - integer (nrow) length of wave in l
!   lnpdis   - integer (nrow) length of wave for l and n
!   jcaprm   - integer remain of unbalance in l
! 
! subprograms called:
!   equdiv   - to compute about equal number of subgroup 
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  jcap,nrow,jcaprm,                                  &
                            jcaphf,nr,n,l,m,lw1pnt,lw2pnt
   integer              ::  lwvdis(nrow),lnpdis(nrow)
   integer, allocatable ::  lw1str(:),lw1end(:)
   integer, allocatable ::  lw2str(:),lw2end(:)
!
   allocate (lw1str(nrow))
   allocate (lw1end(nrow))
   allocate (lw2str(nrow))
   allocate (lw2end(nrow))
!
   jcaphf=(jcap+1)/2
   jcaprm=mod(jcap+1,2)
   call equdiv(jcaphf,nrow,lwvdis)
!
   lw1pnt=1
   lw2pnt=jcap+2
!
   do n = 1,nrow
     lw2end(n)=lw2pnt-1
     lw2pnt=lw2pnt-lwvdis(n)
     lw2str(n)=lw2pnt
     lw1str(n)=lw1pnt
     lw1pnt=lw1pnt+lwvdis(n)
     lw1end(n)=lw1pnt-1
   enddo
   lw1end(nrow)=lw1end(nrow)+jcaprm
!
   do nr = 1,nrow
     m=0
     do l = lw1str(nr),lw1end(nr)
       do n = l,jcap+1
         m=m+1
       enddo
     enddo
     do l = lw2str(nr),lw2end(nr)
       do n = l,jcap+1
         m=m+1
       enddo
     enddo
     lnpdis(nr)=m
   enddo
!
   deallocate (lw1str)
   deallocate (lw1end)
   deallocate (lw2str)
   deallocate (lw2end)
!
   return
   end subroutine equlw12
!-------------------------------------------------------------------------------
