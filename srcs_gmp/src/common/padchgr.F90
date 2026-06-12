!
   module padchgr
!-------------------------------------------------------------------------------
   use paramter, only : mwave,idim,kdim
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ngwri                                             ,&
                            ntrns                                             ,&
                            nfl22                                             ,&
                            npln3                                             ,&
                            nsums                                             ,&
                            nsyms                                             ,&
                            ndjm1                                             ,&
                            nscrch
!
   contains
!-------------------------------------------------------------------------------
   subroutine padchgr_init
!-------------------------------------------------------------------------------
   ngwri=2*((mwave+1)*(mwave+2)) 
   ntrns=(mwave+1)*(mwave+2) 
   nfl22=(mwave+1)*2*1024+(mwave+1)*(mwave+2) 
   npln3=(mwave+1)*(mwave+2)*2+(mwave+1)*(mwave+1)*2+ &
                   (mwave+1)*2*2 + (((mwave+1)*(mwave+2))/2)*2*2 
   nsums=(mwave+1)*2*(mwave+1) 
   nsyms=idim*2*kdim 
   ndjm1=(nfl22/npln3*nfl22+npln3/nfl22*npln3) / (nfl22/npln3+npln3/nfl22)
   nscrch=ndjm1+1
!
   end subroutine padchgr_init
!-------------------------------------------------------------------------------
   end module padchgr
