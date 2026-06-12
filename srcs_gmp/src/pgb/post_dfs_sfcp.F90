!
   subroutine post_dfs_sfcp(sf1,psl1,mt,jl,topo,xtop,ytop,                     &
                            ps,xps,yps,ib,jbw,coslat)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  mt,jl,ib,jbw
   real     ::  psl1(-mt:mt,0:jl)
   real     ::  sf1 (-mt:mt,0:jl)
   real     ::  topo(ib,jbw),ps  (ib,jbw)
   real     ::  xtop(ib,jbw),ytop(ib,jbw)
   real     ::  xps (ib,jbw),yps (ib,jbw)
   real     ::  coslat(jbw,3)
!
! local variables 
!
   real     ::  SANMV(-MT:MT,0:(JL+1)),sanm1(-mt:mt,0:jl)
   integer  ::  i,j,jlg,mtg
!-------------------------------------------------------------------------------
   jlg=jl+1
   mtg=2*mt+1
   sanmv(:,0)=0.
   sanm1(:,0)=0.
!
   CALL dfs_psi2wind( PSL1, SANMV, sanm1 ,mtg, JL ,1)
   CALL dfs_fft_driver( -1,  ps,ib,jbw,1,psl1 ,mtg,jlg  ,1,                    &
                        jlg,1,1,0,coslat,+1 )
   CALL dfs_fft_driver( -1, xps,ib,jbw,1,SANM1 ,mtg,jlg  ,1,                   &
                        jlg,1,1,0,coslat,+1 )
   CALL dfs_fft_driver( -1, yps,ib,jbw,1,SANMV ,mtg,jlg+1,1,                   &
                        jlg+1,1,1,0,coslat,+1 )
!
   CALL dfs_psi2wind( SF1 , SANMV, sanm1 ,mtg, JL ,1)
   CALL dfs_fft_driver( -1,topo,ib,jbw,1, sf1 ,mtg,jlg  ,1,                    &
                        jlg,1,1,0,coslat,+1 )
   CALL dfs_fft_driver( -1,xtop,ib,jbw,1,SANM1 ,mtg,jlg  ,1,                   &
                        jlg,1,1,0,coslat,+1 )
   CALL dfs_fft_driver( -1,ytop,ib,jbw,1,SANMV ,mtg,jlg+1,1,                   &
                        jlg+1,1,1,0,coslat,+1 )
!
   do j = 1,jbw
     do i = 1,ib
       yps (i,j)=-yps (i,j)
       ytop(i,j)=-ytop(i,j)
       ps(i,j)=exp(ps(i,j))
     enddo
   enddo
!     
   return
   end subroutine post_dfs_sfcp
