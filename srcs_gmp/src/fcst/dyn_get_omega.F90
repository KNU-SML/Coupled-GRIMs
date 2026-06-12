#include <define.h>
   subroutine dyn_get_omega(ims2,imx2,kmx,dphi,dlam,cg,ug,vg,                  &
                       dg,del,rcl,vvel,prsi1,prsl)
!-------------------------------------------------------------------------------
!
! abstract : code lifted from post (mcp1840) jun 88--computes vvel (cb/sec)
!            input ps in cb,output vvel in cb/sec
!            do loops altered for better vectorization possibilities.
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only : ILOTS=>LONF2S,levs_
!-------------------------------------------------------------------------------
   integer              ::  ims2,imx2,kmx
   real                 ::  dphi(imx2),dlam(imx2)
   real                 ::  cg(imx2,kmx),ug(imx2,kmx),vg(imx2,kmx)
   real                 ::  dg(imx2,kmx),del(imx2,kmx)
   real                 ::  prsi1(imx2),prsl(imx2,kmx)
!
! vvel contains omega in layers on return from subroutine...
!
   real                 ::  vvel(imx2,kmx)
!
! local array
!
   real                 ::  db(ILOTS,levs_)
   real                 ::  cb(ILOTS,levs_)
   real                 ::  dot(ILOTS,levs_+1)
!
!-------------------------------------------------------------------------------
   nx =ims2
   nz =kmx
   do k = 1,nz+1
     do lo = 1,nx
       dot(lo,k) = 0.e0
     enddo
   enddo
!
!...  compute c=v(true)*del(ln(ps)).divide by cos for del cos for v
!
   do lo = 1,nx
     dphi(lo)=dphi(lo)*rcl
     dlam(lo)=dlam(lo)*rcl
   enddo
   do le = 1,nz
     do lo = 1,nx
       cg(lo,le)=ug(lo,le)*dlam(lo)+vg(lo,le)*dphi(lo)
     enddo
   enddo
   do lo = 1,nx
     db(lo,1)=del(lo,1)/prsi1(lo)*dg(lo,1)
     cb(lo,1)=del(lo,1)/prsi1(lo)*cg(lo,1)
   enddo
   do le = 1,nz-1
     do lo=1,nx
       db(lo,le+1)=db(lo,le)+del(lo,le+1)/prsi1(lo)*dg(lo,le+1)
       cb(lo,le+1)=cb(lo,le)+del(lo,le+1)/prsi1(lo)*cg(lo,le+1)
     enddo
   enddo
!
!...    sigma dot computed only at interior interfaces
!
   do k = 1,nz-1
     do lo = 1,nx
       dot(lo,k+1)=dot(lo,k)+del(lo,k)/prsi1(lo)*                              &
                              (db(lo,nz)+cb(lo,nz)-dg(lo,k)-cg(lo,k))
     enddo
   enddo
   do k = 1,nz
     do lo = 1,nx
       vvel(lo,k)=prsl(lo,k)/prsi1(lo)*(cg(lo,k)-cb(lo,nz)-db(lo,nz))-         &
                 0.5*(dot(lo,k+1)+dot(lo,k))
       vvel(lo,k)=vvel(lo,k)*prsi1(lo)
!
!     vvel(lo,k)=vvel(lo,k)*prsi1(lo)*10.
!
     enddo
   enddo
!
   return
   end subroutine dyn_get_omega
