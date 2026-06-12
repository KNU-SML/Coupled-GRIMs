!
   subroutine mpadjlm(par,lnp,levs,tot,lmp,lng,levsp,nvar)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! add aditional meridional wave components to tot after JLG index
!
! usage:        call mpadjl(par,lnp,levs,tot,lmp,lng,levsp,nvar)
!
!    input argument lists:
!       par     - real (lmp,lnp,ntotal) sub partial field
!       lmp     - integer total spectral grid
!       lnp     - integer total spectral grid
!       ntotal  - integer total set of fields
!
!    output argument list:
!       tot       - real (lmg,lng,ntotal) total field
!
!-------------------------------------------------------------------------------
   use dfsvar, only  :  jlg,kls,kle
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  lmp,lng,lnp,nvar,levsp,levs
   real     ::  par(lmp,lnp,nvar*levs),tot(lmp,lng,nvar*levsp)
   integer  ::  kd,ks,n,k,nv
!
   if (levs.eq.levsp) then
     do k = 1,nvar*levsp
       do n = 1,lnp
         tot(1:lmp,jlg+n,k)=par(1:lmp,n,k)
       enddo
     enddo
   else
     do nv = 1,nvar
       do k = kls,kle
         kd=(nv-1)*levsp+k-kls+1
         ks=(nv-1)*levs+k
         do n = 1,lnp
           tot(1:lmp,jlg+n,kd)=par(1:lmp,n,ks)
         enddo
       enddo
     enddo
   endif
!
   return
   end subroutine mpadjlm
!-------------------------------------------------------------------------------
