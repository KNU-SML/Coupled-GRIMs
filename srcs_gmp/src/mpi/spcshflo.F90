!
   subroutine spcshflo(a,lnt,ntotal,jcap,lwvdef)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    spcshflo
!            
! abstract: shafle spectral coefficients from balancing for output
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call spcshflo(a,lnt,ntotal,jcap,lwvdef)
!
!    input argument lists:
!   a   - real (lnt,ntotal) total spectral field 
!   lnt   - integer total spectral grid
!   jcap   - integer total wavenumber of truncation
!   ntotal   - integer total set of fields
!   lwvdef   - integer (jcap+1) distribution of the index
!
!    output argument list:
!   a   - real (lnt,ntotal) total field
! 
! subprograms called:
!
!-------------------------------------------------------------------------------
   integer         ::  lnt,lntpp,ntotal,jcap,offset,l,n,m,k,j,ll
   real            ::  a(lnt,ntotal)
   integer         ::  lwvdef(jcap+1)
!
   offset(n,l)=(jcap+1)*(jcap+2)-(jcap-l+1)*(jcap-l+2)+2*(n-l)
!
   real,allocatable::  tmp(:,:)
!
   allocate(tmp(lnt,ntotal))
!
   m=0
   do ll = 1,jcap+1
     l=lwvdef(ll)
     do n = l,jcap
       m=m+1
       j=offset(n,l)
       do k = 1,ntotal
         tmp(j+1,k)=a(2*m-1,k)
         tmp(j+2,k)=a(2*m  ,k)
       enddo
     enddo
   enddo
!
   do m = 1,lnt
     do k = 1,ntotal
       a(m,k)=tmp(m,k)
     enddo
   enddo
!
   deallocate(tmp)
!
   return
   end subroutine spcshflo
