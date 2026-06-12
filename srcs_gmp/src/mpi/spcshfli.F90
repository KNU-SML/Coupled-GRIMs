!
   subroutine spcshfli(a,lnt,ntotal,jcap,lwvdef)
!-------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    spcshfli
!            
! abstract: shafle spectral coefficients for balancing from input
!
! program history log:
!    99-06-27  henry juang    finish entire test for gsm
!
! usage:   call spcshfli(a,lnt,ntotal,jcap,lwvdef)
!
!    input argument lists:
!   a   - real (lnt,ntotal) total spectral field 
!   lnt   - integer total spectral grid
!   jcap   - integer total wavenumber of truncation
!   ntotal   - integer total set of fields
!   lwvdef   - integer (jcap+1) distribution index
!
!    output argument list:
!   a   - real (lnt,ntotal) total field
! 
! subprograms called:
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer     ::  lnt,lntpp,ntotal,jcap,                                      &
                   offset,l,n,m,k,j,ll
   real        ::  a(lnt,ntotal)
   integer     ::  lwvdef(jcap+1)
!
   offset(n,l)=(jcap+1)*(jcap+2)-(jcap-l+1)*(jcap-l+2)+2*(n-l)
!
   real,allocatable::tmp(:,:)
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
         tmp(2*m-1,k)=a(j+1,k)
         tmp(2*m  ,k)=a(j+2,k)
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
   end subroutine spcshfli
!-------------------------------------------------------------------------------
