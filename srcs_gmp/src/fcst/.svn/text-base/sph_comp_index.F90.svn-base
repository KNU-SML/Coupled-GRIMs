#include <define.h>
   subroutine sph_comp_index
!-------------------------------------------------------------------------------
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only : mwave  => jcap_                                      ,&
                         mwavep => jcap1_
   use comfcst, only   : indxnn, indxmm
!-------------------------------------------------------------------------------
   l=0
   do m = 1,mwavep
     nend=mwavep-m+1
     do nn = 1,nend
       n=nn+m-1
       l=l+2
       indx=(mwavep*(n-m)-(n-m)*(n-m-1)/2+m)*2-1
       indxnn(l-1)=indx
       indxnn(l  )=indx+1
     enddo
   enddo
!                                                                               
   l=0
   do nn = 1,mwavep                                                         
     lln=mwavep-nn+1
     do ll = 1,lln
       n=ll+nn-1
       m=ll
       indx=(m*mwavep-(mwavep-n)-(m-1)*m/2)*2-1
       l=l+2
       indxmm(l-1)=indx
       indxmm(l  )=indx+1
     enddo
   enddo
!
   return                                                                    
   end subroutine sph_comp_index
