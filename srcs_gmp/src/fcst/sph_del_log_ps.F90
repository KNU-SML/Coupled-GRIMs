#include <define.h>
   subroutine sph_del_log_ps(q,dpdphs,dpdtop,dpdla,llstr,llens,lwvdef)
!-------------------------------------------------------------------------------
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only   :  jcap_,jcap1_
#ifndef RMP
#ifdef MP
   use paramodel, only   :  LLN2S,LLN22S
#else
   use paramodel, only   :  LLN2S=>lnt2_,LLN22S=>lnt22_
#endif
#else
   use paramodel, only   :  LLN2S=>lnt2_,LLN22S=>lnt22_
#endif
   use constant, only    :  rerth_
!-------------------------------------------------------------------------------
   integer              ::  lnep
!
! input q is in ibm triang. order
! output  is in ibm triang. order
!
   integer,intent(in)   ::  llstr,llens
   real,intent(in)      ::  q(LLN22S)
   real,intent(out)     ::  dpdphs(LLN22S)
   real,intent(out)     ::  dpdtop(2,jcap1_)
   real,intent(out)     ::  dpdla(LLN22S)
   integer,intent(in)   ::  lwvdef(jcap1_)
!
! local
!
   real,save,allocatable::  e(:)
   integer,save         ::  ifirst
!-------------------------------------------------------------------------------
!
! function
!
   je(n,l) =((jcap_+2)*(jcap_+3)-(jcap_+2-l)*(jcap_+3-l))/2+n-l
!
   jc(n,l) = (jcap_+1)*(jcap_+2)-(jcap_+1-l)*(jcap_+2-l)+2*(n-l)
!
   data ifirst/1/
!
! --------------------- initial setting---------------
!
   lnep=(jcap_+2)*(jcap_+3)/2
   if(.not.allocated(e)) allocate(e(lnep))
!
   if(ifirst.eq.1) then
     do l = 0,jcap_
       n=l
       ie=je(n,l)+1
       e(ie)=0.
     enddo
     do l = 0,jcap_
       do n=l+1,jcap_+1
         rn=n
         rl=l
         a=(rn*rn-rl*rl)/(4.*rn*rn-1.)
         ie=je(n,l)+1
         e(ie)=sqrt(a)
       enddo
     enddo
     ifirst=0
   endif
!
! ------------------ end of initial setting -----------
!
   aa=1./rerth_
!
   lntx=0
!
   do ll = 1,llens
     l=lwvdef(llstr+ll)
     lnt0=jc(l,l)-lntx
     lntx=lntx+jc(l+1,l+1)-jc(l,l)
!
! ----- l = 0, jcap_
!
     rl=l
     do n = l,jcap_
       icr=jc(n,l)+1 - lnt0
       ici=jc(n,l)+2 - lnt0
       dpdla(ici)= rl*q(icr)
       dpdla(icr)=-rl*q(ici)
     enddo
!
! ------ l = 0, jcap_-1
!
     if( l.le.jcap_-1 ) then
       do n = l,jcap_-1
         ie=je(n+1,l)+1
         icr=jc(n,l)+1 - lnt0
         ici=jc(n,l)+2 - lnt0
         rn=n
         dpdphs(icr)=(rn+2.)*e(ie)*q(icr+2)
         dpdphs(ici)=(rn+2.)*e(ie)*q(ici+2)
       enddo
     endif
!
! ------ l = 0, jcap_
!
     n=jcap_
     icr=jc(n,l)+1 - lnt0
     ici=jc(n,l)+2 - lnt0
     dpdphs(icr)=0.0
     dpdphs(ici)=0.0
!
! ------ l = 0, jcap_-1
!
     if( l.le.jcap_-1 ) then
       do n = l+1,jcap_
         ie=je(n,l)+1
         icr=jc(n,l)+1 - lnt0
         ici=jc(n,l)+2 - lnt0
         rn=n
         dpdphs(icr)=dpdphs(icr)+(1.-rn)*e(ie)*q(icr-2)
         dpdphs(ici)=dpdphs(ici)+(1.-rn)*e(ie)*q(ici-2)
       enddo
     endif
!
! ------ l = 0, jcap_
!
     n=jcap_+1
     rn=n
     ie=je(n,l)+1
     icr=jc(n,l)+1 - lnt0
     ici=jc(n,l)+2 - lnt0
     dpdtop(1,ll)=(1.-rn)*e(ie)*q(icr-2)
     dpdtop(2,ll)=(1.-rn)*e(ie)*q(ici-2)
     dpdtop(1,ll)=dpdtop(1,ll)*aa
     dpdtop(2,ll)=dpdtop(2,ll)*aa
   enddo
!
   do j = 1,LLN2S
     dpdla(j)= dpdla(j)*aa
     dpdphs(j)=dpdphs(j)*aa
   enddo
!
   return
   end subroutine sph_del_log_ps
