#include <define.h>
   module comfspec_vr
!-------------------------------------------------------------------------------
   use paramodel, only  : levs_
#ifndef DFS
   use paramodel, only  : LNT22S,levh_,lnuv_,LEVSS,LEVHS
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                       , dimension(4)    ::  idate
   real   , allocatable          , dimension(:)    ::  relvor                 ,&
                                                       absvor
#ifndef DFS
!
! global spectral coefficients
!
   real   , allocatable          , dimension(:)    ::  eps,epsi,gz,qm
!
! n-1 time step
!
   real   , allocatable, target  , dimension(:,:)  ::  zem
   real                , pointer , dimension(:,:)  ::  dim,tem
!
! n time step
!
   real   , allocatable, target  , dimension(:,:)  ::  ze
   real                , pointer , dimension(:,:)  ::  di,te,uln,vln
   real                , pointer , dimension(:)    ::  dpdphi,dpdlam,q,qlap
!
! n+1 time step
!
   real   , allocatable, target  , dimension(:,:)  ::  z
   real                , pointer , dimension(:,:)  ::  uu,vv,y,x,w
#ifdef NISLQ
   real   , allocatable          , dimension(:,:)  ::  rm,rq,rt
#else
   real                , pointer , dimension(:,:)  ::  rm,rq,rt
#endif
#endif /* ~DFS end */
!
   contains
!-------------------------------------------------------------------------------
   subroutine comfspec_vr_init
!-------------------------------------------------------------------------------
   use comkindex
!-------------------------------------------------------------------------------
   allocate(   relvor   (levs_)                                               ,&
               absvor   (levs_)                                                )
#ifndef DFS
   allocate(   eps      (lnuv_)                                               ,&
               epsi     (lnuv_)                                               ,&
               gz       (LNT22S)                                              ,&
               qm       (LNT22S)                                               )
!
! K level index
!
#ifndef NISLQ
   lots  =5*levs_+levh_+4   ; lotss  =5*LEVSS+LEVHS+4
#else
   lots  =5*levs_+4         ; lotss  =5*LEVSS+4
#endif
   lotst =2*levs_+1         ; lotsts =2*LEVSS+1
   ksz   =1                 ; kszs   =1 
   ksd   =ksz+levs_         ; ksds   =kszs+LEVSS
   kst   =ksd+levs_         ; ksts   =ksds+LEVSS
#ifndef NISLQ
   ksr   =kst+levs_         ; ksrs   =ksts+LEVSS
   ksu   =ksr+levh_         ; ksus   =ksrs+LEVHS
#else
   ksu   =kst+levs_         ; ksus   =ksts+LEVSS
#endif
   ksv   =ksu+levs_         ; ksvs   =ksus+LEVSS
   kspphi=ksv+levs_         ; kspphis=ksvs+LEVSS
   ksplam=kspphi+1          ; ksplams=kspphis+1
   ksp   =ksplam+1          ; ksps   =ksplams+1
   ksplap=ksp   +1          ; ksplaps=ksps   +1
!
#ifndef NISLQ
   lott=5*levs_+levh_+1     ; lotts  =5*LEVSS+LEVHS+1
#else
   lott=5*levs_+1           ; lotts  =5*LEVSS+1
#endif
   kz=1                     ; kzs    =1
   ku=kz+1                  ; kus    =kzs+1
   kv=ku+levs_              ; kvs    =kus+LEVSS
   ky=kv+levs_              ; kys    =kvs+LEVSS
#ifndef NISLQ
   kr=ky+levs_              ; krs    =kys+LEVSS
   kx=kr+levh_              ; kxs    =krs+LEVHS
#else
   kx=ky+levs_              ; kxs    =kys+LEVSS
#endif
   kw=kx+levs_              ; kws    =kxs+LEVSS
!
! n-1 time step
!
#ifndef NISLQ
   allocate(   zem            (1:LNT22S,ksz:3*levs_+levh_)                     )
#else
   allocate(   zem            (1:LNT22S,ksz:3*levs_    )                       )
#endif
               dim    =>zem   (1:LNT22S,ksd:ksd+levs_-1)
               tem    =>zem   (1:LNT22S,kst:kst+levs_-1)
#ifndef NISLQ
               rm     =>zem   (1:LNT22S,ksr:ksr+levh_-1)
#endif
!
! n time step
!
   allocate(   ze             (1:LNT22S,ksz:lots       )                       )
               di     =>ze    (1:LNT22S,ksd:ksd+levs_-1)
               te     =>ze    (1:LNT22S,kst:kst+levs_-1)
#ifndef NISLQ
               rq     =>ze    (1:LNT22S,ksr:ksr+levh_-1)
#endif
               uln    =>ze    (1:LNT22S,ksu:ksu+levs_-1)
               vln    =>ze    (1:LNT22S,ksv:ksv+levs_-1)
               dpdphi =>ze    (1:LNT22S,kspphi)
               dpdlam =>ze    (1:LNT22S,ksplam)
               q      =>ze    (1:LNT22S,ksp   )
               qlap   =>ze    (1:LNT22S,ksplap)
!
! n+1 time step
!
   allocate(   z              (1:LNT22S,kz :lott       )                       )
               uu     =>z     (1:LNT22S,ku :ku+levs_-1 )
               vv     =>z     (1:LNT22S,kv :kv+levs_-1 )
               y      =>z     (1:LNT22S,ky :ky+levs_-1 )
#ifndef NISLQ
               rt     =>z     (1:LNT22S,kr :kr+levh_-1 )
#endif
               x      =>z     (1:LNT22S,kx :kx+levs_-1 )
               w      =>z     (1:LNT22S,kw :kw+levs_-1 )
#ifdef NISLQ
!
! for nislq
!
   allocate(   rm             (1:LNT22S,1:levh_)                              ,&
               rq             (1:LNT22S,1:levh_)                              ,&
               rt             (1:LNT22S,1:levh_)                               )
#endif
#endif /* ~DFS end */
   end subroutine comfspec_vr_init
!-------------------------------------------------------------------------------
   end module comfspec_vr
