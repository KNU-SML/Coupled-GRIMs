#include <machine.h>
#include <define.h>
   module comfver
!-------------------------------------------------------------------------------
   use paramodel, only   : levs_,levm1_,levp1_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   private              :: levs_,levm1_,levp1_
!
#ifndef DFS
   real, allocatable  ::  am(:,:),hm(:,:),tm(:,:)                             ,&
                          bm(:,:),cm(:,:),spdmax(:)                           ,&
                          si(:),sl(:),del(:),tov(:)
#endif
   real, allocatable  ::  rdel2(:),rmsdot(:)                                  ,&
                          ci(:),cl(:),gv(:),sv(:),rpi(:)                      ,&
                          p1(:),p2(:),h1(:),h2(:),rpirec(:)
   real               ::  thour,deltim,sl1,z00,fhour,shour,dtcvav,dtpost      ,&
                          filta,filtb,dk,tk,percut,dtswav,dtlwav,avprs0       ,&
                          cowave,dtwave
   logical            ::  stepone
   integer            ::  kdt,jdt,inistp,limlow,maxstp,numsum,nummax
#ifdef HYBRID
   real, allocatable  ::  ak5(:),bk5(:)
#endif
#ifdef DYNAMIC_ALLOC
   integer            ::  ncldb1,ncpus,ncpus1
#endif
!
   contains
!-------------------------------------------------------------------------------
   subroutine comfver_init
!-------------------------------------------------------------------------------
#ifndef DFS
   allocate(am(levs_,levs_),hm(levs_,levs_),tm(levs_,levs_)                   ,&
            bm(levs_,levs_),cm(levs_,levs_),spdmax(levs_)                     ,&
            si(levp1_),sl(levs_),del(levs_),tov(levs_)                         )
#endif
   allocate(rdel2(levs_),rmsdot(levm1_)                                       ,&
            ci(levp1_),cl(levs_),gv(levs_),sv(levs_),rpi(levm1_)              ,&
            p1(levs_),p2(levs_),h1(levs_),h2(levs_),rpirec(levm1_)             )
#ifdef HYBRID
   allocate(ak5(levp1_),bk5(levp1_))
#endif
!
   end subroutine comfver_init
!-------------------------------------------------------------------------------
   end module comfver
