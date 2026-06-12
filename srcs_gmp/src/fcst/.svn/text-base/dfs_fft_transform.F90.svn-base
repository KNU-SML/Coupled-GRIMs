#include <define.h>
   subroutine dfs_fft_transform
!-------------------------------------------------------------------------------
!
!  ::: structure ::: This file contains ...
!
!      [dfs_fft_transform]
!           |
!           |-- [dfs_fft_driver] *
!           |-- [dfs_fft_setup] *
!           |-- [dfs_sine_transform] *
!           |-- [dfs_pole_cut] *
!           |-- [dfs_fft_2] *
!           |-- [dfs_fft_3] *
!           |-- [dfs_fft_5] *
!           |-- [dfs_fft_all] *
!
! program history log:
!   2000-03-01  hyeong-bin cheong development
!   2006-03-01  hoon park         implementation
!   2008-03-01  hoon park         mpi development
!   2010-03-01  myung-seo koo     dimension allocatable
!   2011-03-01  jung-eun kim      grims structure
!
!-------------------------------------------------------------------------------
   end subroutine dfs_fft_transform
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_fft_driver(isgn,PSI,ib,jbw,lotg,ANM,mt,jt,lota,              &
                             jtg,levsp,levs,nvar,coslat,ipole)
!-------------------------------------------------------------------------------
!
! abstract: Double Fourier Transform on the Sphere : vectorized
!
!     isgn           = 1 : grid to wave
!                    =-1 : wave to grid
!     PSI(ib,jbw,lotg)   : grid data
!     ANM(MT,jtg,lota)   : wave
!     ipole          = 1 : pole cut
!                    = 0 : no pole cut
!     JT                 : local num. of meridional wave
!     JTG                : global num. of meridional wave
!                          For MPI JTG should equal to JLG
!                          For single JTG can be JLG or JLG+1
!     NVAR               : number of multi-level variables
!     PSI                :  (local i, local j, global k) grid
!     ANM                :  (local m, local/global n, global/local k)
!
! program history log:
!   2000-01-01  heongbin cheong        initial development
!   2005-01-01  hoon park              grims implementation
!   2006-04-01  hoon park              mpi
!
!-------------------------------------------------------------------------------
   use dfsvar, only : igs,jgs,ige,jge,iba,jbwa,jbw1a,jba,jls,                  &
                        mls,mle,mgs,mge,mth,mta,mtha,mtg,ngs,nge,kls,jlg,      &
#ifdef ALIASED
                        nls,nle,                                               &
#endif
                        iope,latdef
#ifdef MP
   use commpi, only : mype
#endif
   use constant, only : pi_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                                    ::  isgn,ib,jbw,mt,JT,jtg,nvar
   integer                                    ::  ipole,lotg,lota,levsp,levs
   real                                       ::  PSI(ib,jbw,lotg)
   real, target                               ::  ANM(MT,jt,lota)
                                                  ! cos & 1/cos
   real                                       ::  coslat(jgs:jge,3)
!
! local variables
!
   real                                       ::  data_SIN(MT,jbwa)
   real                                       ::  ALON(0:jbw/2,iba)
   real                                       ::  BLON(0:jbw/2,iba)
!
!  to avoid Bank Conflict
!
   real                                       ::  COSNM(jtg)
   real                                       ::  COSJ(jbwa)

   real,allocatable                           :: FWAVEJ(:,:,:)
#ifdef MP
   real,allocatable,dimension(:,:,:)          ::  psia,FWAVE2
   real,allocatable,dimension(:,:,:),target   ::  anmt
   real,pointer    ,dimension(:,:,:)          ::  anma
#else
   integer,parameter                          ::  mype=0
#endif
   integer                                    ::  jb,k,j,mx,i,mj,ibmx,j2n,MLAST
   integer                                    ::  mjJB,mjJB1,je,ii,kd,ks
   integer                                    ::  n1lv,n1stg,n1stl,lotl,n,kk
   integer                                    ::  mzr,menp,mm,mp,m
   integer                                    ::  jg,mst,jst,ierr,jjdim,jj
   real                                       ::  bb,sum,aa,pin,rib,rinv,winv
!
! save variables
!
   integer,allocatable,save                   ::  jhs(:),jhn(:)
   integer,save                               ::  jend,jodd
   CHARACTER(len=4),save                      ::   FFTOPT
   real,allocatable,save                      ::  TRIGSIB(:),TRIGSJBW(:)
!
! Extra array space
!
   real,allocatable,save                      ::  COSFUN1(:,:),COSFUN2(:,:)
!
! mx= all
!
   real,allocatable,save                      ::  SINHF(:),COSPK(:),SINPK(:)
   LOGICAL                                    ::  JUSDOIT
   data      JUSDOIT/.true./
!-------------------------------------------------------------------------------
#include "abort.h"
!
! initialize
!
   alon=0. ; blon=0. ; cosnm=0. ; cosj=0.
!
   jb=jbw/2
   rib=1./iba
   MLAST= MTA                         ! A part transform for Y-dir
   if(MTA.eq.(IBA/2)) MLAST= MTA-1       ! Full range transform for Y-dir
!
! single level variables
!
   n1lv=lotg-nvar*levs
   lotl=levsp*nvar+n1lv
   n1stg=nvar*levs+1
   n1stl=nvar*levsp+1
#ifdef MP
   mzr=0
   menp=0
   if (mls.eq.0) mzr=mth+1
!
   if (mle.eq.MTA) then
     menp=mt
   endif
!
#else
   mzr=mta+1
   menp=mt
#endif
   IF(JUSDOIT) THEN
#ifdef DBG
     if(iope) then
       write(6,*)'in FFT :ipole,jtg=',ipole,jtg
       write(6,*)'mzr,menp,n1lv,lotg,nvar,lotl,n1stl,mt,mls=',                 &
                  mzr,menp,n1lv,lotg,nvar,lotl,n1stl,mt,mls
     endif
#endif
     JUSDOIT=.FALSE.
     allocate(TRIGSIB(iba*2),TRIGSJBW(jbwa*2))
     allocate(COSFUN1(jgs:jge,0:jbwa),COSFUN2(0:jbwa,jgs:jge))
     allocate(SINHF(0:jba),COSPK(0:jba),SINPK(0:jba))

     FFTOPT='P235'
     J2N= 1
     do k = 1,20
       J2N= J2N*2
       if(IBA.eq.J2N) then
!        FFTOPT='P222'
         EXIT
       endif
     end do        ! 2004.04.02      
#ifdef DBG
     if (iope) then
       write(6,*) ' IBA, FFTSOPT=', IBA, trim(FFTOPT)
     endif
#endif
     CALL dfs_fft_setup(TRIGSIB ,IBA )
     CALL dfs_fft_setup(TRIGSJBW,JBWA)  ! for 2,3,5 prime factor 
!
     PIN= PI_/JBWA
!
     do j = 0,JBA                         ! For SINE and COSINE transform 
       aa=PIN*(j+0.5D-00)
#ifdef IBMSP
       SINHF(j)= SIN(dble(aa))
#else
       SINHF(j)= DSIN(dble(aa))
#endif
     end do
!
     do k = 0,JBA
       bb=PIN*dble(k)
#ifdef IBMSP
       COSPK(k)= COS(dble(bb))
       SINPK(k)= SIN(dble(bb))
#else
       COSPK(k)= DCOS(dble(bb))
       SINPK(k)= DSIN(dble(bb))
#endif
     end do
!
     do j = jgs,jge
       sum= (j-jgs+0.5d0)*PIN          ! For direct COS transform
       do mj = 0,jbw1a
#ifdef IBMSP
         COSFUN1(j,mj)= COS(dble(mj*sum))
         COSFUN2(mj,j)= COS(dble(mj*sum))
#else
         COSFUN1(j,mj)= DCOS(dble(mj*sum))
         COSFUN2(mj,j)= DCOS(dble(mj*sum))
#endif
       end do
     end do
!
     jend=(jbw-1)/2
     allocate(jhs(0:jend))
     allocate(jhn(0:jend))
     if (mod(jbw,2).eq.0) then
       jodd=0
       Do j = 0,jend
         jhs(j)=jend+1-j
         jhn(j)=jend+2+j
#ifdef DBG_DFS
         if(iope) then
           write(6,'(A,4I4)')'even j,jhs,jhn=',j,jhs(j),jhn(j),mype
         endif
#endif
       enddo
     else
       jodd=1
       Do j = 0,jend
         jhs(j)=jend+1-j
         jhn(j)=jend+1+j
#ifdef DBG_DFS
         if(iope) then
           write(6,'(A,4I4)')'odd j,jhs,jhn=',j,jhs(j),jhn(j),mype
         endif
#endif
       enddo
     endif
   ENDIF         ! first
!
!----------------------------------------------------------------------
!
   IF(JTG-1.gt.(IBA/2))  then
     if (iope) then
       write(6,*) ' ==   Invalid range of Y-transform   =='
       write(6,*) ' == JT must not be greater than IB/2 =='
       write(6,*) 'JTG-1,IBA/2=',jtg-1,IBA/2
     endif
     call MPABORT
   ENDIF
!
   PIN= PI_/JBWA
!=============================
   IF(isgn.EQ.+1) THEN
!=============================
!
!-----(I,J)->(mx,J)
!
#ifdef MP
     allocate(psia(iba,jbw,lotl))
     if (nvar.ne.0) then
       call mpxy2yz(psi,ib,lotg,psia,iba,lotl,jbw,levs,levsp,nvar)
     endif
     if (n1lv.ne.0) then
       call mpxy2y(psi(1,1,n1stg),ib,psia(1,1,n1stl),iba,jbw,n1lv)
     endif
#else
#define PSIA PSI
#endif
     allocate(FWAVEJ(mtg,jbw+1,lotl))
!
     do k = 1,lotl
       FWAVEJ(1:mtg,1:jbw+1,k)=0.0
       do I = 1,IBA
         Do J = 0,jend
           BLON(J,I)= PSIA(I,jhn(j),k)
           ALON(J,I)= PSIA(I,jhs(j),k)
         end do
         if (jodd.eq.1) then   
           BLON(0,I)= 0.d0
         endif
       END Do
!   
       CALL dfs_fft_all( ALON,BLON,TRIGSIB, 1,IBA,IBA,1,isgn,JB,jend)
!
       do mx = 2,mlast+1      ! -MTA,-MTA+1,..,-1,1,2,..,MTA
         IBmx=IBA-mx+2
         mm=mta+2-mx          ! -1,-2,..,-MTA
         mp=mta+mx            !  1, 2,.., MTA
         Do J = 0,jend
!
!          img part
!
           FWAVEJ(mm,jhn(j),k)= ( BLON(J,mx)+BLON(J,IBmx))*rib ! FWAVEJ(*,JB)= 0.
           FWAVEJ(mm,jhs(j),k)= ( ALON(J,mx)+ALON(J,IBmx))*rib
!
!          real part
!
           FWAVEJ(mp,jhn(j),k)= (-ALON(J,mx)+ALON(J,IBmx))*rib
           FWAVEJ(mp,jhs(j),k)= ( BLON(J,mx)-BLON(J,IBmx))*rib
         END Do ! par 
       end do
!
       mx=mta+1                     ! For m is 0
       Do J = 0,jend
         FWAVEJ(mx,jhn(j),k)=    BLON(J,1 )*rib
         FWAVEJ(mx,jhs(j),k)=    ALON(J,1 )*rib      ! always A is the last
         if(MTA.eq.(IBA/2)) then
           FWAVEJ(1,jhn(j),k)=    BLON(J,mtg)*rib
           FWAVEJ(1,jhs(j),k)=    ALON(J,mtg)*rib
         endif
       END Do ! par 
!
       Do J = 1,jb
         jg=latdef(jls+j-jgs)
         jj=jbw-j+1
         do mx = mta+3,mtg,2            ! 2, 4,.., MTA
           FWAVEJ( mx,J ,k)= FWAVEJ( mx,J ,k)*COSLAT(jg,2)
           FWAVEJ( mx,JJ,k)= FWAVEJ( mx,jj,k)*COSLAT(jg,2)
         enddo
         do mx = mta-1,1,-2            ! -2,-4,..,-MTA
           FWAVEJ( mx,J ,k)= FWAVEJ( mx,J ,k)*COSLAT(jg,2)
           FWAVEJ( mx,jj,k)= FWAVEJ( mx,jj,k)*COSLAT(jg,2)
         end do
       END Do ! par
     enddo      ! lotl loop
!
     IF(ipole.EQ.+1)                                                           &
       CALL dfs_pole_cut(FWAVEJ,mtg,jbw+1,coslat,lotl)

!-----(mx,J)->(mx,mj)
#ifdef MP
       deallocate(PSIA)
       allocate(FWAVE2(mt,jbwa,lotl))
       call mpyz2mz(FWAVEJ,mtg,jbw,jbw+1,FWAVE2,mt,jbwa,lotl)
       deallocate(FWAVEJ)
#define FWAVEJ FWAVE2
       if (jt.ne.jtg) then
         allocate(anmt(mt,jtg,lotl))
         anma=>anmt(1:mt,1:jtg,1:lotl)
         kk=kls
       else
         anma=>anm(1:MT,1:jt,1:lota)
         kk=1
       endif
#else
#undef PSIA
#endif /* MP end */
       rinv=1./JBA
       winv=1./JBWA
       do k = 1,lotl
         Do J = 1,jbwa
           do mx = 1,mt
             data_SIN(mx,J)= FWAVEJ(mx,j,k)
           end do
           if (mzr.ne.0) then
             data_SIN(mzr,J)= 0.
             COSJ(J)  = FWAVEJ(mzr,J,k)
           endif
         END Do
         CALL dfs_sine_transform(data_SIN, SINHF,COSPK,SINPK,                  &
                    TRIGSJBW, JBWA, JBWA-1, JBA, MT, isgn )
!
!------ This FT replaces COS_FTN caus'f poor vec efficiency
!
         if (mzr.ne.0) then
           do mj=1,jtg
             sum= 0.
#ifdef LINUX_PGI
!pgi$ novector
#endif
             do j=jgs,jge
               sum= sum + COSFUN1(j,mj-1)*COSJ(J-jgs+1)
             end do
             COSNM(mj)= sum*rinv
             if(mj.eq.1) COSNM(mj)= sum*winv
           end do
         endif
!
!------ This FT replaces COS_FTN caus'f poor vec efficiency
!
#ifdef MP
#define ANM ANMA
#endif
         Do mj = 2,jtg
           do mx = 1,mt
             ANM(mx,mj,k)= data_SIN(mx,mj-1)  ! In SUB, mj=0 -> SIN(1)
           END Do
         end do
!
         if (mzr.ne.0) then
           do mj=1,jtg
             ANM(mzr,mj,k)= COSNM(mj)
           end do
           ANM(1:mzr-1,1,k)=0.
           ANM(mzr+1:mt,1,k)=0.
           mj=jtg
           if(jtg-1.eq.(iba/2)) then
             ANM(mzr,jtg,k)= 0.
           endif
         else
           ANM(1:mt,1,k)=0.
         endif
!
!------
!
         if ((MTA.eq.(IBA/2)).and.(menp.ne.0)) then
           do mj=1,jtg
             ANM(menp,mj,k)= 0.
           end do
         endif
       enddo      ! lotl loop
#ifdef ALIASED
      !call dfs_cut_alias(anm,mt,0,jtg-1,lotl)
#endif
!
     deallocate(FWAVEJ)
#ifdef MP
#undef FWAVEJ
#undef ANM
     if (jt.ne.jtg)  then
       if (nvar.ne.0) then
         call mpmz2mn(ANMA,jtg,levsp,ANM,jt,mt,levs,nvar)
       endif
       if (n1lv.ne.0) then
         call mpm2mn(ANMA(1,1,n1stl),jtg,ANM(1,1,n1stg),jt,mt,n1lv)
       endif
       deallocate(anmt)
     endif
     nullify(anma)
#endif
!=============================
   ELSEIF(isgn.EQ.-1) THEN
!=============================
!
!-----(mx,mj)->(mx,J)
!
#ifdef MP
     if (jt.ne.jtg) then
       allocate(ANMT(mt,jtg,lotl))
       anma=>anmt(1:mt,1:jtg,1:lotl)
       if (nvar.ne.0) then
         call mpmn2mz(ANM,jt,levs,ANMA,jlg,mt,levsp,nvar)
       endif
       if (n1lv.ne.0) then
         call mpmn2m(ANM(1,1,n1stg),jt,ANMA(1,1,n1stl),jlg,mt,n1lv)
       endif
     else                  ! The input data was transposed to m,k direction
       ANMA=>ANM(1:MT,1:jt,1:lota)
     endif
#define ANM ANMA
#endif
     allocate(FWAVEJ(mt,jbwa,lotl),stat=ierr)
     if (ierr /= 0) then
       write(6,*)'allocate FWAVEJ error'
       call MPABORT
     endif
!
     do k = 1,lotl
       if (mzr.ne.0) then
         if((jtg-1).eq.(IBA/2)) ANM(mzr,jtg,k)= 0.
         ANM(1:mzr-1,1,k)=0.
         ANM(mzr+1:mt,1,k)=0.
       else
         ANM(1:mt,1,k)=0.
       endif
!
       if((MTA.eq.(IBA/2)).and.(menp.ne.0)) then
         do mj = 1,jtg
           ANM(menp,mj,k)= 0.
         end do
       endif
!
       Do mj = 1,jbwa
         data_SIN(1:mt,mj)= 0.
       END Do
!      
       Do mj = 2,jtg
         do mx = 1,mt
           data_SIN(mx,mj-1)= ANM(mx,mj,k)
         END Do
         if (mzr.ne.0) data_SIN(mzr,mj-1)= 0.
       end do
!
       if (mzr.ne.0) then
         Do mj = 1,jtg
           COSNM(mj) = ANM( mzr,mj,k)
         end do
       endif
!
       CALL dfs_sine_transform( data_SIN, SINHF, COSPK, SINPK,                 &
                     TRIGSJBW, JBWA, JBWA-1,JBA, MT, isgn )
!
!------ This FT replaces COS_FTN caus'f poor vec efficiency
!
       if (mzr.ne.0) then
         do j=jgs,jge
           sum= 0.
#ifdef LINUX_PGI
!pgi$ novector
#endif
           do mj=0,jtg-1
             sum= sum + COSFUN2(mj,j)*COSNM(mj+1)
           end do
           COSJ(j-jgs+1)= sum
         end do
       endif
!
       FWAVEJ(1:mt,1:jbwa,k)=0.0
!
       Do J=1,jbwa
         do mx=1,mt
           FWAVEJ(mx,J,k)= data_SIN(mx,J)
         end do
         if (mzr.ne.0) FWAVEJ(mzr,j,k)= COSJ(j)
       END Do
     enddo
!
!-----(mx,J)->(I,J)
!
#ifdef MP
#undef ANM
     nullify(anma)
     if (allocated(anmt)) deallocate(anmt)
     jjdim=jbw+1
     allocate(FWAVE2(mtg,jjdim,lotl))
     call mpmz2yz(FWAVEJ,mt,jbwa,FWAVE2,mtg,jjdim,lotl)
     deallocate(FWAVEJ)
     allocate(psia(iba,jbw,lotl),stat=ierr)
     if (ierr /= 0) then
       write(6,*)'allocate psia error'
       call MPABORT
     endif
#define FWAVEJ FWAVE2
#define PSI PSIA
#else
     jjdim=jbwa
#endif
!
!    for all even zonal wave, except 0
!
     do k = 1,lotl
       Do J = 1,jb
         jg=latdef(jls+j-jgs)
         jj=jbw-j+1
         do mx = mta+3,mtg,2            ! 2, 4,.., MTA
           FWAVEJ( mx,J ,k)= FWAVEJ( mx,J ,k)*COSLAT(jg,1)
           FWAVEJ( mx,jj,k)= FWAVEJ( mx,jj,k)*COSLAT(jg,1)
         enddo
         do mx = mta-1,1,-2            ! -2,-4,..,-MTA
           FWAVEJ( mx,J ,k)= FWAVEJ( mx,J ,k)*COSLAT(jg,1)
           FWAVEJ( mx,jj,k)= FWAVEJ( mx,jj,k)*COSLAT(jg,1)
         end do
       END Do ! par
     enddo
!
     IF(ipole.EQ.+1)CALL dfs_pole_cut(FWAVEJ,mtg,jjdim,coslat,lotl)
!
     do k=1,lotl
       do I=1,iba
         dO J=0,JB
           BLON(J,I)= 0.
           ALON(J,I)= 0.
         end do
       end do
!
       jst=0
       if (jodd.eq.1) then
         jst=1
       endif
       do mx = 2,mlast+1        ! -MTA,-MTA+1,..,-1,1,2,..,MTA
         IBmx=IBA-mx+2
         mm=mta+2-mx          ! mta+1 is 0 wave number
         mp=mta+mx
         Do J = jst,jend
           BLON(J,  mx)= ( FWAVEJ(mp,jhs(j),k)+FWAVEJ(mm,jhn(j),k))*0.5
           BLON(J,IBmx)= (-FWAVEJ(mp,jhs(j),k)+FWAVEJ(mm,jhn(j),k))*0.5 
           ALON(J,  mx)= ( FWAVEJ(mm,jhs(j),k)-FWAVEJ(mp,jhn(j),k))*0.5 
           ALON(J,IBmx)= ( FWAVEJ(mm,jhs(j),k)+FWAVEJ(mp,jhn(j),k))*0.5 
         END Do
         if (jodd.eq.1) then
           BLON(0,  mx)= ( FWAVEJ(mp,jhs(0),k)               )*0.5 
           BLON(0,IBmx)= (-FWAVEJ(mp,jhs(0),k)               )*0.5 
           ALON(0,  mx)= ( FWAVEJ(mm,jhs(0),k)               )*0.5 
           ALON(0,IBmx)= ( FWAVEJ(mm,jhs(0),k)               )*0.5 
         endif
       end do ! mx=

       Do J = 0,jend
         BLON(J,1)=   FWAVEJ( mta+1,jhn(j),k)
         ALON(J,1)=   FWAVEJ( mta+1,jhs(j),k)
       END Do
       if (jodd.eq.1) BLON(0,  1)=   0.  ! to avoid the same input
       if(MTA.eq.(IBA/2)) then
         Do J = 0,jend
           BLON(J,mta)=   FWAVEJ(1,jhn(j),k)
           ALON(J,mta)=   FWAVEJ(1,jhs(j),k)
         END Do
         if (jodd.eq.1) BLON(0,  mta)=   0.  ! to avoid the same input
       endif
!
       CALL dfs_fft_all( ALON,BLON,TRIGSIB,1,IBA,IBA,1,isgn,JB,JEND)
!
       do I = 1,iba
         DO J = 0,jend
           PSI(I,jhn(j),k )= BLON(J ,I)
           PSI(I,jhs(j),k )= ALON(J ,I)
         end do
       END DO
     enddo
!
     deallocate(FWAVEJ)
#ifdef MP
#undef PSI
#undef FWAVEJ
     if (nvar.ne.0) then
       call mpyz2xy(psia,iba,levsp,PSI,ib,levs,jbw,nvar)
     endif
     if (n1lv.ne.0) then
       call mpy2xy(psia(1,1,n1stl),iba,psi(1,1,n1stg),ib,jbw,n1lv)
     endif
     deallocate(psia)
#endif
!=============================
   ENDIF                  
!=============================
   RETURN
   END SUBROUTINE dfs_fft_driver
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_fft_setup(TRIGS,N)
!-------------------------------------------------------------------------------
!
!        SUBROUTINE 'SETGPFA_8_VEC' : REAL*8 VERSION 
!                                   ; Vectorized version (02 FEB 2005) 
!        = the same as SETGPFA_8.FOR
!        SUBROUTINE 'SETGPFA8' : REAL*8 VERSION
!     <- SUBROUTINE 'SETGPFA'
!        SETUP ROUTINE FOR SELF-SORTING IN-PLACE
!            GENERALIZED PRIME FACTOR (COMPLEX) FFT [GPFA]
!
!        CALL SETGPFA(TRIGS,N)
!
!        INPUT :
!        -----
!        N IS THE LENGTH OF THE TRANSFORMS. N MUST BE OF THE FORM:
!          -----------------------------------
!            N = (2**IP) * (3**IQ) * (5**IR)
!          -----------------------------------
!
!        OUTPUT:
!        ------
!        TRIGS IS A TABLE OF TWIDDLE FACTORS,
!          OF LENGTH 2*IPQR (REAL) WORDS, WHERE:
!          --------------------------------------
!            IPQR = (2**IP) + (3**IQ) + (5**IR)
!          --------------------------------------
!
!        WRITTEN BY CLIVE TEMPERTON 1990
!
!-------------------------------------------------------------------------------
   integer  ::  N
   real     ::  TRIGS(N*2)
   integer  ::  NJ(3)
!-------------------------------------------------------------------------------
!
!     DECOMPOSE N INTO FACTORS 2,3,5
!     ------------------------------
!
   NN = N
   IFAC = 2
!
   DO LL = 1,3
     KK = 0
!
     10 CONTINUE
!
     IF (MOD(NN,IFAC).NE.0) GO TO 20
     KK = KK + 1
     NN = NN / IFAC
     GO TO 10
!
     20 CONTINUE
!
     NJ(LL) = KK
     IFAC = IFAC + LL
   enddo
!
   IF (NN.NE.1) THEN
     WRITE(6,40) N
 40  FORMAT(' *** WARNING!!!',I10,' IS NOT A LEGAL VALUE OF N ***')
     RETURN
   ENDIF
!
   IP = NJ(1)
   IQ = NJ(2)
   IR = NJ(3)
!
!     COMPUTE LIST OF ROTATED TWIDDLE FACTORS
!     ---------------------------------------
!
   NJ(1) = 2**IP
   NJ(2) = 3**IQ
   NJ(3) = 5**IR
!
#ifdef IBMSP
   TWOPI = 4.D0 * ASIN(DBLE(1.0D0)) ! REAL*8
#else
   TWOPI = 4.D0 * DASIN(1.0D0) ! REAL*8
#endif
   I = 1
!
   DO 60 LL = 1 , 3
     NI = NJ(LL)
     IF (NI.EQ.1) GO TO 60
!
#ifdef IBMSP
     DEL = TWOPI / DBLE(NI) ! REAL*8
#else
     DEL = TWOPI / DFLOAT(NI) ! REAL*8
#endif
     IROT = N / NI
     KINK = MOD(IROT,NI)
     KK = 0
!
     DO K = 1,NI
#ifdef IBMSP
       ANGLE = DBLE(KK) * DEL ! REAL*8
       TRIGS(I) = COS(DBLE(ANGLE)) ! REAL*8
       TRIGS(I+1) = SIN(DBLE(ANGLE)) ! REAL*8
#else
       ANGLE = dfloat(KK) * DEL ! REAL*8
       TRIGS(I) = DCOS(dble(ANGLE)) ! REAL*8
       TRIGS(I+1) = DSIN(dble(ANGLE)) ! REAL*8
#endif
       I = I + 2
       KK = KK + KINK
       IF (KK.GT.NI) KK = KK - NI
     enddo
!
   60 CONTINUE
!
   RETURN
   END SUBROUTINE dfs_fft_setup
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_sine_transform                                               &
   ( data, SINHF, COSPK, SINPK,TRIGS, N, N1, ND2, MT, isign )          ! ND2=N/2
!-------------------------------------------------------------------------------
#ifdef MP
   use commpi, only : mype
#endif
!-------------------------------------------------------------------------------
!     In this DFS_global Model,
!     X(j)= SUM_mj ANM(mj)*HALFSINCOSINE is used instead of USUAL
!     form [ X(j)= (1/2) ANM(0) + SUM_mj=1^N ANM(mj)*HALFSINCOSINE ].
!
!     At exit 'CALL dfs_fft_all', the coefficients are all normalized.
!
!     Due to the difference in (1/2)ANM(0), the normalization factors
!     are different for wavenum=0 and Nyquist from USUAL routines.
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer  ::  N, N1, ND2, MT, isign
   real data(mt,0:N1),TRIGS(N*2)
   real     ::  SINHF(0:ND2),COSPK(0:ND2),SINPK(0:ND2)
!
! local variables
!
   real                     ::  ALAT(0:MT/2,0:N1),BLAT(0:MT/2,0:N1)
   real                     ::  partRE(MT,0:ND2),partIM(MT,0:ND2),temp(MT)
   integer                  ::  k,m,j,nk,kom,i,kev,kod,kop,nj
   real                     ::  aa,bb,ty,tx,sca1N,sca2N
   real                     ::  rn
#ifndef MP
   integer, parameter       :: mype=0
#endif
   integer                  ::  iunit
   integer,allocatable,save ::  idl(:),idr(:)
   integer,save             ::  mdim,mend,modd,mst
   logical                  ::  first
   data      first/.true./
!
! save      idl,idr,mend,mdim,modd
!-------------------------------------------------------------------------------
!
!     Calculates the Cosine Fourier transform of a set of n-real valued
!     data points. OUTPUT is NORMALIZED.
!-----
!     f_j= F_1 SIN(pi*(j+0.5)*1/N) + ...... + F_N SIN(pi*(j+0.5)*N/N)
!     F_k= [ sum_j=1^j=N-1 f_j SIN(pi*(j+0.5)*k/N) ] * (2/N)
!     F_N= [ sum_j=0^j=N-1 f_j SIN(pi*(j+0.5)*N/N) ] * (1/N)
!-----
!     FOUR1,REALFT : ARRAY is bounded as 1:N
!     COS_FT       : ARRAY is bounded as 0:N-1
!     SIN_FT       : ARRAY is bounded as 0:N-1
!     REALFT = REALFTN except the NORMALIZATION in REALFTN
!-----
!---------------------------------------------------------------------
!     ISIGN=+1     : INPUT=GRIDS OUTPUT=WAVES
!                                1st data : SIN_1
!                                2nd data : SIN_2 ......
!     ISIGN=-1     : INPUT=WAVES OUTPUT=GRIDS
!                 
!---------------------------------------------------------------------
!     PI = DATAN(1.D-00)*4.D-00
!     PIN= PI/dble(N)
!---------------------------------------------------------------------
!        qCALL= qCALL+1.0D-00
!     IF(qCALL.eq.1.D-00) THEN         ! For a Variable transform
!        do j=0,ND2                    ! Calculates once only
!           aa= PIN*(j+0.5D-00)
!           SINHF(j)= DSIN(aa)
!        end do
!        do k=0,ND2
!           bb= PIN*dble(k)
!           COSPK(k)= DCOS(bb)
!           SINPK(k)= DSIN(bb)
!        end do
!     ENDIF
!---------------------------------------------------------------------
   iunit=100+mype
   rn=1./n
   if ( first) then
     first=.false.
     mdim=mt/2
     mend=(mt-1)/2
     allocate(idl(0:mend),idr(0:mend))
     if (mod(mt,2).eq.0) then
       modd=0
       Do m=0,mend
         idl(m)=mend+1-m
         idr(m)=mend+2+m
#ifdef DBG_DFS
         write(iunit,'(A,3I5)')'even m,idl,idr=',m,idl(m),idr(m)
         call flush(iunit)
#endif
       enddo
     else
       modd=1
       Do m=0,mend
         idl(m)=mend+1-m
         idr(m)=mend+1+m
#ifdef DBG_DFS
         write(iunit,'(A,3I5)')'odd m,idl,idr=',m,idl(m),idr(m)
         call flush(iunit)
#endif
       enddo
     endif
   endif
!=============================
   IF(isign.EQ.+1) THEN         ! GRIDS->WAVES
!=============================
!-----                                 ! variable transform
!               g_j=1/2(f_j-f_(N-j))+(f_j+f_(N-j))*sin(j*pi/N)
!
     do j = 0,ND2-1
       nj = N1-j
       DO m = 1,mt
         aa= (data(m,j )-data(m,nj))*0.5+(data(m,j )+                          &
              data(m,nj))*SINHF(j)
         bb= (data(m,nj)-data(m,j ))*0.5+(data(m,nj)+                          &
              data(m,j ))*SINHF(j)
         data(m, j)= aa
         data(m,nj)= bb
       END DO
     end do
!-----
!-beg------------------------------------------! 2pi periodic FFT ---
!     CALL REALFT( data,N, isign )     ! CALL non-NORMALIZING ROUTINE
!

     do I = 0,N1
       Do m = 0,mend
         BLAT(m,I)= data(idr(m),I)
         ALAT(m,I)= data(idl(m),I)
       end do
       if (modd.eq.1) BLAT(0,I)= 0.d0      
     END Do
!
     CALL dfs_fft_all( ALAT,BLAT,TRIGS,1,N,N,1,isign,mdim,mend)
!
!    Output with isign=+1 : BLAT(Nyquist)= SIN or COS (Nyquist)
!
     do k = 1,ND2-1
       Nk = N-k
       Do m=0,mend
         partRE(idr(m), k )= ( BLAT(m, k )+BLAT(m,Nk))*rn
         partRE(idl(m), k )= ( ALAT(m, k )+ALAT(m,Nk))*rn
         partIM(idr(m), k )= (-ALAT(m, k )+ALAT(m,Nk))*rn
         partIM(idl(m), k )= ( BLAT(m, k )-BLAT(m,Nk))*rn
       END Do
     end do
     Do m = 0,mend
       partRE(idr(m), 0 )=   BLAT(m, 0 )*rn
       partRE(idl(m), 0 )=   ALAT(m, 0 )*rn
       partIM(idr(m), 0 )=   0.d0
       partIM(idl(m), 0 )=   0.d0
       partRE(idr(m),ND2)=   BLAT(m,ND2)*rn
       partRE(idl(m),ND2)=   ALAT(m,ND2)*rn
       partIM(idr(m),ND2)=   0.d0
       partIM(idl(m),ND2)=   0.d0
     END Do
!
!-end------------------------------------------! 2pi periodic FFT ---
!
     do k = 1,ND2-1
       kev = k*2
       DO m = 1,mt
         data(m,kev)= partRE(m,k)*SINPK(k)+partIM(m,k)*COSPK(k)
       END DO
     end do
     DO m = 1,mt
       data(m,0)= partRE(m,ND2)        ! SIN_N (=last) compo.
       data(m,1)= partRE(m, 0 )        ! Factor ()/2 deleted because it is
     END DO
     do k = 1,ND2-1                     ! recurrence
       kod = k*2+1                    ! SIN_k <- COS_2k and SIN_2k (REALFT)
       DO m = 1,mt
         data(m,kod)= (partRE(m,k)*COSPK(k)-partIM(m,k)*SINPK(k))              &
                      +  data(m,kod-2)
       END DO
     end do
!
!-----
!
     do m = 1,mt
       temp(m)   = data(m, 0 )
     END DO
     do j = 0,N-2
       do m = 1,mt
         data(m,j )= data(m,j+1)
       END DO
     end do
     do m = 1,mt
       data(m,N1)= temp(m)
     END DO
!
!=============================
   ELSEIF(isign.EQ.-1) THEN         ! WAVES->GRIDS
!=============================
!
     DO m = 1,mt
       temp(m) = data(m,N1 )          ! To make 1st=SIN_N, 2nd=SIN_1 ...
     END DO 
     do j = N1,1,-1                   ! For here, this is convenient.
       DO m = 1,mt
         data(m,j) = data(m,j-1)
       END DO 
     end do
     DO m = 1,mt
       data(m,0) = temp(m)
     END DO 
!
!-----
!
     DO m = 1,mt
       partRE(m, 0 )= data(m, 1 )      ! to use 'CALL dfs_fft_all'
       partIM(m, 0 )= 0.  
       partRE(m,ND2)= data(m, 0 )
     END DO 
     do k = 1,ND2-1                    ! SIN_k -> COS_2k and SIN_2k  
       kev= k*2                        !          in REALFT
       kop= kev+1
       kom= kev-1
       DO m = 1,mt
         partRE(m,k)= data(m,kev)*SINPK(k)                                     &
                    +(data(m,kop)-data(m,kom))*COSPK(k)  
         partIM(m,k)= data(m,kev)*COSPK(k)                                     &
                    -(data(m,kop)-data(m,kom))*SINPK(k)
       END DO 
     end do
     DO m = 1,mt
       data(m, 0 )= partRE(m, 0 )
       data(m, 1 )= partRE(m,ND2)
     END DO 
     do k = 1,ND2-1
       kev= k*2
       kod= kev+1
       DO m=1,mt
         data(m,kev)= partRE(m, k )     ! COS compo in REALFT
         data(m,kod)= partIM(m, k )     ! SIN compo in REALFT
       END DO 
     end do
!-----
!-beg------------------------------------------! 2pi periodic FFT ---
!     CALL REALFT( data,N, isign )     ! CALL non-NORMALIZING ROUTINE
!
     do k = 0,N1
       do m = 0,mdim
         BLAT(m,k)= 0.
         ALAT(m,k)= 0.
       end do
     end do
!
     mst=0
     if (modd.eq.1) mst=1
!
     do k = 1,ND2-1
       kev= k*2
       kod= kev+1
       Nk= N-k
       Do m = mst,mend
         BLAT(m, k)= ( data(idl(m),kod)+data(idr(m),kev))*0.5
         BLAT(m,Nk)= (-data(idl(m),kod)+data(idr(m),kev))*0.5
         ALAT(m, k)= ( data(idl(m),kev)-data(idr(m),kod))*0.5
         ALAT(m,Nk)= ( data(idl(m),kev)+data(idr(m),kod))*0.5
       end do
       if (modd.eq.1) then
         BLAT(0, k)= ( data(idl(0),kod)             )*0.5
         BLAT(0,Nk)= (-data(idl(0),kod)             )*0.5
         ALAT(0, k)= ( data(idl(0),kev)             )*0.5
         ALAT(0,Nk)= ( data(idl(0),kev)             )*0.5
       endif
     end do
     do m = 0,mend
       BLAT(m, 0)=   data(idr(m),0  )
       ALAT(m, 0)=   data(idl(m),0  )
     enddo
     if (modd.eq.1) BLAT(0, 0)=   0. ! to avoid the same input
     do m = 0,mend
       BLAT(m, ND2)=   data(idr(m),1  )
       ALAT(m, ND2)=   data(idl(m),1  )
     enddo
     if (modd.eq.1) BLAT(0, ND2)=   0. ! to avoid the same input
     !
     CALL dfs_fft_all( ALAT,BLAT,TRIGS, 1, N, N, 1, isign,mdim,mend )

     do I = 0,N1
       Do m = 0,mend
         data(idr(m),I)= BLAT( m,I)
         data(idl(m),I)= ALAT( m,I)
       end do
     END Do
!
!-end------------------------------------------! 2pi periodic FFT ---
!
     do j = 0,ND2-1
       nj=N1-j                       ! GRIDS in REALFT -> GRIDS in COS_FT
#ifdef LINUX_PGI
!pgi$ novector
#endif
       DO m = 1,mt
         ty= data(m,j )-data(m,nj)
         tx=(data(m,nj)+data(m,j ))/(SINHF(j)*2.)
         data(m, j)= (tx+ty)*0.5
         data(m,nj)= (tx-ty)*0.5
       END DO 
     end do
!=============================         ! isign = +-1
   ENDIF
!=============================         ! isign = +-1
!------------------------ NORMALIZATION
! NMZ is Done in 'dfs_fft_all'
!------------------------ NORMALIZATION
!===========
   RETURN
   END SUBROUTINE dfs_sine_transform
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_pole_cut(FWAVEJ,mtg,jdim,COSLAT,lot)
!-------------------------------------------------------------------------------
!
!     To give the UP-LIMIT truncation to escape from POLE PROBLEM    
! To be free from POLE PROBLEM                                        
! ms and js are global zonal wave and global latitude grid index     
!
!-------------------------------------------------------------------------------
   use dfsvar, only : jls,jgs,jge,jba,jbwa,jbw,jb,mgs,mge,iope,latdef
#ifdef MP
   use commpi, only : mype
#endif
!-------------------------------------------------------------------------------
!
   implicit none
   integer              ::  mtg,jdim,lot
   real                 ::  FWAVEJ(mgs:mge,jdim,lot),COSLAT(jgs:jge,3)
!
!-----local-variable-------------------------------------------------!
!
   integer              ::  mx,j,jj,MAXMT,k,iunit
   real                 ::  csinv
   integer,allocatable  ::  M_0_AMP(:)
   LOGICAL   JUSDOIT
   data      JUSDOIT/.TRUE./
   save      M_0_AMP
!-------------------------------------------------------------------------------
!
   IF(JUSDOIT) THEN
     JUSDOIT=.FALSE.
     allocate(M_0_AMP(jb))
!
     MAXMT = jbwa
     csinv=coslat(jbwa/3,2)
     DO J = 1,jb
       jj=latdef(j+jls-1)
       mx= MAXMT * COSLAT(jj,1)*csinv + 1.d0
       M_0_AMP(j)= mx
     END DO
   ENDIF
!
   DO k = 1,lot
     DO j = 1,jb
       mx=M_0_AMP(j)+1
       jj=jbw-j+1
       if (mx.le.mge) then
         FWAVEJ(mx :mge,j ,k)= 0.         ! +
         FWAVEJ(mgs:-mx,j ,k)= 0.         ! -
         FWAVEJ(mx :mge,jj,k)= 0.         ! +
         FWAVEJ(mgs:-mx,jj,k)= 0.         ! -
       endif
     END DO
   enddo
!
   RETURN
   END SUBROUTINE dfs_pole_cut
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine dfs_fft_2                                                        &
               (a,b,trigs,inc,jump,n,mm,lot,isign, JDIM,JEND)
!-------------------------------------------------------------------------------
   IMPLICIT  REAL*8 (A-H, O-Z)
!-------------------------------------------------------------------------------
!
!     dfs_fft_2.FOR : REAL*8 VERSION ; Vectorized version (02 FEB 2005)
!     GPFA2F_8.FOR : REAL*8 VERSION 
!     fortran version of *gpfa2* -
!     radix-2 section of self-sorting, in-place, generalized pfa
!     central radix-2 and radix-8 passes included
!      so that transform length can be any power of 2
!
!     N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. 
!     RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE 
!     (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      
!
!-------------------------------------------------------------------------------
!     dimension a(*), b(*), trigs(*)
   real  ::   a(0:JDIM,n), b(0:JDIM,n), trigs(*) ! 02 FEB 2005 to vectorize
   data lvr/64/
!-------------------------------------------------------------------------------
   n2 = 2**mm
   inq = n/n2
   jstepx = (n2-n) * inc
   ninc = n * inc
   ink = inc * inq
!
   m2 = 0
   m8 = 0
!
   if (mod(mm,2).eq.0) then
     m = mm/2
   else if (mod(mm,4).eq.1) then
     m = (mm-1)/2
     m2 = 1
   else if (mod(mm,4).eq.3) then
     m = (mm-3)/2
     m8 = 1
   endif
!
   mh = (m+1)/2
!
   nblox = 1 + (lot-1)/lvr
   left = lot
   s = Dfloat(isign) ! REAL*8
   istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
!
   do nb = 1,nblox
!
     if (left.le.lvr) then
       nvex = left
     else if (left.lt.(2*lvr)) then
       nvex = left/2
       nvex = nvex + mod(nvex,2)
     else
       nvex = lvr
     endif
!
     left = left - nvex
!
     la = 1
!
!   loop on type I radix-4 passes
!  -----------------------------
!
     mu = mod(inq,4)
     if (isign.eq.-1) mu = 4 - mu
     ss = 1.0
     if (mu.eq.3) ss = -1.0
! 
     if (mh.eq.0) go to 200
!
     do ipass = 1 , mh
       jstep = (n*inc) / (4*la)
       jstepl = jstep - ninc
!
!  k = 0 loop (no twiddle factors)
!  -------------------------------
!
       do jjj = 0 , (n-1)*inc , 4*jstep
         ja = istart + jjj
!
!     "transverse" loop
!     -----------------
!
         do nu = 1 , inq
           jb = ja + jstepl
           if (jb.lt.istart) jb = jb + ninc
           jc = jb + jstepl
           if (jc.lt.istart) jc = jc + ninc
           jd = jc + jstepl
           if (jd.lt.istart) jd = jd + ninc
           j = 0
!
!  loop across transforms
!  ----------------------

!dir$ ivdep, shortloop
           do l = 1 , nvex
             DO jnew=0,JEND ! 02 FEB 2005 to vectorize
               aja = a(jnew,ja+j)
               ajc = a(jnew,jc+j)
               t0 = aja + ajc
               t2 = aja - ajc
               ajb = a(jnew,jb+j)
               ajd = a(jnew,jd+j)
               t1 = ajb + ajd
               t3 = ss * ( ajb - ajd )
               bja = b(jnew,ja+j)
               bjc = b(jnew,jc+j)
               u0 = bja + bjc
               u2 = bja - bjc
               bjb = b(jnew,jb+j)
               bjd = b(jnew,jd+j)
               u1 = bjb + bjd
               u3 = ss * ( bjb - bjd )
               a(jnew,ja+j) = t0 + t1
               a(jnew,jc+j) = t0 - t1
               b(jnew,ja+j) = u0 + u1
               b(jnew,jc+j) = u0 - u1
               a(jnew,jb+j) = t2 - u3
               a(jnew,jd+j) = t2 + u3
               b(jnew,jb+j) = u2 + t3
               b(jnew,jd+j) = u2 - t3
             END DO ! 02 FEB 2005 to vectorize
           j = j + jump
         enddo
         ja = ja + jstepx
         if (ja.lt.istart) ja = ja + ninc
       enddo
     enddo
!
!  finished if n2 = 4
!  ------------------
!
     if (n2.eq.4) go to 490
       kk = 2 * la
!
!  loop on nonzero k
!  -----------------
!
       do k = ink , jstep-ink , ink
         co1 = trigs(kk+1)
         si1 = s*trigs(kk+2)
         co2 = trigs(2*kk+1)
         si2 = s*trigs(2*kk+2)
         co3 = trigs(3*kk+1)
         si3 = s*trigs(3*kk+2)
!
!  loop along transform
!  --------------------
!
         do jjj = k , (n-1)*inc , 4*jstep
           ja = istart + jjj
!
!     "transverse" loop
!     -----------------
!
           do nu = 1 , inq
             jb = ja + jstepl
             if (jb.lt.istart) jb = jb + ninc
             jc = jb + jstepl
             if (jc.lt.istart) jc = jc + ninc
             jd = jc + jstepl
             if (jd.lt.istart) jd = jd + ninc
             j = 0
!
!  loop across transforms
!  ----------------------

!dir$ ivdep,shortloop
             do l = 1 , nvex
               DO jnew=0,JEND ! 02 FEB 2005 to vectorize
                 aja = a(jnew,ja+j)
                 ajc = a(jnew,jc+j)
                 t0 = aja + ajc
                 t2 = aja - ajc
                 ajb = a(jnew,jb+j)
                 ajd = a(jnew,jd+j)
                 t1 = ajb + ajd
                 t3 = ss * ( ajb - ajd )
                 bja = b(jnew,ja+j)
                 bjc = b(jnew,jc+j)
                 u0 = bja + bjc
                 u2 = bja - bjc
                 bjb = b(jnew,jb+j)
                 bjd = b(jnew,jd+j)
                 u1 = bjb + bjd
                 u3 = ss * ( bjb - bjd )
                 a(jnew,ja+j) = t0 + t1
                 b(jnew,ja+j) = u0 + u1
                 a(jnew,jb+j) = co1*(t2-u3) - si1*(u2+t3)
                 b(jnew,jb+j) = si1*(t2-u3) + co1*(u2+t3)
                 a(jnew,jc+j) = co2*(t0-t1) - si2*(u0-u1)
                 b(jnew,jc+j) = si2*(t0-t1) + co2*(u0-u1)
                 a(jnew,jd+j) = co3*(t2+u3) - si3*(u2-t3)
                 b(jnew,jd+j) = si3*(t2+u3) + co3*(u2-t3)
               END DO 
               j = j + jump
             enddo
!
!-----( end of loop across transforms )
!
             ja = ja + jstepx
             if (ja.lt.istart) ja = ja + ninc
           enddo
         enddo
!
!-----( end of loop along transforms )
!
         kk = kk + 2*la
       enddo
!
!-----( end of loop on nonzero k )
!
       la = 4*la
     enddo
!
!-----( end of loop on type I radix-4 passes)
!
!    central radix-2 pass
!  --------------------
     200 continue
     if (m2.eq.0) go to 300
! 
     jstep = (n*inc) / (2*la)
     jstepl = jstep - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
     do jjj = 0 , (n-1)*inc , 2*jstep
       ja = istart + jjj
!
!     "transverse" loop
!     -----------------
!
       do nu = 1 , inq
         jb = ja + jstepl
         if (jb.lt.istart) jb = jb + ninc
         j = 0
!
!  loop across transforms
!  ----------------------

!dir$ ivdep, shortloop
         do l = 1 , nvex
           DO jnew=0,JEND ! 02 FEB 2005 to vectorize
             aja = a(jnew,ja+j)
             ajb = a(jnew,jb+j)
             t0 = aja - ajb
             a(jnew,ja+j) = aja + ajb
             a(jnew,jb+j) = t0
             bja = b(jnew,ja+j)
             bjb = b(jnew,jb+j)
             u0 = bja - bjb
             b(jnew,ja+j) = bja + bjb
             b(jnew,jb+j) = u0
           END DO
           j = j + jump
         enddo
!
!-----(end of loop across transforms)
!
         ja = ja + jstepx
         if (ja.lt.istart) ja = ja + ninc
       enddo
     enddo
!
!  finished if n2=2
!  ----------------
     if (n2.eq.2) go to 490
!
     kk = 2 * la
!
!  loop on nonzero k
!  -----------------
     do k = ink , jstep - ink , ink
       co1 = trigs(kk+1)
       si1 = s*trigs(kk+2)
!
!  loop along transforms
!  ---------------------
       do jjj = k , (n-1)*inc , 2*jstep
         ja = istart + jjj
!
!     "transverse" loop
!     -----------------
         do nu = 1 , inq
           jb = ja + jstepl
           if (jb.lt.istart) jb = jb + ninc
           j = 0
!
!  loop across transforms
!  ----------------------
           if (kk.eq.n2/2) then
!dir$ ivdep, shortloop
             do l = 1 , nvex
               DO jnew=0,JEND ! 02 FEB 2005 to vectorize
                 aja = a(jnew,ja+j)
                 ajb = a(jnew,jb+j)
                 t0 = ss * ( aja - ajb )
                 a(jnew,ja+j) = aja + ajb
                 bjb = b(jnew,jb+j)
                 bja = b(jnew,ja+j)
                 a(jnew,jb+j) = ss * ( bjb - bja )
                 b(jnew,ja+j) = bja + bjb
                 b(jnew,jb+j) = t0
               ENDDO
               j = j + jump
             enddo
!
           else
!dir$ ivdep, shortloop
             do l = 1 , nvex
               DO jnew=0,JEND ! 02 FEB 2005 to vectorize
                 aja = a(jnew,ja+j)
                 ajb = a(jnew,jb+j)
                 t0 = aja - ajb
                 a(jnew,ja+j) = aja + ajb
                 bja = b(jnew,ja+j)
                 bjb = b(jnew,jb+j)
                 u0 = bja - bjb
                 b(jnew,ja+j) = bja + bjb
                 a(jnew,jb+j) = co1*t0 - si1*u0
                 b(jnew,jb+j) = si1*t0 + co1*u0
               END DO
               j = j + jump
             enddo
!
           endif
!
!-----(end of loop across transforms)
!
           ja = ja + jstepx
           if (ja.lt.istart) ja = ja + ninc
         enddo
       enddo
!
!-----(end of loop along transforms)
!
       kk = kk + 2 * la
     enddo
!
!-----(end of loop on nonzero k)
!-----(end of radix-2 pass)
!
     la = 2 * la
     go to 400
!
!  central radix-8 pass
!  --------------------
     300 continue
     if (m8.eq.0) go to 400
     jstep = (n*inc) / (8*la)
     jstepl = jstep - ninc
     mu = mod(inq,8)
     if (isign.eq.-1) mu = 8 - mu
     c1 = 1.0
     if (mu.eq.3.or.mu.eq.7) c1 = -1.0
     c2 = Dsqrt(0.5D0) ! REAL*8
     if (mu.eq.3.or.mu.eq.5) c2 = -c2
     c3 = c1 * c2
!
!  stage 1
!  -------
     do k = 0 , jstep - ink , ink
       do jjj = k , (n-1)*inc , 8*jstep
         ja = istart + jjj
!
!     "transverse" loop
!     -----------------
         do nu = 1 , inq
           jb = ja + jstepl
           if (jb.lt.istart) jb = jb + ninc
           jc = jb + jstepl
           if (jc.lt.istart) jc = jc + ninc
           jd = jc + jstepl
           if (jd.lt.istart) jd = jd + ninc
           je = jd + jstepl
           if (je.lt.istart) je = je + ninc
           jf = je + jstepl
           if (jf.lt.istart) jf = jf + ninc
           jg = jf + jstepl
           if (jg.lt.istart) jg = jg + ninc
           jh = jg + jstepl
           if (jh.lt.istart) jh = jh + ninc
           j = 0
!dir$ ivdep, shortloop
           do l = 1 , nvex
             DO jnew=0,JEND ! 02 FEB 2005 to vectorize
               aja = a(jnew,ja+j)
               aje = a(jnew,je+j)
               t0 = aja - aje
               a(jnew,ja+j) = aja + aje
               ajc = a(jnew,jc+j)
               ajg = a(jnew,jg+j)
               t1 = c1 * ( ajc - ajg )
               a(jnew,je+j) = ajc + ajg
               ajb = a(jnew,jb+j)
               ajf = a(jnew,jf+j)
               t2 = ajb - ajf
               a(jnew,jc+j) = ajb + ajf
               ajd = a(jnew,jd+j)
               ajh = a(jnew,jh+j)
               t3 = ajd - ajh
               a(jnew,jg+j) = ajd + ajh
               a(jnew,jb+j) = t0
               a(jnew,jf+j) = t1
               a(jnew,jd+j) = c2 * ( t2 - t3 )
               a(jnew,jh+j) = c3 * ( t2 + t3 )
               bja = b(jnew,ja+j)
               bje = b(jnew,je+j)
               u0 = bja - bje
               b(jnew,ja+j) = bja + bje
               bjc = b(jnew,jc+j)
               bjg = b(jnew,jg+j)
               u1 = c1 * ( bjc - bjg )
               b(jnew,je+j) = bjc + bjg
               bjb = b(jnew,jb+j)
               bjf = b(jnew,jf+j)
               u2 = bjb - bjf
               b(jnew,jc+j) = bjb + bjf
               bjd = b(jnew,jd+j)
               bjh = b(jnew,jh+j)
               u3 = bjd - bjh
               b(jnew,jg+j) = bjd + bjh
               b(jnew,jb+j) = u0
               b(jnew,jf+j) = u1
               b(jnew,jd+j) = c2 * ( u2 - u3 )
               b(jnew,jh+j) = c3 * ( u2 + u3 )
             END DO
             j = j + jump
           enddo
           ja = ja + jstepx
           if (ja.lt.istart) ja = ja + ninc
         enddo
       enddo
     enddo
!
!  stage 2
!  -------
!
!  k=0 (no twiddle factors)
!  ------------------------
     do jjj = 0 , (n-1)*inc , 8*jstep
       ja = istart + jjj
!
!     "transverse" loop
!     -----------------
       do nu = 1 , inq
         jb = ja + jstepl
         if (jb.lt.istart) jb = jb + ninc
         jc = jb + jstepl
         if (jc.lt.istart) jc = jc + ninc
         jd = jc + jstepl
         if (jd.lt.istart) jd = jd + ninc
         je = jd + jstepl
         if (je.lt.istart) je = je + ninc
         jf = je + jstepl
         if (jf.lt.istart) jf = jf + ninc
         jg = jf + jstepl
         if (jg.lt.istart) jg = jg + ninc
         jh = jg + jstepl
         if (jh.lt.istart) jh = jh + ninc
         j = 0
!dir$ ivdep, shortloop
         do l = 1 , nvex
           DO jnew=0,JEND ! 02 FEB 2005 to vectorize
             aja = a(jnew,ja+j)
             aje = a(jnew,je+j)
             t0 = aja + aje
             t2 = aja - aje
             ajc = a(jnew,jc+j)
             ajg = a(jnew,jg+j)
             t1 = ajc + ajg
             t3 = c1 * ( ajc - ajg )
             bja = b(jnew,ja+j)
             bje = b(jnew,je+j)
             u0 = bja + bje
             u2 = bja - bje
             bjc = b(jnew,jc+j)
             bjg = b(jnew,jg+j)
             u1 = bjc + bjg
             u3 = c1 * ( bjc - bjg )
             a(jnew,ja+j) = t0 + t1
             a(jnew,je+j) = t0 - t1
             b(jnew,ja+j) = u0 + u1
             b(jnew,je+j) = u0 - u1
             a(jnew,jc+j) = t2 - u3
             a(jnew,jg+j) = t2 + u3
             b(jnew,jc+j) = u2 + t3
             b(jnew,jg+j) = u2 - t3
             ajb = a(jnew,jb+j)
             ajd = a(jnew,jd+j)
             t0 = ajb + ajd
             t2 = ajb - ajd
             ajf = a(jnew,jf+j)
             ajh = a(jnew,jh+j)
             t1 = ajf - ajh
             t3 = ajf + ajh
             bjb = b(jnew,jb+j)
             bjd = b(jnew,jd+j)
             u0 = bjb + bjd
             u2 = bjb - bjd
             bjf = b(jnew,jf+j)
             bjh = b(jnew,jh+j)
             u1 = bjf - bjh
             u3 = bjf + bjh
             a(jnew,jb+j) = t0 - u3
             a(jnew,jh+j) = t0 + u3
             b(jnew,jb+j) = u0 + t3
             b(jnew,jh+j) = u0 - t3
             a(jnew,jd+j) = t2 + u1
             a(jnew,jf+j) = t2 - u1
             b(jnew,jd+j) = u2 - t1
             b(jnew,jf+j) = u2 + t1
           END DO
           j = j + jump
         enddo
         ja = ja + jstepx
         if (ja.lt.istart) ja = ja + ninc
       enddo
     enddo
!
     if (n2.eq.8) go to 490
!
!  loop on nonzero k
!  -----------------
     kk = 2 * la
!
     do k = ink , jstep - ink , ink
!
       co1 = trigs(kk+1)
       si1 = s * trigs(kk+2)
       co2 = trigs(2*kk+1)
       si2 = s * trigs(2*kk+2)
       co3 = trigs(3*kk+1)
       si3 = s * trigs(3*kk+2)
       co4 = trigs(4*kk+1)
       si4 = s * trigs(4*kk+2)
       co5 = trigs(5*kk+1)
       si5 = s * trigs(5*kk+2)
       co6 = trigs(6*kk+1)
       si6 = s * trigs(6*kk+2)
       co7 = trigs(7*kk+1)
       si7 = s * trigs(7*kk+2)
!
       do jjj = k , (n-1)*inc , 8*jstep
         ja = istart + jjj
!
!     "transverse" loop
!     -----------------
           do nu = 1 , inq
             jb = ja + jstepl
             if (jb.lt.istart) jb = jb + ninc
             jc = jb + jstepl
             if (jc.lt.istart) jc = jc + ninc
             jd = jc + jstepl
             if (jd.lt.istart) jd = jd + ninc
             je = jd + jstepl
             if (je.lt.istart) je = je + ninc
             jf = je + jstepl
             if (jf.lt.istart) jf = jf + ninc
             jg = jf + jstepl
             if (jg.lt.istart) jg = jg + ninc
             jh = jg + jstepl
             if (jh.lt.istart) jh = jh + ninc
             j = 0
!dir$ ivdep, shortloop
             do l = 1 , nvex
               DO jnew=0,JEND ! 02 FEB 2005 to vectorize
                 aja = a(jnew,ja+j)
                 aje = a(jnew,je+j)
                 t0 = aja + aje
                 t2 = aja - aje
                 ajc = a(jnew,jc+j)
                 ajg = a(jnew,jg+j)
                 t1 = ajc + ajg
                 t3 = c1 * ( ajc - ajg )
                 bja = b(jnew,ja+j)
                 bje = b(jnew,je+j)
                 u0 = bja + bje
                 u2 = bja - bje
                 bjc = b(jnew,jc+j)
                 bjg = b(jnew,jg+j)
                 u1 = bjc + bjg
                 u3 = c1 * ( bjc - bjg )
                 a(jnew,ja+j) = t0 + t1
                 b(jnew,ja+j) = u0 + u1
                 a(jnew,je+j) = co4*(t0-t1) - si4*(u0-u1)
                 b(jnew,je+j) = si4*(t0-t1) + co4*(u0-u1)
                 a(jnew,jc+j) = co2*(t2-u3) - si2*(u2+t3)
                 b(jnew,jc+j) = si2*(t2-u3) + co2*(u2+t3)
                 a(jnew,jg+j) = co6*(t2+u3) - si6*(u2-t3)
                 b(jnew,jg+j) = si6*(t2+u3) + co6*(u2-t3)
                 ajb = a(jnew,jb+j)
                 ajd = a(jnew,jd+j)
                 t0 = ajb + ajd
                 t2 = ajb - ajd
                 ajf = a(jnew,jf+j)
                 ajh = a(jnew,jh+j)
                 t1 = ajf - ajh
                 t3 = ajf + ajh
                 bjb = b(jnew,jb+j)
                 bjd = b(jnew,jd+j)
                 u0 = bjb + bjd
                 u2 = bjb - bjd
                 bjf = b(jnew,jf+j)
                 bjh = b(jnew,jh+j)
                 u1 = bjf - bjh
                 u3 = bjf + bjh
                 a(jnew,jb+j) = co1*(t0-u3) - si1*(u0+t3)
                 b(jnew,jb+j) = si1*(t0-u3) + co1*(u0+t3)
                 a(jnew,jh+j) = co7*(t0+u3) - si7*(u0-t3)
                 b(jnew,jh+j) = si7*(t0+u3) + co7*(u0-t3)
                 a(jnew,jd+j) = co3*(t2+u1) - si3*(u2-t1)
                 b(jnew,jd+j) = si3*(t2+u1) + co3*(u2-t1)
                 a(jnew,jf+j) = co5*(t2-u1) - si5*(u2+t1)
                 b(jnew,jf+j) = si5*(t2-u1) + co5*(u2+t1)
               END DO
               j = j + jump
             enddo
             ja = ja + jstepx
             if (ja.lt.istart) ja = ja + ninc
           enddo
         enddo
         kk = kk + 2 * la
       enddo
!
       la = 8 * la
!
!  loop on type II radix-4 passes
!  ------------------------------
       400 continue
       mu = mod(inq,4)
       if (isign.eq.-1) mu = 4 - mu
       ss = 1.0
       if (mu.eq.3) ss = -1.0
!
       do ipass = mh+1 , m
         jstep = (n*inc) / (4*la)
         jstepl = jstep - ninc
         laincl = la * ink - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
         do ll = 0 , (la-1)*ink , 4*jstep
!
           do jjj = ll , (n-1)*inc , 4*la*ink
             ja = istart + jjj
!
!     "transverse" loop
!     -----------------
             do nu = 1 , inq
               jb = ja + jstepl
               if (jb.lt.istart) jb = jb + ninc
               jc = jb + jstepl
               if (jc.lt.istart) jc = jc + ninc
               jd = jc + jstepl
               if (jd.lt.istart) jd = jd + ninc
               je = ja + laincl
               if (je.lt.istart) je = je + ninc
               jf = je + jstepl
               if (jf.lt.istart) jf = jf + ninc
               jg = jf + jstepl
               if (jg.lt.istart) jg = jg + ninc
               jh = jg + jstepl
               if (jh.lt.istart) jh = jh + ninc
               ji = je + laincl
               if (ji.lt.istart) ji = ji + ninc
               jj = ji + jstepl
               if (jj.lt.istart) jj = jj + ninc
               jk = jj + jstepl
               if (jk.lt.istart) jk = jk + ninc
               jl = jk + jstepl
               if (jl.lt.istart) jl = jl + ninc
               jm = ji + laincl
               if (jm.lt.istart) jm = jm + ninc
               jn = jm + jstepl
               if (jn.lt.istart) jn = jn + ninc
               jo = jn + jstepl
               if (jo.lt.istart) jo = jo + ninc
               jp = jo + jstepl
               if (jp.lt.istart) jp = jp + ninc
               j = 0
!
!  loop across transforms
!  ----------------------

!dir$ ivdep, shortloop
               do l = 1 , nvex
                 DO jnew=0,JEND ! 02 FEB 2005 to vectorize
                   aja = a(jnew,ja+j)
                   ajc = a(jnew,jc+j)
                   t0 = aja + ajc
                   t2 = aja - ajc
                   ajb = a(jnew,jb+j)
                   ajd = a(jnew,jd+j)
                   t1 = ajb + ajd
                   t3 = ss * ( ajb - ajd )
                   aji = a(jnew,ji+j)
                   ajc =  aji
                   bja = b(jnew,ja+j)
                   bjc = b(jnew,jc+j)
                   u0 = bja + bjc
                   u2 = bja - bjc
                   bjb = b(jnew,jb+j)
                   bjd = b(jnew,jd+j)
                   u1 = bjb + bjd
                   u3 = ss * ( bjb - bjd )
                   aje = a(jnew,je+j)
                   ajb =  aje
                   a(jnew,ja+j) = t0 + t1
                   a(jnew,ji+j) = t0 - t1
                   b(jnew,ja+j) = u0 + u1
                   bjc =  u0 - u1
                   bjm = b(jnew,jm+j)
                   bjd =  bjm
                   a(jnew,je+j) = t2 - u3
                   ajd =  t2 + u3
                   bjb =  u2 + t3
                   b(jnew,jm+j) = u2 - t3
!
!----------------------
!
                   ajg = a(jnew,jg+j)
                   t0 = ajb + ajg
                   t2 = ajb - ajg
                   ajf = a(jnew,jf+j)
                   ajh = a(jnew,jh+j)
                   t1 = ajf + ajh
                   t3 = ss * ( ajf - ajh )
                   ajj = a(jnew,jj+j)
                   ajg =  ajj
                   bje = b(jnew,je+j)
                   bjg = b(jnew,jg+j)
                   u0 = bje + bjg
                   u2 = bje - bjg
                   bjf = b(jnew,jf+j)
                   bjh = b(jnew,jh+j)
                   u1 = bjf + bjh
                   u3 = ss * ( bjf - bjh )
                   b(jnew,je+j) = bjb
                   a(jnew,jb+j) = t0 + t1
                   a(jnew,jj+j) = t0 - t1
                   bjj = b(jnew,jj+j)
                   bjg =  bjj
                   b(jnew,jb+j) = u0 + u1
                   b(jnew,jj+j) = u0 - u1
                   a(jnew,jf+j) = t2 - u3
                   ajh =  t2 + u3
                   b(jnew,jf+j) = u2 + t3
                   bjh =  u2 - t3
!
!----------------------
!
                   ajk = a(jnew,jk+j)
                   t0 = ajc + ajk
                   t2 = ajc - ajk
                   ajl = a(jnew,jl+j)
                   t1 = ajg + ajl
                   t3 = ss * ( ajg - ajl )
                   bji = b(jnew,ji+j)
                   bjk = b(jnew,jk+j)
                   u0 = bji + bjk
                   u2 = bji - bjk
                   ajo = a(jnew,jo+j)
                   ajl =  ajo
                   bjl = b(jnew,jl+j)
                   u1 = bjg + bjl
                   u3 = ss * ( bjg - bjl )
                   b(jnew,ji+j) = bjc
                   a(jnew,jc+j) = t0 + t1
                   a(jnew,jk+j) = t0 - t1
                   bjo = b(jnew,jo+j)
                   bjl =  bjo
                   b(jnew,jc+j) = u0 + u1
                   b(jnew,jk+j) = u0 - u1
                   a(jnew,jg+j) = t2 - u3
                   a(jnew,jo+j) = t2 + u3
                   b(jnew,jg+j) = u2 + t3
                   b(jnew,jo+j) = u2 - t3
!
!----------------------
!
                   ajm = a(jnew,jm+j)
                   t0 = ajm + ajl
                   t2 = ajm - ajl
                   ajn = a(jnew,jn+j)
                   ajp = a(jnew,jp+j)
                   t1 = ajn + ajp
                   t3 = ss * ( ajn - ajp )
                   a(jnew,jm+j) = ajd
                   u0 = bjd + bjl
                   u2 = bjd - bjl
                   bjn = b(jnew,jn+j)
                   bjp = b(jnew,jp+j)
                   u1 = bjn + bjp
                   u3 = ss * ( bjn - bjp )
                   a(jnew,jn+j) = ajh
                   a(jnew,jd+j) = t0 + t1
                   a(jnew,jl+j) = t0 - t1
                   b(jnew,jd+j) = u0 + u1
                   b(jnew,jl+j) = u0 - u1
                   b(jnew,jn+j) = bjh
                   a(jnew,jh+j) = t2 - u3
                   a(jnew,jp+j) = t2 + u3
                   b(jnew,jh+j) = u2 + t3
                   b(jnew,jp+j) = u2 - t3
                 END DO
                 j = j + jump
               enddo
!
!-----( end of loop across transforms )
!
               ja = ja + jstepx
               if (ja.lt.istart) ja = ja + ninc
             enddo
           enddo
         enddo
!
!-----( end of double loop for k=0 )
!
!  finished if last pass
!  ---------------------
         if (ipass.eq.m) go to 490
!
         kk = 2*la
!
!     loop on nonzero k
!     -----------------
         do k = ink , jstep-ink , ink
           co1 = trigs(kk+1)
           si1 = s*trigs(kk+2)
           co2 = trigs(2*kk+1)
           si2 = s*trigs(2*kk+2)
           co3 = trigs(3*kk+1)
           si3 = s*trigs(3*kk+2)
!
!  double loop along first transform in block
!  ------------------------------------------
           do ll = k , (la-1)*ink , 4*jstep
!
             do jjj = ll , (n-1)*inc , 4*la*ink
               ja = istart + jjj
!
!     "transverse" loop
!     -----------------
               do nu = 1 , inq
                 jb = ja + jstepl
                 if (jb.lt.istart) jb = jb + ninc
                 jc = jb + jstepl
                 if (jc.lt.istart) jc = jc + ninc
                 jd = jc + jstepl
                 if (jd.lt.istart) jd = jd + ninc
                 je = ja + laincl
                 if (je.lt.istart) je = je + ninc
                 jf = je + jstepl
                 if (jf.lt.istart) jf = jf + ninc
                 jg = jf + jstepl
                 if (jg.lt.istart) jg = jg + ninc
                 jh = jg + jstepl
                 if (jh.lt.istart) jh = jh + ninc
                 ji = je + laincl
                 if (ji.lt.istart) ji = ji + ninc
                 jj = ji + jstepl
                 if (jj.lt.istart) jj = jj + ninc
                 jk = jj + jstepl
                 if (jk.lt.istart) jk = jk + ninc
                 jl = jk + jstepl
                 if (jl.lt.istart) jl = jl + ninc
                 jm = ji + laincl
                 if (jm.lt.istart) jm = jm + ninc
                 jn = jm + jstepl
                 if (jn.lt.istart) jn = jn + ninc
                 jo = jn + jstepl
                 if (jo.lt.istart) jo = jo + ninc
                 jp = jo + jstepl
                 if (jp.lt.istart) jp = jp + ninc
                 j = 0
!
!  loop across transforms
!  ----------------------

!dir$ ivdep, shortloop
                 do l = 1 , nvex
                   DO jnew=0,JEND ! 02 FEB 2005 to vectorize
                     aja = a(jnew,ja+j)
                     ajc = a(jnew,jc+j)
                     t0 = aja + ajc
                     t2 = aja - ajc
                     ajb = a(jnew,jb+j)
                     ajd = a(jnew,jd+j)
                     t1 = ajb + ajd
                     t3 = ss * ( ajb - ajd )
                     aji = a(jnew,ji+j)
                     ajc =  aji
                     bja = b(jnew,ja+j)
                     bjc = b(jnew,jc+j)
                     u0 = bja + bjc
                     u2 = bja - bjc
                     bjb = b(jnew,jb+j)
                     bjd = b(jnew,jd+j)
                     u1 = bjb + bjd
                     u3 = ss * ( bjb - bjd )
                     aje = a(jnew,je+j)
                     ajb =  aje
                     a(jnew,ja+j) = t0 + t1
                     b(jnew,ja+j) = u0 + u1
                     a(jnew,je+j) = co1*(t2-u3) - si1*(u2+t3)
                     bjb =  si1*(t2-u3) + co1*(u2+t3)
                     bjm = b(jnew,jm+j)
                     bjd =  bjm
                     a(jnew,ji+j) = co2*(t0-t1) - si2*(u0-u1)
                     bjc =  si2*(t0-t1) + co2*(u0-u1)
                     ajd =  co3*(t2+u3) - si3*(u2-t3)
                     b(jnew,jm+j) = si3*(t2+u3) + co3*(u2-t3)
!
!----------------------------------------
!
                     ajg = a(jnew,jg+j)
                     t0 = ajb + ajg
                     t2 = ajb - ajg
                     ajf = a(jnew,jf+j)
                     ajh = a(jnew,jh+j)
                     t1 = ajf + ajh
                     t3 = ss * ( ajf - ajh )
                     ajj = a(jnew,jj+j)
                     ajg =  ajj
                     bje = b(jnew,je+j)
                     bjg = b(jnew,jg+j)
                     u0 = bje + bjg
                     u2 = bje - bjg
                     bjf = b(jnew,jf+j)
                     bjh = b(jnew,jh+j)
                     u1 = bjf + bjh
                     u3 = ss * ( bjf - bjh )
                     b(jnew,je+j) = bjb
                     a(jnew,jb+j) = t0 + t1
                     b(jnew,jb+j) = u0 + u1
                     bjj = b(jnew,jj+j)
                     bjg =  bjj
                     a(jnew,jf+j) = co1*(t2-u3) - si1*(u2+t3)
                     b(jnew,jf+j) = si1*(t2-u3) + co1*(u2+t3)
                     a(jnew,jj+j) = co2*(t0-t1) - si2*(u0-u1)
                     b(jnew,jj+j) = si2*(t0-t1) + co2*(u0-u1)
                     ajh =  co3*(t2+u3) - si3*(u2-t3)
                     bjh =  si3*(t2+u3) + co3*(u2-t3)
!
!----------------------------------------
!
                     ajk = a(jnew,jk+j)
                     t0 = ajc + ajk
                     t2 = ajc - ajk
                     ajl = a(jnew,jl+j)
                     t1 = ajg + ajl
                     t3 = ss * ( ajg - ajl )
                     bji = b(jnew,ji+j)
                     bjk = b(jnew,jk+j)
                     u0 = bji + bjk
                     u2 = bji - bjk
                     ajo = a(jnew,jo+j)
                     ajl =  ajo
                     bjl = b(jnew,jl+j)
                     u1 = bjg + bjl
                     u3 = ss * ( bjg - bjl )
                     b(jnew,ji+j) = bjc
                     a(jnew,jc+j) = t0 + t1
                     b(jnew,jc+j) = u0 + u1
                     bjo = b(jnew,jo+j)
                     bjl =  bjo
                     a(jnew,jg+j) = co1*(t2-u3) - si1*(u2+t3)
                     b(jnew,jg+j) = si1*(t2-u3) + co1*(u2+t3)
                     a(jnew,jk+j) = co2*(t0-t1) - si2*(u0-u1)
                     b(jnew,jk+j) = si2*(t0-t1) + co2*(u0-u1)
                     a(jnew,jo+j) = co3*(t2+u3) - si3*(u2-t3)
                     b(jnew,jo+j) = si3*(t2+u3) + co3*(u2-t3)
!
!----------------------------------------
!
                     ajm = a(jnew,jm+j)
                     t0 = ajm + ajl
                     t2 = ajm - ajl
                     ajn = a(jnew,jn+j)
                     ajp = a(jnew,jp+j)
                     t1 = ajn + ajp
                     t3 = ss * ( ajn - ajp )
                     a(jnew,jm+j) = ajd
                     u0 = bjd + bjl
                     u2 = bjd - bjl
                     a(jnew,jn+j) = ajh
                     bjn = b(jnew,jn+j)
                     bjp = b(jnew,jp+j)
                     u1 = bjn + bjp
                     u3 = ss * ( bjn - bjp )
                     b(jnew,jn+j) = bjh
                     a(jnew,jd+j) = t0 + t1
                     b(jnew,jd+j) = u0 + u1
                     a(jnew,jh+j) = co1*(t2-u3) - si1*(u2+t3)
                     b(jnew,jh+j) = si1*(t2-u3) + co1*(u2+t3)
                     a(jnew,jl+j) = co2*(t0-t1) - si2*(u0-u1)
                     b(jnew,jl+j) = si2*(t0-t1) + co2*(u0-u1)
                     a(jnew,jp+j) = co3*(t2+u3) - si3*(u2-t3)
                     b(jnew,jp+j) = si3*(t2+u3) + co3*(u2-t3)
                   enddo
                   j = j + jump
                 enddo
!
!-----(end of loop across transforms)
!
                 ja = ja + jstepx
                 if (ja.lt.istart) ja = ja + ninc
               enddo
             enddo
           enddo
!
!-----( end of double loop for this k )
!
           kk = kk + 2*la
         enddo
!
!-----( end of loop over values of k )
!
         la = 4*la
       enddo
!
!-----( end of loop on type II radix-4 passes )
!-----( nvex transforms completed)
!
       490 continue
       istart = istart + nvex * jump
     enddo
!
!-----( end of loop on blocks of transforms )
!
   return
   end subroutine dfs_fft_2
!
!-------------------------------------------------------------------------------
   subroutine dfs_fft_3                                                        &
              (a,b,trigs,inc,jump,n,mm,lot,isign, JDIM,JEND)
!-------------------------------------------------------------------------------
   IMPLICIT  REAL*8 (A-H, O-Z)
!-------------------------------------------------------------------------------
!
!     dfs_fft_3.FOR : REAL*8 VERSION ; Vectorized version (02 FEB 2005) 
!     GPFA3F_8.FOR : REAL*8 VERSION 
!     fortran version of *gpfa3* -
!     radix-3 section of self-sorting, in-place
!        generalized PFA
!
!-------------------------------------------------------------------------------
!
!     dimension a(*), b(*), trigs(*)
   real  ::  a(0:JDIM,n), b(0:JDIM,n), trigs(*) ! 02 FEB 2005
   data sin60/0.866025403784437D0/
   data lvr/64/
!-------------------------------------------------------------------------------
!
!     N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. 
!     RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE 
!     (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      
!
!-------------------------------------------------------------------------------
   n3 = 3**mm
   inq = n/n3
   jstepx = (n3-n) * inc
   ninc = n * inc
   ink = inc * inq
   mu = mod(inq,3)
   if (isign.eq.-1) mu = 3-mu
   m = mm
   mh = (m+1)/2
   s = Dfloat(isign) ! REAL*8
   c1 = sin60
   if (mu.eq.2) c1 = -c1
!
   nblox = 1 + (lot-1)/lvr
   left = lot
   s = Dfloat(isign) ! REAL*8
   istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
   do nb = 1 , nblox
!
     if (left.le.lvr) then
       nvex = left
     else if (left.lt.(2*lvr)) then
       nvex = left/2
       nvex = nvex + mod(nvex,2)
     else
       nvex = lvr
     endif
     left = left - nvex
!
     la = 1
!
!  loop on type I radix-3 passes
!  -----------------------------
     do ipass = 1 , mh
       jstep = (n*inc) / (3*la)
       jstepl = jstep - ninc
!
!  k = 0 loop (no twiddle factors)
!  -------------------------------
       do jjj = 0 , (n-1)*inc , 3*jstep
         ja = istart + jjj
!
!  "transverse" loop
!  -----------------
         do nu = 1 , inq
           jb = ja + jstepl
           if (jb.lt.istart) jb = jb + ninc
           jc = jb + jstepl
           if (jc.lt.istart) jc = jc + ninc
           j = 0
!
!  loop across transforms
!  ----------------------

!dir$ ivdep, shortloop
           do l = 1 , nvex
             DO jnew=0,JEND ! 02 FEB 2005 to vectorize
               ajb = a(jnew,jb+j)
               ajc = a(jnew,jc+j)
               t1 = ajb + ajc
               aja = a(jnew,ja+j)
               t2 = aja - 0.5D0 * t1 ! REAL*8
               t3 = c1 * ( ajb - ajc )
               bjb = b(jnew,jb+j)
               bjc = b(jnew,jc+j)
               u1 = bjb + bjc
               bja = b(jnew,ja+j)
               u2 = bja - 0.5D0 * u1 ! REAL*8
               u3 = c1 * ( bjb - bjc )
               a(jnew,ja+j) = aja + t1
               b(jnew,ja+j) = bja + u1
               a(jnew,jb+j) = t2 - u3
               b(jnew,jb+j) = u2 + t3
               a(jnew,jc+j) = t2 + u3
               b(jnew,jc+j) = u2 - t3
             END DO ! 02 FEB 2005 to vectorize
             j = j + jump
           enddo
           ja = ja + jstepx
           if (ja.lt.istart) ja = ja + ninc
         enddo
       enddo
!
!  finished if n3 = 3
!  ------------------
       if (n3.eq.3) go to 490
       kk = 2 * la
!
!  loop on nonzero k
!  -----------------
       do k = ink , jstep-ink , ink
         co1 = trigs(kk+1)
         si1 = s*trigs(kk+2)
         co2 = trigs(2*kk+1)
         si2 = s*trigs(2*kk+2)
!
!  loop along transform
!  --------------------
         do jjj = k , (n-1)*inc , 3*jstep
           ja = istart + jjj
!
!  "transverse" loop
!  -----------------
           do nu = 1 , inq
             jb = ja + jstepl
             if (jb.lt.istart) jb = jb + ninc
             jc = jb + jstepl
             if (jc.lt.istart) jc = jc + ninc
             j = 0
!
!  loop across transforms
!  ----------------------

!dir$ ivdep,shortloop
             do l = 1 , nvex
               DO jnew=0,JEND ! 02 FEB 2005 to vectorize
                 ajb = a(jnew,jb+j)
                 ajc = a(jnew,jc+j)
                 t1 = ajb + ajc
                 aja = a(jnew,ja+j)
                 t2 = aja - 0.5D0 * t1 ! REAL*8
                 t3 = c1 * ( ajb - ajc )
                 bjb = b(jnew,jb+j)
                 bjc = b(jnew,jc+j)
                 u1 = bjb + bjc
                 bja = b(jnew,ja+j)
                 u2 = bja - 0.5D0 * u1 ! REAL*8
                 u3 = c1 * ( bjb - bjc )
                 a(jnew,ja+j) = aja + t1
                 b(jnew,ja+j) = bja + u1
                 a(jnew,jb+j) = co1*(t2-u3) - si1*(u2+t3)
                 b(jnew,jb+j) = si1*(t2-u3) + co1*(u2+t3)
                 a(jnew,jc+j) = co2*(t2+u3) - si2*(u2-t3)
                 b(jnew,jc+j) = si2*(t2+u3) + co2*(u2-t3)
               END DO
               j = j + jump
             enddo
!
!-----( end of loop across transforms )
!
             ja = ja + jstepx
             if (ja.lt.istart) ja = ja + ninc
           enddo
         enddo
!
!-----( end of loop along transforms )
!
         kk = kk + 2*la
       enddo
!
!-----( end of loop on nonzero k )
!
       la = 3*la
     enddo
!
!-----( end of loop on type I radix-3 passes)
!
!  loop on type II radix-3 passes
!  ------------------------------
     400 continue
!
   do ipass = mh+1 , m
     jstep = (n*inc) / (3*la)
     jstepl = jstep - ninc
     laincl = la*ink - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
     do ll = 0 , (la-1)*ink , 3*jstep
!
       do jjj = ll , (n-1)*inc , 3*la*ink
         ja = istart + jjj
!
!  "transverse" loop
!  -----------------
         do nu = 1 , inq
           jb = ja + jstepl
           if (jb.lt.istart) jb = jb + ninc
           jc = jb + jstepl
           if (jc.lt.istart) jc = jc + ninc
           jd = ja + laincl
           if (jd.lt.istart) jd = jd + ninc
           je = jd + jstepl
           if (je.lt.istart) je = je + ninc
           jf = je + jstepl
           if (jf.lt.istart) jf = jf + ninc
           jg = jd + laincl
           if (jg.lt.istart) jg = jg + ninc
           jh = jg + jstepl
           if (jh.lt.istart) jh = jh + ninc
           ji = jh + jstepl
           if (ji.lt.istart) ji = ji + ninc
           j = 0
!
!  loop across transforms
!  ----------------------

!dir$ ivdep, shortloop
           do l = 1 , nvex
             DO jnew=0,JEND ! 02 FEB 2005 to vectorize
               ajb = a(jnew,jb+j)
               ajc = a(jnew,jc+j)
               t1 = ajb + ajc
               aja = a(jnew,ja+j)
               t2 = aja - 0.5D0 * t1 ! REAL*8
               t3 = c1 * ( ajb - ajc )
               ajd = a(jnew,jd+j)
               ajb =  ajd
               bjb = b(jnew,jb+j)
               bjc = b(jnew,jc+j)
               u1 = bjb + bjc
               bja = b(jnew,ja+j)
               u2 = bja - 0.5D0 * u1 ! REAL*8
               u3 = c1 * ( bjb - bjc )
               bjd = b(jnew,jd+j)
               bjb =  bjd
               a(jnew,ja+j) = aja + t1
               b(jnew,ja+j) = bja + u1
               a(jnew,jd+j) = t2 - u3
               b(jnew,jd+j) = u2 + t3
               ajc =  t2 + u3
               bjc =  u2 - t3
!
!----------------------
!
               aje = a(jnew,je+j)
               ajf = a(jnew,jf+j)
               t1 = aje + ajf
               t2 = ajb - 0.5D0 * t1 ! REAL*8
               t3 = c1 * ( aje - ajf )
               ajh = a(jnew,jh+j)
               ajf =  ajh
               bje = b(jnew,je+j)
               bjf = b(jnew,jf+j)
               u1 = bje + bjf
               u2 = bjb - 0.5D0 * u1 ! REAL*8
               u3 = c1 * ( bje - bjf )
               bjh = b(jnew,jh+j)
               bjf =  bjh
               a(jnew,jb+j) = ajb + t1
               b(jnew,jb+j) = bjb + u1
               a(jnew,je+j) = t2 - u3
               b(jnew,je+j) = u2 + t3
               a(jnew,jh+j) = t2 + u3
               b(jnew,jh+j) = u2 - t3
!
!----------------------
!
               aji = a(jnew,ji+j)
               t1 = ajf + aji
               ajg = a(jnew,jg+j)
               t2 = ajg - 0.5D0 * t1 ! REAL*8
               t3 = c1 * ( ajf - aji )
               t1 = ajg + t1
               a(jnew,jg+j) = ajc
               bji = b(jnew,ji+j)
               u1 = bjf + bji
               bjg = b(jnew,jg+j)
               u2 = bjg - 0.5D0 * u1 ! REAL*8
               u3 = c1 * ( bjf - bji )
               u1 = bjg + u1
               b(jnew,jg+j) = bjc
               a(jnew,jc+j) = t1
               b(jnew,jc+j) = u1
               a(jnew,jf+j) = t2 - u3
               b(jnew,jf+j) = u2 + t3
               a(jnew,ji+j) = t2 + u3
               b(jnew,ji+j) = u2 - t3
             END DO
             j = j + jump
           enddo
!
!-----( end of loop across transforms )
!
           ja = ja + jstepx
           if (ja.lt.istart) ja = ja + ninc
         enddo
       enddo
     enddo
!
!-----( end of double loop for k=0 )
!
!  finished if last pass
!  ---------------------
     if (ipass.eq.m) go to 490
!
     kk = 2*la
!
!     loop on nonzero k
!     -----------------
     do k = ink , jstep-ink , ink
       co1 = trigs(kk+1)
       si1 = s*trigs(kk+2)
       co2 = trigs(2*kk+1)
       si2 = s*trigs(2*kk+2)
!
!  double loop along first transform in block
!  ------------------------------------------
       do ll = k , (la-1)*ink , 3*jstep
!
         do jjj = ll , (n-1)*inc , 3*la*ink
           ja = istart + jjj
!
!  "transverse" loop
!  -----------------
           do nu = 1 , inq
             jb = ja + jstepl
             if (jb.lt.istart) jb = jb + ninc
             jc = jb + jstepl
             if (jc.lt.istart) jc = jc + ninc
             jd = ja + laincl
             if (jd.lt.istart) jd = jd + ninc
             je = jd + jstepl
             if (je.lt.istart) je = je + ninc
             jf = je + jstepl
             if (jf.lt.istart) jf = jf + ninc
             jg = jd + laincl
             if (jg.lt.istart) jg = jg + ninc
             jh = jg + jstepl
             if (jh.lt.istart) jh = jh + ninc
             ji = jh + jstepl
             if (ji.lt.istart) ji = ji + ninc
             j = 0
!
!  loop across transforms
!  ----------------------

!dir$ ivdep, shortloop
             do l = 1 , nvex
               DO jnew=0,JEND ! 02 FEB 2005 to vectorize
                 ajb = a(jnew,jb+j)
                 ajc = a(jnew,jc+j)
                 t1 = ajb + ajc
                 aja = a(jnew,ja+j)
                 t2 = aja - 0.5D0 * t1 ! REAL*8
                 t3 = c1 * ( ajb - ajc )
                 ajd = a(jnew,jd+j)
                 ajb =  ajd
                 bjb = b(jnew,jb+j)
                 bjc = b(jnew,jc+j)
                 u1 = bjb + bjc
                 bja = b(jnew,ja+j)
                 u2 = bja - 0.5D0 * u1 ! REAL*8
                 u3 = c1 * ( bjb - bjc )
                 bjd = b(jnew,jd+j)
                 bjb =  bjd
                 a(jnew,ja+j) = aja + t1
                 b(jnew,ja+j) = bja + u1
                 a(jnew,jd+j) = co1*(t2-u3) - si1*(u2+t3)
                 b(jnew,jd+j) = si1*(t2-u3) + co1*(u2+t3)
                 ajc =  co2*(t2+u3) - si2*(u2-t3)
                 bjc =  si2*(t2+u3) + co2*(u2-t3)
!
!----------------------
!
                 aje = a(jnew,je+j)
                 ajf = a(jnew,jf+j)
                 t1 = aje + ajf
                 t2 = ajb - 0.5D0 * t1 ! REAL*8
                 t3 = c1 * ( aje - ajf )
                 ajh = a(jnew,jh+j)
                 ajf =  ajh
                 bje = b(jnew,je+j)
                 bjf = b(jnew,jf+j)
                 u1 = bje + bjf
                 u2 = bjb - 0.5D0 * u1 ! REAL*8
                 u3 = c1 * ( bje - bjf )
                 bjh = b(jnew,jh+j)
                 bjf =  bjh
                 a(jnew,jb+j) = ajb + t1
                 b(jnew,jb+j) = bjb + u1
                 a(jnew,je+j) = co1*(t2-u3) - si1*(u2+t3)
                 b(jnew,je+j) = si1*(t2-u3) + co1*(u2+t3)
                 a(jnew,jh+j) = co2*(t2+u3) - si2*(u2-t3)
                 b(jnew,jh+j) = si2*(t2+u3) + co2*(u2-t3)
!
!----------------------
!
                 aji = a(jnew,ji+j)
                 t1 = ajf + aji
                 ajg = a(jnew,jg+j)
                 t2 = ajg - 0.5D0 * t1 ! REAL*8
                 t3 = c1 * ( ajf - aji )
                 t1 = ajg + t1
                 a(jnew,jg+j) = ajc
                 bji = b(jnew,ji+j)
                 u1 = bjf + bji
                 bjg = b(jnew,jg+j)
                 u2 = bjg - 0.5D0 * u1 ! REAL*8
                 u3 = c1 * ( bjf - bji )
                 u1 = bjg + u1
                 b(jnew,jg+j) = bjc
                 a(jnew,jc+j) = t1
                 b(jnew,jc+j) = u1
                 a(jnew,jf+j) = co1*(t2-u3) - si1*(u2+t3)
                 b(jnew,jf+j) = si1*(t2-u3) + co1*(u2+t3)
                 a(jnew,ji+j) = co2*(t2+u3) - si2*(u2-t3)
                 b(jnew,ji+j) = si2*(t2+u3) + co2*(u2-t3)
               END DO
               j = j + jump
             enddo
!
!-----(end of loop across transforms)
!
             ja = ja + jstepx
             if (ja.lt.istart) ja = ja + ninc
           enddo
         enddo
       enddo
!
!-----( end of double loop for this k )
!
       kk = kk + 2*la
     enddo
!
!-----( end of loop over values of k )
!
     la = 3*la
   enddo
!
!-----( end of loop on type II radix-3 passes )
!-----( nvex transforms completed)
!
   490 continue
   istart = istart + nvex * jump
   enddo
!
!-----( end of loop on blocks of transforms )
!
   return
   end subroutine dfs_fft_3
!
!-------------------------------------------------------------------------------
   subroutine dfs_fft_5                                                        &
              (a,b,trigs,inc,jump,n,mm,lot,isign, JDIM,JEND)
!-------------------------------------------------------------------------------
   IMPLICIT  REAL*8 (A-H, O-Z)
!-------------------------------------------------------------------------------
!
!     dfs_fft_5.FOR : REAL*8 VERSION ; Vectorized version (02 FEB 2005) 
!     GPFA5F_8.FOR : REAL*8 VERSION 
!     fortran version of *gpfa5* -
!     radix-5 section of self-sorting, in-place,
!        generalized pfa
!
!     N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. 
!     RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE 
!     (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      
!
!-------------------------------------------------------------------------------
!     dimension a(*), b(*), trigs(*)
   real  ::   a(0:JDIM,n), b(0:JDIM,n), trigs(*) ! 02 FEB 2005 to vectorize
   data sin36/0.587785252292473D0/, sin72/0.951056516295154D0/,                &
         qrt5/0.559016994374947D0/
   data lvr/64/
!-------------------------------------------------------------------------------
   n5 = 5 ** mm
   inq = n / n5
   jstepx = (n5-n) * inc
   ninc = n * inc
   ink = inc * inq
   mu = mod(inq,5)
   if (isign.eq.-1) mu = 5 - mu
!
   m = mm
   mh = (m+1)/2
   s = Dfloat(isign) ! REAL*8
   c1 = qrt5
   c2 = sin72
   c3 = sin36
   if (mu.eq.2.or.mu.eq.3) then
     c1 = -c1
     c2 = sin36
     c3 = sin72
   endif
   if (mu.eq.3.or.mu.eq.4) c2 = -c2
   if (mu.eq.2.or.mu.eq.4) c3 = -c3
!
   nblox = 1 + (lot-1)/lvr
   left = lot
   s = Dfloat(isign) ! REAL*8
   istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
   do nb = 1 , nblox
!
     if (left.le.lvr) then
       nvex = left
     else if (left.lt.(2*lvr)) then
       nvex = left/2
       nvex = nvex + mod(nvex,2)
     else
       nvex = lvr
     endif
     left = left - nvex
!
     la = 1
!
!  loop on type I radix-5 passes
!  -----------------------------
     do ipass = 1 , mh
       jstep = (n*inc) / (5*la)
       jstepl = jstep - ninc
       kk = 0
!
!  loop on k
!  ---------
       do k = 0 , jstep-ink , ink
!
       if (k.gt.0) then
         co1 = trigs(kk+1)
         si1 = s*trigs(kk+2)
         co2 = trigs(2*kk+1)
         si2 = s*trigs(2*kk+2)
         co3 = trigs(3*kk+1)
         si3 = s*trigs(3*kk+2)
         co4 = trigs(4*kk+1)
         si4 = s*trigs(4*kk+2)
       endif
!
!  loop along transform
!  --------------------
         do jjj = k , (n-1)*inc , 5*jstep
           ja = istart + jjj
!
!     "transverse" loop
!     -----------------
           do nu = 1 , inq
             jb = ja + jstepl
             if (jb.lt.istart) jb = jb + ninc
             jc = jb + jstepl
             if (jc.lt.istart) jc = jc + ninc
             jd = jc + jstepl
             if (jd.lt.istart) jd = jd + ninc
             je = jd + jstepl
             if (je.lt.istart) je = je + ninc
             j = 0
!
!  loop across transforms
!  ----------------------
             if (k.eq.0) then
!dir$ ivdep, shortloop
               do l = 1 , nvex
                 DO jnew=0,JEND ! 02 FEB 2005 to vectorize
                   ajb = a(jnew,jb+j)
                   aje = a(jnew,je+j)
                   t1 = ajb + aje
                   ajc = a(jnew,jc+j)
                   ajd = a(jnew,jd+j)
                   t2 = ajc + ajd
                   t3 = ajb - aje
                   t4 = ajc - ajd
                   t5 = t1 + t2
                   t6 = c1 * ( t1 - t2 )
                   aja = a(jnew,ja+j)
                   t7 = aja - 0.25D0 * t5 ! REAL*8
                   a(jnew,ja+j) = aja + t5
                   t8 = t7 + t6
                   t9 = t7 - t6
                   t10 = c3 * t3 - c2 * t4
                   t11 = c2 * t3 + c3 * t4
                   bjb = b(jnew,jb+j)
                   bje = b(jnew,je+j)
                   u1 = bjb + bje
                   bjc = b(jnew,jc+j)
                   bjd = b(jnew,jd+j)
                   u2 = bjc + bjd
                   u3 = bjb - bje
                   u4 = bjc - bjd
                   u5 = u1 + u2
                   u6 = c1 * ( u1 - u2 )
                   bja = b(jnew,ja+j)
                   u7 = bja - 0.25D0 * u5 ! REAL*8
                   b(jnew,ja+j) = bja + u5
                   u8 = u7 + u6
                   u9 = u7 - u6
                   u10 = c3 * u3 - c2 * u4
                   u11 = c2 * u3 + c3 * u4
                   a(jnew,jb+j) = t8 - u11
                   b(jnew,jb+j) = u8 + t11
                   a(jnew,je+j) = t8 + u11
                   b(jnew,je+j) = u8 - t11
                   a(jnew,jc+j) = t9 - u10
                   b(jnew,jc+j) = u9 + t10
                   a(jnew,jd+j) = t9 + u10
                   b(jnew,jd+j) = u9 - t10
                 END DO ! 02 FEB 2005
                 j = j + jump
               enddo

             else
!dir$ ivdep,shortloop
               do l = 1 , nvex
                 DO jnew=0,JEND ! 02 FEB 2005 to vectorize
                   ajb = a(jnew,jb+j)
                   aje = a(jnew,je+j)
                   t1 = ajb + aje
                   ajc = a(jnew,jc+j)
                   ajd = a(jnew,jd+j)
                   t2 = ajc + ajd
                   t3 = ajb - aje
                   t4 = ajc - ajd
                   t5 = t1 + t2
                   t6 = c1 * ( t1 - t2 )
                   aja = a(jnew,ja+j)
                   t7 = aja - 0.25D0 * t5 ! REAL*8
                   a(jnew,ja+j) = aja + t5
                   t8 = t7 + t6
                   t9 = t7 - t6
                   t10 = c3 * t3 - c2 * t4
                   t11 = c2 * t3 + c3 * t4
                   bjb = b(jnew,jb+j)
                   bje = b(jnew,je+j)
                   u1 = bjb + bje
                   bjc = b(jnew,jc+j)
                   bjd = b(jnew,jd+j)
                   u2 = bjc + bjd
                   u3 = bjb - bje
                   u4 = bjc - bjd
                   u5 = u1 + u2
                   u6 = c1 * ( u1 - u2 )
                   bja = b(jnew,ja+j)
                   u7 = bja - 0.25D0 * u5 ! REAL*8
                   b(jnew,ja+j) = bja + u5
                   u8 = u7 + u6
                   u9 = u7 - u6
                   u10 = c3 * u3 - c2 * u4
                   u11 = c2 * u3 + c3 * u4
                   a(jnew,jb+j) = co1*(t8-u11) - si1*(u8+t11)
                   b(jnew,jb+j) = si1*(t8-u11) + co1*(u8+t11)
                   a(jnew,je+j) = co4*(t8+u11) - si4*(u8-t11)
                   b(jnew,je+j) = si4*(t8+u11) + co4*(u8-t11)
                   a(jnew,jc+j) = co2*(t9-u10) - si2*(u9+t10)
                   b(jnew,jc+j) = si2*(t9-u10) + co2*(u9+t10)
                   a(jnew,jd+j) = co3*(t9+u10) - si3*(u9-t10)
                   b(jnew,jd+j) = si3*(t9+u10) + co3*(u9-t10)
                 END DO
                 j = j + jump
               enddo
! 
             endif
!
!-----( end of loop across transforms )
!
             ja = ja + jstepx
             if (ja.lt.istart) ja = ja + ninc
           enddo
         enddo
!
!-----( end of loop along transforms )
!
         kk = kk + 2*la
       enddo
!
!-----( end of loop on nonzero k )
!
       la = 5*la
     enddo
!
!-----( end of loop on type I radix-5 passes)
!
     if (n.eq.5) go to 490
!
!  loop on type II radix-5 passes
!  ------------------------------
  400 continue
!
     do ipass = mh+1 , m
       jstep = (n*inc) / (5*la)
       jstepl = jstep - ninc
       laincl = la * ink - ninc
       kk = 0
!
!     loop on k
!     ---------
       do k = 0 , jstep-ink , ink
!
         if (k.gt.0) then
           co1 = trigs(kk+1)
           si1 = s*trigs(kk+2)
           co2 = trigs(2*kk+1)
           si2 = s*trigs(2*kk+2)
           co3 = trigs(3*kk+1)
           si3 = s*trigs(3*kk+2)
           co4 = trigs(4*kk+1)
           si4 = s*trigs(4*kk+2)
         endif
!
!  double loop along first transform in block
!  ------------------------------------------
         do ll = k , (la-1)*ink , 5*jstep
!
           do jjj = ll , (n-1)*inc , 5*la*ink
   ja = istart + jjj
!
!     "transverse" loop
!     -----------------
             do nu = 1 , inq
               jb = ja + jstepl
               if (jb.lt.istart) jb = jb + ninc
               jc = jb + jstepl
               if (jc.lt.istart) jc = jc + ninc
               jd = jc + jstepl
               if (jd.lt.istart) jd = jd + ninc
               je = jd + jstepl
               if (je.lt.istart) je = je + ninc
               jf = ja + laincl
               if (jf.lt.istart) jf = jf + ninc
               jg = jf + jstepl
               if (jg.lt.istart) jg = jg + ninc
               jh = jg + jstepl
               if (jh.lt.istart) jh = jh + ninc
               ji = jh + jstepl
               if (ji.lt.istart) ji = ji + ninc
               jj = ji + jstepl
               if (jj.lt.istart) jj = jj + ninc
               jk = jf + laincl
               if (jk.lt.istart) jk = jk + ninc
               jl = jk + jstepl
               if (jl.lt.istart) jl = jl + ninc
               jm = jl + jstepl
               if (jm.lt.istart) jm = jm + ninc
               jn = jm + jstepl
               if (jn.lt.istart) jn = jn + ninc
               jo = jn + jstepl
               if (jo.lt.istart) jo = jo + ninc
               jp = jk + laincl
               if (jp.lt.istart) jp = jp + ninc
               jq = jp + jstepl
               if (jq.lt.istart) jq = jq + ninc
               jr = jq + jstepl
               if (jr.lt.istart) jr = jr + ninc
               js = jr + jstepl
               if (js.lt.istart) js = js + ninc
               jt = js + jstepl
               if (jt.lt.istart) jt = jt + ninc
               ju = jp + laincl
               if (ju.lt.istart) ju = ju + ninc
               jv = ju + jstepl
               if (jv.lt.istart) jv = jv + ninc
               jw = jv + jstepl
               if (jw.lt.istart) jw = jw + ninc
               jx = jw + jstepl
               if (jx.lt.istart) jx = jx + ninc
               jy = jx + jstepl
               if (jy.lt.istart) jy = jy + ninc
               j = 0
!
!  loop across transforms
!  ----------------------
               if (k.eq.0) then
!dir$ ivdep, shortloop
                 do l = 1 , nvex
                   DO jnew=0,JEND ! 02 FEB 2005 to vectorize
                     ajb = a(jnew,jb+j)
                     aje = a(jnew,je+j)
                     t1 = ajb + aje
                     ajc = a(jnew,jc+j)
                     ajd = a(jnew,jd+j)
                     t2 = ajc + ajd
                     t3 = ajb - aje
                     t4 = ajc - ajd
                     ajf = a(jnew,jf+j)
                     ajb =  ajf
                     t5 = t1 + t2
                     t6 = c1 * ( t1 - t2 )
                     aja = a(jnew,ja+j)
                     t7 = aja - 0.25D0 * t5 ! REAL*8
                     a(jnew,ja+j) = aja + t5
                     t8 = t7 + t6
                     t9 = t7 - t6
                     ajk = a(jnew,jk+j)
                     ajc =  ajk
                     t10 = c3 * t3 - c2 * t4
                     t11 = c2 * t3 + c3 * t4
                     bjb = b(jnew,jb+j)
                     bje = b(jnew,je+j)
                     u1 = bjb + bje
                     bjc = b(jnew,jc+j)
                     bjd = b(jnew,jd+j)
                     u2 = bjc + bjd
                     u3 = bjb - bje
                     u4 = bjc - bjd
                     bjf = b(jnew,jf+j)
                     bjb =  bjf
                     u5 = u1 + u2
                     u6 = c1 * ( u1 - u2 )
                     bja = b(jnew,ja+j)
                     u7 = bja - 0.25D0 * u5 ! REAL*8
                     b(jnew,ja+j) = bja + u5
                     u8 = u7 + u6
                     u9 = u7 - u6
                     bjk = b(jnew,jk+j)
                     bjc =  bjk
                     u10 = c3 * u3 - c2 * u4
                     u11 = c2 * u3 + c3 * u4
                     a(jnew,jf+j) = t8 - u11
                     b(jnew,jf+j) = u8 + t11
                     aje =  t8 + u11
                     bje =  u8 - t11
                     a(jnew,jk+j) = t9 - u10
                     b(jnew,jk+j) = u9 + t10
                     ajd =  t9 + u10
                     bjd =  u9 - t10
!
!----------------------
!
                     ajg = a(jnew,jg+j)
                     ajj = a(jnew,jj+j)
                     t1 = ajg + ajj
                     ajh = a(jnew,jh+j)
                     aji = a(jnew,ji+j)
                     t2 = ajh + aji
                     t3 = ajg - ajj
                     t4 = ajh - aji
                     ajl = a(jnew,jl+j)
                     ajh =  ajl
                     t5 = t1 + t2
                     t6 = c1 * ( t1 - t2 )
                     t7 = ajb - 0.25D0 * t5 ! REAL*8
                     a(jnew,jb+j) = ajb + t5
                     t8 = t7 + t6
                     t9 = t7 - t6
                     ajq = a(jnew,jq+j)
                     aji =  ajq
                     t10 = c3 * t3 - c2 * t4
                     t11 = c2 * t3 + c3 * t4
                     bjg = b(jnew,jg+j)
                     bjj = b(jnew,jj+j)
                     u1 = bjg + bjj
                     bjh = b(jnew,jh+j)
                     bji = b(jnew,ji+j)
                     u2 = bjh + bji
                     u3 = bjg - bjj
                     u4 = bjh - bji
                     bjl = b(jnew,jl+j)
                     bjh =  bjl
                     u5 = u1 + u2
                     u6 = c1 * ( u1 - u2 )
                     u7 = bjb - 0.25D0 * u5 ! REAL*8
                     b(jnew,jb+j) = bjb + u5
                     u8 = u7 + u6
                     u9 = u7 - u6
                     bjq = b(jnew,jq+j)
                     bji =  bjq
                     u10 = c3 * u3 - c2 * u4
                     u11 = c2 * u3 + c3 * u4
                     a(jnew,jg+j) = t8 - u11
                     b(jnew,jg+j) = u8 + t11
                     ajj =  t8 + u11
                     bjj =  u8 - t11
                     a(jnew,jl+j) = t9 - u10
                     b(jnew,jl+j) = u9 + t10
                     a(jnew,jq+j) = t9 + u10
                     b(jnew,jq+j) = u9 - t10
!
!----------------------
!
                     ajo = a(jnew,jo+j)
                     t1 = ajh + ajo
                     ajm = a(jnew,jm+j)
                     ajn = a(jnew,jn+j)
                     t2 = ajm + ajn
                     t3 = ajh - ajo
                     t4 = ajm - ajn
                     ajr = a(jnew,jr+j)
                     ajn =  ajr
                     t5 = t1 + t2
                     t6 = c1 * ( t1 - t2 )
                     t7 = ajc - 0.25D0 * t5 ! REAL*8
                     a(jnew,jc+j) = ajc + t5
                     t8 = t7 + t6
                     t9 = t7 - t6
                     ajw = a(jnew,jw+j)
                     ajo =  ajw
                     t10 = c3 * t3 - c2 * t4
                     t11 = c2 * t3 + c3 * t4
                     bjo = b(jnew,jo+j)
                     u1 = bjh + bjo
                     bjm = b(jnew,jm+j)
                     bjn = b(jnew,jn+j)
                     u2 = bjm + bjn
                     u3 = bjh - bjo
                     u4 = bjm - bjn
                     bjr = b(jnew,jr+j)
                     bjn =  bjr
                     u5 = u1 + u2
                     u6 = c1 * ( u1 - u2 )
                     u7 = bjc - 0.25D0 * u5 ! REAL*8
                     b(jnew,jc+j) = bjc + u5
                     u8 = u7 + u6
                     u9 = u7 - u6
                     bjw = b(jnew,jw+j)
                     bjo =  bjw
                     u10 = c3 * u3 - c2 * u4
                     u11 = c2 * u3 + c3 * u4
                     a(jnew,jh+j) = t8 - u11
                     b(jnew,jh+j) = u8 + t11
                     a(jnew,jw+j) = t8 + u11
                     b(jnew,jw+j) = u8 - t11
                     a(jnew,jm+j) = t9 - u10
                     b(jnew,jm+j) = u9 + t10
                     a(jnew,jr+j) = t9 + u10
                     b(jnew,jr+j) = u9 - t10
!
!----------------------
!
                     ajt = a(jnew,jt+j)
                     t1 = aji + ajt
                     ajs = a(jnew,js+j)
                     t2 = ajn + ajs
                     t3 = aji - ajt
                     t4 = ajn - ajs
                     ajx = a(jnew,jx+j)
                     ajt =  ajx
                     t5 = t1 + t2
                     t6 = c1 * ( t1 - t2 )
                     ajp = a(jnew,jp+j)
                     t7 = ajp - 0.25D0 * t5 ! REAL*8
                     ax = ajp + t5
                     t8 = t7 + t6
                     t9 = t7 - t6
                     a(jnew,jp+j) = ajd
                     t10 = c3 * t3 - c2 * t4
                     t11 = c2 * t3 + c3 * t4
                     a(jnew,jd+j) = ax
                     bjt = b(jnew,jt+j)
                     u1 = bji + bjt
                     bjs = b(jnew,js+j)
                     u2 = bjn + bjs
                     u3 = bji - bjt
                     u4 = bjn - bjs
                     bjx = b(jnew,jx+j)
                     bjt =  bjx
                     u5 = u1 + u2
                     u6 = c1 * ( u1 - u2 )
                     bjp = b(jnew,jp+j)
                     u7 = bjp - 0.25D0 * u5 ! REAL*8
                     bx = bjp + u5
                     u8 = u7 + u6
                     u9 = u7 - u6
                     b(jnew,jp+j) = bjd
                     u10 = c3 * u3 - c2 * u4
                     u11 = c2 * u3 + c3 * u4
                     b(jnew,jd+j) = bx
                     a(jnew,ji+j) = t8 - u11
                     b(jnew,ji+j) = u8 + t11
                     a(jnew,jx+j) = t8 + u11
                     b(jnew,jx+j) = u8 - t11
                     a(jnew,jn+j) = t9 - u10
                     b(jnew,jn+j) = u9 + t10
                     a(jnew,js+j) = t9 + u10
                     b(jnew,js+j) = u9 - t10
!
!----------------------
!
                     ajv = a(jnew,jv+j)
                     ajy = a(jnew,jy+j)
                     t1 = ajv + ajy
                     t2 = ajo + ajt
                     t3 = ajv - ajy
                     t4 = ajo - ajt
                     a(jnew,jv+j) = ajj
                     t5 = t1 + t2
                     t6 = c1 * ( t1 - t2 )
                     aju = a(jnew,ju+j)
                     t7 = aju - 0.25D0 * t5 ! REAL*8
                     ax = aju + t5
                     t8 = t7 + t6
                     t9 = t7 - t6
                     a(jnew,ju+j) = aje
                     t10 = c3 * t3 - c2 * t4
                     t11 = c2 * t3 + c3 * t4
                     a(jnew,je+j) = ax
                     bjv = b(jnew,jv+j)
                     bjy = b(jnew,jy+j)
                     u1 = bjv + bjy
                     u2 = bjo + bjt
                     u3 = bjv - bjy
                     u4 = bjo - bjt
                     b(jnew,jv+j) = bjj
                     u5 = u1 + u2
                     u6 = c1 * ( u1 - u2 )
                     bju = b(jnew,ju+j)
                     u7 = bju - 0.25D0 * u5 ! REAL*8
                     bx = bju + u5
                     u8 = u7 + u6
                     u9 = u7 - u6
                     b(jnew,ju+j) = bje
                     u10 = c3 * u3 - c2 * u4
                     u11 = c2 * u3 + c3 * u4
                     b(jnew,je+j) = bx
                     a(jnew,jj+j) = t8 - u11
                     b(jnew,jj+j) = u8 + t11
                     a(jnew,jy+j) = t8 + u11
                     b(jnew,jy+j) = u8 - t11
                     a(jnew,jo+j) = t9 - u10
                     b(jnew,jo+j) = u9 + t10
                     a(jnew,jt+j) = t9 + u10
                     b(jnew,jt+j) = u9 - t10
                   END DO
                   j = j + jump
                 enddo
!
               else
!dir$ ivdep, shortloop
                 do l = 1 , nvex
                   DO jnew=0,JEND ! 02 FEB 2005 to vectorize
                     ajb = a(jnew,jb+j)
                     aje = a(jnew,je+j)
                     t1 = ajb + aje
                     ajc = a(jnew,jc+j)
                     ajd = a(jnew,jd+j)
                     t2 = ajc + ajd
                     t3 = ajb - aje
                     t4 = ajc - ajd
                     ajf = a(jnew,jf+j)
                     ajb =  ajf
                     t5 = t1 + t2
                     t6 = c1 * ( t1 - t2 )
                     aja = a(jnew,ja+j)
                     t7 = aja - 0.25D0 * t5 ! REAL*8
                     a(jnew,ja+j) = aja + t5
                     t8 = t7 + t6
                     t9 = t7 - t6
                     ajk = a(jnew,jk+j)
                     ajc =  ajk
                     t10 = c3 * t3 - c2 * t4
                     t11 = c2 * t3 + c3 * t4
                     bjb = b(jnew,jb+j)
                     bje = b(jnew,je+j)
                     u1 = bjb + bje
                     bjc = b(jnew,jc+j)
                     bjd = b(jnew,jd+j)
                     u2 = bjc + bjd
                     u3 = bjb - bje
                     u4 = bjc - bjd
                     bjf = b(jnew,jf+j)
                     bjb =  bjf
                     u5 = u1 + u2
                     u6 = c1 * ( u1 - u2 )
                     bja = b(jnew,ja+j)
                     u7 = bja - 0.25D0 * u5 ! REAL*8
                     b(jnew,ja+j) = bja + u5
                     u8 = u7 + u6
                     u9 = u7 - u6
                     bjk = b(jnew,jk+j)
                     bjc =  bjk
                     u10 = c3 * u3 - c2 * u4
                     u11 = c2 * u3 + c3 * u4
                     a(jnew,jf+j) = co1*(t8-u11) - si1*(u8+t11)
                     b(jnew,jf+j) = si1*(t8-u11) + co1*(u8+t11)
                     aje =  co4*(t8+u11) - si4*(u8-t11)
                     bje =  si4*(t8+u11) + co4*(u8-t11)
                     a(jnew,jk+j) = co2*(t9-u10) - si2*(u9+t10)
                     b(jnew,jk+j) = si2*(t9-u10) + co2*(u9+t10)
                     ajd =  co3*(t9+u10) - si3*(u9-t10)
                     bjd =  si3*(t9+u10) + co3*(u9-t10)
!
!----------------------
!
                     ajg = a(jnew,jg+j)
                     ajj = a(jnew,jj+j)
                     t1 = ajg + ajj
                     ajh = a(jnew,jh+j)
                     aji = a(jnew,ji+j)
                     t2 = ajh + aji
                     t3 = ajg - ajj
                     t4 = ajh - aji
                     ajl = a(jnew,jl+j)
                     ajh =  ajl
                     t5 = t1 + t2
                     t6 = c1 * ( t1 - t2 )
                     t7 = ajb - 0.25D0 * t5 ! REAL*8
                     a(jnew,jb+j) = ajb + t5
                     t8 = t7 + t6
                     t9 = t7 - t6
                     ajq = a(jnew,jq+j)
                     aji =  ajq
                     t10 = c3 * t3 - c2 * t4
                     t11 = c2 * t3 + c3 * t4
                     bjg = b(jnew,jg+j)
                     bjj = b(jnew,jj+j)
                     u1 = bjg + bjj
                     bjh = b(jnew,jh+j)
                     bji = b(jnew,ji+j)
                     u2 = bjh + bji
                     u3 = bjg - bjj
                     u4 = bjh - bji
                     bjl = b(jnew,jl+j)
                     bjh =  bjl
                     u5 = u1 + u2
                     u6 = c1 * ( u1 - u2 )
                     u7 = bjb - 0.25D0 * u5 ! REAL*8
                     b(jnew,jb+j) = bjb + u5
                     u8 = u7 + u6
                     u9 = u7 - u6
                     bjq = b(jnew,jq+j)
                     bji =  bjq
                     u10 = c3 * u3 - c2 * u4
                     u11 = c2 * u3 + c3 * u4
                     a(jnew,jg+j) = co1*(t8-u11) - si1*(u8+t11)
                     b(jnew,jg+j) = si1*(t8-u11) + co1*(u8+t11)
                     ajj =  co4*(t8+u11) - si4*(u8-t11)
                     bjj =  si4*(t8+u11) + co4*(u8-t11)
                     a(jnew,jl+j) = co2*(t9-u10) - si2*(u9+t10)
                     b(jnew,jl+j) = si2*(t9-u10) + co2*(u9+t10)
                     a(jnew,jq+j) = co3*(t9+u10) - si3*(u9-t10)
                     b(jnew,jq+j) = si3*(t9+u10) + co3*(u9-t10)
!
!----------------------
!
                     ajo = a(jnew,jo+j)
                     t1 = ajh + ajo
                     ajm = a(jnew,jm+j)
                     ajn = a(jnew,jn+j)
                     t2 = ajm + ajn
                     t3 = ajh - ajo
                     t4 = ajm - ajn
                     ajr = a(jnew,jr+j)
                     ajn =  ajr
                     t5 = t1 + t2
                     t6 = c1 * ( t1 - t2 )
                     t7 = ajc - 0.25D0 * t5 ! REAL*8
                     a(jnew,jc+j) = ajc + t5
                     t8 = t7 + t6
                     t9 = t7 - t6
                     ajw = a(jnew,jw+j)
                     ajo =  ajw
                     t10 = c3 * t3 - c2 * t4
                     t11 = c2 * t3 + c3 * t4
                     bjo = b(jnew,jo+j)
                     u1 = bjh + bjo
                     bjm = b(jnew,jm+j)
                     bjn = b(jnew,jn+j)
                     u2 = bjm + bjn
                     u3 = bjh - bjo
                     u4 = bjm - bjn
                     bjr = b(jnew,jr+j)
                     bjn =  bjr
                     u5 = u1 + u2
                     u6 = c1 * ( u1 - u2 )
                     u7 = bjc - 0.25D0 * u5 ! REAL*8
                     b(jnew,jc+j) = bjc + u5
                     u8 = u7 + u6
                     u9 = u7 - u6
                     bjw = b(jnew,jw+j)
                     bjo =  bjw
                     u10 = c3 * u3 - c2 * u4
                     u11 = c2 * u3 + c3 * u4
                     a(jnew,jh+j) = co1*(t8-u11) - si1*(u8+t11)
                     b(jnew,jh+j) = si1*(t8-u11) + co1*(u8+t11)
                     a(jnew,jw+j) = co4*(t8+u11) - si4*(u8-t11)
                     b(jnew,jw+j) = si4*(t8+u11) + co4*(u8-t11)
                     a(jnew,jm+j) = co2*(t9-u10) - si2*(u9+t10)
                     b(jnew,jm+j) = si2*(t9-u10) + co2*(u9+t10)
                     a(jnew,jr+j) = co3*(t9+u10) - si3*(u9-t10)
                     b(jnew,jr+j) = si3*(t9+u10) + co3*(u9-t10)
!
!----------------------
!
                     ajt = a(jnew,jt+j)
                     t1 = aji + ajt
                     ajs = a(jnew,js+j)
                     t2 = ajn + ajs
                     t3 = aji - ajt
                     t4 = ajn - ajs
                     ajx = a(jnew,jx+j)
                     ajt =  ajx
                     t5 = t1 + t2
                     t6 = c1 * ( t1 - t2 )
                     ajp = a(jnew,jp+j)
                     t7 = ajp - 0.25D0 * t5 ! REAL*8
                     ax = ajp + t5
                     t8 = t7 + t6
                     t9 = t7 - t6
                     a(jnew,jp+j) = ajd
                     t10 = c3 * t3 - c2 * t4
                     t11 = c2 * t3 + c3 * t4
                     a(jnew,jd+j) = ax
                     bjt = b(jnew,jt+j)
                     u1 = bji + bjt
                     bjs = b(jnew,js+j)
                     u2 = bjn + bjs
                     u3 = bji - bjt
                     u4 = bjn - bjs
                     bjx = b(jnew,jx+j)
                     bjt =  bjx
                     u5 = u1 + u2
                     u6 = c1 * ( u1 - u2 )
                     bjp = b(jnew,jp+j)
                     u7 = bjp - 0.25D0 * u5 ! REAL*8
                     bx = bjp + u5
                     u8 = u7 + u6
                     u9 = u7 - u6
                     b(jnew,jp+j) = bjd
                     u10 = c3 * u3 - c2 * u4
                     u11 = c2 * u3 + c3 * u4
                     b(jnew,jd+j) = bx
                     a(jnew,ji+j) = co1*(t8-u11) - si1*(u8+t11)
                     b(jnew,ji+j) = si1*(t8-u11) + co1*(u8+t11)
                     a(jnew,jx+j) = co4*(t8+u11) - si4*(u8-t11)
                     b(jnew,jx+j) = si4*(t8+u11) + co4*(u8-t11)
                     a(jnew,jn+j) = co2*(t9-u10) - si2*(u9+t10)
                     b(jnew,jn+j) = si2*(t9-u10) + co2*(u9+t10)
                     a(jnew,js+j) = co3*(t9+u10) - si3*(u9-t10)
                     b(jnew,js+j) = si3*(t9+u10) + co3*(u9-t10)
!
!----------------------
!
                     ajv = a(jnew,jv+j)
                     ajy = a(jnew,jy+j)
                     t1 = ajv + ajy
                     t2 = ajo + ajt
                     t3 = ajv - ajy
                     t4 = ajo - ajt
                     a(jnew,jv+j) = ajj
                     t5 = t1 + t2
                     t6 = c1 * ( t1 - t2 )
                     aju = a(jnew,ju+j)
                     t7 = aju - 0.25D0 * t5 ! REAL*8
                     ax = aju + t5
                     t8 = t7 + t6
                     t9 = t7 - t6
                     a(jnew,ju+j) = aje
                     t10 = c3 * t3 - c2 * t4
                     t11 = c2 * t3 + c3 * t4
                     a(jnew,je+j) = ax
                     bjv = b(jnew,jv+j)
                     bjy = b(jnew,jy+j)
                     u1 = bjv + bjy
                     u2 = bjo + bjt
                     u3 = bjv - bjy
                     u4 = bjo - bjt
                     b(jnew,jv+j) = bjj
                     u5 = u1 + u2
                     u6 = c1 * ( u1 - u2 )
                     bju = b(jnew,ju+j)
                     u7 = bju - 0.25D0 * u5 ! REAL*8
                     bx = bju + u5
                     u8 = u7 + u6
                     u9 = u7 - u6
                     b(jnew,ju+j) = bje
                     u10 = c3 * u3 - c2 * u4
                     u11 = c2 * u3 + c3 * u4
                     b(jnew,je+j) = bx
                     a(jnew,jj+j) = co1*(t8-u11) - si1*(u8+t11)
                     b(jnew,jj+j) = si1*(t8-u11) + co1*(u8+t11)
                     a(jnew,jy+j) = co4*(t8+u11) - si4*(u8-t11)
                     b(jnew,jy+j) = si4*(t8+u11) + co4*(u8-t11)
                     a(jnew,jo+j) = co2*(t9-u10) - si2*(u9+t10)
                     b(jnew,jo+j) = si2*(t9-u10) + co2*(u9+t10)
                     a(jnew,jt+j) = co3*(t9+u10) - si3*(u9-t10)
                     b(jnew,jt+j) = si3*(t9+u10) + co3*(u9-t10)
                   END DO
                   j = j + jump
                 enddo
!
               endif
!
!-----(end of loop across transforms)
!
               ja = ja + jstepx
               if (ja.lt.istart) ja = ja + ninc
             enddo
           enddo
         enddo
!
!-----( end of double loop for this k )
!
         kk = kk + 2*la
       enddo
!
!-----( end of loop over values of k )
!
       la = 5*la
     enddo
!
!-----( end of loop on type II radix-5 passes )
!-----( nvex transforms completed)
!
     490 continue
     istart = istart + nvex * jump
   enddo
!
!-----( end of loop on blocks of transforms )
!
   return
   end subroutine dfs_fft_5
!
!-------------------------------------------------------------------------------
   SUBROUTINE dfs_fft_all(A,B,TRIGS,INC,JUMP,N,LOT,isgn, JDIM,JEND) ! 02 FEB 2005
!-------------------------------------------------------------------------------
!
!        SUBROUTINE 'dfs_fft_all' : REAL*8 VERSION ; Vectorized version
!        SUBROUTINE 'GPFA_8' : REAL*8 VERSION
!     <- SUBROUTINE 'GPFA'
!        SELF-SORTING IN-PLACE GENERALIZED PRIME FACTOR (COMPLEX) FFT
!
!        *** THIS IS THE ALL-FORTRAN VERSION **!
!            -------------------------------
!
!        CALL GPFA(A,B,TRIGS,INC,JUMP,N,LOT,isgn)
!
!        A IS FIRST REAL INPUT/OUTPUT VECTOR
!        B IS FIRST IMAGINARY INPUT/OUTPUT VECTOR
!        TRIGS IS A TABLE OF TWIDDLE FACTORS, PRECALCULATED
!              BY CALLING SUBROUTINE 'SETGPFA'
!        INC IS THE INCREMENT WITHIN EACH DATA VECTOR
!        JUMP IS THE INCREMENT BETWEEN DATA VECTORS
!        N IS THE LENGTH OF THE TRANSFORMS:
!          -----------------------------------
!            N = (2**IP) * (3**IQ) * (5**IR)
!          -----------------------------------
!        LOT IS THE NUMBER OF TRANSFORMS
!        isgn = +1 FOR FORWARD TRANSFORM
!              = -1 FOR INVERSE TRANSFORM
!
!        WRITTEN BY CLIVE TEMPERTON
!        RECHERCHE EN PREVISION NUMERIQUE
!        ATMOSPHERIC ENVIRONMENT SERVICE, CANADA
!
!----------------------------------------------------------------------
!
!        DEFINITION OF TRANSFORM
!        -----------------------
!
!        X(J) = SUM(K=0,...,N-1)(C(K)*EXP(isgn*2*I*J*K*PI/N))
!
!---------------------------------------------------------------------
!
!        FOR A MATHEMATICAL DEVELOPMENT OF THE ALGORITHM USED,
!        SEE:
!
!        C TEMPERTON : "A GENERALIZED PRIME FACTOR FFT ALGORITHM
!          FOR ANY N = (2**P)(3**Q)(5**R)",
!          SIAM J. SCI. STAT. COMP., MAY 1992.
!
!----------------------------------------------------------------------
!
!     SUBROUTINE GPFA_8(A,B,TRIGS,INC,JUMP,N,LOT,isgn)
!
!-------------------------------------------------------------------------------
!
!IMPLICIT  REAL*8 (A-H, O-Z)
!     DIMENSION A(*), B(*), TRIGS(*)
   real                 ::  A(0:JDIM,N), B(0:JDIM,N), TRIGS(*) ! 02 FEB 2005
   integer              ::  NJ(3)
!-------------------------------------------------------------------------------
!
!     DECOMPOSE N INTO FACTORS 2,3,5
!     ------------------------------
   NN = N
   IFAC = 2
!
   DO LL = 1,3
     KK = 0
  10 CONTINUE
     IF (MOD(NN,IFAC).NE.0) GO TO 20
     KK = KK + 1
     NN = NN / IFAC
     GO TO 10
  20 CONTINUE
     NJ(LL) = KK
     IFAC = IFAC + LL
   ENDDO
!
   IF (NN.NE.1) THEN
     WRITE(6,40) N
 40  FORMAT(' *** WARNING!!!',I10,' IS NOT A LEGAL VALUE OF N ***')
     RETURN
   ENDIF
!
   IP = NJ(1)
   IQ = NJ(2)
   IR = NJ(3)
!
!     COMPUTE THE TRANSFORM
!     ---------------------
   I = 1
   IF (IP.GT.0) THEN                                  
     CALL dfs_fft_2                                                            &
          (A,B,TRIGS,INC,JUMP,N,IP,LOT,isgn, JDIM,JEND)
     I = I + 2 * ( 2**IP)
   ENDIF
   IF (IQ.GT.0) THEN                                 
     CALL dfs_fft_3                                                            &
          (A,B,TRIGS(I),INC,JUMP,N,IQ,LOT,isgn, JDIM,JEND)
     I = I + 2 * (3**IQ)
   ENDIF
   IF (IR.GT.0) THEN
!
!--- 02 FEB 2005
!
     CALL dfs_fft_5                                                            &
         (A,B,TRIGS(I),INC,JUMP,N,IR,LOT,isgn, JDIM,JEND)
   ENDIF
!
   RETURN
   END SUBROUTINE dfs_fft_all
!-------------------------------------------------------------------------------
