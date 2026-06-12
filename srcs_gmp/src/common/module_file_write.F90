#include "define.h"
   module module_file_write
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [module_file_write]
!      |
!      |--- [file_write_bin] *
!      |--- [file_write_bin2] *
!      |--- [file_write_cps_bin] *
!      |--- [file_write_point] *
!      |--- [file_write_lfm] *
!
!-------------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------------
   subroutine file_write_bin(nn,var,il,jl,lot,ifl)
!-------------------------------------------------------------------------------
#ifdef MP
   use commpi
#endif
#ifdef DFS
   use dfsvar   , only : iope,iba,jbwa
#else
   use paramodel, only : LONF2A=>LONF2F,LATG2A=>LATG2F,LONF2S,LATG2S
   use comio    , only : iope
#ifdef RMP
   use module_trans, only  :  rmp_trans2output_grid 
#else
   use module_trans, only  :  dyn_trans2output_grid
#endif
#endif /* DFS end */
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#ifdef DFS
#define LONF2S il
#define LATG2S jl
#define LONF2A iba
#define LATG2A jbwa
#define MPGP2F mpgp2fd
#endif
#ifdef RMP
#define MPGP2F rmpgp2f
#endif
   integer,intent(in)   ::  il,jl,nn,lot,ifl
   real,intent(in)      ::  var(il*jl,lot)
!
! local variables
!
   integer              :: k
#ifdef MP
   real*4               ::  wor4(LONF2A*LATG2A)
   real                 ::  worka(LONF2A*LATG2A,lot)
#else
   real*4                   wor4(il*jl)
   real                 ::  work(il*jl,lot)
#ifndef DFS
   iope=.true.
#endif
#endif /* MP end */
!
#ifdef MP
   call MPGP2F(var,il,jl,worka,LONF2A,LATG2A,lot)
#define WORK worka
#else
   work(:,:)=var(:,:)
#endif /* MP end */
   if (iope) then
     do k = 1,lot
#ifndef DFS
       if (ifl.eq.1) then
#ifndef RMP
         call dyn_trans2output_grid(WORK(1,k),1)
#else
         call rmp_trans2output_grid(WORK(1,k),1)
#endif
       endif
#endif
       wor4(:)=WORK(:,k)
       write(nn)wor4
     enddo
   endif
#undef LONF2S
#undef LATG2S
#undef LONF2A
#undef LATG2A
#undef MPGP2F
#undef WORK
!
   return
   end subroutine file_write_bin
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine file_write_bin2(nn,var,il,jl,lot,ifl)
!-------------------------------------------------------------------------------
#if defined(DFS) && defined(MP)
   use dfsvar, only        :  iope,iba,jbwa
#endif
   use paramodel, only     :  LONF2S,LATG2S,LONF2A=>LONF2F,LATG2A=>LATG2F
   use comio
#ifdef MP
   use commpi
#endif
#ifdef RMP
   use module_trans, only  :  rmp_trans2output_grid 
#else
   use module_trans, only  :  dyn_trans2output_grid
#endif
!-------------------------------------------------------------------------------
#ifdef DFS
#define LONF2S il
#define LATG2S jl
#define LONF2A iba
#define LATG2A jbwa
#define MPGP2F mpgp2fd
#endif
#ifdef RMP
#define MPGP2F rmpgp2f
#endif
   integer,intent(in)   ::  il,jl,nn,lot
   real,intent(in)      ::  var(il*jl,lot)
!
! local variables
!
#ifdef MP
   real*4               ::  wor4(LONF2A*LATG2A)
   real                 ::  worka(LONF2A*LATG2A,lot)
#else
   real*4               ::  wor4(il*jl)
   real                 ::  work(il*jl,lot)
#ifdef DFS
   logical              ::  iope
#endif
!
   iope=.true.
#endif
!
#ifdef MP
   call MPGP2F(var,LONF2S,LATG2S,worka,LONF2A,LATG2A,lot)
#define WORK worka
#else
   work(:,:)=var(:,:)
#endif
   if (iope) then
     if (ifl.eq.1) then
#ifndef RMP
       call dyn_trans2output_grid(WORK,lot)
#else
       call rmp_trans2output_grid(WORK,lot)
#endif
     endif
     do k = 1,lot
       wor4(:)=WORK(:,k)
       write(nn)wor4
     enddo
   endif
!
#undef LONF2S
#undef LATG2S
#undef LONF2A
#undef LATG2A
#undef MPGP2F
#undef WORK
   return
   end subroutine file_write_bin2
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine file_write_cps_bin(thour,ig,jg,kg,                               &
                       dcu,dcv,dct,dcq,dch,fcu,fcd,                            &
                       deltb,delqb,delhb,cbmf,cbp,ctp,                         &
                       dlt,dlq,dlh)
!-------------------------------------------------------------------------------
   use paramodel, only : levs_
!-------------------------------------------------------------------------------
!
!     to write diagnostic data for various experiments
!     write format (ieee, direct access for grads plot) 
!
!     y.-h. byun               15 Jun 2004
!
!-------------------------------------------------------------------------------
   integer,intent(in)   ::  ig,jg,kg
   real,intent(in)      ::  thour
   real,intent(in)      ::  dcu(ig,kg,jg),dcv(ig,kg,jg),dct(ig,kg,jg)
   real,intent(in)      ::  dcq(ig,kg,jg),dch(ig,kg,jg)
   real,intent(in)      ::  fcu(ig,kg,jg),fcd(ig,kg,jg)
   real,intent(in)      ::  dlt(ig,kg,jg),dlq(ig,kg,jg),dlh(ig,kg,jg)
   real,intent(in)      ::  cbp(ig,jg),   ctp(ig,jg),   cbmf(ig,jg)
   real,intent(in)      ::  deltb(ig,jg), delqb(ig,jg), delhb(ig,jg)
!
   character(len=80)    ::  fno
   real                 ::  dum(ig,jg,kg)
!
   call file_name('sasdiag',7,thour,fno,ncho)
   close(91)
   open(91,file=fno(1:ncho),form='unformatted')
!
   do k = 1,kg
     do j = 1,jg
       dum(:,j,k)=dcu(:,k,j)
     enddo
   enddo
   call file_write_bin(91,dum,ig,jg,levs_,1)
!
   do k = 1,kg
     do j = 1,jg
       dum(:,j,k)=dcv(:,k,j)
     enddo
   enddo
   call file_write_bin(91,dum,ig,jg,levs_,1)
!
   do k = 1,kg
     do j = 1,jg
       dum(:,j,k)=dct(:,k,j)
     enddo
   enddo
   call file_write_bin(91,dum,ig,jg,levs_,1)
!
   do k = 1,kg
     do j = 1,jg
       dum(:,j,k)=dcq(:,k,j)
     enddo
   enddo
   call file_write_bin(91,dum,ig,jg,levs_,1)
!
   do k = 1,kg
     do j = 1,jg
       dum(:,j,k)=dch(:,k,j)
     enddo
   enddo
   call file_write_bin(91,dum,ig,jg,levs_,1)
!
   do k = 1,kg
     do j = 1,jg
       dum(:,j,k)=fcu(:,k,j)
     enddo
   enddo
   call file_write_bin(91,dum,ig,jg,levs_,1)
!
   do k= 1 ,kg
     do j = 1,jg
       dum(:,j,k)=fcd(:,k,j)
     enddo
   enddo
   call file_write_bin(91,dum,ig,jg,levs_,1)
!
   do k = 1,kg
     do j = 1,jg
       dum(:,j,k)=dlt(:,k,j)
     enddo
   enddo
   call file_write_bin(91,dum,ig,jg,levs_,1)
!
   do k = 1,kg
     do j = 1,jg
       dum(:,j,k)=dlq(:,k,j)
     enddo
   enddo
   call file_write_bin(91,dum,ig,jg,levs_,1)
!
   do k = 1,kg
     do j = 1,jg
       dum(:,j,k)=dlh(:,k,j)
     enddo
   enddo
   call file_write_bin(91,dum,ig,jg,levs_,1)
!
   call file_write_bin(91,deltb,ig,jg,1,1)
   call file_write_bin(91,delqb,ig,jg,1,1)
   call file_write_bin(91,delhb,ig,jg,1,1)
   call file_write_bin(91,cbp,ig,jg,1,1)
   call file_write_bin(91,ctp,ig,jg,1,1)
   call file_write_bin(91,cbmf,ig,jg,1,1)
   close(91)
!
   return
   end subroutine file_write_cps_bin
!-------------------------------------------------------------------------------
!
#define LEVP1 levp1
#define LEVS levs
#define MLVARK mlvark
#define SLVARK slvark
!-------------------------------------------------------------------------------
   subroutine file_write_point(npoint,ikfreq,imodk,itnum,svdata,               &
                     lab,fhour,idate,si,sl,nvrken,nptken,nstken,nn)
!-------------------------------------------------------------------------------
   use paramodel, only : levp1_,levs_
!-------------------------------------------------------------------------------
#ifdef DGP
   character(len=8)    ::  lab
   real                ::  svdata(nvrken,nptken,nstken)
   real                ::  si(levp1_),sl(levs_)
   integer             ::  idate(4)
!
   character(len=128)   ::  fno
   if(npoint.le.0) return
!
   call file_name('DGP',nchi,fhour,fno,ncho)
#ifdef ASSIGN
   call assign('assign -R')
#endif
   open(unit=nn,file=fno(1:ncho),form='unformatted',err=900)
!
   go to 901
!
900 continue
   write(6,*) ' error in opening file ',fno(1:ncho)
#ifdef MP
#ifdef RMP
   call rmpabort
#else
   call mpabort
#endif
#else
   call abort
#endif
901 continue
#ifndef NOPRINT
   write(6,*) ' file ',fno(1:ncho),' opened. unit=',nn
#endif
!
!  note, that in several scenerios, itnum=itnum+1 at the
!  bottom of the 10000 loop, so undo it
!
   if (ikfreq.gt.1) then
     if (imodk.le.0) then
       itnum = itnum - 1
     endif
   else
     itnum = itnum - 1
   endif
#ifndef NOPRINT
   print 1047,itnum,npoint
1047 format(1h0,i6,' steps of DGP gridpt data saved for ',            &
               i5,' points')
#endif
   do j = 1,npoint
     do k = 1,itnum
       if (svdata(44,j,k).le.0.) then
         do i = 25,27
           svdata(i,j,k) = svdata(i,j,k-1)
         enddo
         do i = 41,46
           svdata(i,j,k) = svdata(i,j,k-1)
         enddo
         do i = 48,49
           svdata(i,j,k) = svdata(i,j,k-1)
         enddo
         do i = 51,58
           svdata(i,j,k) = svdata(i,j,k-1)
         enddo
         i=38
         svdata(i,j,k) = svdata(i,j,k-1)
         do lv = 1,levs_
           svdata(lv+slvark_+(mlvark_-1)*levs_,j,k)=                           &
                  svdata(lv+slvark_+(mlvark_-1)*levs_,j,k-1)
         enddo
       endif
       if (svdata(50,j,k).le.0.) then
         i = 47
         svdata(i,j,k) = svdata(i,j,k-1)
         i = 50
         svdata(i,j,k) = svdata(i,j,k-1)
       endif
       if (svdata(70,j,k).le.0.) then
         do i = 70, 74
           svdata(i,j,k) = svdata(i,j,k-1)
         enddo
       endif
     enddo
   enddo
!
   rewind nn
   write(nn) lab
   write(nn) fhour,idate,si,sl
   write(nn) nvrken,nptken,nstken,npoint,itnum
!
   do j = 1,npoint
     write(nn) ((svdata(i,j,k),k=1,itnum),i=1,nvrken)
   enddo
   close(nn)
#endif
#undef LEVP1
#undef LEVS
#undef MLVARK
#undef SLVARK
!
   return
   end subroutine file_write_point
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine file_write_lfm(ifstep,thour,nz)
!-------------------------------------------------------------------------------
   use paramodel, only     :  levp1_,levs_,lnt2_,lonf2_,latg2_ 
   use varsfc, only        :  lsoil_,lalbd_
   use comio
   use comlfm
   use module_trans, only  :  dyn_trans2output_grid
!-------------------------------------------------------------------------------
   real                   ::  si(levp1_),sl(levs_)
   integer                ::  idate(4)
   real                   ::  gz(lnt2_)
!  real                   ::  fslmsk(lonf2_,latg2_)
!  real                   ::  fcvb(lonf2_,latg2_),fcvt(lonf2_,latg2_)
   real                   ::  work(lonf2_*latg2_)
!-------------------------------------------------------------------------------
   rewind nz
   read(nz)
   read(nz) dummy,idate,(si(k),k=1,levp1_),(sl(k),k=1,levs_ )
#ifndef NOPRINT
   write(6,*) 'file_write_lfm ifstep=',ifstep,' klenp=',klenp
#endif
   if(ifstep.eq.klenp) then
     vhour=nint(thour-filtwin*0.5)
#ifndef NOPRINT
     write(6,*) 'vhour,idate of filtered output sig=',vhour,idate
#endif
     write(nlfmsgo) lab
     write(nlfmsgo) vhour,idate,(si(k),k=1,levp1_),(sl(k),k=1,levs_)
     read(nz)(gz(i),i=1,lnt2_ )
     write(nlfmsgo)(gz(i),i=1,lnt2_ )
     write(nlfmsgo)(fq(i),i=1,lnt2_ )
     do k = 1,levs_
       write(nlfmsgo) (fte(i,k),i=1,lnt2_ )
     enddo
     do k = 1,levs_
       write(nlfmsgo) (fdi(i,k),i=1,lnt2_ )
       write(nlfmsgo) (fze(i,k),i=1,lnt2_ )
     enddo
     do k = 1,levs_
       write(nlfmsgo) (frq(i,k),i=1,lnt2_ )
     enddo
!
     write(nlfmsfo) lab
     write(nlfmsfo) vhour,idate
!
     call dyn_trans2output_grid(ftsea,1)
!
     write(nlfmsfo) ftsea
     do k = 1,lsoil_
       do j = 1, latg2_
         do i = 1, lonf2_
           ij = (j-1) * lonf2_  + i
           work(ij) = fsmc(i,j,k)
         enddo
       enddo
!
       call dyn_trans2output_grid(work,1)
!
       do j = 1, latg2_
         do i = 1, lonf2_
           ij = (j-1) * lonf2_  + i
           fsmc(i,j,k) = work(ij)
         enddo
       enddo
     enddo
     write(nlfmsfo) fsmc
!
     call dyn_trans2output_grid(fsnoweq,1)
!
     write(nlfmsfo) fsnoweq
     do k = 1,lsoil_
       do j = 1,latg2_
         do i = 1,lonf2_
           ij = (j-1) * lonf2_  + i
           work(ij) = fstc(i,j,k)
         enddo
       enddo
!
       call dyn_trans2output_grid(work,1)
!
       do j = 1,latg2_
         do i = 1,lonf2_
           ij = (j-1) * lonf2_  + i
           fstc(i,j,k) = work(ij)
         enddo
       enddo
     enddo
     write(nlfmsfo) fstc
!
     call dyn_trans2output_grid(ftg3,1)
!
     write(nlfmsfo) ftg3
!
     call dyn_trans2output_grid(fz0cm,1)
!
     write(nlfmsfo) fz0cm
!
     call dyn_trans2output_grid(fcv,1)
!
     write(nlfmsfo) fcv
     do j = 1,latg2_
       do i = 1,lonf2_
         if(wcvb(i,j).ge.0.5) then
           fcvb(i,j)=fcvb(i,j)/wcvb(i,j)
         endif
       enddo
     enddo
!
     call dyn_trans2output_grid(fcvb,1)
!
     write(nlfmsfo) fcvb
     do j = 1,latg2_
       do i = 1,lonf2_
         if(wcvt(i,j).ge.0.5) then
           fcvt(i,j)=fcvt(i,j)/wcvt(i,j)
         endif
       enddo
     enddo
!
     call dyn_trans2output_grid(fcvt,1)
!
     write(nlfmsfo) fcvt
!
     call dyn_trans2output_grid(falbedo,1)
!
     write(nlfmsfo) falbedo
     do j = 1,latg2_
       do i = 1,lonf2_
         if(islmsk(i,j,3).ge.klenp/2) then
           fslmsk(i,j)=2.
         elseif(islmsk(i,j,1).gt.islmsk(i,j,2)) then
           fslmsk(i,j)=0.
         else
           fslmsk(i,j)=1.
         endif
       enddo
     enddo
!
     call dyn_trans2output_grid(fslmsk,1)
!
     write(nlfmsfo) fslmsk
!
     call dyn_trans2output_grid(fplantr,1)
!
     write(nlfmsfo) fplantr
!
     call dyn_trans2output_grid(fcanopy,1)
!
     write(nlfmsfo) fcanopy
!
     call dyn_trans2output_grid(ff10m,1)
!
     write(nlfmsfo) ff10m
   else
     write(nlfmsgo) lab
     write(nlfmsgo) ifstep,idate,(si(k),k=1,levp1_),(sl(k),k=1,levs_)
#ifndef NOPRINT
     write(6,*) 'ifstep,idate of filtered sig=',ifstep,idate
#endif
     write(nlfmsgo)(gz(i),i=1,lnt2_ )
     write(nlfmsgo)(fq(i),i=1,lnt2_ )
     do k = 1,levs_
       write(nlfmsgo) (fte(i,k),i=1,lnt2_ )
     enddo
     do k = 1,levs_
       write(nlfmsgo) (fdi(i,k),i=1,lnt2_ )
       write(nlfmsgo) (fze(i,k),i=1,lnt2_ )
     enddo
     do k = 1,levs_
       write(nlfmsgo) (frq(i,k),i=1,lnt2_ )
     enddo
!
!  surface file
!
     write(nlfmsfo) lab
     write(nlfmsfo) ifstep,idate
#ifndef NOPRINT
     write(6,*) 'ifstep,idate of filtered sfc=',ifstep,idate
#endif
     write(nlfmsfo) ftsea
     write(nlfmsfo) fsmc
     write(nlfmsfo) fsnoweq
     write(nlfmsfo) fstc
     write(nlfmsfo) ftg3
     write(nlfmsfo) fz0cm
     write(nlfmsfo) fcv
     write(nlfmsfo) fcvb,wcvb
     write(nlfmsfo) fcvt,wcvt
     write(nlfmsfo) falbedo
     write(nlfmsfo) islmsk
     write(nlfmsfo) fplantr
     write(nlfmsfo) fcanopy
     write(nlfmsfo) ff10m
   endif
!
   return
   end subroutine file_write_lfm
!-------------------------------------------------------------------------------
   end module module_file_write
