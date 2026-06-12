#include <define.h>
   subroutine sfc_read_file(nsfc,svar,kdim,fld,                                &
#ifdef SMP_RA2SFC
                       idim,jdim,ioflag)
#else
                       ioflag)
#endif
!-------------------------------------------------------------------------------
!
! subroutine: sfc_read_file
!
!  fix field record i/o for non-mpi and mpi
!
!  nsfc .... integer unit number of surface file
!  svar .... character identifying the field (for print only)
!  kdim .... integer number of soil layers (or classifications, e.g. albedo)
!  fld  .... real array to write or read (full or partial depending on the MP)
!  ioflag .. =0 read, =1 write
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use paramodel, only    : LONF2S,LATG2S,LONF2F,LATG2F
#if defined(MP)
   use commpi, only       : mype,master
#endif
#ifndef RMP
   use module_trans, only : dyn_trans2model_grid, dyn_trans2output_grid
#else
   use module_trans, only : rmp_trans2model_grid, rmp_trans2output_grid
#endif
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#include <abort.h>
#if defined(MP)
#ifndef RMP
#define MPGF2P mpgf2p
#define MPGP2F mpgp2f
#else
#define MPGF2P rmpgf2p
#define MPGP2F rmpgp2f
#endif
#endif
!
   character(len=8)     ::  svar  
   integer              ::  nsfc,kdim,ioflag
#ifdef SMP_RA2SFC
   integer              ::  idim,jdim
   real                 ::  fld(idim,jdim,kdim)
#else
   real                 ::  fld(LONF2S,LATG2S,kdim)
#endif
   integer              ::  i,j,k
#ifdef MP
   real, allocatable    :: grid(:,:,:)
#endif
!
#ifdef MP
   allocate (grid(LONF2F,LATG2F,kdim))
#endif
!
   if(ioflag.eq.1) go to 1000
!
! read record
!
#ifdef MP
   if( mype.eq.master ) then
     read(nsfc,err=101) (((grid(i,j,k),i=1,LONF2F),j=1,LATG2F),k=1,kdim)
     !
     go to 100
     !
101  continue
     print *,'sfc file read error for var=',svar
     call MPABORT
100  continue
#ifdef DBG
!
     do k = 1,kdim
       call print_maxmin_six(grid(1,1,k),LONF2F*LATG2F,1,1,1,svar)
     enddo
!
#endif
#ifndef RMP
     call dyn_trans2model_grid(grid,kdim)
#else
     call rmp_trans2model_grid(grid,kdim)
#endif
   endif
   call MPGF2P(grid,LONF2F,LATG2F,fld,LONF2S,LATG2S,kdim)
#else         /* not MP */
   read(nsfc,err=111)                                                          &
#ifdef SMP_RA2SFC
           (((fld(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)
#else
           (((fld(i,j,k),i=1,LONF2F),j=1,LATG2F),k=1,kdim)
#endif
   !
   go to 110
   !
111 continue
   print *,'sfc file read error for var=',svar
   call MPABORT
110 continue
#ifdef DBG
!
   do k = 1,kdim
     call print_maxmin_six(fld(1,1,k),LONF2F*LATG2F,1,1,1,svar)
   enddo
!
#endif
#ifndef RMP
   call dyn_trans2model_grid(fld,kdim)
#else
   call rmp_trans2model_grid(fld,kdim)
#endif
#ifndef NOPRINT
   print *,' fixrdrec completed for ',svar
#endif
#endif
#ifdef MP
   deallocate (grid)
#endif
   return
!--------------------------------------------------------------------
!
! write field
!
1000 continue
#ifdef MP
   call MPGP2F(fld,LONF2S,LATG2S,grid,LONF2F,LATG2F,kdim)
   if( mype.eq.master ) then
#ifndef RMP
     call dyn_trans2output_grid(grid,kdim)
#else
     call rmp_trans2output_grid(grid,kdim)
#endif
     write(nsfc) (((grid(i,j,k),i=1,LONF2F),j=1,LATG2F),k=1,kdim)
   endif
#else
#ifndef RMP
   call dyn_trans2output_grid(fld,kdim)
#else
   call rmp_trans2output_grid(fld,kdim)
#endif
   write(nsfc) (((fld(i,j,k),i=1,LONF2F),j=1,LATG2F),k=1,kdim)
#ifndef RMP
   call dyn_trans2model_grid(fld,kdim)
#else
   call rmp_trans2model_grid(fld,kdim)
#endif
#endif
!
#ifdef DBG
#ifdef MP
!
   if( mype.eq.master ) then
#endif
     write(6,*)' sfc_read_file completed for ',svar
     do k = 1,kdim
#ifdef MP
       call print_maxmin_six(grid(1,1,k),LONF2F*LATG2F,1,1,1,svar)
#else
       call print_maxmin_six(fld(1,1,k),LONF2F*LATG2F,1,1,1,svar)
#endif
     enddo
#ifdef MP
   endif
!
#endif
#endif
#ifdef MP
   deallocate (grid)
#endif
!
   return
   end subroutine sfc_read_file
!-------------------------------------------------------------------------------
