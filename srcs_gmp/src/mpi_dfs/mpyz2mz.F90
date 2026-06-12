#include "define.h"
   subroutine mpyz2mz(a,mtg,latp,latpd,b,mt,latg,ntotal)
!-------------------------------------------------------------------------------
!
! Purpose : Rearrange variables after zonal FFT to gather global 
!           latitudal grid for Meridional Wave to Grid transform.
!           (mtg,jbw,k)->(mt,jbwa,k)
! Input variables
!
!-------------------------------------------------------------------------------
   use commpi, only  :  mype,real_type,comm_col,nrow,ncol,lmstr,lmlen,         &
                        latstr,latlen
   use dfsvar, only  :  jgs,mgs,mge,midxr,jbwa,latdef
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer  ::  mtg,latpd,latp,mt,latg,ntotal
   real     ::  a(mgs:mge,latpd,ntotal),b(mt,latg,ntotal)
!
   integer  ::  i,j,nr,jpe,jpe0,n,mn,len,ierr,llens,m,jj,jp,jbw,js
   integer, pointer  ::  mp(:)
!
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer, allocatable  ::  lensnd(:),lenrcv(:)
   integer, allocatable  ::  locsnd(:),locrcv(:)
!
   jp=latlen(mype)
   allocate(tmpsnd(mtg*jp*2*ntotal))
   allocate(tmprcv(mt*latg*ntotal))
   allocate(lensnd(nrow))
   allocate(lenrcv(nrow))
   allocate(locsnd(nrow))
   allocate(locrcv(nrow))
!
   jpe0=mod(mype,ncol)
!
! cut in l
!
   mn=0
   do nr = 1,nrow
     locsnd(nr)=mn
     jpe=jpe0+(nr-1)*ncol
     llens=lmlen(jpe)
     mp=>midxr(1:llens,nr-1)
     do n = 1,ntotal
       do j = 1,jp
         jj=jp*2-j+1
         do m = 1,llens
           tmpsnd(mn+m)=a(mp(m),j ,n)
         enddo
         mn=mn+llens
         do m = 1,llens
           tmpsnd(mn+m)=a(mp(m),jj,n)
         enddo
         mn=mn+llens
       enddo
     enddo
     lensnd(nr)=mn-locsnd(nr)
   enddo
!
   mn=0
   do nr = 1,nrow
     locrcv(nr)=mn
     jpe=jpe0+(nr-1)*ncol
     lenrcv(nr)=latlen(jpe)*2*lmlen(mype)*ntotal
     mn=mn+lenrcv(nr)
   enddo
!
   call mpi_alltoallv(tmpsnd,lensnd,locsnd,real_type,tmprcv,lenrcv,locrcv,     &
                      real_type,comm_col,ierr)
!
! put to y
!
   mn=0
   do nr = 0,nrow-1
     jpe=jpe0+nr*ncol
     do n = 1,ntotal
       do j = 1,latlen(jpe)
         jj=latdef(latstr(jpe)+j-jgs)
         js=latg-jj+1
         b(1:mt,jj,n)=tmprcv(mn+1:mn+mt)
         mn=mn+mt
         b(1:mt,js,n)=tmprcv(mn+1:mn+mt)
         mn=mn+mt
       enddo
     enddo
   enddo
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(lensnd)
   deallocate(lenrcv)
   deallocate(locsnd)
   deallocate(locrcv)
!
   return
   end subroutine mpyz2mz
!-------------------------------------------------------------------------------
