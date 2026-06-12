#include "define.h"
   subroutine mpmz2yz(a,mt,latg,b,mtg,latp,ntotal)
!------------------------------------------------------------------------------
!
! subprogram documentation block
!
! subprogram:    mpnl2ny
!            
! usage:   call mpnl2ny(a,mt,latg,b,mtg,latp,ntotal)
!
!    input argument lists:
!   a   - real (mt,ntotal,latg) partial field in i k
!   mt   - integer partial longitude grid x 2
!   latg   - integer total latitude grid / 2
!   mtg   - integer total longitude grid x 2
!   latp   - integer partial latitude grid / 2
!   ntotal   - integer total set of fields
!
!    output argument list:
!   b   - real (mtg,ntotal,latp) partial field in j k
! 
! subprograms called:
!   mpi_alltoallv  - send and receive from all pe in the same comm
!
!------------------------------------------------------------------------------
   use commpi, only  :  mype,real_type,comm_col,nrow,ncol,npes,                &      
                        lmstr,lmlen,latstr,latlen,myrow
   use dfsvar, only  :  jgs,mgs,mge,latdef
!------------------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------------------
   integer  ::  mtg,latp,ntotal,mt,latg
   real     ::  a(mt,latg,ntotal),b(mgs:mge,latp,ntotal)
!
   integer  ::  i,j,k,jj,jpe,jpe0,n,mn,ierr,nr,js
   integer  ::  ms,me,lm,mem,msm,lmm,jst,jp,jp2
   integer  ::  mm(mgs:mge),jsp,jnp(mgs:mge),mp,m
!
   real(_mpi_real_), allocatable  ::  tmpsnd(:),tmprcv(:)
   integer, allocatable  ::  lensnd(:),lenrcv(:)
   integer, allocatable  ::  locsnd(:),locrcv(:)
!
   jp=latlen(mype)
   jp2=jp*2
   allocate(tmpsnd(mt*latg*ntotal))
   allocate(tmprcv(mtg*jp2*ntotal))
   allocate(lensnd(nrow))
   allocate(lenrcv(nrow))
   allocate(locsnd(nrow))
   allocate(locrcv(nrow))
!
   jpe0=mod(mype,ncol)
!
! cut in y
!
   mn=0
   do nr = 1,nrow
     locsnd(nr)=mn
     jpe=jpe0+(nr-1)*ncol
     do n = 1,ntotal
       do j = 1,latlen(jpe)
         jj=latdef(latstr(jpe)-jgs+j)
         js=latg-jj+1
         tmpsnd(mn+1:mn+mt)=a(1:mt,jj,n)
         mn=mn+mt
         tmpsnd(mn+1:mn+mt)=a(1:mt,js,n)
         mn=mn+mt
       enddo
     enddo
     lensnd(nr)=mn-locsnd(nr)
   enddo
!
   mn=0
   do nr = 1,nrow
     locrcv(nr)=mn
     jpe=jpe0+(nr-1)*ncol
     lenrcv(nr)=latlen(mype)*2*lmlen(jpe)*ntotal
     mn=mn+lenrcv(nr)
   enddo
!
   call mpi_alltoallv(tmpsnd,lensnd,locsnd,real_type,tmprcv,lenrcv,locrcv,     &
                      real_type,comm_col,ierr)
!
! restore l
!
   mn=0
   do nr = 1,nrow
     jpe=jpe0+(nr-1)*ncol
     ms=lmstr(jpe)
     lmm=lmlen(jpe)/2
     lm=lmm
     if (ms.eq.0) lm=lm+1
     me=ms+lm-1
     msm=-me
     mem=msm+lmm-1
     do n = 1,ntotal
       do j = 1,jp
         js=jp2-j+1
         jj=latdef(latstr(mype)+j-jgs)
         b(msm:mem,j ,n)=tmprcv(mn+1:mn+lmm)
         mn=mn+lmm
         b(ms :me ,j ,n)=tmprcv(mn+1:mn+lm )
         mn=mn+lm 
         b(msm:mem,js,n)=tmprcv(mn+1:mn+lmm)
         mn=mn+lmm
         b(ms :me ,js,n)=tmprcv(mn+1:mn+lm )
         mn=mn+lm 
#ifdef AAA
         if (n==1.and.nr.eq.nrow) then
           call rget(b(msm,j,n),mm,jnp,1)
           call rget(b(ms ,j,n),mp,jsp,1)
           write(96,'(A,7I5)')'hoon:j,js,msm,jnp,ms,jsp,jj=',                  &
                                        j,js,mm(1),jnp(1),mp,jsp,jj
           call rget(b(mgs,j,n),mm,jnp,mge-mgs+1)
           write(96,10)'in mpmz2yz',j,(mm(m),jnp(m),m=mgs,mge)
         endif
#endif
       enddo
     enddo
   enddo
!
10 format(A,:/'j=',I3,10(I5,I3),:/5x,10(I5,I3),:/5x,10(I5,I3),                 &
            :/5x     ,10(I5,I3),:/5x,10(I5,I3),:/5x,10(I5,I3),                 &
            :/5x     ,10(I5,I3),:/5x,10(I5,I3),:/5x,10(I5,I3),                 &
            :/5x     ,10(I5,I3),:/5x,10(I5,I3),:/5x,10(I5,I3),                 &
            :/5x     ,10(I5,I3),:/5x,10(I5,I3),:/5x,10(I5,I3))
!
   deallocate(tmpsnd)
   deallocate(tmprcv)
   deallocate(lensnd)
   deallocate(lenrcv)
   deallocate(locsnd)
   deallocate(locrcv)
!
   return
   end subroutine mpmz2yz
!-------------------------------------------------------------------------------
