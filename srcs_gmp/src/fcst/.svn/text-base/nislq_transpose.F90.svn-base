#include <define.h>
   subroutine nislq_transpose_we2ns(a,b,lev,nsize)
#ifdef NISLQ
!-------------------------------------------------------------------------------
!
! mpi transport from full dimension of west-east to full dimension of 
! north-south with latitude shuffl.
!
!-------------------------------------------------------------------------------
   use nislq, only : lonfull , latfull                                        ,&
                     latpart , lonpart                                        ,&
                     lonhalf , lathalf                                        ,&
                     mylonlen, mylatlen                                       ,&
                     lonstr  , lonlen                                         ,&
                     latstr  , latlen                                         ,&
                     truej   , shflj
#ifdef MP
#ifdef DFS
   use commpi, only : COMM_PE=>comm_col   , MPI_REAL=>real_type
#else
   use commpi, only : COMM_PE=>comm_column, MPI_REAL
#endif
#endif /* MP end */
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                                     ::  lev,nsize                  ,&
                                                   i,j,k,n,mn,jj,lat1,lat2,ierr
   real   , dimension(lonfull,lev,latpart)     ::  a
   real   , dimension(latfull,lev,lonpart)     ::  b
   real   , dimension(2*lev*mylatlen*lonhalf)  ::  works
   real   , dimension(2*lev*mylonlen*lathalf)  ::  workr
   integer, dimension(nsize)                   ::  lensend,lenrecv            ,&
                                                   locsend,locrecv
!
   mn=0
   do n = 1,nsize 
     locsend(n)=mn
     do j = 1,mylatlen
       do k = 1,lev
         do i = lonstr(n),lonstr(n)+lonlen(n)-1
           mn = mn + 1
           works(mn)=a(i        ,k,j)
         enddo
         do i = lonstr(n),lonstr(n)+lonlen(n)-1
           mn = mn + 1
           works(mn)=a(i+lonhalf,k,j)
         enddo
       enddo
     enddo
     lensend(n)=mn-locsend(n)
   enddo
!
   mn=0
   do n = 1,nsize
     locrecv(n)=mn
     lenrecv(n)=lev*latlen(n)*mylonlen*2
     mn=mn+lenrecv(n)
   enddo
!
#ifdef MP
   call mpi_barrier (COMM_PE,ierr)
   call mpi_alltoallv(works,lensend,locsend,MPI_REAL                          ,&
                      workr,lenrecv,locrecv,MPI_REAL,COMM_PE,ierr)
#define WORKR workr
#else
#define WORKR works
#endif /* MP end */
   mn=0
   jj=0
   do n = 1,nsize
     do j = 1,latlen(n)
       jj = jj + 1
       lat1=truej(jj)
       lat2=latfull+1-lat1
       do k = 1,lev
         do i = 1,mylonlen
           mn = mn + 1
           b(lat1,k,i) = WORKR(mn)
         enddo
         do i = 1,mylonlen
           mn = mn + 1
           b(lat2,k,i) = WORKR(mn)
         enddo
       enddo
     enddo
   enddo
!
#undef WORKR
   return
#endif /* NISLQ end */
   end subroutine nislq_transpose_we2ns
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine nislq_transpose_ns2we(a,b,lev,nsize)
#ifdef NISLQ
!-------------------------------------------------------------------------------
!
! mpi transport from full dimension of west-east to full dimension of 
! north-south with latitude shuffl.
!
!-------------------------------------------------------------------------------
   use nislq, only : latfull, lonfull                                         ,&
                     lonpart, latpart                                         ,&
                     lathalf, lonhalf                                         ,&
                     mylonlen,mylatlen                                        ,&
                     lonstr, lonlen                                           ,&
                     truej, latlen
#ifdef MP
#ifdef DFS
   use commpi, only : COMM_PE=>comm_col   , MPI_REAL=>real_type
#else
   use commpi, only : COMM_PE=>comm_column, MPI_REAL
#endif
#endif /* MP end */
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                                     ::  lev,nsize                  ,&
                                                   i,j,k,n,mn,jj,lat1,lat2,ierr
   real   , dimension(latfull,lev,lonpart)     ::  a
   real   , dimension(lonfull,lev,latpart)     ::  b

   real   , dimension(2*lev*mylonlen*lathalf)  ::  works
   real   , dimension(2*lev*mylatlen*lonhalf)  ::  workr
   integer, dimension(nsize)                   ::  lensend,lenrecv            ,&
                                                   locsend,locrecv
!
   mn=0
   jj=0
   do n = 1,nsize
     locsend(n)=mn
     do j = 1,latlen(n)
       jj=jj+1
       lat1=truej(jj)
       lat2=latfull+1-lat1
       do k = 1,lev
         do i = 1,mylonlen
           mn = mn + 1
           works(mn) = a(lat1,k,i)
         enddo
         do i = 1,mylonlen
           mn = mn + 1
           works(mn) = a(lat2,k,i)
         enddo
       enddo
     enddo
     lensend(n)=mn-locsend(n)
   enddo
!
   mn=0
   do n = 1,nsize
     locrecv(n)=mn
     lenrecv(n)=lev*mylatlen*lonlen(n)*2
     mn=mn+lenrecv(n)
   enddo
!
#ifdef MP
   call mpi_barrier (COMM_PE,ierr)
   call mpi_alltoallv(works,lensend,locsend,MPI_REAL                          ,&
                      workr,lenrecv,locrecv,MPI_REAL,COMM_PE,ierr)
#define WORKR workr
#else
#define WORKR works
#endif /* MP end */
!
   mn=0
   do n = 1,nsize 
     do j = 1,mylatlen
       do k = 1,lev
         do i = lonstr(n),lonstr(n)+lonlen(n)-1
           mn = mn + 1
           b(i        ,k,j) = WORKR(mn)
         enddo
         do i = lonstr(n),lonstr(n)+lonlen(n)-1
           mn = mn + 1
           b(i+lonhalf,k,j) = WORKR(mn)
         enddo
       enddo
     enddo
   enddo
#undef WORKR
!
   return
#endif /* NISLQ end */
   end subroutine nislq_transpose_ns2we
!-------------------------------------------------------------------------------
