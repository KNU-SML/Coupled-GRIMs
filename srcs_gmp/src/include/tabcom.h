!
!     common block tabcom contains quantities precomputed in subroutine
!     table for use in the longwave radiation program:
!          em1     =  e1 function, evaluated over the 0-560 and
!                     1200-2200 cm-1 intervals
!          em1wde  =  e1 function, evaluated over the 160-560 cm-1
!                     interval
!          table1  =  e2 function, evaluated over the 0-560 and
!                     1200-2200 cm-1 intervals
!          table2  =  temperature derivative of table1
!          table3  =  mass derivative of table1
!          em3     =  e3 function, evaluated over the 0-560 and
!                     1200-2200 cm-1 intervals
!          source  =  planck function, evaluated at specified temps. for
!                     bands used in cts calculations
!          dsrce   =  temperature derivative of source
!          ind     =  index, with value ind(i)=i. used in fst88
!          indx2   =  index values used in obtaining "lower triangle"
!                     elements of avephi,etc.,in fst88
!          kmaxv   =  index values used in obtaining "upper triangle"
!                     elements of avephi,etc.,in fst88
!          kmaxvm  =  kmaxv(l),used for do loop indices
!
#if !defined(MP) && defined (RMP)
   use comfcst_ser, only : ind,indx2,kmaxv,kmaxvm,idummy2,em1,em1wde,&
                       table1,table2,table3,em3,source,dsrce
#else
   use comfcst, only : ind,indx2,kmaxv,kmaxvm,idummy2,em1,em1wde,&
                       table1,table2,table3,em3,source,dsrce
#endif
