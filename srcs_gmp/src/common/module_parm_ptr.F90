!
   module module_parm_ptr
!-------------------------------------------------------------------------------
   integer  ::  kszs,kszsx,kszsy,ksps,kspsx,kspsy,ksd,ksz,ksu,ksv,kst
   integer  ::  ksq,kso3,ksqc,ksqr,ksqi,ksqs
   integer  ::  kpz,kpu,kpv,kpr,kpq,kpt,kpo,kpa,kpo3,kpqc,kpqr,kpqi
   integer  ::  kpqs,kpqg,kpsun
!
   contains
!-------------------------------------------------------------------------------
   subroutine ini_parm_ptr(levs,ko)
!-------------------------------------------------------------------------------
!
! subprogram: ini_parm_ptr     set both input and output indices
!
!-------------------------------------------------------------------------------
   kszs=1
   kszsx=2
   kszsy=3
   ksps=4
   kspsx=5
   kspsy=6
   ksd=7
   ksz=levs+7
   ksu=2*levs+7
   ksv=3*levs+7
   kst=4*levs+7
   ksq=5*levs+7
!
   kpz=1
   kpu=ko+1
   kpv=2*ko+1
   kpr=3*ko+1
   kpq=4*ko+1
   kpt=5*ko+1
   kpo=6*ko+1
   kpa=7*ko+1
   kpo3=8*ko+1
   kpqc=9*ko+1
   kpqr=10*ko+1
   kpqi=11*ko+1
   kpqs=12*ko+1
   kpqg=13*ko+1
   kpsun=14*ko+1
!
   return
   end subroutine ini_parm_ptr
!-------------------------------------------------------------------------------
   end module module_parm_ptr
   
