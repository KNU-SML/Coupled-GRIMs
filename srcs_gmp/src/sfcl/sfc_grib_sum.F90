!
   subroutine sfc_grib_sum(x,ijdim,kdim,var)
!-------------------------------------------------------------------------------
   real x(ijdim,kdim)
   character*32 var
!
   call sfc_grib_numchar(var,nchar)
   nchar=min(nchar,16)
!
   do k = 1,kdim
     sum=0.
     do ij = 1,ijdim
       sum=sum+x(ij,k)
     enddo
     print *,var(1:nchar),' sum=',sum,' at k=',k
   enddo
!
   return
   end subroutine sfc_grib_sum
!-------------------------------------------------------------------------------
