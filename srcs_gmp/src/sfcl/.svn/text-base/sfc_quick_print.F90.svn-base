!
   subroutine sfc_quick_print(data,imax,jmax)
!-------------------------------------------------------------------------------
   real        ::  data(imax*jmax)
!
   rmax=data(1)
   rmin=data(1)
   do ij = 1,imax*jmax
     rmax=max(rmax,data(ij))
     rmin=min(rmin,data(ij))
   enddo
!
   rmxmn=rmax-rmin
   if(rmax.ne.0.) then
     rn=log10(abs(rmax))
     if(rn.lt.0.) then
       n=rn-1.0
     else
       n=rn
     endif
   else
     n=-9999
   endif
!
   if(rmin.ne.0.) then
     rm=log10(abs(rmin))
     if(rm.lt.0.) then
       m=rm-1.0
     else
       m=rm
     endif
   else
     m=-9999
   endif
!
   n=max(n,m)
   fact=10.**(-n)
   print *,'sfc_quick_print:rmax=',rmax,' rmin=',rmin,' fact=',fact
!
   ilast=0
   i1=1
   i2=80
   if(i2.ge.imax) then
     ilast=1
     i2=imax
   endif
!
1112 continue
   write(6,*) ' '
!
   do j = 1,jmax
     write(6,1111) (nint(data(imax*(j-1)+i)*fact),i=i1,i2)
   enddo
!
   if(ilast.eq.1) return
!
   i1=i1+80
   i2=i1+79
!
   if(i2.ge.imax) then
     ilast=1
     i2=imax
   endif
   !
   go to 1112
   !
1111 format(80i1)
!
!  return
   end subroutine sfc_quick_print
!-------------------------------------------------------------------------------
