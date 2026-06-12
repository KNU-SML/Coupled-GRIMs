!
   subroutine sfc_file_subst(data,imax,jmax,dlon,dlat,ijordr)
!-------------------------------------------------------------------------------
   logical              ::  ijordr
   real                 ::  data(imax*jmax)
   real, allocatable    ::  work(:)
!
   allocate (work(imax*jmax))
!
   if(.not.ijordr.or.                                                          &
      (ijordr.and.(dlat.gt.0..or.dlon.lt.0.))) then
     if(.not.ijordr) then
       ij=0
       do j = 1,jmax
         do i = 1,imax
           ij=(j-1)*imax+i
           ji=(i-1)*jmax+j
           work(ij)=data(ji)
         enddo
       enddo
     else
       do j = 1,jmax
         do i = 1,imax
           ij=(j-1)*imax+i
           work(ij)=data(ij)
         enddo
       enddo
     endif
     do j = 1,jmax
       do i = 1,imax
         if(dlat.gt.0..and.dlon.gt.0.) then
           ij=imax*jmax-imax*j+i
         elseif(dlat.gt.0..and.dlon.lt.0.) then
           ij=imax*jmax-(j-1)*imax-imax+i-1
         elseif(dlat.lt.0..and.dlon.lt.0.) then
           ij=imax*(j-1)+imax-i+1
         endif
         ijo=(j-1)*imax+i
         data(ij)=work(ijo)
       enddo
     enddo
   endif
!
   deallocate (work)
!
   return
   end subroutine sfc_file_subst
!-------------------------------------------------------------------------------
