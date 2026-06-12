!
   subroutine sfc_tsfc_correc(tsfc,orog,slmask,umask,ijdim)                    
!-------------------------------------------------------------------------------
   real, parameter       ::  rlapse=0.65e-2
!                                                                               
   real                  ::   tsfc(ijdim),orog(ijdim),slmask(ijdim)
!
   do ij = 1,ijdim
     if(slmask(ij).eq.umask) then
       tsfc(ij) = tsfc(ij)-orog(ij)*rlapse
     endif
   enddo
!
   return
   end subroutine sfc_tsfc_correc
!-------------------------------------------------------------------------------
