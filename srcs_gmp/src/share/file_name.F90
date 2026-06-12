!
   subroutine file_name(fnin,nchin,fhour,fnout,nchout)
!-------------------------------------------------------------------------------
   character(len=10)   ::   cfhour
   character(len=80)   ::   fnin,fnout
   character(len=4 )   ::   fmt
!
   do n = 1,80
     fnout(n:n)=' '
   enddo
!
   ifh=fhour
!
   if (ifh.eq.0) then
     ndig=1
   else
     ndig=log10(fhour)+1
   endif
!
   if (ndig.eq.1) then
     write(cfhour,'(I2.2)') ifh
     ndig=2
   else
     write(fmt,'(A2,I1,A1)')'(I',ndig,')'
     write(cfhour,fmt) ifh
   endif
!
!  fnout=fnin(1:nchin)//'.ft'//cfhour(1:ndig)
   fnout(1:nchin)=fnin(1:nchin)
   fnout(nchin+1:nchin+3)='.ft'
   fnout(1+nchin+3:ndig+nchin+3)=cfhour(1:ndig)
   nchout=nchin+3+ndig
   !
   return 
   end subroutine file_name
!-------------------------------------------------------------------------------
