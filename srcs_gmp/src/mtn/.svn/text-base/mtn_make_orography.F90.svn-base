!
   subroutine mtn_make_orography(zavg,zvar,zslm,oro,slm,var,var4,              &
                     glat,ist,ien,jst,jen,im,jm,imn,jmn,xlat,numi)
!-------------------------------------------------------------------------------
   real       ::  glat(jmn),xlat(jm)
   integer    ::  zavg(imn,jmn),zvar(imn,jmn),zslm(imn,jmn)
   real       ::  oro(im,jm),slm(im,jm),var(im,jm),var4(im,jm)
   integer    ::  ist(im,jm),ien(im,jm),jst(jm),jen(jm),numi(jm)
   logical        flag
!
!---- global xlat and xlon ( degree )
!
   jm1 = jm - 1
   delxn = 360./imn      ! mountain data resolution
!
   do j = 1,jmn
     glat(j) = -90. + (j-1) * delxn + delxn * 0.5
   enddo
!
!---- find the average of the modes in a grid box
!
   do j = 1,jm
     do i = 1,numi(j)
       im1 = numi(j) - 1
       delx  = 360./numi(j)       ! gaussian grid resolution
       faclon  = delx / delxn
       ist(i,j) = faclon * float(i-1) - faclon * 0.5
       ist(i,j) = ist(i,j) + 1
       ien(i,j) = faclon * float(i) - faclon * 0.5
!x         print*, ' i j ist ien ',i,j,ist(i,j),ien(i,j)
     enddo
   enddo
!
   do j = 1,jm-1
     flag=.true.
     do j1 = 1,jmn
       xxlat = (xlat(j)+xlat(j+1))/2.
       if(flag.and.glat(j1).gt.xxlat) then
         jst(j) = j1
         jen(j+1) = j1 - 1
         flag = .false.
       endif
     enddo
!x         print*, ' j jst jen ',j,jst(j),jen(j),xlat(j),glat(j1)
   enddo
!
   jst(jm) = max(jst(jm-1) - (jen(jm-1)-jst(jm-1)),1)
   jen(1)  = min(jen(2) + (jen(2)-jst(2)),jmn)      
!
!...first, averaged height
!
   do j = 1,jm
     do i = 1,numi(j)
       oro(i,j)  = 0.0
       var(i,j)  = 0.0
       var4(i,j) = 0.0
       xnsum = 0.0
       xland = 0.0
       xwatr = 0.0
       xl1 = 0.0
       xs1 = 0.0
       xw1 = 0.0
       xw2 = 0.0
       xv2 = 0.0
       xw4 = 0.0
       do ii1 = ist(i,j),ien(i,j)
         i1 = ii1
         if(ii1.le.0.)  i1 = i1 + imn
         if(ii1.gt.imn) i1 = i1 - imn
         do j1 = jst(j),jen(j)
           xland = xland + float(zslm(i1,j1))
           xwatr = xwatr + float(1-zslm(i1,j1))
           xnsum = xnsum + 1.
           height = float(zavg(i1,j1)) 
           if(height.lt.-990.) height = 0.0
           xl1 = xl1 + height * float(zslm(i1,j1))
           xs1 = xs1 + height * float(1-zslm(i1,j1))
           xw1 = xw1 + height
           xw2 = xw2 + height ** 2
           xv2 = xv2 + float(zvar(i1,j1)) ** 2
         enddo
       enddo
       if(xnsum.gt.1.) then
         slm(i,j) = float(nint(xland/xnsum))
         if(slm(i,j).ne.0.) then
           oro(i,j)= xl1 / xland
         else
           oro(i,j)= xs1 / xwatr
         endif
         var(i,j)=sqrt(max((xv2+xw2)/xnsum-(xw1/xnsum)**2,0.))
         do ii1 = ist(i,j),ien(i,j)
           i1 = ii1
           if(ii1.lt.0.) i1 = i1 + imn
           if(ii1.gt.imn) i1 = i1 - imn
           do j1 = jst(j),jen(j)
             height = float(zavg(i1,j1)) 
             if(height.lt.-990.) height = 0.0
             xw4 = xw4 + (height-oro(i,j)) ** 4
           enddo
         enddo
         if(var(i,j).gt.1.) then
           var4(i,j) = xw4/xnsum/var(i,j) **4
         endif
       endif
     enddo
   enddo
!
   write(6,*) "! make oro slm var var4 done"
!
   return
   end subroutine mtn_make_orography
!-------------------------------------------------------------------------------
