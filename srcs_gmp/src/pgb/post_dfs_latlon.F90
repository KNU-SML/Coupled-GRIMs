!
   subroutine post_dfs_latlon(din,imxin,jmxin,dout,imxout,jmxout,lot)
!-------------------------------------------------------------------------------
!
!  subprgram:    post_dfs_latlon  
!               
!  abstract:  interpolation from DFS offset grid to lat/lon grid
!
!-------------------------------------------------------------------------------
#undef DBG
   integer imxin,jmxin,imxout,jmxout,lot
   real     din (imxin,jmxin,lot)
   real     dout(imxout,jmxout,lot)
   real, allocatable, dimension(:)   :: phai,phao,ddx,ddy
   integer,allocatable, dimension(:)   :: iindx1,iindx2,jindx1,jindx2
   logical first
   data first/.true./
   save phai,phao,ddx,ddy,iindx1,iindx2,jindx1,jindx2
!-------------------------------------------------------------------------------
!
   if (first) then
     first=.false.
     allocate( phai(jmxin),phao(jmxout))
     allocate( iindx1(imxout))
     allocate( iindx2(imxout))
     allocate( jindx1(jmxout))
     allocate( jindx2(jmxout))
     allocate( ddx(imxout))
     allocate( ddy(jmxout))
     !
     ! DFS grid
     ! 
     dlati=180./jmxin
     do j = 1,jmxin
       phai(j)=90.-(j-0.5)*dlati
     enddo
     !
     ! LAT/LON grid
     ! 
     dlato=180./float(jmxout-1)
     do j = 1,jmxout
       phao(j)=90.-float(j-1)*dlato
     enddo
!
     dxin =360./float(imxin )
     dxout=360./float(imxout)
!
     do i = 1,imxout
       alamd=float(i-1)*dxout
       if(alamd.lt.0.) alamd=alamd+360.
       if(alamd.gt.360.) alamd=alamd-360.
       i1=alamd/dxin+1.001
       iindx1(i)=i1
       i2=i1+1
       if(i2.gt.imxin) i2=1
       iindx2(i)=i2
       ddx(i)=(alamd-float(i1-1)*dxin)/dxin
     enddo
!
     do j = 1,jmxout
       apho=phao(j)
       do jj = 1,jmxin
         if(apho.lt.phai(jj)) cycle
         j2=jj
         !
         go to 42
         !
       enddo
       j2=jmxin
       !
       42 continue
       !
       if(j2.gt.2) go to 43
       j1=1
       j2=2
       !
       go to 44
       !
       43 continue
       if(j2.le.jmxin) go to 45
       j1=jmxin-1
       j2=jmxin
       !
       go to 44
       !
       45 continue
       j1=j2-1
       44 continue
       jindx1(j)=j1
       jindx2(j)=j2
#ifdef DBG
       write(6,'(A,3F10.5,2I4)')'outlat,inlat,j1,j2=',                         &
                             apho,phai(j1),phai(j2),j1,j2
#endif
       ddy(j)=(apho-phai(j1))/(phai(j2)-phai(j1))
     enddo
!
#ifdef DBG
     write(6,*) 'post_dfs_latlon'
     write(6,*) 'iindx1'
     write(6,*) (iindx1(n),n=1,imxout)
     write(6,*) 'iindx2'
     write(6,*) (iindx2(n),n=1,imxout)
     write(6,*) 'jindx1'
     write(6,*) (jindx1(n),n=1,jmxout)
     write(6,*) 'jindx2'
     write(6,*) (jindx2(n),n=1,jmxout)
     write(6,*) 'ddy'
     write(6,*) (ddy(n),n=1,jmxout)
     write(6,*) 'ddx'
     write(6,*) (ddx(n),n=1,jmxout)
#endif
   endif         ! first
!
   do k = 1,lot
!
     do j = 1,jmxout
       y=ddy(j)
       j1=jindx1(j)
       j2=jindx2(j)
       do i = 1,imxout
         x=ddx(i)
         i1=iindx1(i)
         i2=iindx2(i)
         dout(i,j,k)=(1.-x)*(1.-y)*din(i1,j1,k)+(1.-y)*x*din(i2,j1,k)+         &
                        (1.-x)*y*din(i1,j2,k)+x*y*din(i2,j2,k)
       enddo
     enddo
!
     sum1=0.
     sum2=0.
     do i = 1,imxin
       sum1=sum1+din(i,1,k)
       sum2=sum2+din(i,jmxin,k)
     enddo
     sum1=sum1/float(imxin)
     sum2=sum2/float(imxin)
!
     do i = 1,imxout
       dout(i,     1,k)=sum1
       dout(i,jmxout,k)=sum2
     enddo
!
   enddo      ! lot
!
   return
   end
