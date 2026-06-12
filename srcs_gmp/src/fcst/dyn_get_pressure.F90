   subroutine dyn_get_pressure(im,ix,levs,rkap,cp,fv,t,q,psexp,si,             &
                        prsi,prki,prsl,prkl,phii,phil,del)
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    [dyn_get_pressure] --- [dyn_get_height]
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer               ::  im, ix, levs
   real                  ::  rkap, cp, fv
   real                  ::  prsi(ix,levs+1)
   real                  ::  prki(ix,levs+1)
   real                  ::  phii(ix,levs+1)
   real                  ::  phil(ix,levs)
   real                  ::  prsl(ix,levs)
   real                  ::  prkl(ix,levs)
   real                  ::  del(ix,levs)
   real                  ::  T(ix,levs)
   real                  ::  q(ix,levs)
   real                  ::  psexp(ix)
   real                  ::  si(levs)
   real                  ::  rkapi, rkapp1, tem, dphib, dphit
   integer               ::  i, k
   real,parameter        ::  cons_0 = 0.d0
!-------------------------------------------------------------------------------
   rkapi  = 1.0 / rkap
   rkapp1 = 1.0 + rkap
!
   do k = 1,levs
     do i = 1,im
       del(i,k) = prsi(i,k) - prsi(i,k+1)
     enddo
   enddo
!
   if (prki(1,1) .le. 0.0) then
     do i = 1,im
       prki(i,1) = (prsi(i,1)*0.01) ** rkap
     enddo
!
     do k = 1,levs
       do i = 1,im
         prki(i,k+1) = (prsi(i,k+1)*0.01) ** rkap
         tem         = rkapp1 * del(i,k)
         prkl(i,k)   = (prki(i,k)*PRSI(i,k)-prki(i,k+1)*PRSI(i,k+1))/tem
       enddo
     enddo
   elseif (prkl(1,1) .le. 0.0) then
     do k = 1,levs
       do i = 1,im
         tem         = rkapp1 * del(i,k)
         prkl(i,k)   = (prki(i,k)*prsi(i,k)-prki(i,k+1)*prsi(i,k+1))/tem
       enddo
     enddo
   endif
!
   if (prsl(1,1) .le. 0.0) then
     do k = 1,levs
       do i = 1,im
         prsl(i,k)   = 100.0 * prkl(i,k) ** rkapi
       enddo
     enddo
   endif
!
   if (phil(1,levs) .le. 0.0) then
     do i = 1,im
       phii(i,1)   = 0.0           ! Ignoring topography height here
     enddo
     do k = 1,levs
       do i = 1,im
         tem         = cp * t(i,k) * (1.0 + fv * max(q(i,k),cons_0))/prkl(i,k)
         dphib       = (prki(i,k) - prkl(i,k)) * tem
         dphit       = (prkl(i,k) - prki(i,k+1)) * tem
         phil(i,k)   = phii(i,k) + dphib
         phii(i,k+1) = phil(i,k) + dphit
       enddo
     enddo
   endif
!
   return
   end
!
!===============================================================================
   subroutine dyn_get_height(im,ix,levs,cp,fv,t,q,prki,prkl,phii,phil)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer               ::  im, ix, levs
   real                  ::  rkap, cp, fv
   real                  ::  prsi(ix,levs+1)
   real                  ::  prki(ix,levs+1)
   real                  ::  phii(ix,levs+1)
   real                  ::  phil(ix,levs)
   real                  ::  prsl(ix,levs)
   real                  ::  prkl(ix,levs)
   real                  ::  del(ix,levs)
   real                  ::  T(ix,levs)
   real                  ::  q(ix,levs)
   real                  ::  psexp(ix)
   real                  ::  si(levs)
   real                  ::  rkapi, rkapp1, tem, dphib, dphit
   integer               :: i, k
   real,parameter        :: cons_0 = 0.d0
!
   do i = 1,im
     phii(i,1)   = 0.0           ! Ignoring topography height here
   enddo
!
   do k = 1,levs
     do i = 1,im
       tem         = cp * t(i,k) * (1.0 + fv * max(q(i,k),cons_0)) / prkl(i,k)
       dphib       = (prki(i,k) - prkl(i,k)) * tem
       dphit       = (prkl(i,k) - prki(i,k+1)) * tem
       phil(i,k)   = phii(i,k) + dphib
       phii(i,k+1) = phil(i,k) + dphit
     enddo
   enddo
!
   return
   end 
