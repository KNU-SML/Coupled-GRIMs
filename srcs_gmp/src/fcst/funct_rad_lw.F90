#include <define.h>
   module funct_rad_lw
!-------------------------------------------------------------------------------
!
!  ::: structure :::
!
!    for rad_lw_gfdl
!      --- [rad_lw_ele290]
!      --- [rad_lw_e290]
!      --- [rad_lw_e2spec]
!      --- [rad_lw_e3v88]
!      --- [rad_lw_spa88]
!
!    for rad_lw_nasa
!      --- [rad_lw_nasa_b10exps]
!      --- [rad_lw_nasa_b10kdis]
!      --- [rad_lw_nasa_cfcexps]
!      --- [rad_lw_nasa_cfckdis]
!      --- [rad_lw_nasa_ch4exps]
!      --- [rad_lw_nasa_ch4kdis]
!      --- [rad_lw_nasa_co2exps]
!      --- [rad_lw_nasa_co2kdis]
!      --- [rad_lw_nasa_column]
!      --- [rad_lw_nasa_comexps]
!      --- [rad_lw_nasa_comkdis]
!      --- [rad_lw_nasa_conexps]
!      --- [rad_lw_nasa_h2oexps]
!      --- [rad_lw_nasa_h2okdis]
!      --- [rad_lw_nasa_n2oexps]
!      --- [rad_lw_nasa_n2okdis]
!      --- [rad_lw_nasa_tablup]
!
!-------------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_ele290(ipts,g1,g2,g3,g4,g5,emiss,fxoe1,dte1,fxoe2,dte2,   &
                            avephi,temp,t)
!-------------------------------------------------------------------------------
!
!     subroutine rad_lw_ele290 computes the exchange terms in the flux equation
!  for longwave radiation for all terms except the exchange with the
!  top of the atmosphere. the method is a table lookup on a pre-
!  computed e2 function (defined in ref. (4)).
!      the e1 function  calculations (formerly done in subroutine
!  e1v88 compute the flux resulting from the exchange of photons
!  between a layer and the top of the atmosphere.  the method is a
!  table lookup on a pre-computed e1 function.
!     calculations are done in two frequency ranges:
!       1) 0-560,1200-2200 cm-1   for q(approx)
!       2) 160-560 cm-1           for q(approx,cts).
!  motivation for these calculations is in references (1) and (4).
!       inputs:                    (common blocks)
!     table1,table2,table3,em1,em1wde  tabcom
!     avephi                           tfcom
!     temp                             radisw
!     t                                kdacom
!     fxoe1,dte1                argument list
!     fxoe2,dte2                argument list
!       outputs:
!     emiss                            tfcom
!     g1,g2,g3                  argument list,for 1st freq. range
!     g4,g5                     argument list,for 2nd freq. range
!
!        called by :     rad_lw
!        calls     :
!
!-------------------------------------------------------------------------------
   use paramodel
   use hcon
   use rdparm
!-------------------------------------------------------------------------------
#include <tabcom.h>
   real                 ::  temp(imbx,lp1),t(imbx,lp1)
   real                 ::  avephi(imbx,lp1),emiss(imbx,lp1)
!
   integer              ::  it1(imbx,ll3p),ival(imbx,lp1)
   real                 ::  fyo(imbx,lp1),du(imbx,lp1)
   real                 ::  ww1(imbx,lp1),ww2(imbx,lp1)
   real                 ::  tmp3(imbx,lp1),tmp5(imax),tmp9(imax)
!
!---variables equivalenced to common block variables
!
   real                 ::  t1(5040),t2(5040),t4(5040)
   real                 ::  em1v(5040),em1vw(5040)
!
!---variables in the argument list
!
   real                 ::  fxoe1(imbx,lp1),dte1(imbx,lp1)
   real                 ::  fxoe2(imbx,lp1),dte2(imbx,lp1)
   real                 ::  g1(imbx,lp1),g2(imbx,l),g3(imbx,lp1)
   real                 ::  g4(imbx,lp1),g5(imbx,l)
!
   n=0
   do k = 1,180
     do i = 1,28
       n=n+1
       em1v(n)=em1(i,k)
       em1vw(n)=em1wde(i,k)
       t1(n)=table1(i,k)
       t2(n)=table2(i,k)
       t4(n)=table3(i,k)
     enddo
   enddo
!
!---first we obtain the emissivities as a function of temperature
!   (index fxo) and water amount (index fyo). this part of the code
!   thus generates the e2 function. the fxo indices have been
!   obtained in rad_lw, for convenience.
!
!---this subroutine evaluates the k=1 case only--
!
!---this loop replaces loops going fromi=1,imax and kp=2,lp1 plus
!   the special case for the lp1th layer.
!
   do kp = 1,lp1
     do i=1,ipts
       tmp3(i,kp)=log10(avephi(i,kp))+h16e1
       fyo(i,kp)=aint(tmp3(i,kp)*ten)
       du(i,kp)=tmp3(i,kp)-hp1*fyo(i,kp)
       fyo(i,kp)=h28e1*fyo(i,kp)
       ival(i,kp)=fyo(i,kp)+fxoe2(i,kp)
       emiss(i,kp)=t1(ival(i,kp))+du(i,kp)*t2(ival(i,kp))                      &
                           +dte2(i,kp)*t4(ival(i,kp))
     enddo
   enddo
!
!---the special case emiss(i,l) (layer kp) is obtained now
!   by averaging the values for l and lp1:
!
   do i = 1,ipts
     emiss(i,l)=haf*(emiss(i,l)+emiss(i,lp1))
   enddo
!
!   calculations for the kp=1 layer are not performed, as
!   the radiation code assumes that the top flux layer (above the
!   top data level) is isothermal, and hence contributes nothing
!   to the fluxes at other levels.
!
!***the following is the calculation for the e1 function, formerly
!    done in subroutine e1v88. the move to e1e288 is due to the
!    savings in obtaining index values (the temp. indices have
!    been obtained in rad_lw, while the u-indices are obtained
!    in the e2 calcs.,with k=1).
!
!   for terms involving top layer, du is not known; in fact, we
!   use index 2 to repersent index 1 in prev. code. this means that
!    the it1 index 1 and llp1 has to be calculated separately. the
!   index llp2 gives the same value as 1; it can be omitted.
!
   do i = 1,ipts
     it1(i,1)=fxoe1(i,1)
     ww1(i,1)=ten-dte1(i,1)
     ww2(i,1)=hp1
   enddo
!
   do kp = 1,l
     do i = 1,ipts
       it1(i,kp+1)=fyo(i,kp)+fxoe1(i,kp+1)
       it1(i,kp+lp1)=fyo(i,kp)+fxoe1(i,kp)
       ww1(i,kp+1)=ten-dte1(i,kp+1)
       ww2(i,kp+1)=hp1-du(i,kp)
     enddo
   enddo
!
   do kp = 1,l
     do i = 1,ipts
       it1(i,kp+llp1)=fyo(i,kp)+fxoe1(i,1)
     enddo
   enddo
!
!  g3(i,1) has the same values as g1 (and did all along)
!
   do i = 1,ipts
     g1(i,1)=ww1(i,1)*ww2(i,1)*em1v(it1(i,1))+                                 &
           ww2(i,1)*dte1(i,1)*em1v(it1(i,1)+1)
     g3(i,1)=g1(i,1)
   enddo
!
   do kp = 1,l
     do i = 1,ipts
       g1(i,kp+1)=ww1(i,kp+1)*ww2(i,kp+1)*em1v(it1(i,kp+1))+                   &
           ww2(i,kp+1)*dte1(i,kp+1)*em1v(it1(i,kp+1)+1)+                       &
           ww1(i,kp+1)*du(i,kp)*em1v(it1(i,kp+1)+28)+                          &
           dte1(i,kp+1)*du(i,kp)*em1v(it1(i,kp+1)+29)
       g2(i,kp)=ww1(i,kp)*ww2(i,kp+1)*em1v(it1(i,kp+lp1))+                     &
           ww2(i,kp+1)*dte1(i,kp)*em1v(it1(i,kp+lp1)+1)+                       &
           ww1(i,kp)*du(i,kp)*em1v(it1(i,kp+lp1)+28)+                          &
           dte1(i,kp)*du(i,kp)*em1v(it1(i,kp+lp1)+29)
     enddo
   enddo
!
   do kp = 2,lp1
     do i = 1,ipts
       g3(i,kp)=ww1(i,1)*ww2(i,kp)*em1v(it1(i,ll+kp))+                         &
           ww2(i,kp)*dte1(i,1)*em1v(it1(i,ll+kp)+1)+                           &
           ww1(i,1)*du(i,kp-1)*em1v(it1(i,ll+kp)+28)+                          &
           dte1(i,1)*du(i,kp-1)*em1v(it1(i,ll+kp)+29)
     enddo
   enddo
!
   do i = 1,ipts
     g4(i,1)=ww1(i,1)*ww2(i,1)*em1vw(it1(i,1))+                                &
           ww2(i,1)*dte1(i,1)*em1vw(it1(i,1)+1)
   enddo
!
   do kp = 1,l
     do i = 1,ipts
       g4(i,kp+1)=ww1(i,kp+1)*ww2(i,kp+1)*em1vw(it1(i,kp+1))+                  &
           ww2(i,kp+1)*dte1(i,kp+1)*em1vw(it1(i,kp+1)+1)+                      &
           ww1(i,kp+1)*du(i,kp)*em1vw(it1(i,kp+1)+28)+                         &
           dte1(i,kp+1)*du(i,kp)*em1vw(it1(i,kp+1)+29)
       g5(i,kp)=ww1(i,kp)*ww2(i,kp+1)*em1vw(it1(i,kp+lp1))+                    &
           ww2(i,kp+1)*dte1(i,kp)*em1vw(it1(i,kp+lp1)+1)+                      &
           ww1(i,kp)*du(i,kp)*em1vw(it1(i,kp+lp1)+28)+                         &
           dte1(i,kp)*du(i,kp)*em1vw(it1(i,kp+lp1)+29)
     enddo
   enddo
!
   return
   end subroutine rad_lw_ele290
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_e290(ipts,emissb,emiss,avephi,klen,fxoe2,dte2)
!-------------------------------------------------------------------------------
!
!     subroutine rad_lw_e290 computes the exchange terms in the flux equation
!  for longwave radiation for all terms except the exchange with the
!  top of the atmosphere. the method is a table lookup on a pre-
!  computed e2 function (defined in ref. (4)).
!     calculations are done in the frequency range:
!       1) 0-560,1200-2200 cm-1   for q(approx)
!  motivation for these calculations is in references (1) and (4).
!       inputs:                    (common blocks)
!     table1,table2,table3,            tabcom
!     avephi                           tfcom
!     fxoe2,dte2,klen           argument list
!       outputs:
!     emiss,emissb                     tfcom
!
!        called by :     rad_lw
!        calls     :
!
!-------------------------------------------------------------------------------
   use paramodel
   use hcon
   use rdparm
!-------------------------------------------------------------------------------
#include <tabcom.h>
   real                 ::  emissb(imbx,lp1),emiss(imbx,lp1),avephi(imbx,lp1)
   integer              ::  ival(imbx,lp1)
   real, pointer        ::  dt(:,:)
   real                 ::  fyo(imbx,lp1),du(imbx,lp1)
!---tmp3 may be equivalenced to dt in vtemp
   real, target         ::  tmp3(imbx,lp1)
!---variables equivalenced to common block variables
   real                 ::  t1(5040),t2(5040),t4(5040)
!---variables in the argument list
   real                 ::  fxoe2(imbx,lp1),dte2(imbx,lp1)
!
   dt=>tmp3
   n=0
   do k = 1,180
     do i = 1,28
       n=n+1
       t1(n)=table1(i,k)
       t2(n)=table2(i,k)
       t4(n)=table3(i,k)
     enddo
   enddo
!
!---first we obtain the emissivities as a function of temperature
!   (index fxo) and water amount (index fyo). this part of the code
!   thus generates the e2 function.
!
!---calculations for varying kp (from kp=k+1 to lp1, including special
!   case: results are in emiss
!
   do k = 1,lp2-klen
     do i = 1,ipts
       tmp3(i,k)=log10(avephi(i,k+klen-1))+h16e1
       fyo(i,k)=aint(tmp3(i,k)*ten)
       du(i,k)=tmp3(i,k)-hp1*fyo(i,k)
       fyo(i,k)=h28e1*fyo(i,k)
       ival(i,k)=fyo(i,k)+fxoe2(i,k+klen-1)
       emiss(i,k+klen-1)=t1(ival(i,k))+du(i,k)*t2(ival(i,k))                   &
                              +dte2(i,k+klen-1)*t4(ival(i,k))
     enddo
   enddo
!
!---the special case emiss(i,l) (layer kp) is obtained now
!   by averaging the values for l and lp1:
!
   do i = 1,ipts
     emiss(i,l)=haf*(emiss(i,l)+emiss(i,lp1))
   enddo
!
!---note that emiss(i,lp1) is not useful after this point.
!
!---calculations for kp=klen and varying k; results are in emissb.
!  in this case, the temperature index is unchanged, always being
!  fxo(i,klen-1); the water index changes, but is symmetrical with
!  that for the varying kp case.note that the special case is not
!  involved here.
!     (fixed level) k varies from (klen+1) to lp1; results are in
!   emissb(i,(klen) to l)
!
   do k = 1,lp1-klen
     do i = 1,ipts
       dt(i,k)=dte2(i,klen-1)
       ival(i,k)=fyo(i,k)+fxoe2(i,klen-1)
     enddo
   enddo
!
   do k = 1,lp1-klen
     do i = 1,ipts
       emissb(i,k+klen-1)=t1(ival(i,k))+du(i,k)*t2(ival(i,k))                  &
                                   +dt(i,k)*t4(ival(i,k))
     enddo
   enddo
!
   return
   end subroutine rad_lw_e290
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_e2spec(ipts,emiss,avephi,fxosp,dtsp)
!-------------------------------------------------------------------------------
!
!     subroutine rad_lw_e2spec computes the exchange terms in the flux equation
!  for longwave radiation for 2 terms used for nearby layer compu-
!  tations. the method is a table lookup on a pre-
!  computed e2 function (defined in ref. (4)).
!     calculations are done in the frequency range:
!        0-560,1200-2200 cm-1
!  motivation for these calculations is in references (1) and (4).
!       inputs:                    (common blocks)
!     table1,table2,table3,            tabcom
!     avephi                           tfcom
!     fxosp,dtsp                argument list
!       outputs:
!     emiss                            tfcom
!
!        called by :     rad_lw
!        calls     :
!
!-------------------------------------------------------------------------------
   use paramodel
   use hcon
   use rdparm
!-------------------------------------------------------------------------------
#include <tabcom.h>
   real                 ::  avephi(imbx,lp1),emiss(imbx,lp1)
   real                 ::  fyo(imbx,lp1),du(imbx,lp1)
   real                 ::  tmp3(imbx,lp1)
   integer              ::  ival(imbx,lp1)
!
!---variables equivalenced to common block variables
!
   real                 ::  t1(5040),t2(5040),t4(5040)
!
!---variables in the argument list
!
   real                 ::  fxosp(imbx,2),dtsp(imbx,2)
!
   n=0
   do k = 1,180
     do i = 1,28
       n=n+1
       t1(n)=table1(i,k)
       t2(n)=table2(i,k)
       t4(n)=table3(i,k)
     enddo
   enddo
!
!---first we obtain the emissivities as a function of temperature
!   (index fxo) and water amount (index fyo). this part of the code
!   thus generates the e2 function.
!
   do k = 1,2
     do i = 1,ipts
       tmp3(i,k)=log10(avephi(i,k))+h16e1
       fyo(i,k)=aint(tmp3(i,k)*ten)
       du(i,k)=tmp3(i,k)-hp1*fyo(i,k)
       ival(i,k)=h28e1*fyo(i,k)+fxosp(i,k)
       emiss(i,k)=t1(ival(i,k))+du(i,k)*t2(ival(i,k))+                         &
                               dtsp(i,k)*t4(ival(i,k))
     enddo
   enddo
!
   return
   end subroutine rad_lw_e2spec
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_e3v88(ipts,emv,tv,av)
!-------------------------------------------------------------------------------
   use paramodel
   use hcon
   use rdparm
!-------------------------------------------------------------------------------
#include <tabcom.h>
   integer              ::  it(imbx,llp1)
   real, pointer        ::  ww1(:,:), ww2(:,:)
   real                 ::  dt(imbx,llp1)
   real                 ::  du(imbx,llp1)
!
!   the following arrays are equivalenced to vtemp arrays
!
   real, target         ::  fxo(imbx,llp1),fyo(imbx,llp1)
   real                 ::  tmp3(imbx,llp1)
!
!    dimensions of arrays in argument list
!
   real                 ::  emv(imbx,llp1),tv(imbx,llp1),av(imbx,llp1)
!
!   the following array is equivalenced to an array in tabcom.h
!
   real                 ::  em3v(5040)
!
   n=0
   do k = 1,180
     do i = 1,28
       n=n+1
       em3v(n)=em3(i,k)
     enddo
   enddo
   ww1=>fxo
   ww2=>fyo
!
!---the following loop replaces a double loop over i (1-imax) and
!   k (1-llp1)
!
   do k = 1,llp1
     do i = 1,ipts
       fxo(i,k)=aint(tv(i,k)*hp1)
       tmp3(i,k)=log10(av(i,k))+h16e1
       dt(i,k)=tv(i,k)-ten*fxo(i,k)
       fyo(i,k)=aint(tmp3(i,k)*ten)
       du(i,k)=tmp3(i,k)-hp1*fyo(i,k)
!
!---obtain index for table lookup; this value will have to be
!   decremented by 9 to account for table temps starting at 100k.
!
       it(i,k)=fxo(i,k)+fyo(i,k)*h28e1
       ww1(i,k)=ten-dt(i,k)
       ww2(i,k)=hp1-du(i,k)
       emv(i,k)=ww1(i,k)*ww2(i,k)*em3v(it(i,k)-9)+                             &
              ww2(i,k)*dt(i,k)*em3v(it(i,k)-8)+                                &
              ww1(i,k)*du(i,k)*em3v(it(i,k)+19)+                               &
              dt(i,k)*du(i,k)*em3v(it(i,k)+20)
     enddo
   enddo
!
   return
   end subroutine rad_lw_e3v88
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_spa88(ipts,excts,ctso3,gxcts,sorc,csour,                  &
#ifdef CLR
                    excts0,ctso30,gxcts0,                                      &
#endif
                    cldfac,temp,press,var1,var2,                               &
                    p,delp,delp2,totvo2,to3sp,to3spc,                          &
                    co2sp1,co2sp2,co2sp)
!-------------------------------------------------------------------------------
   use paramodel
   use hcon
   use rdparm
   use rnddta
!-------------------------------------------------------------------------------
   real                 ::  sorc(imbx,lp1,nbly),csour(imbx,lp1)
   real                 ::  cldfac(imbx,lp1,lp1)
   real                 ::  temp(imbx,lp1),press(imbx,lp1)
   real                 ::  var1(imbx,l),var2(imbx,l)
   real                 ::  p(imbx,lp1),delp(imbx,l),delp2(imbx,l)
   real                 ::  totvo2(imbx,lp1),to3spc(imbx,l),to3sp(imbx,lp1)
   real                 ::  co2sp1(imbx,lp1),co2sp2(imbx,lp1),co2sp(imbx,lp1)
   real                 ::  excts(imbx,l),ctso3(imbx,l),gxcts(imax)
#ifdef CLR
   real                 ::  excts0(imbx,l),ctso30(imbx,l),gxcts0(imax)
   real                 ::  ctmp0(imbx,lp1),ctmp20(imbx,lp1),ctmp30(imbx,lp1)
#endif
!
   real                 ::  phitmp(imbx,l),psitmp(imbx,l)
   real                 ::  tt(imbx,l)
   real                 ::  fac1(imbx,l),fac2(imbx,l)
   real                 ::  ctmp(imbx,lp1),x(imbx,l),y(imbx,l)
   real                 ::  topm(imbx,l),topphi(imbx,l)
   real                 ::  ctmp3(imbx,lp1),ctmp2(imbx,lp1)
   real                 ::  f(imbx,l),ff(imbx,l),ag(imbx,l),agg(imbx,l)
!
!---compute temperature quantities for use in program
!
   do k = 1,l
     do i = 1,imax
       x(i,k)=temp(i,k)-h25e2
       y(i,k)=x(i,k)*x(i,k)
     enddo
   enddo
!
!---initialize ctmp(i,1),ctmp2(i,1),ctmp3(i,1) to unity; these are
!   transmission fctns at the top.
!
   do i = 1,ipts
     ctmp(i,1)=one
     ctmp2(i,1)=1.
     ctmp3(i,1)=1.
   enddo
!
!...   ditto for the clear sky calculation...
!
#ifdef CLR
   do i = 1,ipts
     ctmp0(i,1)=one
     ctmp20(i,1)=1.
     ctmp30(i,1)=1.
   enddo
#endif
!
!***begin loop on frequency bands (1)**!
!
!---calculation for band 1 (combined band 1)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(1)*x(i,k)+bpcm(1)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(1)*x(i,k)+btpcm(1)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(1)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(1)*topphi(i,k))
       tt(i,k)=exp(hm1ez*fac1(i,k)/sqrt(1.+fac2(i,k)))
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=sorc(i,k,1)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=sorc(i,k,1)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=                 tt(i,l)*sorc(i,l,1)+                            &
         (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +                     &
         tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                          &
         (sorc(i,lp1,1)-sorc(i,l,1))
#ifdef CLR
     gxcts0(i)=gxcts(i)
#endif
     gxcts(i)=cldfac(i,lp1,1)*gxcts(i)
   enddo
!
!-----calculation for band 2 (combined band 2)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(2)*x(i,k)+bpcm(2)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(2)*x(i,k)+btpcm(2)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(2)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(2)*topphi(i,k))
       tt(i,k)=exp(hm1ez*fac1(i,k)/sqrt(1.+fac2(i,k)))
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,2)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,2)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,2)+                   &
         (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +                     &
         tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                          &
         (sorc(i,lp1,2)-sorc(i,l,2)))
   enddo
#ifdef CLR
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+               tt(i,l)*sorc(i,l,2)+                   &
         (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +                     &
         tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                          &
         (sorc(i,lp1,2)-sorc(i,l,2))
   enddo
#endif
!
!-----calculation for band 3 (combined band 3)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(3)*x(i,k)+bpcm(3)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(3)*x(i,k)+btpcm(3)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(3)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(3)*topphi(i,k))
       tt(i,k)=exp(hm1ez*fac1(i,k)/sqrt(1.+fac2(i,k)))
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,3)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,3)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,3)+                   &
            (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +                  &
            tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                       &
            (sorc(i,lp1,3)-sorc(i,l,3)))
   enddo
#ifdef CLR
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,3)+                                  &
            (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +                  &
            tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                       &
            (sorc(i,lp1,3)-sorc(i,l,3))
   enddo
#endif
!
!-----calculation for band 4 (combined band 4)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(4)*x(i,k)+bpcm(4)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(4)*x(i,k)+btpcm(4)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i=1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(4)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(4)*topphi(i,k))
       tt(i,k)=exp(hm1ez*fac1(i,k)/sqrt(1.+fac2(i,k)))
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,4)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,4)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,4)+                   &
            (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +                  &
            tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                       &
            (sorc(i,lp1,4)-sorc(i,l,4)))
   enddo
#ifdef CLR
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,4)+                                  &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,4)-sorc(i,l,4))
   enddo
#endif
!
!-----calculation for band 5 (combined band 5)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(5)*x(i,k)+bpcm(5)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(5)*x(i,k)+btpcm(5)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i=1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(5)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(5)*topphi(i,k))
       tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+                       &
              betacm(5)*totvo2(i,k+1)*sko2d))
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,5)* (ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,5)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,5)+                   &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,5)-sorc(i,l,5)))
   enddo
#ifdef CLR
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,5)+                                  &
            (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +                  &
            tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                       &
            (sorc(i,lp1,5)-sorc(i,l,5))
   enddo
#endif
!
!-----calculation for band 6 (combined band 6)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(6)*x(i,k)+bpcm(6)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(6)*x(i,k)+btpcm(6)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(6)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(6)*topphi(i,k))
       tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+                       &
              betacm(6)*totvo2(i,k+1)*sko2d))
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,6)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,6)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,6)+                   &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,6)-sorc(i,l,6)))
   enddo
#ifdef CLR
!
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+ tt(i,l)*sorc(i,l,6)+                                 &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,6)-sorc(i,l,6))
   enddo
#endif
!
!-----calculation for band 7 (combined band 7)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(7)*x(i,k)+bpcm(7)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(7)*x(i,k)+btpcm(7)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(7)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(7)*topphi(i,k))
       tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+                       &
              betacm(7)*totvo2(i,k+1)*sko2d))
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,7)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,7)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,7)+                   &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,7)-sorc(i,l,7)))
   enddo
#ifdef CLR
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+(tt(i,l)*sorc(i,l,7)+                                 &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,7)-sorc(i,l,7)))
   enddo
#endif
!
!-----calculation for band 8 (combined band 8)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(8)*x(i,k)+bpcm(8)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(8)*x(i,k)+btpcm(8)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(8)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(8)*topphi(i,k))
       tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+                       &
              betacm(8)*totvo2(i,k+1)*sko2d))
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,8)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,8)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,8)+                   &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,8)-sorc(i,l,8)))
   enddo
#ifdef CLR
!
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,8)+                                  &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,8)-sorc(i,l,8))
   enddo
#endif
!
!-----calculation for band 9 ( 560-670 cm-1; includes co2)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(9)*x(i,k)+bpcm(9)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(9)*x(i,k)+btpcm(9)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(9)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(9)*topphi(i,k))
       tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+                       &
              betacm(9)*totvo2(i,k+1)*sko2d))*co2sp1(i,k+1)
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,9)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,9)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,9)+                   &
            (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +                  &
            tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                       &
            (sorc(i,lp1,9)-sorc(i,l,9)))
   enddo
#ifdef CLR
!
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,9)+                                  &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,9)-sorc(i,l,9))
   enddo
#endif
!
!-----calculation for band 10 (670-800 cm-1; includes co2)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(10)*x(i,k)+bpcm(10)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(10)*x(i,k)+btpcm(10)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(10)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(10)*topphi(i,k))
       tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+                       &
              betacm(10)*totvo2(i,k+1)*sko2d))*co2sp2(i,k+1)
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,10)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,10)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,10)+                  &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,10)-sorc(i,l,10)))
   enddo
#ifdef CLR
!
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,10)+                                 &
            (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +                  &
            tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                       &
            (sorc(i,lp1,10)-sorc(i,l,10))
   enddo
#endif
!
!-----calculation for band 11 (800-900 cm-1)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(11)*x(i,k)+bpcm(11)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(11)*x(i,k)+btpcm(11)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(11)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(11)*topphi(i,k))
       tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+                       &
              betacm(11)*totvo2(i,k+1)*sko2d))
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,11)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,11)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,11)+                  &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,11)-sorc(i,l,11)))
   enddo
#ifdef CLR
!
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+(tt(i,l)*sorc(i,l,11)+                                &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,11)-sorc(i,l,11)))
   enddo
#endif
!
!-----calculation for band 12 (900-990 cm-1)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(12)*x(i,k)+bpcm(12)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(12)*x(i,k)+btpcm(12)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(12)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(12)*topphi(i,k))
       tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+                       &
              betacm(12)*totvo2(i,k+1)*sko2d))
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,12)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,12)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,12)+                  &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,12)-sorc(i,l,12)))
   enddo
#ifdef CLR
!
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,12)+                                 &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,12)-sorc(i,l,12))
   enddo
#endif
!
!-----calculation for band 13 (990-1070 cm-1; includes o3))
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(13)*x(i,k)+bpcm(13)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(13)*x(i,k)+btpcm(13)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(13)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(13)*topphi(i,k))
       tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+                       &
              betacm(13)*totvo2(i,k+1)*sko2d +to3spc(i,k)))
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,13)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,13)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,13)+                  &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,13)-sorc(i,l,13)))
   enddo
#ifdef CLR
!
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,13)+                                 &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,13)-sorc(i,l,13))
   enddo
#endif
!
!-----calculation for band 14 (1070-1200 cm-1)
!
!---obtain temperature correction (capphi,cappsi),then multiply
!   by optical path (var1,var2) to compute temperature-corrected
!   optical path and mean pressure for a layer (phitmp,psitmp)
!
   do k = 1,l
     do i = 1,ipts
       f(i,k)=h44194m2*(apcm(14)*x(i,k)+bpcm(14)*y(i,k))
       ff(i,k)=h44194m2*(atpcm(14)*x(i,k)+btpcm(14)*y(i,k))
       ag(i,k)=(h1p41819+f(i,k))*f(i,k)+one
       agg(i,k)=(h1p41819+ff(i,k))*ff(i,k)+one
       phitmp(i,k)=var1(i,k)*(((( ag(i,k)*ag(i,k))**2)**2)**2)
       psitmp(i,k)=var2(i,k)*(((( agg(i,k)*agg(i,k))**2)**2)**2)
     enddo
   enddo
!
!---obtain optical path,mean pressure from the top to the pressure
!   p(k) (topm,topphi)
!
   do i = 1,ipts
     topm(i,1)=phitmp(i,1)
     topphi(i,1)=psitmp(i,1)
   enddo
!
   do k = 2,l
     do i = 1,ipts
       topm(i,k)=topm(i,k-1)+phitmp(i,k)
       topphi(i,k)=topphi(i,k-1)+psitmp(i,k)
     enddo
   enddo
!
!---tt is the cloud-free cts transmission function
!
   do k = 1,l
     do i = 1,ipts
       fac1(i,k)=acomb(14)*topm(i,k)
       fac2(i,k)=fac1(i,k)*topm(i,k)/(bcomb(14)*topphi(i,k))
       tt(i,k)=exp(hm1ez*(fac1(i,k)/sqrt(one+fac2(i,k))+                       &
              betacm(14)*totvo2(i,k+1)*sko2d))
       ctmp(i,k+1)=tt(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp0(i,k+1)=tt(i,k)
#endif
     enddo
   enddo
!
!---excts is the cts cooling rate accumulated over frequency bands
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)+sorc(i,k,14)*(ctmp(i,k+1)-ctmp(i,k))
#ifdef CLR
       excts0(i,k)=excts0(i,k)+sorc(i,k,14)*(ctmp0(i,k+1)-ctmp0(i,k))
#endif
     enddo
   enddo
!
!---gxcts is the exact cts top flux accumulated over frequency bands
!
   do i = 1,ipts
     gxcts(i)=gxcts(i)+cldfac(i,lp1,1)*(tt(i,l)*sorc(i,l,14)+                  &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,14)-sorc(i,l,14)))
   enddo
#ifdef CLR
!
   do i = 1,ipts
     gxcts0(i)=gxcts0(i)+tt(i,l)*sorc(i,l,14)+                                 &
               (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) +               &
               tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *                    &
               (sorc(i,lp1,14)-sorc(i,l,14))
   enddo
#endif
!
!   obtain cts flux at the top by integration of heating rates and
!   using cts flux at the bottom (current value of gxcts). note
!   that the pressure quantities and conversion factors have not
!   been included either in excts or in gxcts. these cancel out, thus
!   reducing computations!
!
   do k = 1,l
     do i = 1,ipts
       gxcts(i)=gxcts(i)-excts(i,k)
#ifdef CLR
       gxcts0(i)=gxcts0(i)-excts0(i,k)
#endif
     enddo
   enddo
!
!   now scale the cooling rate (excts) by including the pressure
!   factor (delp) and the conversion factor (radcon)
!
   do k = 1,l
     do i = 1,ipts
       excts(i,k)=excts(i,k)*radcon*delp(i,k)
#ifdef CLR
       excts0(i,k)=excts0(i,k)*radcon*delp(i,k)
#endif
     enddo
   enddo
!
!---this is the end of the exact cts computations; at this point
!   excts has its appropriate value.
!
!*** compute approximate cts heating rates for 15um and 9.6 um bands
!     (ctso3)
!
   do k = 1,l
     do i = 1,ipts
       ctmp2(i,k+1)=co2sp(i,k+1)*cldfac(i,k+1,1)
       ctmp3(i,k+1)=to3sp(i,k)*cldfac(i,k+1,1)
#ifdef CLR
       ctmp20(i,k+1)=co2sp(i,k+1)
       ctmp30(i,k+1)=to3sp(i,k)
#endif
     enddo
   enddo
!
   do k = 1,l
     do i = 1,ipts
       ctso3(i,k)=radcon*delp(i,k)*                                            &
               (csour(i,k)*(ctmp2(i,k+1)-ctmp2(i,k)) +                         &
               sorc(i,k,13)*(ctmp3(i,k+1)-ctmp3(i,k)))
     enddo
   enddo
#ifdef CLR
!
   do k = 1,l
     do i = 1,ipts
       ctso30(i,k)=radcon*delp(i,k)*                                           &
               (csour(i,k)*(ctmp20(i,k+1)-ctmp20(i,k)) +                       &
               sorc(i,k,13)*(ctmp30(i,k+1)-ctmp30(i,k)))
     enddo
   enddo
#endif
!
   return
   end subroutine rad_lw_spa88
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_b10exps(ipts,m,np,dh2o,dcont,dco2,dn2o,pa,dt         &
                                 ,h2oexp,conexp,co2exp,n2oexp)
!-------------------------------------------------------------------------------
!
!   compute band3a exponentials for individual layers.
!
!---- input parameters
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of layers (np)
!  layer h2o amount for line absorption (dh2o)
!  layer h2o amount for continuum absorption (dcont)
!  layer co2 amount (dco2)
!  layer n2o amount (dn2o)
!  layer pressure (pa)
!  layer temperature minus 250k (dt)
!
!---- output parameters
!
!  exponentials for each layer (h2oexp,conexp,co2exp,n2oexp)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  m,n,np,i,k,ipts
!
!---- input parameters -----
!
   real                 ::  dh2o(m,np),dcont(m,np),dn2o(m,np)
   real                 ::  dco2(m,np),pa(m,np),dt(m,np)
!
!---- output parameters -----
!
   real                 ::  h2oexp(m,np,6),conexp(m,np,3),co2exp(m,np,6,2)
   real                 ::  n2oexp(m,np,4)
!
!---- temporary arrays -----
!
   real                 ::  xx,xx1,xx2,xx3
!
   do k = 1,np
     do i = 1,ipts   
!
!-----compute the scaled h2o-line amount for the sub-band 3a
!     table 1, chou et al. (jas, 1993)

       xx=dh2o(i,k)*(pa(i,k)/500.0)*(1.+(0.0149+6.20e-5*dt(i,k))*dt(i,k))
!
!-----six exponentials by powers of 8
!
       h2oexp(i,k,1)=exp(-xx*0.10624)
!
       xx=h2oexp(i,k,1)*h2oexp(i,k,1)
       xx=xx*xx
       h2oexp(i,k,2)=xx*xx
!
       xx=h2oexp(i,k,2)*h2oexp(i,k,2)
       xx=xx*xx
       h2oexp(i,k,3)=xx*xx
!
       xx=h2oexp(i,k,3)*h2oexp(i,k,3)
       xx=xx*xx
       h2oexp(i,k,4)=xx*xx
!
       xx=h2oexp(i,k,4)*h2oexp(i,k,4)
       xx=xx*xx
       h2oexp(i,k,5)=xx*xx
!
       xx=h2oexp(i,k,5)*h2oexp(i,k,5)
       xx=xx*xx
       h2oexp(i,k,6)=xx*xx
!
!-----compute the scaled co2 amount for the sub-band 3a
!     table 1, chou et al. (jas, 1993)

       xx=dco2(i,k)*(pa(i,k)/300.0)**0.5*                                      &
               (1.+(0.0179+1.02e-4*dt(i,k))*dt(i,k))
!
!-----six exponentials by powers of 8
!
       co2exp(i,k,1,1)=exp(-xx*2.656e-5)
!
       xx=co2exp(i,k,1,1)*co2exp(i,k,1,1)
       xx=xx*xx
       co2exp(i,k,2,1)=xx*xx
!
       xx=co2exp(i,k,2,1)*co2exp(i,k,2,1)
       xx=xx*xx
       co2exp(i,k,3,1)=xx*xx
!
       xx=co2exp(i,k,3,1)*co2exp(i,k,3,1)
       xx=xx*xx
       co2exp(i,k,4,1)=xx*xx
!
       xx=co2exp(i,k,4,1)*co2exp(i,k,4,1)
       xx=xx*xx
       co2exp(i,k,5,1)=xx*xx
!
       xx=co2exp(i,k,5,1)*co2exp(i,k,5,1)
       xx=xx*xx
       co2exp(i,k,6,1)=xx*xx
!
!-----one exponential of h2o continuum for sub-band 3a
!
       conexp(i,k,1)=exp(-dcont(i,k)*1.04995e+2)
!
!-----compute the scaled n2o amount for sub-band 3a
!
       xx=dn2o(i,k)*(1.+(1.4476e-3+3.6656e-6*dt(i,k))*dt(i,k))
!
!-----two exponential2 by powers of 58
!
       n2oexp(i,k,1)=exp(-xx*0.22646)
!
       xx=n2oexp(i,k,1)*n2oexp(i,k,1)
       xx1=xx*xx
       xx1=xx1*xx1
       xx2=xx1*xx1
       xx3=xx2*xx2
       n2oexp(i,k,2)=xx*xx1*xx2*xx3
!
     enddo
   enddo
!
   return
   end subroutine rad_lw_nasa_b10exps
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_b10kdis(ipts,m,np,k,                                 &
                                  h2oexp,conexp,co2exp,n2oexp,                 &
                                  th2o,tcon,tco2,tn2o,tran)
!-------------------------------------------------------------------------------
!
!   compute h2o (line and continuum),co2,n2o transmittances between
!   levels k1 and k2 for m x n soundings
!   using the k-distribution method with linear pressure scaling.
!
!---- input parameters
!   number of grid intervals in zonal direction (m)
!   number of grid intervals in meridional direction (n)
!   number of levels (np)
!   current level (k)
!   exponentials for h2o line absorption (h2oexp)
!   exponentials for h2o continuum absorption (conexp)
!   exponentials for co2 absorption (co2exp)
!   exponentials for n2o absorption (n2oexp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to h2o line absorption
!     for the various values of the absorption coefficient (th2o)
!   transmittance between levels k1 and k2 due to h2o continuum absorption
!     for the various values of the absorption coefficient (tcon)
!   transmittance between levels k1 and k2 due to co2 absorption
!     for the various values of the absorption coefficient (tco2)
!   transmittance between levels k1 and k2 due to n2o absorption
!     for the various values of the absorption coefficient (tn2o)
!   total transmittance (tran)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  m,n,np,k,i,ipts
!
!---- input parameters -----
!
   real                 ::  h2oexp(m,np,6),conexp(m,np,3),co2exp(m,np,6,2)
   real                 ::  n2oexp(m,np,4)
!
!---- updated parameters -----
!
   real                 ::  th2o(m,6),tcon(m,3),tco2(m,6,2),tn2o(m,4)
   real                 ::  tran(m)
!
!---- temporary arrays -----
!
   real                 ::  xx
!
!-----initialize tran
!
   do i = 1,ipts   
     tran(i)=1.0
   enddo
!
!-----for h2o line
!
   do i = 1,ipts   
     th2o(i,1)=th2o(i,1)*h2oexp(i,k,1)
     xx=   0.3153*th2o(i,1)
!
     th2o(i,2)=th2o(i,2)*h2oexp(i,k,2)
     xx=xx+0.4604*th2o(i,2)
!
     th2o(i,3)=th2o(i,3)*h2oexp(i,k,3)
     xx=xx+0.1326*th2o(i,3)
!
     th2o(i,4)=th2o(i,4)*h2oexp(i,k,4)
     xx=xx+0.0798*th2o(i,4)
!
     th2o(i,5)=th2o(i,5)*h2oexp(i,k,5)
     xx=xx+0.0119*th2o(i,5)
!
     tran(i)=tran(i)*xx
   enddo
!
!-----for h2o continuum
!
   do i = 1,ipts   
     tcon(i,1)=tcon(i,1)*conexp(i,k,1)
     tran(i)=tran(i)*tcon(i,1)
   enddo
! 
!-----for co2 
!
   do i = 1,ipts   
     tco2(i,1,1)=tco2(i,1,1)*co2exp(i,k,1,1)
     xx=    0.2673*tco2(i,1,1)
!
     tco2(i,2,1)=tco2(i,2,1)*co2exp(i,k,2,1)
     xx=xx+ 0.2201*tco2(i,2,1)
!
     tco2(i,3,1)=tco2(i,3,1)*co2exp(i,k,3,1)
     xx=xx+ 0.2106*tco2(i,3,1)
!
     tco2(i,4,1)=tco2(i,4,1)*co2exp(i,k,4,1)
     xx=xx+ 0.2409*tco2(i,4,1)
!
     tco2(i,5,1)=tco2(i,5,1)*co2exp(i,k,5,1)
     xx=xx+ 0.0196*tco2(i,5,1)
!
     tco2(i,6,1)=tco2(i,6,1)*co2exp(i,k,6,1)
     xx=xx+ 0.0415*tco2(i,6,1)
!
     tran(i)=tran(i)*xx
   enddo
!
!-----for n2o
!
   do i = 1,ipts   
     tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
     xx=   0.970831*tn2o(i,1)
!
     tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
     xx=xx+0.029169*tn2o(i,2)
!
     tran(i)=tran(i)*(xx-1.0)
   enddo
!
   return
   end subroutine rad_lw_nasa_b10kdis
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_cfcexps(ipts,ib,m,np,a1,b1,fk1,                      &
                                  a2,b2,fk2,dcfc,dt,cfcexp)
!-------------------------------------------------------------------------------
!
!   compute cfc(-11, -12, -22) exponentials for individual layers.
!
!---- input parameters
!  spectral band (ib)
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of layers (np)
!  parameters for computing the scaled cfc amounts
!             for temperature scaling (a1,b1,a2,b2)
!  the absorption coefficients for the
!     first k-distribution function due to cfcs (fk1,fk2)
!  layer cfc amounts (dcfc)
!  layer temperature minus 250k (dt)
!
!---- output parameters
!  1 exponential for each layer (cfcexp)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ib,m,n,np,i,k,ipts
!
!---- input parameters -----
!
   real                 ::  dcfc(m,np),dt(m,np)
!
!---- output parameters -----
!
   real                 ::  cfcexp(m,np)
!
!---- static data -----
!
   real                 ::  a1,b1,fk1,a2,b2,fk2
!
!---- temporary arrays -----
!
   real                 ::  xf
!
   do k = 1,np
     do i = 1,ipts   
!
!-----compute the scaled cfc amount (xf) and exponential (cfcexp)
!
       if (ib.eq.4) then
         xf=dcfc(i,k)*(1.+(a1+b1*dt(i,k))*dt(i,k))
         cfcexp(i,k)=exp(-xf*fk1)
       else
         xf=dcfc(i,k)*(1.+(a2+b2*dt(i,k))*dt(i,k))
         cfcexp(i,k)=exp(-xf*fk2)
       endif
!
     enddo
   enddo
!
   return
   end subroutine rad_lw_nasa_cfcexps
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_cfckdis(ipts,m,np,k,cfcexp,tcfc,tran)
!-------------------------------------------------------------------------------
!
!   compute cfc-(11,12,22) transmittances between levels k1 and k2 for m x n
!   soundings using the k-distribution method with linear pressure scaling.
!
!---- input parameters
!   number of grid intervals in zonal direction (m)
!   number of grid intervals in meridional direction (n)
!   number of levels (np)
!   current level (k)
!   exponentials for cfc absorption (cfcexp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to cfc absorption
!     for the various values of the absorption coefficient (tcfc)
!   total transmittance (tran)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  m,n,np,k,i,ipts
!
!---- input parameters -----
!
   real                 ::  cfcexp(m,np)
!
!---- updated parameters -----
!
   real                 ::  tcfc(m),tran(m)
!
!-----tcfc is the exp factors between levels k1 and k2. 
!
   do i = 1,ipts   
     tcfc(i)=tcfc(i)*cfcexp(i,k)
     tran(i)=tran(i)*tcfc(i)
   enddo
!
   return
   end subroutine rad_lw_nasa_cfckdis
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_ch4exps(ipts,ib,m,np,dch4,pa,dt,ch4exp)
!-------------------------------------------------------------------------------
!
!   compute ch4 exponentials for individual layers.
!
!---- input parameters
!  spectral band (ib)
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of layers (np)
!  layer pressure (pa)
!  layer ch4 amount (dch4)
!  layer temperature minus 250k (dt)
!
!---- output parameters
!  1 or 4 exponentials for each layer (ch4exp)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ib,m,n,np,i,k,ipts
!
!---- input parameters -----
!
   real                 ::  dch4(m,np),pa(m,np),dt(m,np)
!
!---- output parameters -----
!
   real                 ::  ch4exp(m,np,4)
!
!---- temporary arrays -----
!
   real                 ::  xc
!
   do k = 1,np
     do i = 1,ipts   
!
!-----four exponentials for band 6
!
       if (ib.eq.6) then
         xc=dch4(i,k)*(1.+(1.7007e-2+1.5826e-4*dt(i,k))*dt(i,k))
         ch4exp(i,k,1)=exp(-xc*4.81988e-3)
!
!-----four exponentials by powers of 12 for band 7
!
       else
         xc=dch4(i,k)*(pa(i,k)/500.0)**0.65                                    &
          *(1.+(5.9590e-4-2.2931e-6*dt(i,k))*dt(i,k))
         ch4exp(i,k,1)=exp(-xc*6.07297e-2)
!
         xc=ch4exp(i,k,1)*ch4exp(i,k,1)*ch4exp(i,k,1)
         xc=xc*xc
         ch4exp(i,k,2)=xc*xc
!   
         xc=ch4exp(i,k,2)*ch4exp(i,k,2)*ch4exp(i,k,2)
         xc=xc*xc
         ch4exp(i,k,3)=xc*xc
!   
         xc=ch4exp(i,k,3)*ch4exp(i,k,3)*ch4exp(i,k,3)
         xc=xc*xc
         ch4exp(i,k,4)=xc*xc
!
       endif

     enddo
   enddo
!
   return
   end subroutine rad_lw_nasa_ch4exps
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_ch4kdis(ipts,ib,m,np,k,ch4exp,tch4,tran)
!-------------------------------------------------------------------------------
!
!   compute ch4 transmittances between levels k1 and k2 for m x n soundings
!   using the k-distribution method with linear pressure scaling.
!
!---- input parameters
!   spectral band (ib)
!   number of grid intervals in zonal direction (m)
!   number of grid intervals in meridional direction (n)
!   number of levels (np)
!   current level (k)
!   exponentials for ch4 absorption (ch4exp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to ch4 absorption
!     for the various values of the absorption coefficient (tch4)
!   total transmittance (tran)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ib,m,n,np,k,i,ipts
!
!---- input parameters -----
!
   real                 ::  ch4exp(m,np,4)
!
!---- updated parameters -----
!
   real                 ::  tch4(m,4),tran(m)
!
!---- temporary arrays -----
!
   real                 ::  xc
!
!-----tch4 is the 2 exp factors between levels k1 and k2. 
!     xc is the total ch4 transmittance

   do i = 1,ipts   
!
!-----band 6
!
     if (ib.eq.6) then
!
       tch4(i,1)=tch4(i,1)*ch4exp(i,k,1)
       xc= tch4(i,1)
!
!-----band 7
!
     else
       tch4(i,1)=tch4(i,1)*ch4exp(i,k,1)
       xc=   0.610650*tch4(i,1)
!
       tch4(i,2)=tch4(i,2)*ch4exp(i,k,2)
       xc=xc+0.280212*tch4(i,2)
!
       tch4(i,3)=tch4(i,3)*ch4exp(i,k,3)
       xc=xc+0.107349*tch4(i,3)
!
       tch4(i,4)=tch4(i,4)*ch4exp(i,k,4)
       xc=xc+0.001789*tch4(i,4)
     endif
!
     tran(i)=tran(i)*xc
   enddo
!
   return
   end subroutine rad_lw_nasa_ch4kdis
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_co2exps(ipts,m,np,dco2,pa,dt,co2exp)
!-------------------------------------------------------------------------------
!
!   compute co2 exponentials for individual layers.
!
!---- input parameters
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of layers (np)
!  layer co2 amount (dco2)
!  layer pressure (pa)
!  layer temperature minus 250k (dt)
!
!---- output parameters
!  6 exponentials for each layer (co2exp)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  m,n,np,i,k,ipts
!
!---- input parameters -----
!
   real                 ::  dco2(m,np),pa(m,np),dt(m,np)
!
!---- output parameters -----
!
   real                 ::  co2exp(m,np,6,2)
!
!---- temporary arrays -----
!
   real                 ::  xc
!
   do k = 1,np
     do i = 1,ipts   
!
!-----compute the scaled co2 amount from eq. (27) for band-wings
!     (sub-bands 3a and 3c).
!
       xc = dco2(i,k)*(pa(i,k)/300.0)**0.5                                     &
                *(1.+(0.0182+1.07e-4*dt(i,k))*dt(i,k))
!
!-----six exponentials by powers of 8 (table 7).
!
       co2exp(i,k,1,1)=exp(-xc*2.656e-5)
!
       xc=co2exp(i,k,1,1)*co2exp(i,k,1,1)
       xc=xc*xc
       co2exp(i,k,2,1)=xc*xc
!
       xc=co2exp(i,k,2,1)*co2exp(i,k,2,1)
       xc=xc*xc
       co2exp(i,k,3,1)=xc*xc
!
       xc=co2exp(i,k,3,1)*co2exp(i,k,3,1)
       xc=xc*xc
       co2exp(i,k,4,1)=xc*xc
!
       xc=co2exp(i,k,4,1)*co2exp(i,k,4,1)
       xc=xc*xc
       co2exp(i,k,5,1)=xc*xc
!
       xc=co2exp(i,k,5,1)*co2exp(i,k,5,1)
       xc=xc*xc
       co2exp(i,k,6,1)=xc*xc
!
!-----compute the scaled co2 amount from eq. (27) for band-center
!     region (sub-band 3b).
!
       xc = dco2(i,k)*(pa(i,k)/30.0)**0.85                                     &
                *(1.+(0.0042+2.00e-5*dt(i,k))*dt(i,k))
!
       co2exp(i,k,1,2)=exp(-xc*2.656e-3)
!
       xc=co2exp(i,k,1,2)*co2exp(i,k,1,2)
       xc=xc*xc
       co2exp(i,k,2,2)=xc*xc
!
       xc=co2exp(i,k,2,2)*co2exp(i,k,2,2)
       xc=xc*xc
       co2exp(i,k,3,2)=xc*xc
!
       xc=co2exp(i,k,3,2)*co2exp(i,k,3,2)
       xc=xc*xc
       co2exp(i,k,4,2)=xc*xc
!
       xc=co2exp(i,k,4,2)*co2exp(i,k,4,2)
       xc=xc*xc
       co2exp(i,k,5,2)=xc*xc
!
       xc=co2exp(i,k,5,2)*co2exp(i,k,5,2)
       xc=xc*xc
       co2exp(i,k,6,2)=xc*xc
!
     enddo
   enddo
!
   return
   end subroutine rad_lw_nasa_co2exps
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_co2kdis(ipts,m,np,k,co2exp,tco2,tran)
!-------------------------------------------------------------------------------
!
!   compute co2 transmittances between levels k1 and k2 for m x n soundings
!   using the k-distribution method with linear pressure scaling.
!
!   computations follow eq. (34).
!
!---- input parameters
!   number of grid intervals in zonal direction (m)
!   number of grid intervals in meridional direction (n)
!   number of levels (np)
!   current level (k)
!   exponentials for co2 absorption (co2exp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to co2 absorption
!     for the various values of the absorption coefficient (tco2)
!   total transmittance (tran)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  m,n,np,k,i,ipts
!
!---- input parameters -----
!
   real                 ::  co2exp(m,np,6,2)
!
!---- updated parameters -----
!
   real                 ::  tco2(m,6,2),tran(m)
!
!---- temporary arrays -----
!
   real                 ::  xc
!
!-----tco2 is the 6 exp factors between levels k1 and k2. 
!     xc is the total co2 transmittance given by eq. (53).
!
   do i = 1,ipts   
!
!-----band-wings
!
     tco2(i,1,1)=tco2(i,1,1)*co2exp(i,k,1,1)
     xc=   0.1395 *tco2(i,1,1)
!
     tco2(i,2,1)=tco2(i,2,1)*co2exp(i,k,2,1)
     xc=xc+0.1407 *tco2(i,2,1)
!
     tco2(i,3,1)=tco2(i,3,1)*co2exp(i,k,3,1)
     xc=xc+0.1549 *tco2(i,3,1)
!
     tco2(i,4,1)=tco2(i,4,1)*co2exp(i,k,4,1)
     xc=xc+0.1357 *tco2(i,4,1)
!
     tco2(i,5,1)=tco2(i,5,1)*co2exp(i,k,5,1)
     xc=xc+0.0182 *tco2(i,5,1)
!
     tco2(i,6,1)=tco2(i,6,1)*co2exp(i,k,6,1)
     xc=xc+0.0220 *tco2(i,6,1)
!
!-----band-center region
!
     tco2(i,1,2)=tco2(i,1,2)*co2exp(i,k,1,2)
     xc=xc+0.0766 *tco2(i,1,2)
!
     tco2(i,2,2)=tco2(i,2,2)*co2exp(i,k,2,2)
     xc=xc+0.1372 *tco2(i,2,2)
!
     tco2(i,3,2)=tco2(i,3,2)*co2exp(i,k,3,2)
     xc=xc+0.1189 *tco2(i,3,2)
!
     tco2(i,4,2)=tco2(i,4,2)*co2exp(i,k,4,2)
     xc=xc+0.0335 *tco2(i,4,2)
!
     tco2(i,5,2)=tco2(i,5,2)*co2exp(i,k,5,2)
     xc=xc+0.0169 *tco2(i,5,2)
!
     tco2(i,6,2)=tco2(i,6,2)*co2exp(i,k,6,2)
     xc=xc+0.0059 *tco2(i,6,2)
!
     tran(i)=tran(i)*xc
   enddo
!
   return
   end subroutine rad_lw_nasa_co2kdis
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_column (ipts,m,np,pa,dt,sabs0,sabs,spre,stem)
!-------------------------------------------------------------------------------
!
! compute column-integrated (from top of the model atmosphere)
!     absorber amount (sabs), absorber-weighted pressure (spre) and
!     temperature (stem).
!     computations of spre and stem follows eqs. (37) and (38).
!
! input parameters
!   number of soundings in zonal direction (m)
!   number of soundings in meridional direction (n)
!   number of atmospheric layers (np)
!   layer pressure (pa)
!   layer temperature minus 250k (dt)
!   layer absorber amount (sabs0)
!
! output parameters
!   column-integrated absorber amount (sabs)
!   column absorber-weighted pressure (spre)
!   column absorber-weighted temperature (stem)
!
! units of pa and dt are mb and k, respectively.
!    units of sabs are g/cm**2 for water vapor and (cm-atm)stp for co2 and o3
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  m,n,np,i,k,ipts
!
! input parameters
!
   real                 ::  pa(m,np),dt(m,np),sabs0(m,np)
!
!  output parameters
!
   real                 ::  sabs(m,np+1),spre(m,np+1),stem(m,np+1)
!
   do i = 1,ipts
     sabs(i,1)=0.0
     spre(i,1)=0.0
     stem(i,1)=0.0
   enddo
!
   do k = 1,np
     do i = 1,ipts
       sabs(i,k+1)=sabs(i,k)+sabs0(i,k)
       spre(i,k+1)=spre(i,k)+pa(i,k)*sabs0(i,k)
       stem(i,k+1)=stem(i,k)+dt(i,k)*sabs0(i,k)
     enddo
   enddo
!
   return
   end subroutine rad_lw_nasa_column
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_comexps(ipts,ib,m,np,dcom,dt,comexp)
!-------------------------------------------------------------------------------
!
!   compute co2-minor exponentials for individual layers.
!
!---- input parameters
!  spectral band (ib)
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of layers (np)
!  layer co2 amount (dcom)
!  layer temperature minus 250k (dt)
!
!---- output parameters
!  2 exponentials for each layer (comexp)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ib,m,n,np,i,k,ipts
!
!---- input parameters -----
!
   real                 ::  dcom(m,np),dt(m,np)
!
!---- output parameters -----
!
   real                 ::  comexp(m,np,2)
!
!---- temporary arrays -----
!
   real                 ::  xc,xc1,xc2
!
   do k = 1,np
     do i = 1,ipts   
!
!-----two exponentials by powers of 60 for band 4
!
       if (ib.eq.4) then
         xc=dcom(i,k)*(1.+(3.5775e-2+4.0447e-4*dt(i,k))*dt(i,k))
         comexp(i,k,1)=exp(-xc*1.95404e-5)
!
         xc=comexp(i,k,1)*comexp(i,k,1)*comexp(i,k,1)
         xc=xc*xc
         xc1=xc*xc
         xc=xc1*xc1
         xc=xc*xc
         comexp(i,k,2)=xc*xc1
!
!-----two exponentials by powers of 44 for band 5
!
       else
         xc=dcom(i,k)*(1.+(3.4268e-2+3.7401e-4*dt(i,k))*dt(i,k))
         comexp(i,k,1)=exp(-xc*4.25830e-5)
!
         xc=comexp(i,k,1)*comexp(i,k,1)
         xc1=xc*xc
         xc2=xc1*xc1
         xc=xc2*xc2
         xc=xc*xc
         comexp(i,k,2)=xc1*xc2*xc
       endif
     enddo
   enddo
!
   return
   end subroutine rad_lw_nasa_comexps
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_comkdis(ipts,ib,m,np,k,comexp,tcom,tran)
!-------------------------------------------------------------------------------
!
!   compute co2-minor transmittances between levels k1 and k2 for m x n
!   soundings using the k-distribution method with linear pressure scaling.
!
!---- input parameters
!   spectral band (ib)
!   number of grid intervals in zonal direction (m)
!   number of grid intervals in meridional direction (n)
!   number of levels (np)
!   current level (k)
!   exponentials for co2-minor absorption (comexp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to co2-minor absorption
!     for the various values of the absorption coefficient (tcom)
!   total transmittance (tran)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ib,m,n,np,k,i,ipts
!
!---- input parameters -----
!
   real                 ::  comexp(m,np,2)
!
!---- updated parameters -----
!
   real                 ::  tcom(m,2),tran(m)
!
!---- temporary arrays -----
!
   real                 ::  xc
!
!-----tcom is the 2 exp factors between levels k1 and k2. 
!     xc is the total co2-minor transmittance
!
   do i = 1,ipts   
!
!-----band 4
!
     if (ib.eq.4) then
       tcom(i,1)=tcom(i,1)*comexp(i,k,1)
       xc=   0.972025*tcom(i,1)
       tcom(i,2)=tcom(i,2)*comexp(i,k,2)
       xc=xc+0.027975*tcom(i,2)
!
!-----band 5
!
     else
       tcom(i,1)=tcom(i,1)*comexp(i,k,1)
       xc=   0.961324*tcom(i,1)
       tcom(i,2)=tcom(i,2)*comexp(i,k,2)
       xc=xc+0.038676*tcom(i,2)
     endif
     tran(i)=tran(i)*xc
   enddo
!
   return
   end subroutine rad_lw_nasa_comkdis
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_conexps(ipts,ib,m,np,dcont,xke,ne,conexp)
!-------------------------------------------------------------------------------
!
!   compute exponentials for continuum absorption in individual layers.
!
!---- input parameters
!  spectral band (ib)
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of layers (np)
!  layer scaled water vapor amount for continuum absorption (dcont)
!  absorption coefficients for the first k-distribution function
!     due to water vapor continuum absorption (xke)
!  number of terms used in each band to compute h2o continuum
!     transmittance (ne)
!
!---- output parameters
!  1 or 3 exponentials for each layer (conexp)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ib,m,n,np,i,k,iq,ipts
!
!---- input parameters ------
!
   real                 ::  dcont(m,np)
!
!---- updated parameters -----
!
   real                 ::  conexp(m,np,3)
!
!---- static data -----
!
   integer              ::  ne(9)
   real                 ::  xke(9)
!
   do k = 1,np
     do i = 1,ipts
       conexp(i,k,1) = exp(-dcont(i,k)*xke(ib))
     enddo
   enddo
!
   if (ib .eq. 3) then
!
!-----the absorption coefficients for sub-bands 3b (iq=2) and 3a (iq=3)
!     are, respectively, double and quadruple that for sub-band 3c (iq=1)
!     (table 6).
!
     do iq = 2,3
       do k = 1,np
         do i = 1,ipts
           conexp(i,k,iq) = conexp(i,k,iq-1) *conexp(i,k,iq-1)
         enddo
       enddo
     enddo
   endif
!
   return
   end subroutine rad_lw_nasa_conexps
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_h2oexps(ipts,ib,m,np,dh2o,pa,                        &
                                  dt,xkw,aw,bw,pm,mw,h2oexp)
!-------------------------------------------------------------------------------
!
!   compute exponentials for water vapor line absorption
!   in individual layers.
!
!---- input parameters
!  spectral band (ib)
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of layers (np)
!  layer water vapor amount for line absorption (dh2o) 
!  layer pressure (pa)
!  layer temperature minus 250k (dt)
!  absorption coefficients for the first k-distribution
!     function due to h2o line absorption (xkw)
!  coefficients for the temperature and pressure scaling (aw,bw,pm)
!  ratios between neighboring absorption coefficients for
!     h2o line absorption (mw)
!
!---- output parameters
!  6 exponentials for each layer  (h2oexp)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ib,m,n,np,i,k,ik,ipts
!
!---- input parameters ------
!
   real                 ::  dh2o(m,np),pa(m,np),dt(m,np)
!
!---- output parameters -----
!
   real                 ::  h2oexp(m,np,6)
!
!---- static data -----
!
   integer              ::  mw(9)
   real                 ::  xkw(9),aw(9),bw(9),pm(9)
!
!---- temporary arrays -----
!
   real                 ::  xh,xh1
!-------------------------------------------------------------------------------
!    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
!    and bw.  therefore, h2oexp for these sub-bands are identical.
!-------------------------------------------------------------------------------
   do k = 1,np
     do i = 1,ipts   
!
!-----xh is   the scaled water vapor amount for line absorption
!     computed from (27).
!
       xh = dh2o(i,k)*(pa(i,k)/500.)**pm(ib)                                   &
           * ( 1.+(aw(ib)+bw(ib)* dt(i,k))*dt(i,k) )
!
!-----h2oexp is the water vapor transmittance of the layer (k2-1)
!     due to line absorption
!
       h2oexp(i,k,1) = exp(-xh*xkw(ib))

     enddo
   enddo
!
   do ik = 2,6
     if (mw(ib).eq.6) then
       do k = 1,np
         do i = 1,ipts   
           xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
           h2oexp(i,k,ik) = xh*xh*xh
         enddo
       enddo
     elseif (mw(ib).eq.8) then
       do k = 1,np
         do i = 1,ipts   
           xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
           xh = xh*xh
           h2oexp(i,k,ik) = xh*xh
         enddo
       enddo
     elseif (mw(ib).eq.9) then
       do k = 1,np
         do i = 1,ipts   
           xh=h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
           xh1 = xh*xh
           h2oexp(i,k,ik) = xh*xh1
         enddo
       enddo
     else
       do k = 1,np
         do i = 1,ipts   
           xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
           xh = xh*xh
           xh = xh*xh
           h2oexp(i,k,ik) = xh*xh
         enddo
       enddo
     endif
   enddo
!
   return
   end subroutine rad_lw_nasa_h2oexps
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_h2okdis(ipts,ib,m,np,k,                              &
                                  fkw,gkw,ne,h2oexp,conexp,                    &
                                  th2o,tcon,tran)
!-------------------------------------------------------------------------------
!
!   compute water vapor transmittance between levels k1 and k2 for
!   m x n soundings using the k-distribution method.
!
!   computations follow eqs. (34), (46), (50) and (52).
!
!---- input parameters
!  spectral band (ib)
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of levels (np)
!  current level (k)
!  planck-weighted k-distribution function due to
!    h2o line absorption (fkw)
!  planck-weighted k-distribution function due to
!    h2o continuum absorption (gkw)
!  number of terms used in each band to compute water vapor
!     continuum transmittance (ne)
!  exponentials for line absorption (h2oexp) 
!  exponentials for continuum absorption (conexp) 
!
!---- updated parameters
!  transmittance between levels k1 and k2 due to
!    water vapor line absorption (th2o)
!  transmittance between levels k1 and k2 due to
!    water vapor continuum absorption (tcon)
!  total transmittance (tran)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ib,m,n,np,k,i,ipts
!
!---- input parameters ------
!
   real                 ::  conexp(m,np,3),h2oexp(m,np,6)
   real                 ::  fkw(6,9),gkw(6,3)
   integer              ::  ne(9)
!
!---- updated parameters -----
!
   real                 ::  th2o(m,6),tcon(m,3),tran(m)
!
!---- temporary arrays -----
!
   real                 ::  trnth2o
!
!-----tco2 are the six exp factors between levels k1 and k2 
!     tran is the updated total transmittance between levels k1 and k2

!-----th2o is the 6 exp factors between levels k1 and k2 due to
!     h2o line absorption. 

!-----tcon is the 3 exp factors between levels k1 and k2 due to
!     h2o continuum absorption.

!-----trnth2o is the total transmittance between levels k1 and k2 due
!     to both line and continuum absorption computed from eq. (52).
!
   do i = 1,ipts   
     th2o(i,1) = th2o(i,1)*h2oexp(i,k,1)
     th2o(i,2) = th2o(i,2)*h2oexp(i,k,2)
     th2o(i,3) = th2o(i,3)*h2oexp(i,k,3)
     th2o(i,4) = th2o(i,4)*h2oexp(i,k,4)
     th2o(i,5) = th2o(i,5)*h2oexp(i,k,5)
     th2o(i,6) = th2o(i,6)*h2oexp(i,k,6)
   enddo
!
   if (ne(ib).eq.0) then
     do i = 1,ipts   
       trnth2o      =(fkw(1,ib)*th2o(i,1)                                      &
                    + fkw(2,ib)*th2o(i,2)                                      &
                    + fkw(3,ib)*th2o(i,3)                                      &
                    + fkw(4,ib)*th2o(i,4)                                      &
                    + fkw(5,ib)*th2o(i,5)                                      &
                    + fkw(6,ib)*th2o(i,6))
       tran(i)=tran(i)*trnth2o
     enddo
   elseif (ne(ib).eq.1) then
     do i = 1,ipts   
       tcon(i,1)= tcon(i,1)*conexp(i,k,1)
       trnth2o  =(fkw(1,ib)*th2o(i,1)                                          &
                + fkw(2,ib)*th2o(i,2)                                          &
                + fkw(3,ib)*th2o(i,3)                                          &
                + fkw(4,ib)*th2o(i,4)                                          &
                + fkw(5,ib)*th2o(i,5)                                          &
                + fkw(6,ib)*th2o(i,6))*tcon(i,1)
       tran(i)=tran(i)*trnth2o
     enddo
   else
     do i = 1,ipts   
       tcon(i,1)= tcon(i,1)*conexp(i,k,1)
       tcon(i,2)= tcon(i,2)*conexp(i,k,2)
       tcon(i,3)= tcon(i,3)*conexp(i,k,3)
       trnth2o  = (  gkw(1,1)*th2o(i,1)                                        &
                   + gkw(2,1)*th2o(i,2)                                        &
                   + gkw(3,1)*th2o(i,3)                                        &
                   + gkw(4,1)*th2o(i,4)                                        &
                   + gkw(5,1)*th2o(i,5)                                        &
                   + gkw(6,1)*th2o(i,6) ) * tcon(i,1)                          &
                + (  gkw(1,2)*th2o(i,1)                                        &
                   + gkw(2,2)*th2o(i,2)                                        &
                   + gkw(3,2)*th2o(i,3)                                        &
                   + gkw(4,2)*th2o(i,4)                                        &
                   + gkw(5,2)*th2o(i,5)                                        &
                   + gkw(6,2)*th2o(i,6) ) * tcon(i,2)                          &
                + (  gkw(1,3)*th2o(i,1)                                        &
                   + gkw(2,3)*th2o(i,2)                                        &
                   + gkw(3,3)*th2o(i,3)                                        &
                   + gkw(4,3)*th2o(i,4)                                        &
                   + gkw(5,3)*th2o(i,5)                                        &
                   + gkw(6,3)*th2o(i,6) ) * tcon(i,3)
       tran(i)=tran(i)*trnth2o
     enddo
   endif
!
   return
   end subroutine rad_lw_nasa_h2okdis
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_n2oexps(ipts,ib,m,np,dn2o,pa,dt,n2oexp)
!-------------------------------------------------------------------------------
!
!   compute n2o exponentials for individual layers.
!
!---- input parameters
!  spectral band (ib)
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of layers (np)
!  layer pressure (pa)
!  layer n2o amount (dn2o)
!  layer temperature minus 250k (dt)
!
!---- output parameters
!  2 or 4 exponentials for each layer (n2oexp)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ib,m,n,np,i,k,ipts
!
!---- input parameters -----
!
   real                 ::  dn2o(m,np),pa(m,np),dt(m,np)
!
!---- output parameters -----
!
   real                 ::  n2oexp(m,np,4)
!
!---- temporary arrays -----
!
   real                 ::  xc,xc1,xc2
!
   do k = 1,np
     do i = 1,ipts
!
!-----four exponential by powers of 21 for band 6
!
       if (ib.eq.6) then
         xc=dn2o(i,k)*(1.+(1.9297e-3+4.3750e-6*dt(i,k))*dt(i,k))
         n2oexp(i,k,1)=exp(-xc*5.51803e-2)
!
         xc=n2oexp(i,k,1)*n2oexp(i,k,1)*n2oexp(i,k,1)
         xc1=xc*xc
         xc2=xc1*xc1
         n2oexp(i,k,2)=xc*xc1*xc2
!
!-----four exponential by powers of 8 for band 7
!
       else
          xc=dn2o(i,k)*(pa(i,k)/500.0)**0.48                                   &
                  *(1.+(1.3804e-3+7.4838e-6*dt(i,k))*dt(i,k))
          n2oexp(i,k,1)=exp(-xc*5.20113e-2)
!
          xc=n2oexp(i,k,1)*n2oexp(i,k,1)
          xc=xc*xc
          n2oexp(i,k,2)=xc*xc
          xc=n2oexp(i,k,2)*n2oexp(i,k,2)
          xc=xc*xc
          n2oexp(i,k,3)=xc*xc
          xc=n2oexp(i,k,3)*n2oexp(i,k,3)
          xc=xc*xc
          n2oexp(i,k,4)=xc*xc
        endif
      enddo
   enddo
!
   return
   end subroutine rad_lw_nasa_n2oexps
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_n2okdis(ipts,ib,m,np,k,n2oexp,tn2o,tran)
!-------------------------------------------------------------------------------
!
!   compute n2o transmittances between levels k1 and k2 for m x n soundings
!   using the k-distribution method with linear pressure scaling.
!
!---- input parameters
!   spectral band (ib)
!   number of grid intervals in zonal direction (m)
!   number of grid intervals in meridional direction (n)
!   number of levels (np)
!   current level (k)
!   exponentials for n2o absorption (n2oexp)
!
!---- updated parameters
!   transmittance between levels k1 and k2 due to n2o absorption
!     for the various values of the absorption coefficient (tn2o)
!   total transmittance (tran)
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  ib,m,n,np,k,i,ipts
!
!---- input parameters -----
!
   real                 ::  n2oexp(m,np,4)
!
!---- updated parameters -----
!
   real                 ::  tn2o(m,4),tran(m)
!
!---- temporary arrays -----
!
   real                 ::  xc
!
!-----tn2o is the 2 exp factors between levels k1 and k2. 
!     xc is the total n2o transmittance
!
   do i = 1,ipts   
!
!-----band 6
!
     if (ib.eq.6) then
       tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
       xc=   0.940414*tn2o(i,1)
!
       tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
       xc=xc+0.059586*tn2o(i,2)
!
!-----band 7
!
     else
       tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
       xc=   0.561961*tn2o(i,1)
!
       tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
       xc=xc+0.138707*tn2o(i,2)
!
       tn2o(i,3)=tn2o(i,3)*n2oexp(i,k,3)
       xc=xc+0.240670*tn2o(i,3)
!
       tn2o(i,4)=tn2o(i,4)*n2oexp(i,k,4)
       xc=xc+0.058662*tn2o(i,4)
     endif
     tran(i)=tran(i)*xc
   enddo
!
   return
   end subroutine rad_lw_nasa_n2okdis
!
!-------------------------------------------------------------------------------
   subroutine rad_lw_nasa_tablup(ipts,k1,k2,m,np,nx,nh,nt,                     &
                                 sabs,spre,stem,w1,p1,                         &
                                 dwe,dpe,coef1,coef2,coef3,tran)
!-------------------------------------------------------------------------------
!
!   compute water vapor, co2, and o3 transmittances between levels k1 and k2
!   using table look-up for m x n soundings.
!
!   calculations follow eq. (40) of chou and suarez (1995)
!
!---- input ---------------------
!  indices for pressure levels (k1 and k2)
!  number of grid intervals in zonal direction (m)
!  number of grid intervals in meridional direction (n)
!  number of atmospheric layers (np)
!  number of pressure intervals in the table (nx)
!  number of absorber amount intervals in the table (nh)
!  number of tables copied (nt)
!  column-integrated absorber amount (sabs)
!  column absorber amount-weighted pressure (spre)
!  column absorber amount-weighted temperature (stem)
!  first value of absorber amount (log10) in the table (w1) 
!  first value of pressure (log10) in the table (p1) 
!  size of the interval of absorber amount (log10) in the table (dwe)
!  size of the interval of pressure (log10) in the table (dpe)
!  pre-computed coefficients (coef1, coef2, and coef3)
!
!---- updated ---------------------
!  transmittance (tran)
!
!  note:
!   (1) units of sabs are g/cm**2 for water vapor and (cm-atm)stp for co2 and o3.
!   (2) units of spre and stem are, respectively, mb and k.
!   (3) there are nt identical copies of the tables (coef1, coef2, and
!       coef3).  the prupose of using the multiple copies of tables is
!       to increase the speed in parallel (vectorized) computations.
!       if such advantage does not exist, nt can be set to 1.
!       [06/19/00] fixed nt=1
!   
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  k1,k2,m,n,np,nx,nh,nt,i,k,ipts
!
!---- input parameters -----
!
   real                 ::  w1,p1,dwe,dpe
   real                 ::  sabs(m,np+1),spre(m,np+1),stem(m,np+1)
   real                 ::  coef1(nx,nh,nt),coef2(nx,nh,nt),coef3(nx,nh,nt)
!
!---- update parameter -----
!
   real                 ::  tran(m)
!
!---- temporary variables -----
!
   real                 ::  x1,x2,x3,we,pe,fw,fp,fm,pa,pb
   real                 ::  pc,ax,ba,bb,t1,ca,cb,t2
   integer              ::  iw,ip,nn
!
   do i = 1,ipts   
     nn=mod(i,nt)+1
!
     x1=sabs(i,k2)-sabs(i,k1)
     x2=(spre(i,k2)-spre(i,k1))/x1
     x3=(stem(i,k2)-stem(i,k1))/x1
!
     we=(log10(x1)-w1)/dwe
     pe=(log10(x2)-p1)/dpe
!
     we=max(we,w1-2.*dwe)
     pe=max(pe,p1)
!
     iw=int(we+1.5)
     ip=int(pe+1.5)
!
     iw=min(iw,nh-1)
     iw=max(iw, 2)
!
     ip=min(ip,nx-1)
     ip=max(ip, 1)
!
     fw=we-float(iw-1)
     fp=pe-float(ip-1)
!
!-----linear interpolation in pressure
!
     pa = coef1(ip,iw-1,nn)*(1.-fp)+coef1(ip+1,iw-1,nn)*fp
     pb = coef1(ip,iw,  nn)*(1.-fp)+coef1(ip+1,iw,  nn)*fp
     pc = coef1(ip,iw+1,nn)*(1.-fp)+coef1(ip+1,iw+1,nn)*fp
!
!-----quadratic interpolation in absorber amount for coef1
!
     ax = (-pa*(1.-fw)+pc*(1.+fw)) *fw*0.5 + pb*(1.-fw*fw)
!
!-----linear interpolation in absorber amount for coef2 and coef3
!
     ba = coef2(ip,iw,  nn)*(1.-fp)+coef2(ip+1,iw,  nn)*fp
     bb = coef2(ip,iw+1,nn)*(1.-fp)+coef2(ip+1,iw+1,nn)*fp
     t1 = ba*(1.-fw) + bb*fw
!
     ca = coef3(ip,iw,  nn)*(1.-fp)+coef3(ip+1,iw,  nn)*fp
     cb = coef3(ip,iw+1,nn)*(1.-fp)+coef3(ip+1,iw+1,nn)*fp
     t2 = ca*(1.-fw) + cb*fw
!
!-----update the total transmittance between levels k1 and k2
!
     tran(i)= (ax + (t1+t2*x3) * x3)*tran(i)
   enddo
!
   return
   end subroutine rad_lw_nasa_tablup
!
!-------------------------------------------------------------------------------
   end module funct_rad_lw
